
source("functions.R")

library(tidybayes)
library(modelr)

# synthetic data for all examples
set.seed(42)
N_total <- 10
beta0_true <- 0
beta1_true <- 2
ssd_true <- 5

x_all <- runif(N_total, -5, 5)
y_all <- beta0_true + beta1_true * x_all + rnorm(N_total, sd = ssd_true)

df <- data.frame(x = x_all, y = y_all)

# prior plots ####

# Plot all priors together

flat_density  <- function(x) dnorm(x, 0, 1e2)  # ~uniform over [-5,5]
ridge_density <- function(x) dnorm(x, 0, 1)  # N(0,1)
laplace_density <- function(x, b = 1) 1 / (2 * b) * exp(-abs(x) / b)  # Laplace(0,1)

x_grid <- seq(-5, 5, length.out = 1000)

prior_df <- tibble(
  beta   = x_grid,
  uninformativ = flat_density(x_grid),
  Ridge  = ridge_density(x_grid),
  Lasso  = laplace_density(x_grid, b = 1)
)

prior_long <- prior_df %>%
  pivot_longer(-beta, names_to = "prior", values_to = "density")

pp_priors <- ggplot(prior_long, aes(x = beta, y = density, color = prior)) +
  geom_line(linewidth = 0.8) +
  scale_color_brewer("Priori", palette = "Set2") +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
  labs(
    title = "",
    x     = expression(theta),
    y     = "Dichte"
  ) +
  theme_minimal() +
  theme(legend.position = "inside",
    legend.position.inside = c(0.95, 0.95),
    legend.justification = c("right","top"),
    legend.background = element_rect(
      fill   = alpha("white", 0.8),
      color  = "black",
      linewidth   = 0.2
    ))

# Plot effect of priors on posteriors

# Likelihood grid
beta_grid <- seq(-5, 5, length.out = 1001)
# log-likelihood up to constant
loglik <- sapply(beta_grid, function(b) {
  sum(dnorm(y_all, mean = b * x_all, sd = ssd_true, log = TRUE))
})

# exponentiate and re-center for numerical stability
lik <- exp(loglik - max(loglik))

# unnormalized posteriors, then normalize by trapezoid rule
step <- diff(beta_grid)[1]

post_unnorm <- tibble(
  beta = beta_grid,
  uninformativ = flat_density(beta_grid) * lik,
  Ridge       = ridge_density(beta_grid)*lik,
  Lasso       = laplace_density(beta_grid)*lik
)

post_norm <- post_unnorm %>%
  mutate_at(vars(-beta), ~ .x / sum(.x * step))

post_long <- post_norm %>%
  pivot_longer(-beta, names_to = "Prior", values_to = "Density")

pp_posteriors <- ggplot(post_long, aes(x = beta, y = Density, color = Prior)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = beta1_true, linetype = "dashed", color = "black") +
  scale_color_brewer("Prior", palette = "Set2") +
  scale_x_continuous(limits = c(-2, 5), breaks = c(-2, 0, 2, 4)) +
  labs(
    title = "Posterior Verteilungen für θ",
    x     = expression(theta),
    y     = "Dichte"
  ) +
  theme_minimal() +
  theme(legend.position = "right")


# build everything together and arrange it for the slides
pp_model <- ggplot() +
  # first line
  annotate(
    "text", x = 0.5, y = 0.70,
    label = "y[i] == theta * x[i] + epsilon[i]",
    parse = TRUE, size = 4
  ) +
  # second line
  annotate(
    "text", x = 0.5, y = 0.55,
    label = "theta == 2",
    parse = TRUE, size = 4
  ) +
  # third line
  annotate(
    "text", x = 0.5, y = 0.40,
    label = "epsilon[i] %~% N(0,5)",
    parse = TRUE, size = 4
  ) +
  # left arrow
  annotate(
    "segment",
    x     = 0.05, xend = 0.8,
    y     = 0.20, yend = 0.20,
    arrow = arrow(length = unit(0.3, "cm"))
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void()


ggsave("../figures/plot_priors.png", pp_priors, width = 4, height = 4)

p_update <- ((pp_priors + scale_y_continuous(limits = c(0, 0.8)) + theme(legend.position = "none")) +
               pp_model + pp_posteriors) +
  plot_layout(widths = c(11, 5, 14))


ggsave("../figures/plot_posteriors.png", pp_posteriors, width = 5.5, height = 3)

ggsave("../figures/plot_update.png", p_update, width = 11, height = 5)


# general bayes plot ####

fit <- brm(
  formula = y ~ x,
  data    = df,
  family  = gaussian(),
  prior   = prior(normal(0,5), class="b"),
  chains  = 2, iter = 2000, refresh = 0, seed = 123
)

# plot prior N(0, 5)

# fit + other possibilities

ex_draws <- df %>%
  add_epred_draws(fit, ndraws = 20) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(aes(y = .epred, group = paste(.draw)), alpha = .3) +
  geom_point(data = df) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "mögliche Modelle",
        caption = "mit Ziehungen aus der Parameter Posterior") +
  theme_minimal()

# fit + CI
ex_model <- df %>%
  add_epred_draws(fit) %>%
  ggplot(aes(x = x, y = y)) +
  stat_lineribbon(aes(y = .epred), alpha = 0.3) +
  geom_point(data = df) +
  scale_fill_brewer(palette = "Greys", name = "Credibility\nInterval") +
  labs(title = "MAP-Modell mit Credibility Intervallen") +
  theme_minimal()

# prediction (of training data) + CIs form PPD

ex_ppd <- df %>%
  add_predicted_draws(fit) %>%
  ggplot(aes(x = x, y = y)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50), alpha = 0.3) +
  geom_point(data = df) +
  scale_fill_brewer(palette = "Greys", name = "Credibility\nInterval") +
  labs(title = "Vorhersagen für Trainingsdaten") +
  theme_minimal()


plot_layout(guides = "collect", axes = "collect")

# parameter prior
prior_draws <- data.frame(b_x = rnorm(50000, mean = 0, sd = 5))

ex_prior <- prior_draws %>%
  ggplot(aes(x = b_x)) +
  geom_density(fill = "grey", alpha = 0.6) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "red") +
  scale_x_continuous(limits = c(-15, 15)) +
  labs(
    title = "Prior Verteilung von θ",
    subtitle ="wahres θ = 2 in rot",
    x = expression(theta),
    y = "Density"
  ) +
  theme_minimal()


# parameter posterior
ex_posterior <- fit %>%
  spread_draws(b_x) %>%
  ggplot(aes(x = b_x)) +
  geom_density(alpha=0.6, fill = "grey") +
  geom_vline(xintercept = 2, linetype="dashed", color = "red") +
  scale_x_continuous(limits = c(-10, 10)) +
  labs(title="Posterior Verteilung von θ",
       x=expression(theta), y="Density"
  ) +
  theme_minimal()

pp_example <- (ex_prior | ex_posterior) / (ex_draws | ex_model) +
  plot_layout(axes = "collect")

ggsave("../figures/plot_example.png", pp_example, width = 10, height = 5)

# plot PPD ####

df_new <- data.frame(x = c(-3, 0, 5, 20))
df_sample <- df[sample(nrow(df), 5),]

ppd_sample <- df_sample %>%
  add_predicted_draws(fit, ndraws = 500) %>%
  ggplot(aes(x = .prediction, y = factor(x), fill = factor(x))) +
  stat_slab(alpha = 0.6) +
  geom_point(data = df_sample,
             aes(x = y, y = factor(x))) +
  scale_y_discrete(
    labels = function(lvls) sprintf("%.2f", as.numeric(lvls))
  ) +
  labs(
    title = "PPD für zufällig gezogene Trainingsdaten",
    x     = expression(hat(y)),
    y     = expression(x)
  ) +
  theme_minimal() +
  theme(legend.position = "none")


ppd_new <- df_new %>%
  add_predicted_draws(fit, ndraws = 500) %>%
  ggplot(aes(x = .prediction, y = factor(x), fill = factor(x))) +
  stat_slab(alpha = 0.6) +
  # geom_point(data = df_new,
  #            aes(x = y, y = factor(x))) +
  labs(
    title = "PPD für neue Daten",
    x = expression(hat(y)),
    y = expression(x)
  ) +
  theme_minimal() +
  theme(legend.position = "none")

pp_ppd <- (ex_ppd + theme(legend.position = "bottom")) + ppd_sample + ppd_new +
  plot_layout(widths = c(3, 3, 4), axes = "collect")

ggsave("../figures/plot_ppd.png", pp_ppd, width = 13, height = 5.5)


