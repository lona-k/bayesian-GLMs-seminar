
source("functions.R")


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
  scale_color_brewer("Prior", palette = "Set2") +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
  labs(
    title = "Prior Verteilungen für θ",
    x     = expression(theta),
    y     = "Dichte"
  ) +
  theme_minimal()

# Plot effect of priors on posteriors
# synthetic data:
set.seed(42)
N_total    <- 10
beta0_true <- 0
beta1_true <- 2
ssd_true   <- 5

x_all <- runif(N_total, -5, 5)
y_all <- beta0_true + beta1_true * x_all + rnorm(N_total, sd = ssd_true)

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


ggsave("../figures/plot_priors.png", pp_priors, width = 8, height = 5)

p_update <- ((pp_priors + scale_y_continuous(limits = c(0, 0.8)) + theme(legend.position = "none")) +
               pp_model + pp_posteriors) +
  plot_layout(widths = c(11, 5, 14))


ggsave("../figures/plot_posteriors.png", pp_posteriors, width = 5.5, height = 3)

ggsave("../figures/plot_update.png", p_update, width = 11, height = 5)

