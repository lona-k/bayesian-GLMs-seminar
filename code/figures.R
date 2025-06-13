
source("functions.R")


# Plot all priors together

flat_density  <- function(x) dnorm(x, 0, 1e3)  # ~uniform over [-5,5]
ridge_density <- function(x) dnorm(x, 0, 1)  # N(0,1)
laplace_density <- function(x, b = 1) 1 / (2 * b) * exp(-abs(x) / b)  # Laplace(0,1)

x_grid <- seq(-5, 5, length.out = 1000)

prior_df <- tibble(
  beta   = x_grid,
  Flat   = flat_density(x_grid),
  Ridge  = ridge_density(x_grid),
  Lasso  = laplace_density(x_grid, b = 1)
)

prior_long <- prior_df %>%
  pivot_longer(-beta, names_to = "prior", values_to = "density")

pp_priors <- ggplot(prior_long, aes(x = beta, y = density, color = prior)) +
  geom_line(size = 1) +
  scale_color_brewer("Prior", palette = "Set2") +
  labs(
    title = "Priors on regression coefficients",
    x     = expression(theta[j]),
    y     = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
