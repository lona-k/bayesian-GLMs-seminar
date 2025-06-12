# Packages ####

library(ggplot2)
library(dplyr)
library(arm)  # bayesglm() for LA
library(MCMCpack)  # MCMC for Bayesian GLMs
library(checkmate)

# Regression and Regularization ####
data_regr <- function() {

}
data_class <- function() {

}


# Approximate inference ####

# function to simulate data for the approximate inference experiment
# input: size of data n, true coefficient beta_true, GLM family
# returns: data frame with simulated data
data_sim <- function(
    n = 1000,
    beta_true = c(0.5, 2), # can be any dimension
    family = "gaussian") {

  assert_count(n)
  assert_numeric(beta_true)
  assert_choice(family, c("gaussian", "binomial"))

  dim <- length(beta_true) - 1
  x <- rmvnorm(n, mean = rep(0, dim), sigma = diag(dim))  # x_i sim N(0, 1)
  eta   <- beta_true[[1]] + x %*% beta_true[2:length(beta_true)]  # linear predictor

  # mean with link function
  mu   <- switch(family,
                 binomial = plogis(eta),
                 gaussian = eta
  )

  # simulate y kind of randomly
  y   <- switch(family,
                binomial = rbinom(n, size=1, prob=mu),
                gaussian = rnorm(n, mean=mu, sd=1)
  )

  dat <- data.frame(x = x, y = y)
  colnames(dat) <- c(paste0("x_", seq_len(ncol(x))), "y")
  dat
}

# function for glm with Laplace Approximation via bayesglm()
# input: data dat, prior for model (always gaussian), model family
# output: list with
  # - time elapsed
  # - MAP estimate for coefficients (MAP of posterior mean)
  # - MAP estimate for covariance (MAP of posterior variance)

la_sim <- function(dat,
                   prior = c(mean = 0, scale = 100), family = "gaussian") {
  assert_choice(family, c("gaussian", "binomial"))
  assert_numeric(prior, length = 2)

  t0  <- proc.time()
  form <- sprintf("y ~ %s", paste0(colnames(dat)[-ncol(dat)], collapse = " + "))
  fit_la <- bayesglm(
    as.formula(form),
    family = family,
    data   = dat,
    prior.mean  = prior["mean"],
    prior.scale = prior["scale"],
    prior.df    = Inf # effectively Gaussian
  )
  la_time <- (proc.time() - t0)["elapsed"]

  theta_map <- coef(fit_la)
  sigma_map <- vcov(fit_la)

  list(time = la_time, mu_post = theta_map, sigma_post = sigma_map)
}

mcmc_sim <- function(dat, mcmc_iters = 10000, burnin = 1000,
                     prior = c(mean = 0, var = 100), family = "gaussian") {

  t0 <- proc.time()
  form <- sprintf("y ~ %s", paste0(colnames(dat)[-ncol(dat)], collapse = " + "))

  mcmc_fit <- switch(family,
     binomial = MCMClogit(as.formula(form),
       data = dat,
       b0 = rep(prior["mean"], ncol(dat)), B0 = diag(1/prior["var"], ncol(dat)),
       mcmc = mcmc_iters,
       # beta.start default: MLE estimate of beta
       burnin = burnin),

     # actually uses gibbs sampling because of the NIG prior
     gaussian = MCMCregress(as.formula(form),
       data = dat,
       # prior.density default: multivariate normal prior
       b0 = rep(prior["mean"], ncol(dat)), B0 = diag(1/prior["var"], ncol(dat)),
       sigma.var = 1e6, # flat variance prior -> effectively only Gaussian prior
       mcmc = mcmc_iters,
       # beta.start default: MLE estimate of beta
       burnin = burnin,
       marginal.likelihood = "Laplace")
  )

  mh_time <- (proc.time() - t0)["elapsed"]

  mh_mean <- summary(mcmc_fit)$statistics[, 1]
  mh_sd <- summary(mcmc_fit)$statistics[, 2]
  # mh_quantiles <- summary(mcmc_fit)$quantiles
  # plot(posterior)  #> plots markov chain
  # summary(posterior)  #> gives model summary

  # acceptance rate is difficult to record, can only be done with hacks,
  # see https://stackoverflow.com/questions/46473147/is-there-a-way-to-save-the-acceptance-rate-of-the-mcmc-algorithm-used-by-mcmclog

  list(time = mh_time, mu_post = mh_mean, sigma_post = mh_sd)
}


# Plots ####
