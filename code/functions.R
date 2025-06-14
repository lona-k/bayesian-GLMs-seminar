# Packages ####

library(ggplot2)
library(tidyverse)
library(patchwork)

library(dplyr)
library(purrr)
library(checkmate)

library(bayesreg)  # easy Baysian Regr/Classif
library(brms)  # more complicated stuff. Necessary for custom prior

library(numDeriv)
library(mvtnorm)
library(MCMCpack)  # MCMC for Bayesian GLMs
library(MCMCglmm)  # more MCMC
library(INLA)  # Laplace Approximation. Installed with install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE) on MAC


# Regression and Regularization ####

## data generation ####
# function to simulate data for the regularization experiment
# generates a data set with *standardized* covariates
# input: size of data n, covariate correlation matrix, true coefficient theta_true,
#       GLM family, error variance (if > 1 -> noisy data)
# returns: data frame with simulated data
data_regularization <- function(
    n = 100,  # data size
    Sigma,  # correlation between covariates
    theta_true,  # true parameter vector
    family = "gaussian",
    sigma2 = 10  # "error" variance, i.e. y ~ N(X * theta, sigma * I)
) {

  assert_choice(family, c("gaussian", "binomial"))
  assert_count(n)
  assert_numeric(theta_true)

  p <- length(theta_true) - 1
  assert_true(all(dim(Sigma) == c(p, p)))

  # draw X
  X <- rmvnorm(n, mean = rep(0, p), sigma = Sigma) # X ~ N(0, Sigma)
  eta <- theta_true[1] + X %*% theta_true[-1]

  # simulate y from X
  y <- switch(family,
              gaussian = mvrnorm(n = 1, mu = eta, Sigma = sigma2 * diag(n)),
              binomial = rbinom(n, size=1, prob = plogis(eta)))

  dat <- data.frame(x = X, y = y)
  colnames(dat) <- c(paste0("x_", seq_len(p)), "y")
  dat
}

## fitting models ####
# function to fit models:
# flat prior, lasso, ridge
fit_models <- function(dat, family = "gaussian",
                       sample = 2000, burnin = 500, thin = 5) {

  assert_choice(family, c("gaussian", "binomial"))
  assert_count(sample)
  assert_count(burnin)
  assert_count(thin)


  model_fam_flat <- if (family == "gaussian") gaussian() else bernoulli(link = "logit")
  model_fam <- ifelse(family == "gaussian", "normal", "logistic")
  form <- as.formula(paste("y ~", paste(names(dat)[-ncol(dat)], collapse = "+")))

  if (family == "binomial") dat$y <- factor(dat$y)

  # priors
  priors_flat <- if (family == "gaussian") {
    c(prior(normal(0, 1e6), class = b)
      + prior(inv_gamma(0.001, 0.001), class = sigma))  # uninformative
  } else {
    prior(normal(0, 1e6), class = b)  # uninformative
  }

  # fit models
  flat <- brm(y ~ ., data = dat, family = model_fam_flat,
              prior = priors_flat, sample_prior = "yes",
              chains = 1,
              iter = sample,
              warmup = burnin,
              thin = thin)

  ridge <- bayesreg(form, data = dat, model = model_fam, prior = "ridge",
                    n.samples = sample, burnin = burnin, thin = thin)

  lasso <- bayesreg(form, data = dat, model = model_fam, prior = "lasso",
                    n.samples = sample, burnin = burnin, thin = thin)

  list(flat = flat, ridge = ridge, lasso = lasso)
}


# Approximate inference ####

# function to simulate data for the approximate inference experiment
# input: size of data n, true coefficient theta_true, GLM family
# returns: data frame with simulated data
data_sim <- function(
    n = 1000,
    theta_true = c(0.5, 2), # can be any dimension
    family = "gaussian") {

  assert_count(n)
  assert_numeric(theta_true)
  assert_choice(family, c("gaussian", "binomial"))

  dim <- length(theta_true) - 1
  x <- rmvnorm(n, mean = rep(0, dim), sigma = diag(dim))  # x_i sim N(0, 1)
  eta   <- theta_true[[1]] + x %*% theta_true[2:length(theta_true)]  # linear predictor

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

mcmc_sim_null <- function(dat, mcmc_iters = 10000, burnin = 1000,
                     prior = c(mean = 0, var = 100), sigma2 = 100,
                     family = "gaussian") {
  assert_choice(family, c("gaussian", "binomial"))
  assert_number(sigma2)
  assert_numeric(prior, len = 2)
  assert_number(mcmc_iters)
  assert_number(burnin)

  t0 <- proc.time()
  form <- sprintf("y ~ %s", paste0(colnames(dat)[-ncol(dat)], collapse = " + "))

  mcmc_fit <- switch(
    family,
    binomial = MCMClogit(
      as.formula(form),
      data = dat,
      b0 = rep(prior["mean"], ncol(dat)), B0 = diag(1/prior["var"], ncol(dat)),
      mcmc = mcmc_iters,
      # beta.start default: MLE estimate of theta
      burnin = burnin),

    # actually uses Gibbs sampling because of the non-informative IG prior
    gaussian = MCMCregress(
      as.formula(form),
      data = dat,
      # prior.density default: multivariate normal prior
      b0 = rep(prior["mean"], ncol(dat)), B0 = diag(1/prior["var"], ncol(dat)),
      sigma.mu = sqrt(sigma2), sigma.var = 0.001,  # fix variance at sigma2
      mcmc = mcmc_iters,
      # beta.start default: MLE estimate of theta
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


# function for glm with MCMC via MCMCglmm
# input: data dat, prior parameters for model (always Gaussian prior), model family, fixed variance
#        mcmc iterations, burin length
# output: list with
# - time elapsed
# - posterior mean and sd (for each coefficient)
# alternative mcmc function with MCMCglmm
mcmc_sim <- function(dat, mcmc_iters = 5000, burnin = 500,
                     prior = c(mean = 0, var = 100),
                     sigma2 = 100,
                     family = "gaussian") {

  assert_choice(family, c("gaussian", "binomial"))
  assert_number(sigma2)
  assert_numeric(prior, len = 2)
  assert_number(mcmc_iters)
  assert_number(burnin)

  form <- sprintf("y ~ %s", paste0(colnames(dat)[-ncol(dat)], collapse = " + "))
  dim <- ncol(dat)
  family <- ifelse(family == "binomial", "categorical", "gaussian")

  prior_lst <- list(
    B = list(mu = rep(prior["mean"], dim),  # theta ~ N(theta_mean, theta_var I)
             V  = diag(prior["var"], dim)),
    R = list(V  = sigma2,  # residual variance fixed at sigma2
             nu = 0)   # zero df -> point‐mass at V
  )

  t0 <- proc.time()

  mcmc_fit <- MCMCglmm(
    as.formula(form),
    data = dat,
    family = family,
    prior = prior_lst,
    nitt  = mcmc_iters,
    burnin = burnin,
    verbose = FALSE
  )

  mh_time <- (proc.time() - t0)["elapsed"]

  mh_mean <- summary(mcmc_fit)$solutions[, 1]
  mh_sd <- apply(mcmc_fit$Sol, 2, sd)

  # ess <- summary(mcmc_fit)$solutions[, 4]

  list(time = mh_time, mu_post = mh_mean, sigma_post = mh_sd)
}


# function for glm with Laplace Approximation via bayesglm()
# input: data dat, prior parameters for model (always Gaussian prior), model family, fixed variance
# output: list with
  # - time elapsed
  # - posterior mean, sd, and mode (for each coefficient)
la_sim <- function(dat,
                   prior = c(mean = 0, var = 100),
                   sigma2 = 100,  # FIXED variance of y | theta ~ N(theta, sigma2)
                   family = "gaussian") {

  assert_choice(family, c("gaussian", "binomial"))
  assert_number(sigma2)
  assert_numeric(prior, len = 2)

  form <- as.formula(paste("y ~", paste(colnames(dat)[-ncol(dat)], collapse = " + ")))
  prior_prec <- 1 / (prior["var"])
  fixed_tau <- 1 / sigma2  # tau = 1/sigma^2

  # FIX noise precision at 'fixed_tau', i.e. remove it from inference for Gaussian LM
  cfam <- if (family == "gaussian") {
    list(hyper = list(prec = list(initial = log(1/sigma2), fixed = TRUE)))
  } else list()

  t0 <- proc.time()

  la_fit <- inla(
    formula = form,
    family = family,
    data = dat,
    control.family = cfam,
    control.fixed = list(
      mean = prior["mean"],
      prec = prior_prec
    ),
    control.inla = list(
      int.strategy = "eb",  # only use the posterior-mode of theta (no grid or CCD)
      strategy     = "laplace"  # full Laplace for the conditional field (not the default “simplified”)
    )
  )

  la_time <- (proc.time() - t0)["elapsed"]

  # mean, sd, and mode for the distribution of theta
  la_mean <- summary(la_fit)$fixed[, "mean"]
  la_sd <- summary(la_fit)$fixed[, "sd"]
  la_mode <- summary(la_fit)$fixed[, "mode"]

  list(time = la_time, mu_post = la_mean, sigma_post = la_sd, mode_post = la_mode)
}



# Plots ####
