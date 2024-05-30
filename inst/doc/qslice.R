## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library("qslice")

## ----target-------------------------------------------------------------------
## Target is a Gamma(shape = alpha, rate = 1) distribution
alpha <- 2.5
ltarget <- function(x) ifelse(x > 0.0, (alpha - 1.0) * log(x) - x, -Inf)

## ----step---------------------------------------------------------------------
## Set up MCMC
n_iter <- 1e3
x_sample <- numeric(n_iter + 1)
x_sample[1] <- 0.5  # initialize
n_eval <- 0  # count target evaluations

## Run step-and-shrink slice sampler
for (i in 2:(n_iter+1)) {
  state <- slice_stepping_out(x = x_sample[i-1], 
                              log_target = ltarget,
                              w = 2.0)
  n_eval <- n_eval + state$nEvaluations
  x_sample[i] <- state$x
}

## Check samples
n_eval / n_iter  # target evaluations per MCMC iteration
hist(x_sample, freq = FALSE, n = 30)
curve(dgamma(x, shape = alpha), col = "blue", lwd = 2, add = TRUE)

## Null hypothesis: Samples are iid Gamma(alpha, 1)
ks.test(x_sample[seq(1, n_iter, by = 5)], pgamma, shape = alpha)

## ----pseudo1------------------------------------------------------------------
pseu <- pseudo_list(family = "t", 
                    params = list(loc = 0, sc = 3, degf = 1), 
                    lb = 0)

## ----vis_pseudo1, fig.width=4.5, fig.height=3---------------------------------
utility_pseudo(pseudo = pseu, log_target = ltarget, 
               type = "function", 
               utility_type = "AUC", plot = TRUE)

## ----qslice1, fig.show='hold'-------------------------------------------------
## Set up MCMC
n_iter <- 1e3
x_sample <- psi_sample <- numeric(n_iter + 1)
x_sample[1] <- psi_sample[1] <- 0.5  # initialize
n_eval <- 0  # count target evaluations

## Run step-and-shrink slice sampler
for (i in 2:(n_iter+1)) {
  state <- slice_quantile(x = x_sample[i-1], 
                          log_target = ltarget,
                          pseudo = pseu)
  n_eval <- n_eval + state$nEvaluations
  x_sample[i] <- state$x
  psi_sample[i] <- state$u
}

## Check samples
n_eval / n_iter  # target evaluations per iteration of MCMC

hist(x_sample, freq = FALSE, n = 30)
curve(dgamma(x, shape = alpha), col = "blue", lwd = 2, add = TRUE)

hist(psi_sample, freq = FALSE, n = 30)

auc(u = psi_sample) # calculate AUC from transformed samples

## Null hypothesis: Samples are iid Gamma(alpha, 1)
ks.test(x_sample[seq(1, n_iter, by = 10)], pgamma, shape = alpha)

## ----pseudo2, fig.width=4.5, fig.height=3-------------------------------------
pseu_opt <- pseudo_opt(log_target = ltarget, 
                       type = "function",
                       family = "t", degf = c(1, 5, 20), 
                       lb = 0,
                       utility_type = "AUC", plot = TRUE)

## ----pseudo3------------------------------------------------------------------
pseu_opt <- pseudo_opt(samples = x_sample, 
                       type = "samples",
                       family = "t", degf = c(1, 5, 20), 
                       lb = 0,
                       utility_type = "AUC", plot = TRUE, nbins = 20)

names(pseu_opt)
names(pseu_opt$pseudo)
pseu_opt$pseudo$txt

## ----qslice2, fig.show='hold'-------------------------------------------------
## Set up MCMC
n_iter <- 1e3
x_sample <- psi_sample <- numeric(n_iter + 1)
x_sample[1] <- psi_sample[1] <- 0.5  # initialize
n_eval <- 0  # count target evaluations

## Run step-and-shrink slice sampler
for (i in 2:(n_iter+1)) {
  state <- slice_quantile(x = x_sample[i-1], 
                          log_target = ltarget,
                          pseudo = pseu_opt$pseudo)
  n_eval <- n_eval + state$nEvaluations
  x_sample[i] <- state$x
  psi_sample[i] <- state$u
}

## Check samples
n_eval / n_iter  # target evaluations per iteration of MCMC

hist(x_sample, freq = FALSE, n = 30)
curve(dgamma(x, shape = alpha), col = "blue", lwd = 2, add = TRUE)

hist(psi_sample, freq = FALSE, n = 30)

auc(u = psi_sample) # calculate AUC from transformed samples

## Null hypothesis: Samples are iid Gamma(alpha, 1)
ks.test(x_sample[seq(1, n_iter, by = 10)], pgamma, shape = alpha)

