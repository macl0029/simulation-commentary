rm(list=ls())
## ------------------------------------------------------------
## Simulation for misclassified binary outcome with Se = 1, Sp < 1
## ------------------------------------------------------------
set.seed(20260424)
## Simulation parameters
R      <- 10^6    # number of Monte Carlo replicates
n      <- 100     # sample size per replicate
p_true <- 0.20    # true prevalence P(X = 1)

## Specificity distribution used in the simulation
a_S   <- 20
b_S   <- 5
sp_fix <- a_S/(a_S+b_S)    # fixed plug-in value for methods 2 and 3

## Beta mean/variance for specificity
mean_sp_beta <- sp_fix
var_sp_beta  <- (a_S * b_S) / (((a_S + b_S)^2) * (a_S + b_S + 1))


## ------------------------------------------------------------
## functions

## Misclassification correction when Se = 1
##  p = (q + sp - 1) / sp
adj_p <- function(qhat, sp) {
  (qhat + sp - 1) / sp
}

## Reported variance for the crude estimator
v_crude <- function(qhat, n) {
  qhat * (1 - qhat) / n
}

## Variance for adjusted estimator; procedure 2; when sp is treated as fixed
v_adj2 <- function(qhat, sp, n) {
  (1 / sp)^2 * v_crude(qhat, n)
}

## Variance for adjusted estimator; procedure 3; when sp is treated as random
v_adj3 <- function(qhat, sp, n, var_sp) {
  var_q <- v_crude(qhat, n)
  d_q   <- 1 / sp
  d_sp  <- (1 - qhat) / (sp^2)
  
  d_q^2 * var_q + d_sp^2 * var_sp
}

## ------------------------------------------------------------
## Storage
## ------------------------------------------------------------

method_names <- c(
  "1_crude",
  "2_adj2_fixed_sp",
  "3_adj3_var_sp"
)

est_ran<-est_fix  <- matrix(NA_real_, nrow = R, ncol = length(method_names))
vhat_ran<-vhat_fix <- matrix(NA_real_, nrow = R, ncol = length(method_names))

colnames(est_ran)  <- method_names
colnames(vhat_ran) <- method_names
colnames(est_fix)  <- method_names
colnames(vhat_fix) <- method_names

sp_draw_store_fix <- numeric(R)
qhat_store_fix    <- numeric(R)
sp_draw_store_ran <- numeric(R)
qhat_store_ran    <- numeric(R)

## ------------------------------------------------------------
## Simulation loop
## ------------------------------------------------------------

for (r in 1:R) {
  sp_r_fix<-mean_sp_beta
  sp_r_ran <- rbeta(1, a_S, b_S)
  
  

  ## Generate W directly given Sp - FIXED
  pw_fix <- p_true+(1-sp_r_fix)*(1-p_true)
  W_fix <- rbinom(n, size=1, prob = pw_fix)
  ## Generate W directly given Sp - Rand
  pw_ran <- p_true+(1-sp_r_ran)*(1-p_true)
  W_ran <- rbinom(n, size=1, prob = pw_ran)
  
  ## Observed proportion of W = 1
  qhat_fix <- mean(W_fix)
  qhat_ran <- mean(W_ran)
  
  ## Store
  sp_draw_store_fix[r] <- sp_r_fix
  qhat_store_fix[r]    <- qhat_fix
  sp_draw_store_ran[r] <- sp_r_ran
  qhat_store_ran[r]    <- qhat_ran

    ## Method 1: crude
  est_fix[r, 1]  <- qhat_fix
  vhat_fix[r, 1] <- v_crude(qhat_fix, n)
  est_ran[r, 1]  <- qhat_ran
  vhat_ran[r, 1] <- v_crude(qhat_ran, n)
  
  ## Method 2: adjusted, plug in fixed sp 
  est_fix[r, 2]  <- adj_p(qhat_fix, sp_fix)
  vhat_fix[r, 2] <- v_adj2(qhat_fix, sp_fix, n)
  est_ran[r, 2]  <- adj_p(qhat_ran, sp_fix)
  vhat_ran[r, 2] <- v_adj2(qhat_ran, sp_fix, n)
  
  ## Method 3: same point estimator as 2, but var for sp
  est_fix[r, 3]  <- est_fix[r, 2]
  vhat_fix[r, 3] <- v_adj3(qhat_fix, sp_fix, n, var_sp_beta)
  est_ran[r, 3]  <- est_ran[r, 2]
  vhat_ran[r, 3] <- v_adj3(qhat_ran, sp_fix, n, var_sp_beta)
}

## ------------------------------------------------------------
## Summaries
## ------------------------------------------------------------

summarize_methods <- function(est_mat, vhat_mat, p_true) {
  
  est_mean <- colMeans(est_mat)
  emp_var  <- apply(est_mat, 2, var)
  mean_v   <- colMeans(vhat_mat)
  
  zcrit <- qnorm(0.975)
  lower <- est_mat - zcrit * sqrt(vhat_mat)
  upper <- est_mat + zcrit * sqrt(vhat_mat)
  
  out <- data.frame(
    method = colnames(est_mat),
    mean_estimate      = est_mean,
    bias               = est_mean - p_true,
    empirical_variance = emp_var,
    mean_reported_var  = mean_v,
    var_ratio          = mean_v / emp_var,
    mse                = colMeans((est_mat - p_true)^2),
    coverage_95        = colMeans(lower <= p_true & upper >= p_true),
    mean_CI_width      = colMeans(upper - lower),
    row.names = NULL
  )
  
  out
}

results_fix <- summarize_methods(est_fix, vhat_fix, p_true)
results_ran <- summarize_methods(est_ran, vhat_ran, p_true)

print(results_fix, digits = 4)
print(results_ran, digits = 4)

## Optional: inspect the simulated specificity draws
summary(sp_draw_store_fix)
summary(sp_draw_store_ran)

## Write the results table to a CSV file in your working directory
write.csv(results_fix, file = "/Users/rmaclehose/Dropbox/Projects - Dropbox/completeness/manuscript/completeness/results/simulation_fix.csv", row.names = FALSE)
write.csv(results_ran, file = "/Users/rmaclehose/Dropbox/Projects - Dropbox/completeness/manuscript/completeness/results/simulation_ran.csv", row.names = FALSE)


