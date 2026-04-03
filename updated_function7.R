##############################################
# BAYESIAN-HYBRID ANALYSIS
# Computes Estimates, Standard Errors, and P-values
# Across Frequentist, Bayesian, and Hybrid Approaches
##############################################

# Define function
bayes_hybrid_regression <- function(Y, X, mu1, std1) {
  # Inputs:
  # Y: dependent variable, n x 1 vector
  # X matrix: independent variables n x p matrix
  # mu1: prior mean for beta1
  # std1: prior standard deviation for beta1
  
  ############################
  # Frequentist estimates
  ############################
  freq_mu <- solve(t(X) %*% X) %*% (t(X) %*% Y)
  sigmasq <- c(sum((Y - X %*% freq_mu)^2) / (nrow(X) - ncol(X)))  # Estimated sigma
  freq_var <- sigmasq * solve(t(X) %*% X)   # Frequentist variance
  
  ############################
  # Bayesian estimates
  ############################
  # Prior information
  bayes_prior_mu <- matrix(rep(0, 5), nrow = 5)
  bayes_prior_var <- diag(10000, 5)
  
  bayes_prior_mu[2] <- mu1
  bayes_prior_var[2, 2] <- std1^2
  
  # Posterior mean and variance (Exact inference w/ informative prior; plug-in sigma^2 = sigmasq)
  
  Sigma_inv <- solve(bayes_prior_var)
  
  bayes_pos_mu <- solve( t(X) %*% X + sigmasq * Sigma_inv ) %*%
    ( t(X) %*% Y + sigmasq * Sigma_inv %*% bayes_prior_mu )
  
  bayes_pos_var <- sigmasq * solve( t(X) %*% X + sigmasq * Sigma_inv )
  
  ############################
  # Hybrid Bayesian-Frequentist estimates
  ############################
  hb_mu <- mu1
  hb_var <- std1^2
  
  # Get initial frequentist estimates
  freq_mu1 <- solve(t(X) %*% X) %*% (t(X) %*% Y)
  
  # Bayesian update for beta1 (hybrid conditional mean; exact informative prior)
  Y_temp <- Y - X[, -2, drop = FALSE] %*% freq_mu1[-2]
  X_temp <- matrix(X[, 2], ncol = 1)
  
  hb_mu1_pos <- solve(t(X_temp) %*% X_temp + sigmasq * solve(hb_var)) %*%
    (t(X_temp) %*% Y_temp + sigmasq * solve(hb_var) * hb_mu)
  
  
  # Frequentist update
  Y_temp <- Y - X[, 2] * as.vector(hb_mu1_pos)
  X_temp <- X[, c(1,3,4,5)]
  hb_mu1_fre <- solve(t(X_temp) %*% X_temp) %*% (t(X_temp) %*% Y_temp)
  
  # Iterate for 100 times
  for (i in 1:100) {
    # Bayesian update for beta1 (conditional on current frequentist params)
    Y_temp <- Y - X[, -2, drop = FALSE] %*% hb_mu1_fre
    X_temp <- matrix(X[, 2], ncol = 1)
    
    hb_mu1_pos_new <- solve(t(X_temp) %*% X_temp + sigmasq * solve(hb_var)) %*%
      (t(X_temp) %*% Y_temp + sigmasq * solve(hb_var) * hb_mu)
    
    
    # Frequentist update for others
    Y_temp <- Y - X[, 2] * as.vector(hb_mu1_pos_new)
    X_temp <- X[, -2]
    hb_mu1_fre_new <- solve(t(X_temp) %*% X_temp) %*% (t(X_temp) %*% Y_temp)
    
    hb_mu1_pos <- hb_mu1_pos_new
    hb_mu1_fre <- hb_mu1_fre_new
  }
  
  ############################
  # Store Hybrid estimates
  ############################
  hb_est <- c(hb_mu1_fre[1], hb_mu1_pos, hb_mu1_fre[2], hb_mu1_fre[3], hb_mu1_fre[4])
  est3 <- cbind(freq_mu, bayes_pos_mu, hb_est, hb_est)
  
  ############################
  # Variance estimates
  ############################
  var_freq <- diag(freq_var)
  var_bayes <- diag(bayes_pos_var)
  
  #-----------------------------------------------------------------
  # --- Hybrid conditional variance --- 
  #-----------------------------------------------------------------
  
  hb_var <- std1^2         # prior variance
  
  # Rearrange so beta_1 is first
  X_re <- X[, c(2,1,3,4,5), drop = FALSE]
  
  # Prior precision matrix
  hb_var <- std1^2
  Sigma_inv <- matrix(0, ncol(X_re), ncol(X_re))
  Sigma_inv[1,1] <- 1 / hb_var
  
  # Equation (5) in manuscript
  XtX <- t(X_re) %*% X_re
  Q <- XtX + sigmasq * Sigma_inv
  Q_inv <- solve(Q)
  V <- sigmasq * Q_inv %*% XtX %*% Q_inv
  
  # Partition V
  V11 <- V[1,1, drop = FALSE]
  V12 <- V[1,-1, drop = FALSE]
  V21 <- V[-1,1, drop = FALSE]
  V22 <- V[-1,-1, drop = FALSE]
  
  # Finite-sample conditional variance
  V_cond_beta1 <- as.numeric(
    V11 - V12 %*% solve(V22) %*% V21
  )
  
  #----------------------------------------------------------
  # --- Hybrid unconditional variance ---
  #----------------------------------------------------------
  bayes_idx <- 2  # bayesian parameter in full design matrix
  
  # x1 (Bayesian covariate) and X_F (frequentist design matrix)
  X_1  <- X[, bayes_idx, drop = FALSE]      # n x 1
  X_F <- X[, -bayes_idx, drop = FALSE]     # n x (d-1)
  
  # I_F: remove the Bayesian row/parameter
  I_F <- diag(ncol(X))[-bayes_idx, , drop = FALSE]   # (d-1) x d
  
  # I: n x n identity
  I_n <- diag(nrow(X))
  
  # A = (x1'x1 + sigma^2 / sigma1^2)^(-1)
  A <- 1 / (as.numeric(t(X_1) %*% X_1) + (sigmasq / hb_var))
  
  # B = I - XF IF (X'X)^(-1) X'
  B <- I_n - X_F %*% I_F %*% solve(t(X) %*% X) %*% t(X)
  
  # Var = sigma^2 A x1' B B' x1 A
  V_uncond_beta1 <- as.numeric(
    sigmasq * A * (t(X_1) %*% B %*% t(B) %*% X_1) * A
  )
  
  # Assemble variances
  #hybrid conditional
  var_hb <- diag(freq_var)
  var_hb[2] <- V_cond_beta1
  
  # Unconditional Hybrid
  var_uncond_hb <- diag(freq_var)
  var_uncond_hb[2] <- V_uncond_beta1
  
  
  VAR <- cbind(var_freq, var_bayes, var_hb, var_uncond_hb)
  se3 <- sqrt(VAR)
  
  # Compute p-values
  p3 <- matrix(0, nrow = 5, ncol = 4)
  for (i in 1:5) {
    for (j in 1:4) {
      z_score <- est3[i, j] / se3[i, j]
      if (pt(z_score, df = length(Y) - 5) < 0.5) {
        p3[i, j] <- 2 * pt(z_score, df = length(Y) - 5)
      } else {
        p3[i, j] <- 2 * (1 - pt(z_score, df = length(Y) - 5))
      }
    }
  }
  return(list(estimates = est3, se = se3, p = p3))
}
