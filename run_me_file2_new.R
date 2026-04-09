# This file sets an example to call the function "bayes_est"
# Inputs include the outcome variable, prior mean, prior std



# Part 1, the function file
hybrid_est<-function(outcome_l,pi_mu,pi_sig){
  ## This function returns estimate and SE of phenotype variable 
  ## in the regression model using the Bayesian inference
  # Inputs: response variable, prior mean, prior STD of the 
  #         regression parameter for phenotype; 
  #                  
  #        Prior of other variables are non-informative
  #        Output variable may be standardized, and transformed in log scale
  #        Example: "outcome_l<-log(camkv$sum_CAMKV+1) #sum_CAMV RNA count, log scale" 
  #        Length of the outcome has to be the same as each independent variable      
  #
  #        The prior mean and std should be decided by pilot data set and AI 
  #        inputs of the prior strength.
  #
  # Output: estimate and standard error of the regression parameter for 
  #         phenotype (case control status) from the Bayesian inference
  
  library(readxl)
  # Load independent variables from excel spreadsheet;
  # Here make sure to specify the location of the data file with covariates
  camkv <- read_excel(".\\camkv2.xlsx")
  
  pheno <- camkv$group_n  #1 AD, 0 control 
  ncell_l <- log(camkv$n_cell_patient+1)  # number of cells, log scale
  apoe <- camkv$apoe  # Apoe genotype, 0: 23, 33; 1: 24, 34, 44
  age <- camkv$age_bl #Age
  sex <- camkv$sex_c
  
  # design matrix from ind_data
  Y <- outcome_l-ncell_l
  X <- matrix(c(rep(1, length(Y)), pheno, age,apoe, sex),ncol=5)
  
  # frequentist estimates
  freq_mu1 <- solve(t(X)%*%X)%*%(t(X)%*%Y)
  sigmasq <- c(sum((Y - X %*% freq_mu1)^2) / (nrow(X) - ncol(X)))  # estimated sigma
  freq_sig1 <- sigmasq * solve(t(X) %*% X)
  
  # SET Prior information
  prior_mu1 <- matrix(c(0,sign(freq_mu1[2])*abs(pi_mu),0,0,0),nrow=5)
  # we use "sign(freq_mu1[2])" so that the effect of the prior is in the 
  # same direction as data
  prior_sig1 <- matrix(0,5,5)
  diag(prior_sig1) <- 1000000000
  prior_sig1[2,2] <- pi_sig^2
  
  # Full Bayesian estimates
  # Posterior mean and variance
  #mu1_pos <- solve(solve(prior_sig1)+t(X)%*%X)%*%(solve(prior_sig1)%*%prior_mu1+t(X)%*%Y)
  mu1_pos <- solve(solve(prior_sig1)*sigmasq +t(X)%*%X) %*% (sigmasq*solve(prior_sig1)%*%prior_mu1+t(X)%*%Y)
  #sig1_pos <- solve(solve(prior_sig1)*sigmasq+t(X)%*%X)*sigmasq
  
  #----------------------------------------------------------
  # --- Hybrid unconditional variance ---
  #----------------------------------------------------------
  bayes_idx <- 2  # bayesian parameter in full design matrix
  hb_var <- pi_sig^2
  
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
  # Unconditional Hybrid
  var_b <- V_uncond_beta1
  
  est_b <- mu1_pos[2]
  se_b <- sqrt(var_b)
  
  z_score <- est_b/se_b
  if (pnorm(z_score)<0.5) 
  {p_value <- 2*pnorm(z_score)} else   
  {p_value <- 2*(1-pnorm(z_score))}
  
  # Initiate outpus
  result=list(est=as.null(),se=as.null(), p_value = as.null())
  result$est <- est_b
  result$se  <- se_b 
  result$p_value <- p_value
  return(result) 
}






##### Part 2, Set up the inputs

# Outcome: log scale of the original sum
camkv <- read_excel('.\\camkv2.xlsx')
outcome_l<-log(camkv$sum_CAMKV+1)


# Prior depends on the ChatGPT outcome of Cohen's D and 
# the pooled variance estimate from pilot data


# Inputs from the AI levels, 0 to 5
ai_level=2

if (ai_level == 0) {
  ai_effect = 0
} else if (ai_level == 1) {
  ai_effect = 0.2  
} else if (ai_level == 2) {
  ai_effect = 0.5  
} else if (ai_level == 3) {
  ai_effect = 0.9  
} else if (ai_level == 4) {
  ai_effect = 1.5
} else if (ai_level == 5) {
  ai_effect = 3.25  
}

# Pooled variance from pilot data
# Find the pooled_variance from prior data set;
n1= 12
n2= 9
s1= 0.393
s2= 0.551
sp= sqrt(((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2))

# inputs about the prior
pi_sig <- sp
pi_mu <- ai_effect*sp

##### Part 3, Run the function (from Part 1) with inputs from Part 2
example_hybridresult = hybrid_est(outcome_l,pi_mu,pi_sig)




