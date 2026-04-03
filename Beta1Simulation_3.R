##############################################
# SIMULATION: β1 ESTIMATION & SE SUMMARY
# 36 Scenarios: n × σ(Y|X) × σ(β1) × True Effect
# Output: Mean(SD) of β1 and Median [Q1,Q3] of SE
##############################################

library(officer)
library(flextable)


n_sim <- 10000 # total simulations

true_effects <- c(0, 0.5,1)
prior_means  <- c(0, 0.2, 0.5, 1)
sample_sizes <- c(50, 100, 300)
data_sds     <- c(0.58, 3)     # σ(Y|X)
prior_sds    <- c(0.1, 0.7, 1)      # σ(β1)

summary_results <- list()
set.seed(2025)
##############################################
# MAIN SIMULATION LOOPS
##############################################
for (n in sample_sizes) {
  for (sigma_y in data_sds) {
    for (std1 in prior_sds) {
      for (true_eff in true_effects) {
        
        cat("\n===== Running scenario:",
            "n =", n, "| sigma_y =", sigma_y,
            "| prior_sd =", std1, "| true_eff =", true_eff, "=====\n")
        
        # Containers for results across prior means
        result_df <- data.frame(
          Prior_Effect = prior_means,
          Frequentist = NA,
          Bayesian = NA,
          HB_Conditional = NA,
          HB_Unconditional = NA
        )
        
        result_se <- data.frame(
          Prior_Effect = prior_means,
          Frequentist_SE = NA,
          Bayesian_SE = NA,
          HB_Conditional_SE = NA,
          HB_Unconditional_SE = NA
        )
        
        ##############################################
        # LOOP OVER PRIOR MEANS
        ##############################################
        for (j in seq_along(prior_means)) {
          mu1 <- prior_means[j]
          cat("  -> Prior mean =", mu1, "\n")
          
          # Storage vectors
          b_f <- b_b <- b_hc <- b_hu <- numeric(n_sim)
          se_f <- se_b <- se_hc <- se_hu <- numeric(n_sim)
          
          ##############################################
          # SIMULATION REPLICATIONS
          ##############################################
          for (i in 1:n_sim) {
            # Covariates
            x1 <- rep(c(0, 1), each = n / 2)
            x2 <- rnorm(n, mean = 7, sd = 3)
            x3 <- rbinom(n, size = 1, prob = 0.5)
            x4 <- rbinom(n, size = 1, prob = 0.5)
            X <- cbind(1, x1, x2, x3, x4)
            
            # Outcome Y with current sigma_y
            Y <- rnorm(n, mean = 7 + true_eff * x1, sd = sigma_y)
            Y <- matrix(Y, ncol = 1)
            
            # Run estimation
            res <- bayes_hybrid_regression(Y, X, mu1, std1)
            
            # β1 estimates and SEs (row 2)
            b_f[i]  <- res$estimates[2, 1]
            b_b[i]  <- res$estimates[2, 2]
            b_hc[i] <- res$estimates[2, 3]
            b_hu[i] <- res$estimates[2, 4]
            
            se_f[i]  <- res$se[2, 1]
            se_b[i]  <- res$se[2, 2]
            se_hc[i] <- res$se[2, 3]
            se_hu[i] <- res$se[2, 4]
          }
          
          # β1 Mean (SD)
          result_df$Frequentist[j]      <- sprintf("%.3f (%.3f)", mean(b_f), sd(b_f))
          result_df$Bayesian[j]         <- sprintf("%.3f (%.3f)", mean(b_b), sd(b_b))
          result_df$HB_Conditional[j]   <- sprintf("%.3f (%.3f)", mean(b_hc), sd(b_hc))
          result_df$HB_Unconditional[j] <- sprintf("%.3f (%.3f)", mean(b_hu), sd(b_hu))
          
          # SE Median [Q1,Q3]
          q_f  <- quantile(se_f,  probs = c(0.25, 0.75))
          q_b  <- quantile(se_b,  probs = c(0.25, 0.75))
          q_hc <- quantile(se_hc, probs = c(0.25, 0.75))
          q_hu <- quantile(se_hu, probs = c(0.25, 0.75))
          
          result_se$Frequentist_SE[j]      <- sprintf("%.3f [%.3f, %.3f]", median(se_f), q_f[1], q_f[2])
          result_se$Bayesian_SE[j]         <- sprintf("%.3f [%.3f, %.3f]", median(se_b), q_b[1], q_b[2])
          result_se$HB_Conditional_SE[j]   <- sprintf("%.3f [%.3f, %.3f]", median(se_hc), q_hc[1], q_hc[2])
          result_se$HB_Unconditional_SE[j] <- sprintf("%.3f [%.3f, %.3f]", median(se_hu), q_hu[1], q_hu[2])
        }
        
        ##############################################
        # STORE RESULTS
        ##############################################
        key <- paste0("n", n, "_sigY", sigma_y,
                      "_priorSD", std1, "_trueEff", true_eff)
        summary_results[[key]] <- list(beta = result_df, se = result_se)
      }
    }
  }
}

##############################################
# DISPLAY RESULTS IN CONSOLE AFTER SIMULATION
##############################################

cat("\n==============================\n")
cat("SUMMARY OF ALL SCENARIOS\n")
cat("==============================\n\n")


for (name in names(summary_results)) {
  scenario <- summary_results[[name]]
  
  cat("\n------------------------------------------\n")
  cat("Scenario:", name, "\n")
  cat("------------------------------------------\n")
  
  cat("\nMean (SD) of β₁ estimates:\n")
  print(scenario$beta)
  
  cat("\nMedian [Q1, Q3] of Standard Errors:\n")
  print(scenario$se)
  
  cat("\n------------------------------------------\n")
}

##############################################
# WORD DOCUMENT OUTPUT
##############################################
doc <- read_docx()

for (name in names(summary_results)) {
  scenario <- summary_results[[name]]
  
  # Parse scenario name for heading
  parts <- unlist(strsplit(name, "_"))
  heading <- paste(
    "n =", sub("n", "", parts[1]),
    "| σ(Y|X) =", sub("sigY", "", parts[2]),
    "| σ(β₁) =", sub("priorSD", "", parts[3]),
    "| True Effect =", sub("trueEff", "", parts[4])
  )
  
  doc <- body_add_par(doc, heading, style = "heading 2")
  doc <- body_add_par(doc, "Mean (SD) of β₁ Estimates", style = "heading 3")
  doc <- body_add_flextable(doc, flextable(scenario$beta))
  doc <- body_add_par(doc, "Median [Q1, Q3] of Standard Errors", style = "heading 3")
  doc <- body_add_flextable(doc, flextable(scenario$se))
  doc <- body_add_par(doc, "")
}

print(doc, target = "Final_Simulation_Results5.docx")
cat("\n Simulation and report generation complete: 'Simulation_Results.docx'\n")
getwd()
