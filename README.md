# Bayesian–Hybrid Inference Simulation

This repository provides code to implement and evaluate a hybrid Bayesian–frequentist regression framework, with a focus on estimating **β₁** (bayesian parameter) under varying data and prior conditions.

---

## Repository Structure

- **`updated_function7.R`**  
  Core estimation function implementing:
  - Frequentist regression
  - Bayesian regression with informative prior
  - Hybrid Bayesian–frequentist estimation
  - Conditional and unconditional variance estimators

- **`Beta1Simulation_3.R`**  
  Simulation script to evaluate estimator performance across multiple scenarios.

---

## Requirements

- R (version ≥ 4.0)

Install required packages:
```r
install.packages(c("officer", "flextable"))
```

---

## Running the Code

### Option 1: Set Working Directory (Recommended)

Set your working directory to the folder containing the repository files:

```r
setwd("path_to_your_folder")
```

Then run:

```r
source("updated_function7.R")
source("Beta1Simulation_3.R")
```

---

### Option 2: Without Setting Working Directory

If you prefer not to set a working directory, ensure both `.R` files are in the same folder and run:

```r
source("full_path/updated_function7.R")
source("full_path/Beta1Simulation_3.R")
```

---

## Output

Running the simulation will generate:

- Console summaries for all simulation scenarios
- A Word document:

```
Final_Simulation_Results5.docx
```

This file contains:
- Mean (SD) of **β₁** estimates  
- Median [Q1, Q3] of standard errors  

---

## Computational Considerations

The simulation is computationally intensive.

- Each scenario involves 10,000 Monte Carlo replications  
- A total of 36 scenarios are evaluated across combinations of sample size, outcome variance, prior variance, and true effect size  
- The hybrid estimation procedure includes iterative updates (100 iterations per replication)

As a result, the full simulation may take several hours to complete depending on system specifications.

For exploratory runs, users may reduce the number of simulations (`n_sim`) in `Beta1Simulation_3.R`; final results should be generated using the full simulation size.  

---

## Reproducibility

A fixed randomization seed is used:

```r
set.seed(2025)
```

Running the code with the same seed will reproduce identical results.

---

## Notes

- The simulation evaluates performance across combinations of:
  - Sample size
  - Outcome variance
  - Prior variance
  - True effect size
- All computations are fully self-contained within the provided scripts.
- No additional data input is required.

---

## Citation

If this code is used, please cite the associated manuscript.

---

## Contact

For questions regarding the implementation, please contact the authors.
