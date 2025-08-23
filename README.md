# Luria-Delbr-ck-Fluctuation-Assays-Maximum-Likelihood-Estimation-R-Tool
Computational Analysis of Luria-Delbrück Fluctuation Assays Using Maximum Likelihood Estimation: A Comprehensive R-Based Analysis Tool
Readme

Fluctuation Assay Analysis in R
Enhanced implementation combining MSS Maximum Likelihood and Two-Parameter Model Analysis
Based on the work by Abdalla, Walker & Ishimori (2024) and Original MATLAB functions: Greg Lang, Harvard University (Lang & Murray, 2008). 
Overview
This repository provides R implementations for analyzing Luria-Delbrück fluctuation assay data using both the Ma-Sandri-Sarkar (MSS) Maximum Likelihood method and two-parameter model analysis. The code combines and extends functionality from published algorithms to provide comprehensive mutation rate analysis.
Features
•	Complete MSS Maximum Likelihood Analysis with bootstrap confidence intervals
•	Two-Parameter Model Analysis for post-plating growth effects
•	Comprehensive Visualization with cumulative distribution plots
•	Mutation Rate Calculations with confidence intervals
•	MATLAB-to-R Conversion of established algorithms
Quick Start
r
# Load your mutation count data
data <- c(0, 0, 0, 1, 1, 2, 3, 4, 7, 15, 27, 94, 145)  # example data

# Analyze (nt = total cells per culture is MANDATORY)
results <- analyze_fluctuation_data(data, nt = 320000)

# Results include:
# - Two-parameter model estimates (m, d)
# - MSS maximum likelihood estimates
# - Bootstrap confidence intervals  
# - Mutation rate calculations
# - Comprehensive visualization
Installation
Simply source the R script:
r
source("fluctuation_assay_analysis.R")
Requirements
•	R (≥ 3.5.0)
•	Base R packages only (no additional dependencies)
Input Data Format
Your data should be a vector where each element represents the number of mutants observed in one culture:
r
data <- c(0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 7, 15, 27, 94, 145)
Parameters
Required Parameters
•	data: Vector of mutation counts per culture
•	nt: Total number of cells tested per culture (required for mutation rate calculation)
Optional Parameters
•	nBootstraps: Number of bootstrap samples (default: 1000)
•	seed: Random seed for reproducibility (default: 4)
Output
The analysis provides:
1. Two-Parameter Model Results
•	m: Expected number of mutations per culture
•	d: Post-plating growth parameter
•	Based on algorithms by Lang & Murray
2. MSS Maximum Likelihood Results
•	m: Maximum likelihood estimate
•	Bootstrap 95% CI: Non-parametric confidence intervals
•	Normal approximation 95% CI: Parametric confidence intervals
•	Mutation rate (μ): Per-cell mutation rate with confidence intervals
3. Visualization
•	Cumulative distribution plot comparing: 
o	Empirical data (red stepped line)
o	One-parameter model (green line)
o	Two-parameter model (purple circles)
Example Analysis
r
# Example with the data from the original paper
data <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 4, 4, 6, 7, 7, 7, 7, 8, 8, 9, 10, 11, 11, 12, 13, 13,
          13, 15, 16, 21, 27, 31, 35, 81, 94, 102, 145)

# Run complete analysis
results <- analyze_fluctuation_data(data, nt = 320000, nBootstraps = 1000, seed = 4)

# Access specific results
cat("Two-parameter model m:", results$matlab_equivalent_analysis$m_estimate, "\n")
cat("Two-parameter model d:", results$matlab_equivalent_analysis$d_estimate, "\n")
cat("MSS ML estimate:", results$mss_analysis$m_estimate, "\n")
cat("Mutation rate:", sprintf("%.6e", results$mss_analysis$mutation_rate), "\n")
Expected output:
Two-parameter model m: 1.0237
Two-parameter model d: 0
MSS ML estimate: 1.106733  
Mutation rate: 3.458540e-06
Algorithm Details
MSS Maximum Likelihood Method
Implements the non-recursive Ma-Sandri-Sarkar algorithm for maximum likelihood estimation of mutation rates in fluctuation assays, as described in:
•	Sarkar S, Ma WT, Sandri GH. On fluctuation analysis: a new, simple and efficient method for computing the expected number of mutants. Genetica. 1992;85(2):173-9.
Two-Parameter Model
Implements the combined Luria-Delbrück/Poisson model to account for post-plating growth effects, as described in:
•	Lang GI, Murray AW. Estimating the per-base-pair mutation rate in the yeast Saccharomyces cerevisiae. Genetics. 2008;178(1):67-82.
Bootstrap Confidence Intervals
Uses non-parametric bootstrapping (resampling with replacement) to generate robust confidence intervals for mutation rate estimates.
File Structure
├── README.md                           # This file
├── fluctuation_assay_analysis.R        # Main R script
├── example_analysis.R                  # Usage examples
├── matlab_original/                    # Original MATLAB implementations
│   ├── findMLmTwoParam.m              # Two-parameter ML estimation
│   ├── generateLD.m                   # Luria-Delbrück distribution
│   ├── generatePO.m                   # Poisson distribution  
│   ├── generateTwoParam.m             # Combined distribution
│   ├── SSDScoreLD.m                   # One-parameter scoring
│   ├── SSDScoreTwoParam.m             # Two-parameter scoring
│   ├── finalCurve.m                   # Main MATLAB script
│   └── other supporting functions...

    
Validation
The R implementations have been validated against the original MATLAB code and produce identical results:
•	Two-parameter estimates match MATLAB output exactly
•	Visualizations reproduce the same cumulative distribution plots
•	All intermediate calculations verified for consistency
Contributing
Contributions are welcome! Please:
1.	Fork the repository
2.	Create a feature branch
3.	Submit a pull request with clear description of changes
License
This code is provided under the MIT License. Please cite appropriately if used in publications.
Citation
If you use this code in your research, please cite:
Primary Citation:
Abdalla, O. M., Walker, C. G., & Ishimori, K. (2024). R-code for calculating fluctuation assay results and 95% confidence intervals based on Ma-Sandri-Sarkar maximum likelihood. Software Impacts, 21, 100661. https://doi.org/10.1016/j.simpa.2024.100661
And cite the original methodological papers:
•	Sarkar S, Ma WT, Sandri GH. On fluctuation analysis: a new, simple and efficient method for computing the expected number of mutants. Genetica. 1992;85(2):173-9.
•	Lang GI, Murray AW. Estimating the per-base-pair mutation rate in the yeast Saccharomyces cerevisiae. Genetics. 2008;178(1):67-82.
Related Resources
•	Code Ocean Repository: Estimation of mutation rate for fluctuation assay via MSS Maximum Likelihood
•	Original Publication: Abdalla, O. M., Walker, C. G., & Ishimori, K. (2024). Software Impacts, 21, 100661.
Contact 
ola_abdallah@science.tanta.edu.eg


Acknowledgments
• Primary R Implementation: Abdalla, O. M., Walker, C. G., & Ishimori, K. (2024)
• Original MATLAB functions: Greg Lang, Harvard University (Lang & Murray, 2008). 
The following functions are based on his publication:
• generateLD() - Generates Luria-Delbruck distribution
• generatePO() - Generates Poisson distribution
• generateTwoParam() - Generates combined distribution
• findMLm() - Finds ML estimate for single parameter
• findMLmTwoParam() - Finds ML estimates for two parameters
•	MSS algorithm: Sarkar, Ma, and Sandri
•	Statistical methods: Based on Foster (2006) and Rosche & Foster (2000)
