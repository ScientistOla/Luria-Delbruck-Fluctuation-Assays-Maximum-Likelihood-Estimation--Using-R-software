# Fluctuation Assay Analysis in R
# Based on Abdalla, O. M., Walker, C. G. (2024)and original MATLAB script by Greg Lang, Harvard University (Lang & Murray, 2008)
# Enhanced with Two-Parameter Model Analysis
#
# SIMPLE USAGE FOR RESEARCHERS:
# 1. Replace 'data' with your mutation counts
# 2. Replace 'nt' with your total cells per culture  
# 3. Run the script
# 4. Get results: 
#    - Parameter m (Two-Parameter Model)
#    - Parameter d (Two-Parameter Model) 
#    - Parameter m (MSS Maximum Likelihood Method)
#    - 95% confidence intervals for parameter m (Normal approximation)
#    - 95% confidence intervals for parameter m (Bootstrap)
#    - Mutation rate with confidence intervals
#    - Comparative plot

# =============================================================================
# USER INPUT SECTION - EDIT ONLY THIS PART
# =============================================================================

# INPUT YOUR DATA HERE:
data <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 4, 4, 6, 7, 7, 7, 7, 8, 8, 9, 10, 11, 11, 12, 13, 13,
          13, 15, 16, 21, 27, 31, 35, 81, 94, 102, 145)

# INPUT YOUR TOTAL CELLS PER CULTURE:
nt <- 320000

# Optional settings (you can leave these as default):
nBootstraps <- 1000    # Number of bootstrap samples for confidence intervals
seed <- 4              # Random seed for reproducible results

# =============================================================================
# ANALYSIS FUNCTIONS - DO NOT EDIT BELOW THIS LINE
# =============================================================================

# Generate Luria-Delbruck distribution
generateLD <- function(m, max_val) {
  ldDist <- exp(-m)
  for (i in 1:max_val) {
    x <- 0:(length(ldDist) - 1)
    ldDist <- c(ldDist, (m / length(x)) * sum(ldDist / (length(x) - x + 1)))
  }
  return(ldDist)
}

# Generate Poisson distribution  
generatePO <- function(lambda, max_val) {
  dist <- numeric(max_val + 1)
  for (x in 0:max_val) {
    dist[x + 1] <- (lambda^x / factorial(x)) * exp(-lambda)
  }
  return(dist)
}

# Generate Two-Parameter distribution
generateTwoParam <- function(m, d, max_val) {
  LD <- generateLD(m, max_val)
  PO <- generatePO(m * d, max_val)
  PO[is.nan(PO)] <- 0
  
  coDist <- numeric(max_val + 1)
  for (i in 0:max_val) {
    if (i == 0) {
      coDist[i + 1] <- LD[1] * PO[1]
    } else {
      coDist[i + 1] <- sum(LD[1:(i + 1)] * rev(PO[1:(i + 1)]))
    }
  }
  return(coDist)
}

# Score functions
scoreData <- function(data, m) {
  tabdata <- numeric(max(data) + 1)
  for (i in 0:max(data)) {
    tabdata[i + 1] <- length(which(data == i))
  }
  dist <- generateLD(m, max(data))
  return(sum(-log(dist^tabdata)))
}

scoreDataTwoParam <- function(data, m, d) {
  tabdata <- numeric(max(data) + 1)
  for (i in 0:max(data)) {
    tabdata[i + 1] <- length(which(data == i))
  }
  dist <- generateTwoParam(m, d, max(data))
  return(sum(-log(dist^tabdata)))
}

# Find ML estimate for single parameter
findMLm <- function(data) {
  m <- median(data) + 0.1
  if (m > 10) m <- 10
  
  for (precision in c(1, 0.1, 0.01, 0.001, 0.0001)) {
    repeat {
      ms <- c(m - precision, m, m + precision)
      ms[ms < 0] <- 0
      scores <- c(scoreData(data, ms[1]), scoreData(data, ms[2]), scoreData(data, ms[3]))
      new_m <- ms[which.min(scores)]
      if (new_m == m) break
      m <- new_m
    }
  }
  return(m)
}

# Find ML estimates for two parameters
findMLmTwoParam <- function(data) {
  m <- median(data) + 0.1
  if (m > 10) m <- 10
  d <- 0
  
  previous_m <- 0
  while (m != previous_m) {
    previous_m <- m
    
    # Optimize m
    for (precision in c(1, 0.1, 0.01, 0.001, 0.0001)) {
      repeat {
        ms <- c(m - precision, m, m + precision)
        scores <- c(scoreDataTwoParam(data, ms[1], d), 
                    scoreDataTwoParam(data, ms[2], d), 
                    scoreDataTwoParam(data, ms[3], d))
        scores[ms < 0] <- max(scores) + 1
        new_m <- ms[which.min(scores)]
        if (new_m == m) break
        m <- new_m
      }
    }
    
    # Optimize d
    for (precision in c(0.1, 0.01)) {
      repeat {
        ds <- c(d - precision, d, d + precision)
        ds[ds < 0] <- precision
        scores <- c(scoreDataTwoParam(data, m, ds[1]), 
                    scoreDataTwoParam(data, m, ds[2]), 
                    scoreDataTwoParam(data, m, ds[3]))
        scores[ds < 0] <- max(scores) + 1
        new_d <- ds[which.min(scores)]
        if (new_d == d) break
        d <- new_d
      }
    }
  }
  return(c(m, d))
}

# MSS Maximum Likelihood functions
conditionalProb_r <- function(m, r) {
  partials <- numeric(r + 1)
  partials[1] <- exp(-1 * m)
  if (r > 0) {
    if (r == 1) {
      return((m / r) * partials[1] / (r + 1))
    }
    for (i in 1:r) {
      partialSum <- 0
      for (j in 1:i) {
        partialSum <- partialSum + partials[j] / (i - (j - 1) + 1)
      }
      partials[i + 1] <- (m / i) * partialSum
    }
    return(partials[r + 1])
  }
  return(partials[1])
}

loglik <- function(m, rVec, rCard) {
  objfn <- 0
  for (i in 1:length(rVec)) {
    objfn <- objfn + rCard[i] * log(conditionalProb_r(m, rVec[i]))
  }
  return(objfn)
}

# =============================================================================
# MAIN ANALYSIS - AUTOMATIC EXECUTION
# =============================================================================

cat("=== FLUCTUATION ASSAY ANALYSIS ===\n")
cat("Based on Abdalla, Walker & Ishimori (2024)\n\n")

# Validate inputs
if (length(data) == 0) stop("ERROR: Please provide mutation count data")
if (nt <= 0) stop("ERROR: Total cells per culture (nt) must be positive")

cat("Data Summary:\n")
cat("- Number of cultures:", length(data), "\n")
cat("- Mutation counts range:", min(data), "to", max(data), "\n") 
cat("- Total cells per culture:", nt, "\n\n")

# =============================================================================
# 1. TWO-PARAMETER MODEL ANALYSIS
# =============================================================================

cat("--- TWO-PARAMETER MODEL RESULTS ---\n")
m_d <- findMLmTwoParam(data)
m_two_param <- m_d[1]
d_two_param <- m_d[2]

cat("m (expected mutations per culture) =", round(m_two_param, 6), "\n")
cat("d (post-plating growth parameter) =", round(d_two_param, 6), "\n\n")

# =============================================================================
# 2. MSS MAXIMUM LIKELIHOOD ANALYSIS  
# =============================================================================

cat("--- MSS MAXIMUM LIKELIHOOD RESULTS ---\n")

# Convert data to frequency format
rDt <- as.numeric(names(table(data)))
cDt <- as.numeric(table(data))

# Find MLE
f <- function(z) {-1 * loglik(z, rDt, cDt)}
m_mss <- optimize(f, c(0, 15))$minimum

cat("Parameter m (MSS Maximum Likelihood Method) =", round(m_mss, 6), "\n")

# Normal approximation confidence intervals
sigma <- 1.225 * (m_mss^-0.315) / sqrt(length(data))
u <- log(m_mss) + 1.96 * sigma * (exp(1)^(1.96 * sigma))^-0.315
l <- log(m_mss) - 1.96 * sigma * (exp(1)^(1.96 * sigma))^-0.315
normal_ci <- exp(c(l, u))

cat("95% confidence intervals for parameter m (Normal approximation): [", round(normal_ci[1], 4), ",", round(normal_ci[2], 4), "]\n")

# Bootstrap confidence intervals
cat("Calculating bootstrap confidence intervals...\n")
set.seed(seed)
bReplicates <- numeric(nBootstraps)
C <- length(data)

for (b in 1:nBootstraps) {
  indices <- sample(1:C, C, replace = TRUE)
  bRep <- data[indices]
  rBRep <- as.numeric(names(table(bRep)))
  cBRep <- as.numeric(table(bRep))
  f <- function(z) {-1 * loglik(z, rBRep, cBRep)}
  bReplicates[b] <- optimize(f, c(0, 15))$minimum
}

boot_ci <- quantile(bReplicates, c(0.025, 0.975))
cat("95% confidence intervals for parameter m (Bootstrap): [", round(boot_ci[1], 4), ",", round(boot_ci[2], 4), "]\n\n")

# Mutation rate
mu <- m_mss / nt
mu_ci_normal <- normal_ci / nt
mu_ci_bootstrap <- boot_ci / nt

cat("--- MUTATION RATE RESULTS ---\n")
cat("Mutation rate (μ) =", sprintf("%.6e", mu), "\n")
cat("Mutation rate 95% CI (Normal approximation): [", sprintf("%.6e", mu_ci_normal[1]), ",", sprintf("%.6e", mu_ci_normal[2]), "]\n")
cat("Mutation rate 95% CI (Bootstrap): [", sprintf("%.6e", mu_ci_bootstrap[1]), ",", sprintf("%.6e", mu_ci_bootstrap[2]), "]\n\n")

# =============================================================================
# 3. GENERATE PLOT
# =============================================================================

cat("--- GENERATING PLOT ---\n")

# Calculate cumulative distributions
maxValue <- max(data)
SSDScoreLDResult <- cumsum(generateLD(findMLm(data), maxValue))
SSDScoreTwoParamResult <- cumsum(generateTwoParam(m_two_param, d_two_param, maxValue))

# Prepare plotting data
dataLength <- length(data)
SSDScoreLDResultForDraw <- numeric(dataLength)
SSDScoreTwoParamResultForDraw <- numeric(dataLength)

for (c in 1:dataLength) {
  SSDScoreLDResultForDraw[c] <- SSDScoreLDResult[data[c] + 1]
  SSDScoreTwoParamResultForDraw[c] <- SSDScoreTwoParamResult[data[c] + 1]
}

# Create plot
plot(data, SSDScoreTwoParamResultForDraw, 
     pch = 1, col = rgb(0.4940, 0.1840, 0.5560), 
     lwd = 2.5, cex = 1.2,
     xlab = "Mutants per culture", 
     ylab = "Cumulative distribution",
     main = "Fluctuation Assay Analysis",
     xlim = c(0, maxValue), ylim = c(0, 1))

# Add empirical data (step function)
data_sorted <- sort(data)
empirical_y <- (1:length(data)) / length(data)
lines(data_sorted, empirical_y, col = "red", lwd = 2, type = "s")
points(data_sorted, empirical_y, pch = 16, cex = 0.5, col = "red")

# Add one-parameter model
lines(data, SSDScoreLDResultForDraw, 
      col = rgb(0.4660, 0.6740, 0.1880), lwd = 2)

# Add legend
legend("bottomright", 
       legend = c("Two-parameter Model", "Data", "One-parameter Model"),
       col = c(rgb(0.4940, 0.1840, 0.5560), "red", rgb(0.4660, 0.6740, 0.1880)),
       pch = c(1, 16, NA),
       lty = c(NA, 1, 1),
       lwd = c(2.5, 2, 2),
       cex = 0.9)

# =============================================================================
# 4. SUMMARY OF ALL RESULTS
# =============================================================================

cat("=== FINAL RESULTS SUMMARY ===\n")
cat("Two-Parameter Model:\n")
cat("  Parameter m =", round(m_two_param, 6), "\n")
cat("  Parameter d =", round(d_two_param, 6), "\n\n")

cat("MSS Maximum Likelihood Method:\n")
cat("  Parameter m =", round(m_mss, 6), "\n")
cat("  95% CL for parameter m (Normal approximation): [", round(normal_ci[1], 4), ",", round(normal_ci[2], 4), "]\n")
cat("  95% CL for parameter m (Bootstrap): [", round(boot_ci[1], 4), ",", round(boot_ci[2], 4), "]\n\n")

cat("Mutation Rate Results:\n")
cat("  Mutation rate (μ) =", sprintf("%.6e", mu), "\n")
cat("  Mutation rate 95% CI (Normal): [", sprintf("%.6e", mu_ci_normal[1]), ",", sprintf("%.6e", mu_ci_normal[2]), "]\n")
cat("  Mutation rate 95% CI (Bootstrap): [", sprintf("%.6e", mu_ci_bootstrap[1]), ",", sprintf("%.6e", mu_ci_bootstrap[2]), "]\n\n")

cat("Plot generated showing comparison of models.\n")
cat("Analysis complete!\n")

# =============================================================================
# INSTRUCTIONS FOR USERS
# =============================================================================

cat("\n=== HOW TO USE THIS SCRIPT ===\n")
cat("1. Replace 'data' vector with your mutation counts\n") 
cat("2. Set 'nt' to your total cells per culture\n")
cat("3. Run the entire script\n")
cat("4. Results will be displayed and plot will be generated\n\n")

cat("Example data format:\n")
cat("data <- c(0, 0, 1, 1, 2, 3, 5, 8, 12)\n")
cat("nt <- 100000\n\n")

cat("Citation:\n")
cat("Abdalla, O. M., Walker, C. G. (2024).\n")
cat("Software Impacts, 21, 100661.\n")