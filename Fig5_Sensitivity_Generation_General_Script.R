# ============================================================
# Sensitivity Analysis of Early Warning Signals
# ============================================================
#
# This script evaluates the sensitivity of an Early Warning
# Signal (EWS) to the choice of Gaussian detrending bandwidth
# and rolling-window size.
#
# The example below calculates the sensitivity of meanAR
# (mean lag-1 autocorrelation) using Kendall's tau.
#
# The resulting sensitivity matrix can subsequently be read
# into MATLAB for visualization using the contourf command.
#
# ============================================================


library(EWSmethods)


# ------------------------------------------------------------
# 1. Load Data
# ------------------------------------------------------------

data <- read.table(
  "D:/Cell fate EWSs/Data/Forward_TS1.csv",
  header = FALSE
)

# Convert the imported data into a numeric matrix
data <- as.matrix(data)


# ------------------------------------------------------------
# 2. Select Time Series
# ------------------------------------------------------------

m <- 1
n <- 910   # Data considered up to the identified transition time

# Time index
t <- seq(m, n)

# Extract Notch, Delta, and NICD time series
N <- as.numeric(data[m:n, 3])
D <- as.numeric(data[m:n, 4])
I <- as.numeric(data[m:n, 5])

# Construct the multivariate time series.
# The first column represents the time index.
combined_data <- cbind(t, N, D, I)


# ------------------------------------------------------------
# 3. Define Sensitivity Parameters
# ------------------------------------------------------------

# Gaussian detrending bandwidths
bandwidth_list <- seq(5, 75, by = 10)

# Rolling-window sizes
winsize_list <- seq(25, 75, by = 10)


# ------------------------------------------------------------
# 4. Initialize Sensitivity Matrix
# ------------------------------------------------------------

# Rows    : Gaussian bandwidths
# Columns : Rolling-window sizes
#
# Each element of the matrix contains Kendall's tau
# measuring the temporal trend of the EWS.

result_matrix <- matrix(
  NA,
  nrow = length(bandwidth_list),
  ncol = length(winsize_list)
)


# ------------------------------------------------------------
# 5. Sensitivity Analysis for meanAR
# ------------------------------------------------------------

for (i in seq_along(bandwidth_list)) {
  
  bandwidthsize <- bandwidth_list[i]
  
  # Gaussian detrending
  combined_ts <- detrend_ts(
    combined_data,
    method = "gaussian",
    bandwidth = bandwidthsize
  )
  
  
  for (j in seq_along(winsize_list)) {
    
    winsize <- winsize_list[j]
    
    
    # --------------------------------------------------------
    # Calculate meanAR using a rolling window
    # --------------------------------------------------------
    
    ews_result <- multiEWS(
      combined_ts,
      metrics = "meanAR",
      method = "rolling",
      winsize = winsize
    )
    
    
    # Extract the raw meanAR time series
    datax <- ews_result$EWS$raw[["meanAR"]]
    
    
    # --------------------------------------------------------
    # Calculate Kendall's tau
    # --------------------------------------------------------
    
    if (!is.null(datax)) {
      
      # Remove non-finite values
      valid <- is.finite(datax)
      
      if (sum(valid) > 2) {
        
        timevec <- seq_along(datax)
        
        KtAR <- cor.test(
          timevec[valid],
          datax[valid],
          method = "kendall",
          alternative = "two.sided"
        )
        
        # Store Kendall's tau
        result_matrix[i, j] <- as.numeric(KtAR$estimate)
      }
    }
    
    
    # Display progress
    cat(
      "Bandwidth =", bandwidthsize,
      "| Window size =", winsize,
      "| Kendall's tau =", result_matrix[i, j],
      "\n"
    )
  }
}


# ------------------------------------------------------------
# 6. Display Sensitivity Matrix
# ------------------------------------------------------------

print(result_matrix)


# ------------------------------------------------------------
# 7. Save Sensitivity Matrix
# ------------------------------------------------------------

# Only numerical values are stored.
# Row names and column names are intentionally omitted.

write.table(
  result_matrix,
  file = "D:/Cell fate EWSs/Sensitivity_Data_meanAR.csv",
  sep = ",",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)


cat("\nSensitivity matrix saved successfully:\n")
cat("Sensitivity_Data_meanAR.csv\n")


# ============================================================
# MATLAB Visualization
# ============================================================
#
# The generated CSV file can be imported into MATLAB and
# visualized using the contourf command.
#
# For example:
#
#   A = readmatrix('Sensitivity_Data_meanAR.csv');
#   contourf(A);
#
# The rows correspond to the selected Gaussian bandwidths,
# while the columns correspond to the selected rolling-window
# sizes.
#
# ============================================================
#
# NOTE:
# For computational efficiency during exploratory sensitivity
# analysis, the bandwidth and rolling-window sizes can be
# varied using a larger increment.
#
# For the final sensitivity analysis, the increment can be
# reduced to 2:
#
#   bandwidth_list <- seq(5, 75, by = 2)
#   winsize_list   <- seq(25, 75, by = 2)
#
# To calculate the sensitivity of another EWS, replace
# "meanAR" in the multiEWS() function and in the extraction
# step with the desired EWS metric.
#
# ============================================================