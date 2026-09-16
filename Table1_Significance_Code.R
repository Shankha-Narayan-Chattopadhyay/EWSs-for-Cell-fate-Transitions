library(EWSmethods)
library(tictoc)
tic()
# Read CSV
data <- read.table(
  "D:/Cell fate EWSs/Data/Forward_TS1.csv",
  header = FALSE
)

# Convert to numeric matrix
data <- as.matrix(data)

m <- 1
n <- 910  #Take upto the identified transition time

# Time index
t <- seq(m, n)

# Extract Notch, Delta and NICD
N <- as.numeric(data[m:n, 3])
D <- as.numeric(data[m:n, 4])
I <- as.numeric(data[m:n, 5])

combined_data <- cbind(t, N, D, I)


# Testing of significance on raw Data 
# Here "arima" for ARIMA surrogate 
# Put "sample" for Random Sampling surrogate
# Put "red.noise" for Stochastic Red Noise surrogate

perm_ews_eg_roll <- perm_rollEWS(
  data = combined_data,
  metrics = c(
    "meanAR", "maxAR",
    "meanSD", "maxSD",
    "eigenMAF", "mafAR",
    "mafSD", "pcaAR", "pcaSD",
    "eigenCOV", "maxCOV", "mutINFO"
  ),
  variate = "multi",
  perm.meth = "arima",
  winsize = 50,
  iter = 50
)

print(perm_ews_eg_roll$EWS$cor) 
toc()