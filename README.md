(1) Download the datasets & the codes marked by the unique Fig number to reproduce the very same of Chattopadhyay & Gupta 'Do multivariate measures presage cell-state transitions in the Notch pathway?' paper.

(2) First run Stocha_Simu_Gillespie.m code to generate the stochastic trajectory.

(3) Then, use the Detrending_Code.R code to detrend that stochastic trajectory using Gaussian kernel smoothing, LOESS filtering, first differencing or linear detrending. This will compute residuals, free from nonstationarities.

(4) Finally, to generate the multivariate early warning signal from the residuals as in Figs. 3, 7,8 and 10 use EWSmethods toolbox (Duncan O'Brien; https://duncanobrien.github.io/EWSmethods/). 

(5) For additional information regarding the necessary hyperparameter configurations follow the manuscript.

Note: The MATLAB version used for the analysis is R2024b, and the R version is 4.5.1.
