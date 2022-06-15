# Robust Forecast and Error Analysis

The files contained in this folder are necessary to estimate the models, make predictions and evaluate which strategy performs best for the HVAC system analyzed.

There are three main files to execute in R:
 
 * the `arimax.R` script estimate the ARIMAX models for the *PMV* and *DGI indices* and evaluate the forecast errors.
 * the `robReg.R` script estimate the multi-interactions Kernel-based additive models for the *PMV* and *DGI indices* and evaluate the forecast errors.
 * the `bench.R` script evaluate the most standard indices of forecast performance and returns useful statistics for prediction-strategy selection.

The final results are stored in `.RData` files.

**ATTENTION: The files in this folder contain sequential code only. Their execution in computers with multicore technology may still require over 6 month of computations.**
 
## Required R packages

 * `stats` for dealing with standard statistical functions
 * `quantreg` for the median-LASSO regression
 * `Rcpp` for interfacing R with Cpp scripts

## Script-Tools required by the main files

 * `arimax.cpp` is a Cpp script file including the functions for a recurrent fitting of the traditional ARIMAX models.
 * `npfit.cpp` is a Cpp script file including the functions for a recurrent median-LASSO fitting of the multi-interactions Kernel-based additive model.
 * `benchforex.cpp` is a Cpp script file containing all the functions to calculate the benchmark indices of forecast performance (such as the RMSE, MAE, MAPE and many more).
 * `merge.R` is an R script file which is used to merge the results containing the forecast errors.
