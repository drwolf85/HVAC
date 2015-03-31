# Smart automation of an HVAC system

The current folder contains the code to optimize the process behind the selection of action of an HVAC system. An empirical approach based on the **decision networks** is followed. This implies that the identification of the structure of the network, the study of the evolution of its edges, the regression on the quantities to optimize, their **forecast and the multi-objective optimization** must be performed.

To facilitate the understanding of the process, the code is organized in the following sub-folders:
 1. Network Identification;
 2. Robust Forecasts;
 3. Latent Variable;
 4. Optimization.

The main non-interpreted languages are C and C++, while the [R software](http://www.r-project.org/) is used to interpret the scripts with `.R` extension.

A `README.md` file is present in each folder. It describes the file included, the R packages required for the analysis, and the steps to use build the shared libraries adopted inside the R-scripts.

## Requirements

 * R software and [GNU compilers](https://gcc.gnu.org/).

 * [R packages](http://cran.r-project.org/):
     * `igraph` for dealing with graphs and network analysis
     * `stats` for dealing with standard statistical functions
     * `quantreg` for the median-LASSO regression
     * `Rcpp` for interfacing R with Cpp scripts
     * `biglm` for iterative regression of big linear models

## Aknoledgments

This software was developed with F.S.E. founding program of Veneto region (Italy).

![Veneto region](https://salute.regione.veneto.it/portal-web-theme/images/custom/logo-regione.png)
