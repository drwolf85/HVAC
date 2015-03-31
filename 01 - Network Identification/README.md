# Network Identification Analysis

The files contained in this folder are necessary to identify the structure of the network. The main file to be executed in R is the `netAnalysis.R` script. To evaluate if the structure of the network change with time the file `netStable.R` must be executed successively.
(The file `finalNet.R` evaluates separately the connections for the *PMV* and the *DGI index*. This is done in order to display better the resulting networks as sub-graphs.)

The identified structure of the networks is displayed as final output of the analysis.

## Required R packages

 * `igraph` for dealing with graphs and network analysis

## Script-Tools required by the main file

 * `readData.R` is a script file including the commands to read the data-sets.

 * `missing.R` is a script file including the instructions to adjust the data-sets and dealing with missing data.

 * `Entropy.R` is a script file containing all the functions to calculate the entropy, mutual information and all related statistics.

