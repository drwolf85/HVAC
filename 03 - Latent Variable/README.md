# Latent Response Regression

The files contained in this folder are necessary to estimate the models, make predictions and inference on the energy consumption.

**ATTENTION: to execute the main R-script the `EMalgo.c` file must be compiled to be loaded as a dynamical library. See the section "Compiling the `EMalgo.c`".**

The main file to execute in R is the `LatentRegression.R` script, which identifies a consistent latent process and estimates the model for *Energy Consumption*.

The final results are stored in `.RData` files.
 
## Required R packages

 * `stats` for dealing with standard statistical functions
 * `quantreg` for the median-LASSO regression
 * `Rcpp` for interfacing R with Cpp scripts
 * `biglm` for iterative regression of big linear models

## Script-Tools required by the main files

 * `EMalgo.c` is a C script file including the functions for identify the hidden values of the response variable via the EM algorithm.
 * `npfit.cpp` is a Cpp script file including the functions for a calculating the design matrix of the multi-interactions Kernel-based additive model.
 * `readData.R` is a script file including the commands to read the data-sets.
 * `missing.R` is a script file including the  instructions to adjust the data-sets and dealing with missing data.

## Compiling the `EMalgo.c`
On **Windows** the `Rtools.exe` must be installed and the environmental variable `PATH` must contain the folders containing the GNU compilers provided.

By opening the command prompt (`cmd.exe`), it is possible to point to this directory with `cd` and then execute the compiler via `Rcmd.exe` with the following command:

```bash
Rcmd.exe SHLIB EMalgo.c --preclean
```

On **Linux** it is necessary to open a terminal and entering the following command:
```bash
R CMD SHLIB EMalgo.c --preclean
```
If no compiler is installed, one of the following options is required, depending on the Linux distribution in use.

On Debian / Ubuntu / Mint:
```bash
sudo apt-get install gcc*
```
On OpenSuSE:
```bash
sudo zypper install gcc*
```
On Mandriva / Mageia
```bash
urpmi gcc*
```
On Fedora / RedHat / CentOS
```bash
yum install gcc*
```
On ArchLinux
```bash
sudo pacman -S gcc*
```

On **MacOSX** it is necessary to install XCode, and from its interface installing the GNU compilers.
By opening a terminal and entering the following command:
```bash
R CMD SHLIB EMalgo.c --preclean
```
it is possible to obtain the compiled shared object.
