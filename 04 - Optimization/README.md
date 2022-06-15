# Stochastic multi-objective optimization

The files contained in this folder are necessary to select the optimal action to perform before the next datum is observed. The algorithm is based on a multi-objective optimization with stochastic objectives.

**ATTENTION: to execute the script `OrganizeData.R` the `EMalgo.c` file must be compiled to be loaded as a dynamical library. See the section "Compiling the `EMalgo.c`".**

The main file to execute in R is the `pareto.R` script, which identifies the optimal action to perform in HVAC system.

The final results is a vector of dominance probabilities, which provides information on the most dominant solution to select.

## Required R packages

 * `stats` for dealing with standard statistical functions
 * `quantreg` for the median-LASSO regression
 * `Rcpp` for interfacing R with Cpp scripts
 * `biglm` for iterative regression of big linear models

## Script-Tools required by the main files

 * `EMalgo.c` is a C script file including the functions for identify the hidden values of the response variable via the EM algorithm.
 * `pareto.cpp` is a Cpp script file including the functions for optimize the HVAC system by reducing the energy consumption and increasing the comfort conditions.
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
