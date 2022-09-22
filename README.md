# Rapster
Rapid population synthesis code for binary black hole mergers in dynamical environments.

Author: Konstantinos Kritos <kkritos1@jhu.edu>

Date: September 22, 2022

### Contents:
1. Overview
2. Dependancies
3. Input parameters
4. Running a simulation
5. Output files
6. Applications of the code
7. Citing this work
8. Reporting bugs
9. Thanks

---

### 1. Overview

The repository provides the source codes, files ./main.py and ./functions.py, and all necessary data files in folder ./MzamsMrem/, for the rapid evolution of dense star cluster environments and the dynamical assembly of binary black hole mergers.

The modeling accounts for the necessary physical processes regarding the formation of binary black holes employing semi-analytic prescriptions as described in Sec. 2 of [K. Kritos, V. Strokov, V. Baibhav, E. Berti, to appear].

##### Note:
For computational efficiency, the folder ./MzamsMrem/ contains 12 files with pre-calculated tables of stellar remnants masses on a grid of zero-age main sequence values up to 340$M_\odot$ and 12 values of absolute metallicity in the range from $10^{-4}$ to $1.7\times10^{-2}$ as calculated with the $\tt SEVN$ code [M. Spera & M. Mapelli, (2017)](https://academic.oup.com/mnras/article/470/4/4739/3883764).

### 2. Depedancies:

The following Python packages are required

- $\tt precession$ (1.0.3)
- $\tt astropy$ (5.0.4)
- $\tt argparse$ (1.1)
- $\tt cmath$
- $\tt numpy$ (1.12.3)
- $\tt scipy$ (1.8.0)
- $\tt math$
- $\tt time$

The code is tested with packages in the versions shown in parentheses above, however it is likely that other versions work too.

##### Note:
It is suggested that the $\tt precession$ package is used in its latest version 1.0.3 [D. Gerosa & M. Kesden, (2016)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.93.124066).

### 3. Input parameters:

Our code accepts parameters with flag options.

For a description of all input parameters, run the following command in the command line interface:

> python main.py --help

or see Table 1 from [K. Kritos, V. Strkov, V. Baibhav, E. Berti, to appear].
