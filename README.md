# Rapster
Rapid population synthesis code for binary black hole mergers in dynamical environments.

Author: Konstantinos Kritos <kkritos1@jhu.edu>

Date: September 22, 2022

### Contents:
1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Input parameters](#inputparameters)
4. [Running a simulation](#runningasimulation)
5. [Output files](#outputfiles)
6. [Applications of the code](#applicationsofthecode)
7. [Citing this work](#citingthiswork)
8. [Reporting bugs](#reportingbugs)
9. [Thanks](#thanks)

---

##### Abbreviations:

- BH: black hole
- BBH: binary black hole
- GW: gravitational wave

<a name="overview"></a>
### 1. Overview

The repository provides the source codes, files ./rapster.py and ./functions.py, and all necessary data files in folder ./MzamsMrem/, for the rapid evolution of dense star cluster environments and the dynamical assembly of binary black hole mergers.

The modeling accounts for the necessary physical processes regarding the formation of binary black holes employing semi-analytic prescriptions as described in Sec. 2 of [K. Kritos et al. (2022)]().

##### Note:
For computational efficiency, the folder ./MzamsMrem/ contains 12 files with pre-calculated tables of stellar remnants masses on a grid of zero-age main sequence values up to $340M_\odot$ and 12 values of absolute metallicity in the range from $10^{-4}$ to $1.7\times10^{-2}$ as calculated with the $\tt SEVN$ code [M. Spera & M. Mapelli (2017)](https://academic.oup.com/mnras/article/470/4/4739/3883764).

<a name="requirements"></a>
### 2. Requirements

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
It is suggested that the $\tt precession$ package is used in its latest version 1.0.3 [D. Gerosa & M. Kesden (2016)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.93.124066).

<a name="inputparameters"></a>
### 3. Input parameters

Our code accepts parameters with flag options.

For a description of all input parameters, run the following command in the command line interface:

> python3 rapster.py --help

or see Table 1 from [K. Kritos et al. (2022)]().

For the user’s convenience we paste the list of optional arguments in the form of a Table here as well:

| Flag | Description | Type | Default |
|:--- |:--- |:--- |:--- |
| -Mcl, --ClusterMass | Initial cluster mass $(M_\odot)$ | float | $10^6M_\odot$ |
| -rh, --HalfMassRadius | Initial half-mass radius (pc) | float | 1pc |
| -rhoC, --CentralDensity | Initial central star density $(M_\odot{\rm pc}^{-3})$ | float | $4\times10^5M_\odot{\rm pc}^{-3}$ |
| -Rgal, --GalactocentricRadius | Initial galactocentric radius (kpc) | float | 8kpc |
| -Z, --Metallicity | Cluster metallicity $(Z_\odot)$ | float | $0.1Z_\odot$ |
| -fb, --BinaryFraction | Initial binary star fraction | float | 10\% |
| -w, --NatalKickParameter | Natal velocity kick parameter of BHs (km/s) | float | 256km/s |
| -chi, --NatalSpinParameter | Natal spin parameter of first generation (1g) BHs | float | 0 |
| -SM, --NatalSpinDistribution | Natal spin distribution (1 for monochromatic, 0 for uniform) | int | 0 |
| -tMax, --SimulationTime | Maximum simulation time (Myr) | float | 13,800Myr |
| -dtMin, --MinimumTimeStep | Minimum simulation timestep (Myr) | float | 0.1Myr |
| -dtMax, --MaximumTimeStep | Maximum simulation timestep (Myr) | float | 50Myr |
| -z, --FormationRedshift | Redshift of cluster formation | float | 3 |
| -aIMF, --HighMassIMFslope | High mass initial star mass function slope | floar | -2.3 |
| -ZAMSmax, --MaximumZAMSmass | Maximum ZAMS star mass $(M_\odot)$ | float | $150M_\odot$ |
| -c, --ConcentrationNFW | Concentration parameter of the NFW profile | float | 10 |
| -Rs, --ScaleRadiusNFW | Scale radius of the NFW profile (kpc) | float | 50kpc |
| -Mh, --DarkMatterHaloMass | Total DM halo mass of the host galaxy $(M_\odot)$ | float | $10^{12}M_\odot$ |
| -s, --Seed | Random number generator seed | int | 123456789 |
| -MF, --MergersFile | Name of output file with BBH merger source parameters | str | ``mergers'' |
| -EF, --EvolutionFile | Name of output file with time-dependent quantities | str | ``evolution'' |
| -BF, --BlackHoleFile | Name of output file containing the masses of all 1g BHs in $M_\odot$ | str | ``blackholes'' |

<a name="runningasimulation"></a>
### 4. Running a simulation

usage: rapster.py [-h] [-Mcl ] [-rh ] [-rhoC ] [-Rgal ] [-Z ] [-fB ] [-w ] [-chi ] [-SM ] [-tMax ] [-dtMin ] [-dtMax ] [-z ] [-aIMF ] [-ZAMSmax ] [-c ] [-Rs ] [-Mh ] [-s ] [-MF ] [-EF ] [-BF ]

As an example we give the commands that produce data used to generate the results in Fig.4 of [K. Kritos, V. Strokov, V. Baibhav, E. Berti, to appear]:

  > python3 rapster.py -Mcl  1.36e5 -rh 1.6 -rhoC  5.6e4 -EF ev_a -MF me_a -BF bh_a -Z 0.08 -z 3 -Rgal  8

  > python3 rapster.py -Mcl  5.40e5 -rh 1.6 -rhoC 26.9e4 -EF ev_b -MF me_b -BF bh_b -Z 0.08 -z 3 -Rgal  8

  > python3 rapster.py -Mcl 10.82e5 -rh 1.6 -rhoC 15.1e4 -EF ev_c -MF me_c -BF bh_c -Z 0.03 -z 3 -Rgal 20

The default values are assumed for other parameters not entered in the commands above.

##### Suggestion:
Taking different values of seed number corresponds to different realizations of the system under the same initial conditions. 
Passing the argument $$\tt\$ RANDOM$$ in the -s flag, simulates the star cluster with a pseudo-randomly generated number.

<a name="outputfiles"></a>
### 5. Output files:

At the end of each simulation the code generates three .txt files, one with the black hole masses of all first generation black holes that are initially retained in the cluster, a second file with information about all dynamical mergers that took place during the simulation, and finally a file that keeps track of time-dependent quantities.

a) Column description of mergers .txt file: 

| Column | Variable | Description |
|:--- |:--- |:--- |
| 1 | channel | Formation mechanism of BBH |
| 2 | $a$ | Semimajor axis (AU) when BBH enters the GW regime |
| 3 | $e$ | Eccentricity of BBH when BBH enters the GW regime |
| 4 | $m_1$ | Primary mass $(M_\odot)$ |
| 5 | $m_2$ | Secondary mass $(M_\odot)$ |
| 6 | $\chi_1$ | Primary dimensionless spin parameter |
| 7 | $\chi_2$ | Secondary dimensionless spin parameter |
| 8 | $g_1$ | Generation of primary |
| 9 | $g_2$ | Generation of secondary |
| 10 | $t_{\rm form}$ | Time (Myr) BBH formed since simulation started |
| 11 | $t_{\rm merge}$ | Time (Myr) BBH merged since simulation started |
| 12 | $z_{\rm form}$ | Redshift BBH formed |
| 13 | $z_{\rm merge}$ | Redshift BBH merged |
| 14 | $N_{\rm har}$ | Number of hardening interactions |
| 15 | $N_{\rm sub}$ | Number of BH exchanges |
| 16 | $q$ | Mass ratio |
| 17 | $\chi_{\rm eff}$ | Effective spin parameter |
| 18 | $\theta_1$ | Polar angle (rad) of primary's spin with orbital angular momentum |
| 19 | $\theta_2$ | Polar angle (rad) of secondary's spin with orbital angular momentum |
| 20 | $\Delta\phi$ | Azimuthal angle (rad) between BH spins in the orbital plane |
| 21 | $m_{\rm rem}$ | Remnant mass $(M_\odot)$ |
| 22 | $\chi_{\rm rem}$ | Dimensionless remnant spin parameter |
| 23 | $g_{\rm rem}$ | Remnant generation |
| 24 | $v_{\rm GW}$ | GW kick (km/s) of remnant BH |
| 25 | $j_1$ | Reference to primary’s progenitor BBH (if hierarchical product) |
| 26 | $j_2$ | Reference to secondary’s progenitor BBH (if hierarchical product) |
| 27 | $M_{\rm cl,0}$ | Initial cluster mass $(M_\odot)$ |
| 28 | $z_{\rm cl,\ form}$ | Cluster formation redshift |

##### Note:
BBH assembly channel (first column of mergers file), the ``-'' sign means BBH was ejected from the cluster:
- (-)1: exchange processes
- (-)2: two-body capture
- (-)3: three-BH binary induced
- (-)4: von Zeipel-Lidov-Kozai merger

b) Column description of evolution .txt file:

| Columnn | Variable | Description |
|:--- |:--- |:--- |
| 1 | $t$                             |  Simulation time (Myr) |
| 2 | $z$                            |  Redshift |
| 3 | $M_{\rm cl}$                   |   Cluster mass $(M_\odot)$  |
| 4 | $r_{\rm h}$                    |         Cluster half-mass radius (pc) |
| 5 | $R_{\rm gal}$                  |        Cluster’s galactocentric radius (kpc) |
| 6 | $N_{\rm BH}$                   |    Total number of BHs in cluster |
| 7 | $N_{\rm BH,\ sin}$            |    Number of single BHs in cluster |             
| 8 | $N_{\rm BH-star}$             |    Number of BH-star pairs in cluster |
| 9 | $N_{\rm BBH}$                 |    Number of BBHs in cluster |
| 10 | $N_{\rm triples}$            |      Number of triple BHs in cluster |
| 11 | $N_{\rm me}$                 |      Number of BBH mergers in total |
| 12 | $N_{\rm me,\ re}$            |        Number of retained mergers |
| 13 | $N_{\rm me,\ ej}$            |         Number of ejected mergers |
| 14 | $N_{\rm me,\ out}$           |     Number of dynamical BBHs that merge in field |
| 15 | $N_{\rm ZLK}$                |      Number of ZLK mergers |
| 16 | $N_{\rm cap}$                |       Number of captured BBHs |
| 17 | $N_{\rm ej}$                 |          Number of BBH ejections |
| 18 | $v_{\rm esc}$                |           Escape velocity (km/s) |
| 19 | $n_{\rm star-star}$          |        Number density of binary stars $({\rm pc}^{-3})$ |
| 20 | $\overline{V_{\rm seg}}$     |            Mean segregation volume of BHs $({\rm pc}^3)$ |
| 21 | $\xi$                       |        Equipartition parameter |
| 22 | $v_{\rm star}$              |            Stellar velocity (km/s) |
| 23 | $v_{\rm BH}$                 |          Mean BH velocity (km/s) |
| 24 | $\tau_{\rm cap}$             | Mean capture timescale (Myr) |
| 25 | $\tau_{\rm ex,1}$            | Mean star-star $\to$ BH-star timescale (Myr) |
| 26 | $\tau_{\rm ex,2}$            | Mean BH-star $\to$ BBH timescale (Myr) |
| 27 | $\tau_{\rm 3bb}$             | Mean 3bb formation timescale (Myr) |
| 28 | $\tau_{bb}$ | Mean binary-binary interaction timescale (Myr) |
| 29 | $\tau_{pp}$ | Mean BH-star $-$ BH-star interaction timescale (Myr) |
| 30 | ${\cal N}_{\rm ex,1}$ | Cumulative number of star-star$\to$BH-star interactions |
| 31 | ${\cal N}_{\rm ex,2}$ | Cumulative number of BH-star$\to$BBH interactions |
| 32 | ${\cal N}_{\rm 3bb}$ | Cumulative number of 3bb interactions |
| 33 | ${\cal N}_{\rm cap}$ | Cumulative number of capture interactions |
| 34 | ${\cal N}_{bb}$ | Cumulative number of BBH$-$BBH interactions |
| 35 | ${\cal N}_{pp}$ | Cumulative number of BH-star$-$BH-star interactions |
| 36 | $\overline{m_{\rm BH}}$ | Mean BH mass $(M_\odot)$ in cluster |

<a name="applicationsofthecode"></a>
### 6. Applications of the code

The code can be useful when executed multiple times, for instance when simulating a set of clusters and generating a population of dynamically formed binary black hole mergers.

Although the program itself is not computationally expensive (we have tested in a laptop that we generate a few binary black hole mergers per second), independent parallelization is still encouraged when simulating a very large number of star clusters for efficiency.

<a name="citingthiswork"></a>
### 7. Citing this work

If you utilize this code in your research, please cite the following reference:

[K. Kritos et al. (2022)]().

<a name="reportingbugs"></a>
### 8. Reporting bugs

If you find a bug in the code, please contact us in kkritos1@jhu.edu with a description of the bug.

Suggestions and pull requests are welcome :)

<a name="thanks"></a>
### 9. Thanks

V. Strokov, V. Baibhav, E. Berti, A. Antonelli, M. Cheung, R. Cotesta, G. Franciolini, T. Helfer, V. Kapil, I. Krommydas, L. Reali, C. Rodriguez.


