# Rapster
Rapid population synthesis code for binary black-hole mergers in dynamical environments.

$\tt Rapster$ stands for $\rm RAPid\ cluSTEP$ evolution. (Thanks to M. Cheung for coming up with the short-hand version!)

Author: Konstantinos Kritos <kkritos1@jhu.edu>

Version: July 7, 2023

![LOGO](./Rapster_LOGO.png)
(Thanks to H. Cruz for digitizing my hand-drawn logo!)

### Contents:
1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Units](#units)
4. [Input parameters](#inputparameters)
5. [Running a simulation](#runningasimulation)
6. [Output files](#outputfiles)
7. [Applications of the code](#applicationsofthecode)
8. [Citing this work](#citingthiswork)
9. [Reporting bugs](#reportingbugs)
10. [Thank you](#thanks)

---

<a name="overview"></a>
### 1. Overview

The repository provides the source codes of the current version, files ``./rapster2.py``, ``./functions2.py``, ``./constants2.py``, ``ThreeBodyBinary2.py``, ``BBHevol2.py``, ``TwoBodyCapture2.py``, ``Exchanges2.py``, ``triples2.py``, ``Planck18_lookup_table.npz``, and all necessary data files in folder ``./MzamsMrem/``, for the rapid evolution of dense star cluster environments and the dynamical assembly of binary black hole mergers.

The modeling accounts for the necessary physical processes regarding the formation of binary black holes employing semi-analytic prescriptions as described in Sec. 2 of [K. Kritos et al. (2022)](https://arxiv.org/abs/2210.10055). This is our code paper we wrote together with V. Strokov, V. Baibhav, and E. Berti.

##### Note:
For computational efficiency, the folder ``./MzamsMrem/`` contains 12 files with pre-calculated look-up tables of stellar remnants masses on a grid of zero-age main sequence values up to $340M_\odot$ and 12 values of absolute metallicity in the range from $10^{-4}$ to $1.7\times10^{-2}$ as calculated with the $\tt SEVN$ code [M. Spera & M. Mapelli (2017)](https://academic.oup.com/mnras/article/470/4/4739/3883764).

In the current version, we have also included look-up tables for the Mandel-Muller, and the Fryer et al. (2012) delayed and rapid remnant-mass prescription models in files ``Mueller_Mandel.txt``, ``MzamsMrem_F12d.txt``, and ``MzamsMrem_F12r.txt``, respectively.
Finally, ``Planck18_lookup_table.npz`` is a look-up table for the redshidt-lookback time and lookback time-redshift relations assuming the Planck 2018 cosmology.

##### Abbreviations:

- BH: black hole
- BBH: binary black hole
- GW: gravitational wave

<a name="requirements"></a>
### 2. Requirements

The following Python packages are required

- $\tt precession$ (1.0.3)
- $\tt astropy$ (5.0.4)
- $\tt argparse$ (1.1)
- $\tt numpy$ (>=1.12.3)
- $\tt scipy$ (1.8.0)
- $\tt pandas$
- $\tt time$

The code is tested with packages in the versions shown in parentheses above, however, it is likely that other versions work too.

##### Note:
It is suggested that the $\tt precession$ package is used in the version 1.0.3 [D. Gerosa & M. Kesden (2016)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.93.124066).

<a name="units"></a>
### 3. Units

The current version of the code uses astrophysical units:
 - Mass: solar mass ($M_\odot$)
 - Distance: parsec ($\rm pc$)
 - Time: $\rm Myr$
 - Velocity: $\rm km\ s^{-1}$
 - Gravitational constant: $G=1/232$

<a name="inputparameters"></a>
### 4. Input parameters

The code accepts parameters with flag options.

For a description of all input parameters, run the following command in the command line interface:

> python3 rapster2.py --help

or see Table 1 from [K. Kritos et al. (2022)](https://arxiv.org/abs/2210.10055).

For the user’s convenience we paste the list of optional arguments in the form of a Table here as well:

| Flag | Description | Type | Default |
|:--- |:--- |:--- |:--- |
| -N, --number | Initial number of stars | int | ``1000000`` |
| -r, --half_mass_radius | Initial half-mass radius [pc] | float | ``1`` |
| -mm, --minimum_star_mass | Smallest ZAMS mass [Msun] | float | ``0.08`` |
| -mM, --maximum_star_mass | Largest ZAMS mass [Msun] | float | ``150`` |
| -Z, --metallicity | Absolute metallicity | float | ``0.001`` |
| -z, --cluster_formation_redshift | Redshift of cluster formation | float | ``3.0`` |
| -n, --central_stellar_density | Central stellar number density [pc^-3] | float | ``1e6`` |
| -fb, --binary_fraction | Initial binary star fraction | float | ``0.1`` |
| -S, --seed | Seed number | int | ``1234567890`` |
| -dtm, --minimum_time_step | Minimum simulation time-step [Myr] | float | ``0.1`` |
| -dtM, --maximum_time_step | Maximum simulation time-step [Myr] | float | ``50.0`` |
| -tM, --maximum_time | Maximum simulation time [Myr] | float | ``140000`` |
| -wK, --supernova_kick_parameter | One-dimensional supernova kick parameter [km/s] | float | ``265.0`` |
| -K, --natal_kick_prescription | Natal kick prescription (0 for fallback, 1 for momentum conservation) | int | ``0`` |
| -R, --galactocentric_radius | Initial galactocentric radius [kpc] | float | ``8.0`` |
| -vg, --galactocentric_velocity | Galactocentric circular velocity [km/s] | float | ``220.0`` |
| -s, --spin_parameter | Natal spin parameter of first generation (1g) BHs | float | ``0.0`` |
| -SD, --spin_distribution | Natal spin distribution model (0 for uniform, 1 for monochromatic) | int | ``0`` |
| -P, --print_information | Print runtime information (0 for no, 1 for yes) | int | ``1`` |
| -Mi, --mergers_file_indicator | Export mergers file (0 for no, 1 for yes) | int | ``1`` |
| -MF, --mergers_file_name | Name of .txt output file with BBH merger source parameters | str | ``mergers`` |
| -Ei, --evolution_file_indicator | Export evolution file (0 for no, 1 for yes) | int | ``1`` |
| -EF, --evolution_file_name | Name of .txt output file with time-dependent quantities | str | ``evolution`` |
| -Hi, --hardening_file_indicator | Export hardening file (0 for no, 1 for yes) | int | ``1`` |
| -HF, --hardening_file_name | Name of .txt output file with BBH time evolution information | str | ``hardening`` |
| -BIi, --blackholes_in_file_indicator | Use external BH file (0 for no, 1 for yes) | int | ``0`` |
| -BIF, --blackholes_in_file_name | Name of .npz input file with initial BH masses | str | ``input_BHs.npz`` |
| -BOi, --blackholes_out_file_indicator | Export BH masses file (0 for no, 1 for yes) | int | ``1`` |
| -BOF, --blackholes_out_file_name | Name of .npz file with the masses of all BHs in solar masses | str | ``output_BHs.npz`` |
| -RP, --remnant_mass_prescription | Remnant mass prescription (0 for SEVN delayed, 1 for Fryer+2012 delayed, 2 for Fryer+2012 rapid) | int | ``1`` |

<a name="runningasimulation"></a>
### 5. Running a simulation

usage: rapster2.py [-h] [-N] [-r] [-mm] [-mM] [-Z] [-z] [-n] [-fb] [-S] [-dtm] [-dtM] [-tM] [-wK] [-K] [-R] [-vg] [-s] [-SD] [-P] [-Mi] [-MF] [-Ei] [-EF] [-Hi] [-HF] [-BIi] [-BIF] [-BOi] [-BOF] [-RP]

As an example, we give the commands that produce data used to generate the results in Fig.4 of [K. Kritos et al. (2022)](https://arxiv.org/abs/2210.10055):

  > python3 rapster2.py -N 200000 -r 1.6 -n 9.5e4 -R 8 -Z 0.001 -MF meA -EF evA -HF haA -BOF bhA

  > python3 rapster2.py -N 800000 -r 1.6 -n 45.6e4 -R 8 -Z 0.001 -MF meB -EF evB -HF haB -BOF bhC

  > python3 rapster2.py -N 1600000 -r 1.6 -n 257e4 -R 20 -Z 0.0005 -MF meC -EF evC -HF haC -BOF bhC

The default values are assumed for other parameters not entered in the commands above.

##### Note:
To test the code, execute the program with all default values:

  > python3 rapster2.py
  
This should create four files ``mergers.txt``, ``evolution.txt``, ``hardening.txt``, and ``output_BHs.npz`` in your current directory. To check and verify whether you have produced these files correctly, we include the corresponding files ``mergers_TEST.txt``, ``evolution_TEST.txt``, ``hardening_TEST.txt``, and ``output_BHs_TEST.npz`` in folder ``./Testing2/`` in this repository with data that should match your output.

##### Suggestion:
Taking different values of seed number corresponds to different realizations of the system under the same initial conditions. 
Passing the argument $$\tt\$ RANDOM$$ in the -s flag simulates the star cluster with a pseudo-randomly generated number. Notice this syntax works only in the bash environment.

<a name="outputfiles"></a>
### 6. Output files:

At the end of each simulation the code generates by default three .txt and one .npz file, one with information about all dynamical mergers that took place during the simulation, a second file that keeps track of time-dependent quantities, a third file that stores information about the hardening evolution of each BBH, and finally a file with the masses of the BHs that are initially retained and at end of the simulation, respectively.

a) Column description of mergers .txt file: 

| Column | Variable | Description |
|:--- |:--- |:--- |
| 1 | $\rm seed$ | seed of the simulation |
| 2 | $\rm ind$ | Unique integer binary identifier |
| 3 | $\rm channel$ | Formation mechanism of BBH |
| 4 | $a$ | Semimajor axis ($\rm pc$) when BBH enters the GW regime |
| 5 | $e$ | Eccentricity when BBH enters the GW regime |
| 6 | $m_1$ | Primary mass ($M_\odot$) |
| 7 | $m_2$ | Secondary mass ($M_\odot$) |
| 8 | $\chi_1$ | Primary dimensionless spin parameter |
| 9 | $\chi_2$ | Secondary dimensionless spin parameter |
| 10 | $g_1$ | Generation of primary |
| 11 | $g_2$ | Generation of secondary |
| 12 | $\theta_1$ | Polar angle ($\rm rad$) of primary's spin with orbital angular momentum |
| 13 | $\theta_2$ | Polar angle ($\rm rad$) of secondary's spin with orbital angular momentum |
| 14 | $\Delta\phi$ | Azimuthal angle ($\rm rad$) between BH spins in the orbital plane |
| 15 | $t_{\rm form}$ | Time ($\rm Myr$) BBH formed since simulation started |
| 16 | $z_{\rm form}$ | Redshift BBH formed |
| 17 | $t_{\rm merge}$ | Time ($\rm Myr$) BBH merged since simulation started |
| 18 | $z_{\rm merge}$ | Redshift BBH merged |
| 19 | $m_{\rm rem}$ | Remnant mass ($M_\odot$) |
| 20 | $\chi_{\rm rem}$ | Remnant dimensionless spin parameter |
| 21 | $g_{\rm rem}$ | Remnant generation |
| 22 | $v_{\rm GW}$ | GW kick (km/s) of remnant BH |
| 23 | $\chi_{\rm eff}$ | Effective dimensionless spin parameter |
| 24 | $q$ | Mass ratio |

##### Note:
BBH assembly channel (first column of mergers file), the ``-`` sign means BBH was ejected and merged outside the cluster:
- (-)1: exchange processes
-    2: two-body capture
- (-)3: three-BH binary induced
- (-)4: von Zeipel-Lidov-Kozai (ZLK) merger
- (-)5: ZLK remnant BBH
-    6: binary-single capture

b) Column description of evolution .txt file:

| Columnn | Variable | Description |
|:--- |:--- |:--- |
| 1 | $\rm seed$ | seed of the simulation |
| 2 | $t$ | Simulation time ($\rm Myr$) |
| 3 | $z$ | Redshift |
| 4 | $dt$ | Timestep ($\rm Myr$) |
| 5 | $\overline{m}$ | Average mass ($M_\odot$) |
| 6 | $M_{\rm cl}$ | Cluster mass ($M_\odot$) |
| 7 | $r_{\rm h}$ | Half-mass radius ($\rm pc$) |
| 8 | $R_{\rm gal}$ | Galactocentric radius ($\rm kpc$) |
| 9 | $v_{\rm gal}$ | Galactocentric velocity ($\rm km\ s^{-1}$) |
| 10 | $\tau_{\rm rlx}$ | Half-mass relaxation timescale ($\rm Myr$) |
| 11 | $\tau_{\rm rlx,BH}$ | BH half-mass relaxation time ($\rm Myr$) |
| 12 | $n_{\rm star}$ | Stellar central number density ($\rm pc^{-3}$) |
| 13 | $N_{\rm BH}$ | Number of BHs |
| 14 | $\overline{m}_{\rm BH}$ | Average BH mass ($M_\odot$) |
| 15 | $m_{\rm BH}^{\rm max}$ | Heaviest BH mass ($M_\odot$) |
| 16 | $r_{\rm h,BH}$ | BH half-mass radius ($\rm pc$) |
| 17 | $r_{\rm c,BH}$ | BH core radius ($\rm pc$) |
| 18 | $S$ | Spitzer parameter |
| 19 | $\xi$ | Equipartition parameter |
| 20 | $\psi$ | Multimass relaxation factor |
| 21 | $\psi_{\rm BH}$ | BH multimass relaxation factor |
| 22 | $\tau_{\rm 3bb}$ | 3bb timescale ($\rm Myr$) |
| 23 | $\tau_{\rm 2,cap}$ | 2-capture timescale ($\rm Myr$) |
| 24 | $k_{\rm 3bb}$ | Number of 3bb in current step |
| 25 | $k_{\rm 2,cap}$ | Number of 2-captures in current step |
| 26 | $N_{\rm me}$ | Cumulative number of mergers |
| 27 | $N_{\rm BBH}$ | Current number of BBHs |
| 28 | $N_{\rm me,Re}$ | Cumulative number of ejected merger remnants |
| 29 | $N_{\rm me,Ej}$ | Cumulative number of retained merger remnants |
| 30 | $v_{\rm star}$ | Stellar velocity dispersion ($\rm km\ s^{-1}$) |
| 31 | $v_{\rm BH}$ | BH velocity dispersion ($\rm km\ s^{-1}$) |
| 32 | $n_{\rm h,BH}$ | Half-mass BH number density ($\rm pc^{-3}$) |
| 33 | $n_{\rm c,BH}$ | Core BH number density ($\rm pc^{-3}$) |
| 34 | $n_{\rm a,BH}$ | Average BH number density ($\rm pc^{-3}$) |
| 35 | $N_{\rm 3bb}$ | Cumulative number of 3bb |
| 36 | $N_{\rm 2,cap}$ | Cumulative number of 2-captures |
| 37 | $N_{\rm 3,cap}$ | Cumulative number of 3-captures |
| 38 | $N_{\rm BH,ej}$ | Cumulative number of single BH ejections |
| 39 | $N_{\rm BBH,ej}$ | Cumulative number of BBH ejections |
| 40 | $N_{\rm dis}$ | Cumulative number of BBH disruptions |
| 41 | $N_{\rm ex}$ | Cumulative number of BBH-BH exchanges |
| 42 | $\tau_{\rm bb}$ | BBH-BBH interaction timescale ($\rm Myr$) |
| 43 | $N_{\rm bb}$ | Cumulative number of BBH-BBH interactions |
| 44 | $N_{\rm me,Fi}$ | Cumulative number of field mergers |
| 45 | $N_{\rm me,2b}$ | Cumulative number of in-cluster 2-body mergers |
| 46 | $\tau_{\rm ex,1}$ | star-star$\to$BH-star timescale ($\rm Myr$) |
| 47 | $\tau_{\rm ex,2}$ | BH-star$\to$BBH timescale ($\rm Myr$) |
| 48 | $k_{\rm ex,1}$ | Number of star-star$\to$BH-star exchanges in this step |
| 49 | $k_{\rm ex,2}$ | Number of BH-star$\to$BBH exchanges in this step |
| 50 | $N_{\rm ex,1}$ | Cumulative number of star-star$\to$BH-star exchanges |
| 51 | $N_{\rm ex,2}$ | Cumulative number of BH-star$\to$BBH exchanges |
| 52 | $N_{\rm BH-star}$ | Current number of BH-star pairs |
| 53 | $\tau_{\rm pp}$ | BH-star--BH-star interaction timescale ($\rm Myr$) |
| 54 | $k_{\rm pp}$ | Number of BH-star--BH-star interactions in this step |
| 55 | $N_{\rm pp}$ | Cumulative number of BH-star--BH-star interactions |
| 56 | $v_{\rm esc}$ | Escape velocity ($\rm km\ s^{-1}$) |
| 57 | $v_{\rm esc,BH}$ | Escape velocity from BH subsystem ($\rm km\ s^{-1}$) |
| 58 | $N_{\rm triples}$ | Cumulative number of BH triples |
| 59 | $N_{\rm ZLK}$ | Cumulative number of ZLK mergers |

c) Column description of hardening .txt file:

| Columnn | Variable | Description |
|:--- |:--- |:--- |
| 1 | $t$ | Global time ($\rm Myr$) |
| 2 | $dt$ | Global timestep ($\rm Myr$) |
| 3 | $t_{\rm local}$ | Local time ($\rm Myr$) |
| 4 | $dt_{\rm local}$ | Local timestep ($\rm Myr$) |
| 5 | $\rm ind$ | Binary's unique integer identifier |
| 6 | $a$ | Semimajor axis ($\rm pc$) |
| 7 | $e$ | Eccentricity |
| 8 | $m_1$ | Primary mass ($M_\odot$) |
| 9 | $m_2$ | Secondary mass ($M_\odot$) |
| 10 | $q$ | Mass ratio |
| 11 | $\rm condition$ | BBH status |
| 12 | $N_{\rm ex}$ | Number of BBH-BH exchanges |

##### Note:
Condition or binary status (last column of hardening file):
- 0: BBH available to evolve (see the flowchart of our algorithm in Fig.3 of [K. Kritos et al. (2022)](https://arxiv.org/abs/2210.10055))
- 1: local time exceeds global time
- 2: 2-body merger (the BBH hardens and merges in the cluster after entering the GW regime)
- 3: binary ionized
- 4: binary-single capture
- 5: binary ejected

Unless ${\rm condition}=0$, the local simulation is terminated.

d) The output_BHs.npz file (if exported) contains two arrays, called ``mBH_ini`` and ``mBH_fin`` which provide in $M_\odot$ the masses of all single BHs that are retained at the start and the end of the simulation.

<a name="applicationsofthecode"></a>
### 7. Applications of the code

The code can be useful when executed multiple times, for instance when simulating a set of clusters and generating a population of dynamically formed BBH mergers.

Although the program itself is not computationally expensive (we have tested in a laptop that we generate a few binary black hole mergers per second), independent parallelization is still encouraged when simulating a very large number of star clusters for efficiency.

<a name="citingthiswork"></a>
### 8. Citing this work

If you utilize this code in your research, please cite the following reference:

[K. Kritos, V. Strokov, V. Baibhav & E. Berti (2022)](https://arxiv.org/abs/2210.10055).

$\tt Rapster$ has been used in the following works:

- [K. Kritos, E. Berti & J. Silk (2022)](https://arxiv.org/abs/2212.06845)

- [A. Antonelli, K. Kritos, K. Ng, R. Cotesta, E. Berti (2023)](https://arxiv.org/abs/2306.11088)

- [K. Ng, K. Kritos, A. Antonelli, R. Cotesta, E. Berti (2023)](https://arxiv.org/abs/2307.03227)

<a name="reportingbugs"></a>
### 9. Reporting bugs

If you find a bug in the code, please contact us in kkritos1@jhu.edu with a description of the bug.

Suggestions and pull requests are welcome :)

<a name="thanks"></a>
### 10. Thank you

Vladimir Strokov, Vishal Baibhav, Emanuele Berti, Andrea Antonelli, Fabio Antonini, Dany Atallah, Muhsin Aljaf, Mark Cheung, Roberto Cotesta, Hector Cruz, Giacomo Fragione, Gabriele Franciolini, Thomas Helfer, Veome Kapil, Kyle Kremer, Iason Krommydas, Miguel Martinez, Luca Reali, Carl Rodriguez, Newlin Weatherford, Xiao-Xiao Kou.


