# Changelog

All notable changes to Rapster will be documented in this file.

## [Unreleased]

## [2.11.4] - 2026-09-15

### Added
- Beta spin distribution option (`-SD 2`). Draws natal BH spins from a Beta(1.4, 3.6) distribution scaled by `s1g_max`. Usage: `python -m rapster.run_cluster -SD 2 -s <max_spin>`.
- This CHANGELOG file to track changes to the project.
- Random mass pairing option (`-RMP 1` / `--random_mass_pairing_2body_3body 1`). When enabled, 3-body binary formation and 2-body capture use uniform random pairing instead of mass-weighted (m^5 and m^2 respectively). Default is mass-weighted (`-RMP 0`).
- Automatic diagnostic plots (`-plot 1`). Generates 10 PNG plots (cluster evolution, radii, cluster mass, BH mass function, merger masses, merger channels, merger spins, eccentricity by channel, TDEs, hardening) in `Results/plots/`. New `plot_cluster.py` module.
- Post-simulation analysis summary via `analyze_cluster()` (`-analyze 1`). Prints merger statistics (total, in-cluster, ejected, per channel, retained), maximum dynamically-formed BH mass, BH generation counts, TDE summary, and final cluster state. New `analyze_cluster.py` module.
- All output is now saved to `Results/log.txt`. With `-P 1` output goes to both screen and log; with `-P 0` output goes only to the log file.
- Initial BH mass distribution options (`-BMD`): `0` for Kroupa+collapse (default), `1` for uniform, `2` for Salpeter power law (m^-2.35), `3` for log-uniform. BH count is always determined from the Kroupa IMF; for BMD>0 only the masses are resampled in [`-mBH1gMin`, `-mBH1gMax`] (default [3, 60] Msun). Momentum-conservation SN kicks are applied for BMD>0 (fallback kicks are unavailable since there is no stellar progenitor).
- Residual gas and gas accretion physics via the new `compact_accretion.py` module. Models Eddington-limited gas accretion onto compact objects (BHs and NSs) during the cluster's residual-gas phase, including spin evolution during accretion from ISCO energetics and NS structure/equation-of-state handling.
- Residual-gas model flags: `-sfe` (star formation efficiency ε ∈ (0,1], default `1.0`; sets the initial residual gas mass, with `1.0` recovering gas-free Rapster), `-fEdd` (Eddington ratio ceiling for gas accretion, default `1.0`), `-fge` (gas expulsion timescale in units of the initial crossing time, default `5.0`), and `-cs` (gas sound speed [km/s], default `10.0`).
- Neutron star equation of state option (`-EoS`, either `APR` or `AU`, default `APR`). Adds NS radius–mass lookup tables `Data/APR.txt`, `Data/APR4.txt`, `Data/AU.txt`, and the `Data/.tovseq_to_eos.py` helper used to generate them.
- GW recoil kick model option (`-RK`): `0` for Gerosa & Kesden (2016, default), `1` for `gwModel_kick_prec_flow` from Islam & Wadekar (2025). Model 1 requires the optional dependency `gwModels[kicks]`.
- TDE sampling flags: `-fA` (fraction of a disrupted star accreted by the compact object, default `0.5`) and `-mb` (mass bias power index p for drawing stars from IMF·m^p for TDEs, default `0.0`).
- Bimodal neutron star mass distribution (Rocha et al. 2023) via `-NS 2`, alongside the existing monochromatic 1.4 Msun option (`-NS 1`).
- BBH–star tidal disruption channels: TDEs occurring during binary–single and binary–binary strong interactions, including hard and soft μTDEs of a star by a BH binary (TDE types 4, 21, 22). Adds an `N_tdeBBHstar` counter to the evolution output and TDE file.
- Micro-TDE (μTDE) treatment during binary–single encounters, with dedicated stellar-mass sampling for disrupted stars.

### Changed
- Refactored `run_cluster.py` from a monolithic script into modular functions: `parse_args()`, `initialize_cluster()`, `compute_cluster_properties()`, `compute_timescales()`, `form_binaries()`, `evolve_interactions()`, `compute_external_params()`, `evolve_tdes()`, `record_evolution()`, `update_cluster()`, `print_status()`, and `write_output()`. All simulation state is bundled in a `state` dictionary. No logic changes.
- Added docstrings to all 12 functions in `run_cluster.py` describing their purpose, arguments, and return values.
- Resolved TODO comments in `run_cluster.py`: added visual delimiter before `__main__` block and inline docstring.
- Replaced all 27 `ADDME` placeholder comments in `run_cluster.py` with descriptive comments explaining each code block.
- Added `.ipynb_checkpoints/` to `.gitignore`.
- Split `run_cluster.py` into `cluster_evolution.py` (11 evolution functions) and `run_cluster.py` (`parse_args()` + `main()` entry point). Command `python -m rapster.run_cluster` still works as before.
- Replaced 36 repetitive `np.load` lines in `stellar_evolution.py` with a `_load_grid()` helper function. Data loading for SEVN delayed/rapid remnant masses and CO core masses is now 3 one-liners.
- Added `#`-prefixed column headers to all output `.txt` files (mergers, evolution, hardening, tdes). Compatible with `np.loadtxt` which skips `#` comment lines.
- Updated README input parameters table with new flags (`-SD 2`, `-RMP`, `-plot`, `-analyze`) and output files section (`log.txt`, `plots/`).
- Added usage examples to README section 5 (Running a simulation).
- Updated README input parameters table with `-BMD`, `-mBH1gMin`, `-mBH1gMax` flags.
- Regenerated `Example/Results_Test/` with current code (includes column headers and log.txt).
- Performance: substantial speedups across the hot paths — vectorized `IMF_kroupa` (hybrid scalar + `np.select` array path); a faster, bit-identical TDE spin integrator (`evolve_spin_during_accretion`, formerly `evolve_v2`); list-accumulation of the evolution and hardening arrays instead of `np.append` (O(N^2) → O(N)); removal of `np.vectorize` from the remnant-mass functions in favor of their vectorized array branch; and `np.clip` replaced with `min`/`max` in the `compact_accretion` hot path.
- BH/NS spins are now clipped to the Thorne limit χ ≤ 0.998 via a new `THORNE_SPIN_LIMIT` constant (default `0.998`) in `constants.py`, applied during accretion spin-up.
- Escape velocity now accounts for the residual gas mass: `M_gas` is included in the cluster potential when computing `v_esc`, and the gas mass at the time of BH formation is used when comparing natal SN kicks against the escape velocity.
- TDE stellar-mass sampling now includes a rate-dependent mass factor in p(m_star | TDE), and pericenter (`r_p`) sampling is unified across single–BBH interactions (with or without a μTDE).
- Poisson rate parameter (μ = dt/t_channel) is now capped for all interaction channels to prevent overflow in `poisson.rvs`.

### Fixed
- Corrected a sign error in the average-mass evolution: the stellar mass-loss factor was `(t/t_sev)^ν_sev` and is now correctly `(t/t_sev)^(−ν_sev)`. This affects the cluster's average-mass evolution and everything downstream of it.
- Corrected pericenter sampling: `r_p` is now sampled with probability ∝ `r_p` in the gravitational-focusing regime (for both single–BBH and `N_IMS−1` cases), rather than uniformly.
- Fixed initialization of NS masses when the bimodal distribution (`-NS 2`) is selected.
- Fixed redshift–lookback interpolation: added range checks when interpolating `z_merge` for 2-body capture, 3-body GW, and ejected mergers, and handled negative lookback times.
- Fixed the unconditional definition of `M_gas_BHform` in the gas/accretion path.
- Fixed numerous index/channel bugs in `cluster_evolution.py`, `binary_evolution.py`, `exchanges.py`, `triples.py`, and `tidal_disruptions.py` — including the pair–pair and ex1/ex2 channels, `k1`/`k2`/`kss`/`kp` index handling when creating pairs, a non-existent third BH in the BH subcluster, `a1`/`a2` semimajor-axis sampling, and correct unpacking of binary properties (m, s, g, h) after a TDE.
- Fixed TDE bookkeeping during star–star–BH and BH–star–BH interactions.
- Fixed IMF integration and TDE plotting (scipy `quad`), and ensured mergers in `triples.py` occur by z = 0.
- Renamed `type` variable to `tde_type` in `tidal_disruptions.py` and `cluster_evolution.py` to avoid shadowing Python's builtin `type()`.
- Populated `__init__.py` with public API exports. Users can now do `from rapster import initialize_cluster, analyze_cluster, generate_all_plots` etc. CLI-specific functions (`parse_args`, `main`) remain in `run_cluster.py`.
- Replaced `np.transpose(array)[:][i]` with `array[:, i]` across `cluster_evolution.py` (15 instances), `binary_evolution.py` (4 instances), and `exchanges.py` (3 instances). Standard NumPy column indexing idiom — more readable and avoids unnecessary transpose.
- Replaced all 10 bare `except:` clauses in `cluster_evolution.py` with `except Exception:`. Bare `except:` catches everything including `KeyboardInterrupt` and `SystemExit`, which means Ctrl-C could not stop the simulation and errors were silently ignored. `except Exception:` still handles expected errors (e.g., `ValueError`, `ZeroDivisionError` from empty arrays or zero rates) but allows Ctrl-C and system exits to propagate normally.
