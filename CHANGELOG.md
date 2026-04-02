# Changelog

All notable changes to Rapster will be documented in this file.

## [Unreleased]

### Added
- Beta spin distribution option (`-SD 2`). Draws natal BH spins from a Beta(1.4, 3.6) distribution scaled by `s1g_max`. Usage: `python -m rapster.run_cluster -SD 2 -s <max_spin>`.
- This CHANGELOG file to track changes to the project.
- Random mass pairing option (`-RMP 1` / `--random_mass_pairing_2body_3body 1`). When enabled, 3-body binary formation and 2-body capture use uniform random pairing instead of mass-weighted (m^5 and m^2 respectively). Default is mass-weighted (`-RMP 0`).
- Automatic diagnostic plots (`-plot 1`). Generates 10 PNG plots (cluster evolution, radii, cluster mass, BH mass function, merger masses, merger channels, merger spins, eccentricity by channel, TDEs, hardening) in `Results/plots/`. New `plot_cluster.py` module.
- Post-simulation analysis summary via `analyze_cluster()` (`-analyze 1`). Prints merger statistics (total, in-cluster, ejected, per channel, retained), maximum dynamically-formed BH mass, BH generation counts, TDE summary, and final cluster state. New `analyze_cluster.py` module.
- All output is now saved to `Results/log.txt`. With `-P 1` output goes to both screen and log; with `-P 0` output goes only to the log file.

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

- Initial BH mass distribution options (`-BMD`): `0` for Kroupa+collapse (default), `1` for uniform, `2` for Salpeter power law (m^-2.35), `3` for log-uniform. BH count is always determined from the Kroupa IMF; for BMD>0 only the masses are resampled in [`-mBH1gMin`, `-mBH1gMax`] (default [3, 60] Msun). Momentum-conservation SN kicks are applied for BMD>0 (fallback kicks are unavailable since there is no stellar progenitor).

### Fixed
- Renamed `type` variable to `tde_type` in `tidal_disruptions.py` and `cluster_evolution.py` to avoid shadowing Python's builtin `type()`.
- Populated `__init__.py` with public API exports. Users can now do `from rapster import initialize_cluster, analyze_cluster, generate_all_plots` etc. CLI-specific functions (`parse_args`, `main`) remain in `run_cluster.py`.
- Replaced `np.transpose(array)[:][i]` with `array[:, i]` across `cluster_evolution.py` (15 instances), `binary_evolution.py` (4 instances), and `exchanges.py` (3 instances). Standard NumPy column indexing idiom — more readable and avoids unnecessary transpose.
- Replaced all 10 bare `except:` clauses in `cluster_evolution.py` with `except Exception:`. Bare `except:` catches everything including `KeyboardInterrupt` and `SystemExit`, which means Ctrl-C could not stop the simulation and errors were silently ignored. `except Exception:` still handles expected errors (e.g., `ValueError`, `ZeroDivisionError` from empty arrays or zero rates) but allows Ctrl-C and system exits to propagate normally.
