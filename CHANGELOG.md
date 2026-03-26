# Changelog

All notable changes to Rapster will be documented in this file.

## [Unreleased]

### Added
- Beta spin distribution option (`-SD 2`). Draws natal BH spins from a Beta(1.4, 3.6) distribution scaled by `s1g_max`. Usage: `python -m rapster.run_cluster -SD 2 -s <max_spin>`.
- This CHANGELOG file to track changes to the project.
- Random mass pairing option (`-RMP 1` / `--random_mass_pairing_2body_3body 1`). When enabled, 3-body binary formation and 2-body capture use uniform random pairing instead of mass-weighted (m^5 and m^2 respectively). Default is mass-weighted (`-RMP 0`).

### Changed
- Refactored `run_cluster.py` from a monolithic script into modular functions: `parse_args()`, `initialize_cluster()`, `compute_cluster_properties()`, `compute_timescales()`, `form_binaries()`, `evolve_interactions()`, `compute_external_params()`, `evolve_tdes()`, `record_evolution()`, `update_cluster()`, `print_status()`, and `write_output()`. All simulation state is bundled in a `state` dictionary. No logic changes.
- Added docstrings to all 12 functions in `run_cluster.py` describing their purpose, arguments, and return values.
- Resolved TODO comments in `run_cluster.py`: added visual delimiter before `__main__` block and inline docstring.
- Replaced all 27 `ADDME` placeholder comments in `run_cluster.py` with descriptive comments explaining each code block.
- Added `.ipynb_checkpoints/` to `.gitignore`.
