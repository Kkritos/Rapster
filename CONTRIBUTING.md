# Contributing to Rapster

Thanks for your interest in Rapster! Bug reports, physics corrections, new channels, and documentation improvements are all welcome.

## Reporting bugs
Open an issue using the **Bug report** template. Please include the full command line (all arguments), the Rapster version or commit hash, and the complete traceback. If the problem only appears for certain cluster parameters, include them.

## Suggesting features or physics changes
Open an issue using the **Feature request** template before writing code. For changes to the physics (prescriptions, rates, fitting formulae), please cite the relevant references.

## Development workflow
1. Fork the repository and create a branch from `main`.
2. Make your changes. Keep each pull request focused on a single change.
3. Run the tests with `python3 -m pytest tests/` (install with `pip install pytest`). The same tests run automatically on GitHub for every pull request and must pass before merging.
   - If your change is *meant* to alter results (new physics, a bug fix), `tests/test_regression.py` will fail. Regenerate its reference values with `python3 tests/test_regression.py --update`, paste them into the test, and explain the change in the pull request.
4. Update the README and docstrings if you add or change user-facing arguments.
5. Open a pull request and fill in the template.

## Working with the code
- **Single BHs** (and retained NSs) are stored in `state['BH']`, a `BHPool` (`rapster/bhpool.py`), not in separate `mBH`, `sBH`, `gBH`, `hBH` arrays. Read them through `BH.m`, `BH.s`, `BH.g`, `BH.h` (views of the first `BH.size` slots) or `BH.get(k)`. Change them only through the pool's methods: `BH.add(m, s, g, h)` (returns the new slot), `BH.remove(k)`, `BH.remove_two(k1, k2)`, `BH.replace(k, m, s, g, h)`, `BH.set_mass(k, m)`, `BH.set_spin(k, s)`, `BH.set_tde_count(k, h)` and `BH.set_all_masses(m)`. Writing into `BH.m` directly would leave the pool's weights and sums out of date.
- **Slots move**: removing a BH moves the last BH into the freed slot, so a slot number is valid only until the next removal.
- **Weighted draws** from the single BHs use `BH.draw(FAM['name'])` and `BH.total(FAM['name'])`, with the weight families listed in `FAMILIES` (`bhpool.py`). A sampler that needs a new weight (a new power of the BH mass) adds a family there.
- **Random numbers** must come from the global NumPy state (`np.random.*`), which `-S` seeds, so that runs are reproducible; do not create separate generators such as `np.random.default_rng()`.
- **Performance changes**: if a change is meant to leave results unchanged, run the old and new code with the same seed on a few configurations (small and compact clusters, with and without NSs, gas, seed BH) and compare all output files byte for byte. If it only keeps the distributions unchanged (e.g. a new sampler), say so, test the sampler against the exact distribution (as in `tests/test_bhpool.py`), and compare ensembles of runs with different seeds.

## Code style
- Follow PEP 8 where practical.
- Use descriptive variable names for physical quantities, and state units in comments or docstrings.
- New command-line options go through `argparse` with a sensible default and help text.

## Questions
For questions about using Rapster, open an issue or contact the maintainer (see README).
