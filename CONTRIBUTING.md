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

## Code style
- Follow PEP 8 where practical.
- Use descriptive variable names for physical quantities, and state units in comments or docstrings.
- New command-line options go through `argparse` with a sensible default and help text.

## Questions
For questions about using Rapster, open an issue or contact the maintainer (see README).
