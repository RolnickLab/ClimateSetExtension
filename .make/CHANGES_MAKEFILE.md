# CHANGES.md

## [Unreleased](https://github.com/RolnickLab/lab-advanced-template/tree/main) (latest)

______________________________________________________________________

<!-- (New changes here in list form) -->

## [1.3.0](https://github.com/RolnickLab/lab-advanced-template/tree/makefile-1.3.0) (2026-02-12)

______________________________________________________________________

- Fix test targets
- Move `docs` and `lab` dependency groups as optional groups `[extras]`
- Remove upper bounds of dependencies
- Add `test-notebooks` target to facilitate running tests using notebooks through the `nbval` library
- Add `self-test` and `self-test-autoinit` targets to Makefile utils
- Add new test for auto-initialization script
- Add `mypy` and `autotyping`

## [1.2.0](https://github.com/RolnickLab/lab-advanced-template/tree/makefile-1.2.0) (2026-01-21)

______________________________________________________________________

- Add auto-initialization script and corresponding makefile targets

## [1.1.0](https://github.com/RolnickLab/lab-advanced-template/tree/makefile-1.1.0) (2026-01-20)

______________________________________________________________________

- Add MkDocs dependencies and skeleton structure for MkDocs pages
- Add `docs` makefile targets
- Refactor base package to `src/core` instead of `src/` and improve package structure to follow current python best practices
- Improve and fix bugs/typos from the conda, poetry and uv targets

## [1.0.0](https://github.com/RolnickLab/lab-advanced-template/tree/makefile-1.0.0) (2025-11-21)

______________________________________________________________________

- BREAKING CHANGE - Make default version of project use `uv`
  - Consists of a rework of the `pyproject.toml` file that no longer works with `poetry<2.0.0`
- Added `poetry python install` functionality to the makefike
- Refactored `conda` installation to use miniforge and micromamba instead of miniconda
- Improve determination of build tool and environment by makefile to make experience simpler
- Refactor target group enablement via `Makefile.variables` file instead of commenting
  out lines in `Makefile`
- Add link checker to `pre-commit`
- Refactor tests to reduce duplication
- Remove target that installed `poetry` inside conda environment
- Update documentation and README.md
- Convert to Google docstring format

## [0.7.1](https://github.com/RolnickLab/lab-advanced-template/tree/makefile-0.7.1) (2025-09-17)

______________________________________________________________________

## [0.7.0](https://github.com/RolnickLab/lab-advanced-template/tree/makefile-0.7.0) (2025-09-10)

______________________________________________________________________

- Add `mdformat` tool for markdown linting
- Modularize the Makefile structure, allowing to choose which tools are available
- Add tests for the Makefile

## [0.6.0](https://github.com/RolnickLab/lab-advanced-template/tree/makefile-0.6.0) (2025-05-29)

______________________________________________________________________

- Improve venv support
- Add venv remove target
- Refactor `poetry-install-auto` to really be automatic
- Add `mamba-install` target as an alternative/complement to conda for local dev
- Add `autoflake`, `autopep8` and `ruff` targets
- Fix some typos that caused targets to fail

## [0.5.0](https://github.com/RolnickLab/lab-advanced-template/tree/makefile-0.5.0) (2025-03-11)

______________________________________________________________________

- Add venv support
- Update and fix `poetry-install-auto` target

## [0.4.0](https://github.com/RolnickLab/lab-advanced-template/tree/makefile-0.4.0) (2024-10-30)

______________________________________________________________________

- Add cyclomatic complexity check target

## [0.3.0](https://github.com/RolnickLab/lab-advanced-template/tree/makefile-0.3.0) (2024-10-30)

______________________________________________________________________

- Add utilities to track makefile versions
- Fix many of the targets relating to standalone `poetry` installation
