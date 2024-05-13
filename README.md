# ClimatesetExtension

## Requirements

This project has only been tested in a Linux (Debian based) environment and assumes
some basic tools for development are already installed.

The project uses a Makefile to automate most operations. If make is available on your 
machine there's a good chance this will work.

The following Makefile files should not be modified, but can be consulted:

* [Makefile](Makefile) : Orchestration of the different files
* [base.make](.make/base.make) : Shared utilities, project agnostic.

The following Makefile files are project or user specific:

* [Makefile.variables](Makefile.variables) : Shared project variables.
* [Makefile.targets](Makefile.targets) : Shared project targets.
* [Makefile.private](Makefile.private.example) : User specific variables and targets.

## Basic Information

The different targets and their description can be examined by executing the command
`make targets`

![](img/make_targets.png)

## Installation

This project assumes environment management will be done with `Conda` or directly through
`Poetry`. 

* [Poetry](https://python-poetry.org/docs/basic-usage/)
* [Conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html)

While it is possible to manage the environment with, for example, pyenv or virtualenv, 
those specific use cases are not supported by the Makefile and require users to set up 
their own environments beforehand.

For detailed information about `Poetry` and `Conda`:

If you want to use something else than `Conda` or `Poetry` to manage environment isolation, 
it is recommended to follow 
[Poetry's guidelines on managing environments](https://python-poetry.org/docs/managing-environments/)

Poetry is not included in the [environment.yml](environment.yml), due to some possible problems
in compute cluster environments, but will be installed automatically if needed
by most `install` targets.

Currently, the project runs on Python version 3.10.

### Environment management choices

Environment management can become quite complicated. Using Conda allows a certain
ease of management since the Poetry installation is contained inside the created Conda 
environment.

However, some computing environments do not permit the use of Conda (like certain SLURM
clusters). This is why the `pipx` option for Poetry is also enabled in this project.

**Unless you really know what you are doing, it is not recommended to install Poetry
as a standalone tool (with pipx) while also using Conda environments.**

### Conda targets

If you need or want to install Conda:
```
make conda-install 
```

To create the conda environment:
```
make conda-create-env
```

To remove the conda environment:
```
make conda-clean-env
```

Make sure to activate the configured environment before installing this package.

### Poetry targets

There are 2 possibilities available in how to manage `Poetry` if it is not already
configured on your system.

**It is not recommended to install `Poetry` in the same environment that will be managed
by `Poetry`; It should be installed by `Conda` (as in, `conda install poetry`) or 
on the system using `pipx`, in order to minimize dependency conflicts.**

The following target will first try to install `Poetry` in the active `Conda` 
environment; if it can't find `Conda`, it will proceed to install via `pipx`

```
make poetry-install-auto
```

To install `Poetry` using `Conda`:

```
make conda-poetry-install
```

Using `pipx` will instead allow environment management directly with `Poetry`. The
following target will also make `Poetry` use a Python 3.10 environment for this
project.

```
make poetry-install-pipx
```

This will also create a virtual environment managed by `Poetry`.

A standalone environment can also be created later using the `make poetry-create-env`
command, and removed with the `make poetry-remove-env` command.

Information about the currently active environment used by Poetry, 
whether Conda or Poetry, can be seen using the `make poetry-env-info` command.

Both install methods can also be cleaned up:

```
make conda-poetry-uninstall
# or
make poetry-uninstall-pipx
```

**Important note!**

If you have an active `Conda` environment  and install `Poetry` using `pipx`,
you will have to use `poetry run python <your_command_or_script_path>` instead of 
`python <your_command_or_script_path>`, (which is normal when using Poetry) as 
`python` will use Conda's active environment.

### Install targets

**All `install` targets will first check if `Poetry` is available and try to install
it with the `make poetry-install-auto` target.**

To install the package, development dependencies and CLI tools (if available):
```
make install
```

To install only the package, without development tools:
```
make install-package
```

## First time user quick setup

### Conda
The easiest and quickest way to get up and running with Conda.

Create Conda environment (will check for Conda and install it if not found):

```
make conda-create-env
```

Activate Conda environment (substitute with your <CONDA_TOOL> if something else 
than `conda`:

```
conda activate climateset
```

Install package:

```
make install
```

### Poetry

The easiest and quickest way to get up and running with Poetry.

Install pipx and Poetry and activate project environment :

```
make poetry-install
```

**Or, if Poetry is already available:**

```
make poetry-create-env
```

Install package:

```
make install
```

## Basic automations

To run linting checks with `flake8`, `pylint`, `black` and `isort`:
```
make check-lint
```

To fix linting with `black`, `flynt` and `isort`:
```
make fix-lint
```

To run a `pre-commit` check before actually committing:
```
make precommit
```

To run tests:
```
make test
```


## Data

## Experiment tracking

Nothing is set up for now, but since Weights and Bias is accessible to MILA and DRAC, it
will probably be the way to go.


## Training

## Contributing to this repository

See [Contributing guidelines](CONTRIBUTING.md)


### Configurations
Configurations are in the [config/](config) folder.


### Tests
