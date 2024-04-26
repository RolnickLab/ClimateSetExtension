# ClimatesetExtension

## Requirements

This project has only been tested in a Linux (Debian based) environment and assumes
some basic tools for development are already installed.

The project uses a Makefile to automate most operations. If you can run make on your 
machine there's a good chance this will work on your machine.

Please consult the [Makefile](Makefile) before trying to use it.

## Installation

This project assumes environment management will be done with Conda. It would, however,
be possible to create you own environment with, for example, poetry, pyenv or virtualenv.

**Do note that this project is configured to use [Poetry](https://python-poetry.org/)
and as such, this tool is required to install it.**

Poetry is included in the [environment.yml](environment.yml) used by the commands below to
create a conda environment.

If you need or want to install Conda:
```
make install-conda 
```

To create the conda environment:
```
make create-env
```

After you have created your Conda environment, or you want to manage your environment 
yourself, make sure to activate it properly before installing this package.

To install the package and its CLI tools:
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

To run a `pre-commit` check before actually commiting:
```
make pre-commit
```

To run tests:
```
make test
```

To run all tests, including the integration test:
```
make test-all
```


## Data

## Experiment tracking

Nothing is set up for now, but since Weights and Bias is accessible to MILA and DRAC, it
will probably be the way to go.


## Training

## Contributing to this repository

### Design patterns
Two main considerations have been made in the structure of this package.

First, a polymorphic approach, using abstract classes and their concrete implementation,
has been but forward in order to increase maintainability and extensibility.

Therefore, new additions should try to follow this design pattern and either implement
new concrete classes or create new abstract classes and their implementations for 
completely new behavior or needs.

Secondly, a dependency-injection approach is to be favored, as well as a composition 
approach when creating new modules or extending existing ones.

### Configurations
Configurations are in the [config/](config) folder.


### Tests
