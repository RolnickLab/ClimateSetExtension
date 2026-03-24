# Environment Setup Guide

This guide covers the necessary steps to configure your environment for the ClimatesetExtension project.

## Requirements


This project is configured to use `Miniforge` or `Micromamba` for environment isolation and non-Python dependencies, while using [Poetry](https://python-poetry.org/) as the primary build tool for managing Python packages. You can choose which one to install using the `make conda-install` command.

The project is configured to use the `mamba` tool instead of `conda` or `micromamba`. See [Makefile.private.example](Makefile.private.example) if you require something other than `mamba` locally.

## Makefile Automation

The project relies on a `Makefile` to automate common tasks, including environment management, package installation, linting, testing, and documentation generation.

For a comprehensive list of targets and detailed explanations of the Makefile structure, please consult the [Makefile Documentation](.make/README.md). 

You can discover all available commands and view your current project configuration directly in your terminal:

```bash
# List all available Makefile targets
make targets

# View current project configuration (Python version, active target groups, etc.)
make info
```

## Quick Setup

Follow these steps to initialize the environment using the configured Conda and Poetry toolchain.

1. **Create the Conda environment**:
   This step handles non-Python dependencies and isolates the environment.
   ```bash
   make conda-create-env
   ```
2. **Activate the environment**:
   ```bash
   eval $(make conda-activate)
   ```
3. **Install the Python package and development tools via Poetry**:
   This target will automatically verify Poetry is available and use it to install the package dependencies.
   ```bash
   make install
   ```

## Experiment Tracking

Weights and Biases (WandB) is the planned solution for experiment tracking, as it is natively accessible on MILA and DRAC clusters.

## Additional Resources

- **Configurations**: Review the [configs/](configs) directory for YAML configurations.
- **Contributing**: Refer to the [Contribution Guidelines](CONTRIBUTING.md) for repository standards.
