## -- Basic variables
PROJECT_PATH := $(dir $(abspath $(firstword $(MAKEFILE_LIST))))
SHELL := /usr/bin/env bash
BUMP_TOOL := bump-my-version
DOCKER_COMPOSE := docker compose

## -- Conda variables
# CONDA_TOOL can be overridden in Makefile.private file
CONDA_TOOL ?= conda
CONDA_ENVIRONMENT := climateset


## -- Private variables import to override variables for local 
-include Makefile.private

#### Install ####

.PHONY: install-conda
install-conda: ## Install conda on your local machine
	@$(CONDA_TOOL) --version; \
	if [ $$? != "0" ]; then \
		echo " "; \
		echo "Your defined Conda tool [$(CONDA_TOOL)] has not been found."; \
		echo " "; \
		echo "If you know you already have [$(CONDA_TOOL)] or some other Conda tool installed,"; \
		echo "Check your [CONDA_TOOL] variable in the Makefile.private for typos."; \
		echo " "; \
		echo "If your conda tool has not been initiated through your .bashrc file,"; \
		echo "consider using the full path to its executable instead when"; \
		echo "defining your [CONDA_TOOL] variable"; \
		echo " "; \
		echo "If in doubt, don't install Conda and manually create and activate"; \
		echo "your own Python environment."; \
		echo " "; \
		echo -n "Would you like to install Miniconda ? [y/N]: "; \
		read ans; \
		case $$ans in \
			[Yy]*) \
				echo "Fetching and installing miniconda"; \
				echo " "; \
				wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh; \
    			bash ~/miniconda.sh -b -p $HOME/.conda; \
    			export PATH=$HOME/.conda/bin:$PATH; \
    			conda init; \
				/usr/bin/rm ~/miniconda.sh; \
				;; \
			*) \
				echo "Skipping installation."; \
				echo " "; \
				;; \
		esac; \
	else \
		echo "Conda tool [$(CONDA_TOOL)] has been found, skipping installation"; \
	fi;


.PHONY: create-env
create-env: install-conda ## Create a local conda environment
	$(CONDA_TOOL) env create -f environment.yml
	@echo ""
	@echo "#####################################################################"
	@echo "#                                                                   #"
	@echo "# Please activate your new environment using the following command: #"
	@echo "#                                                                   #"
	@echo "#####################################################################"
	@echo ""
	@echo "$(CONDA_TOOL) activate $(CONDA_ENVIRONMENT)"
	@echo ""
	@echo ""

.PHONY: install
install: install-precommit ## Install the application and it's dependencies

.PHONY: install-poetry
install-poetry: ## Install the application and it's dependencies
	@poetry --version; \
	if [ $$? != "0" ]; then \
		echo "Poetry not found, proceeding to install Poetry..."; \
		pip install poetry==1.8.2; \
	fi;

.PHONY: install-precommit
install-precommit: install-dev## Install the application and it's dependencies
	@if [ -f .git/hooks/pre-commit ]; then \
		echo "Pre-commit hook found"; \
	else \
	  	echo "Pre-commit hook not found, proceeding to configure it"; \
		poetry run pre-commit install; \
	fi;

.PHONY: install-package
install-package: install-poetry ## Install the application and it's dependencies
	@poetry install

.PHONY: install-dev
install-dev: install-poetry ## Install the application and it's dependencies
	@poetry install --with dev

.PHONY: install-with-lab
install-with-lab: install-poetry ## Install the application and it's dependencies, including Jupyter
	poetry install --with lab

#### Versionning ####

# Use the "dry" target for a dry-run version bump, ex.
# make bump-major dry
BUMP_ARGS ?= --verbose 
ifeq ($(filter dry, $(MAKECMDGOALS)), dry)
	BUMP_ARGS := $(BUMP_ARGS) --dry-run --allow-dirty
endif

.PHONY: bump-major
bump-major: ## Bump application major version  <X.0.0>
	$(BUMP_TOOL) $(BUMP_ARGS) major

.PHONY: bump-minor
bump-minor: ## Bump application major version  <0.X.0>
	$(BUMP_TOOL) $(BUMP_ARGS) minor

.PHONY: bump-patch
bump-patch: ## Bump application major version  <0.0.X>
	$(BUMP_TOOL) $(BUMP_ARGS) patch

#### Docker ####

#### Linting ####

.PHONY: check-lint
check-lint: ## Check code linting
	poetry run tox -e black,isort,flake8,pylint

.PHONY: check-flake8
check-flake8: ## Check code linting with only flake8
	poetry run tox -e flake8

.PHONY: check-pylint
check-pylint: ## Check code linting with only flake8
	poetry run tox -e pylint

.PHONY: pre-commit
pre-commit: ## Fix code linting
	poetry run tox -e precommit

.PHONY: fix-lint
fix-lint: ## Fix code linting
	poetry run tox -e fix

#### Tests ####

.PHONY: test
test: ## Run tests
	poetry run tox -e test

.PHONY: test-all
test-all: ## Run tests
	poetry run tox -e test-all

#### Application ####

