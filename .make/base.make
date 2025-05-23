########################################################################################
# DO NOT MODIFY!!!
# If necessary, override the corresponding variable and/or target, or create new ones
# in one of the following files, depending on the nature of the override :
#
# Makefile.variables, Makefile.targets or Makefile.private,
#
# The only valid reason to modify this file is to fix a bug or to add new
# files to include.
#
# Please report bugs to francis.pelletier@mila.quebec
########################################################################################

# Basic variables
PROJECT_PATH := $(dir $(abspath $(firstword $(MAKEFILE_LIST))))
MAKEFILE_NAME := $(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))
SHELL := /usr/bin/env bash
BUMP_TOOL := bump-my-version
MAKEFILE_VERSION := 0.5.0
DOCKER_COMPOSE ?= docker compose
AUTO_INSTALL ?=

# Conda variables
# CONDA_TOOL can be overridden in Makefile.private file
CONDA_TOOL := conda
CONDA_ENVIRONMENT ?=
CONDA_YES_OPTION ?=

# Default environment to install package
# Can be overridden in Makefile.private file
DEFAULT_INSTALL_ENV ?=
DEFAULT_POETRY_INSTALL_ENV ?=

# Colors
_SECTION := \033[1m\033[34m
_TARGET  := \033[36m
_NORMAL  := \033[0m

.DEFAULT_GOAL := help

# Project and Private variables and targets import to override variables for local
# This is to make sure, sometimes the Makefile includes don't work.
-include Makefile.variables
-include Makefile.private
## -- Informative targets ------------------------------------------------------------------------------------------- ##

.PHONY: all
all: help

# Auto documented help targets & sections from comments
#	detects lines marked by double #, then applies the corresponding target/section markup
#   target comments must be defined after their dependencies (if any)
#	section comments must have at least a double dash (-)
#
# 	Original Reference:
#		https://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
# 	Formats:
#		https://misc.flogisoft.com/bash/tip_colors_and_formatting
#
#	As well as influenced by it's implementation in the Weaver Project
#		https://github.com/crim-ca/weaver/tree/master

.PHONY: help
# note: use "\#\#" to escape results that would self-match in this target's search definition
help: ## print this help message (default)
	@echo ""
	@echo "Please use 'make <target>' where <target> is one of below options."
	@echo ""
	@for makefile in $(MAKEFILE_LIST); do \
        grep -E '\#\#.*$$' "$(PROJECT_PATH)/$${makefile}" | \
            awk 'BEGIN {FS = "(:|\-\-\-)+.*\#\# "}; \
            	/\--/ {printf "$(_SECTION)%s$(_NORMAL)\n", $$1;} \
				/:/  {printf "    $(_TARGET)%-24s$(_NORMAL) %s\n", $$1, $$2} ' 2>/dev/null ; \
    done

.PHONY: targets
targets: help

.PHONY: version
version: ## display current version
	@echo "version: $(APP_VERSION)"

## -- Virtualenv targets -------------------------------------------------------------------------------------------- ##

VENV_PATH := $(PROJECT_PATH).venv
VENV_ACTIVATE := $(VENV_PATH)/bin/activate

.PHONY: venv-create
venv-create: ## Create a virtualenv '.venv' at the root of the project folder
	@virtualenv $(VENV_PATH)
	@make -s venv-activate

.PHONY: venv-activate
venv-activate: ## Print out the shell command to activate the project's virtualenv.
	@echo "source $(VENV_ACTIVATE)"

## -- Poetry targets ------------------------------------------------------------------------------------------------ ##

.PHONY: poetry-install-auto
poetry-install-auto: ## Install Poetry in Conda environment, or with pipx in a virtualenv if Conda not found
	@poetry --version; \
    	if [ $$? != "0" ]; then \
			echo "Poetry not found, proceeding to install Poetry..."; \
			if [ "$(DEFAULT_POETRY_INSTALL_ENV)" == "conda" ]; then \
				ans_where="conda"; \
			elif [ "$(DEFAULT_POETRY_INSTALL_ENV)" == "venv" ]; then \
				ans_where="venv"; \
			else\
				echo -n "Where would you like to install Poetry, in a dedicated virtualenv (venv), or a conda environment? [venv/conda]: "; \
				read ans_where; \
			fi; \
			case $$ans_where in \
				"venv" | "Venv" |"VENV") \
					make AUTO_INSTALL=true -s poetry-install-venv; \
					;; \
				"conda" | "Conda" | "CONDA") \
					echo "Installing poetry with Conda"; \
					make AUTO_INSTALL=true -s conda-poetry-install; \
					;; \
				*) \
					echo ""; \
					echo -e "\e[1;39;41m-- WARNING --\e[0m Option $$ans_how not found, exiting process."; \
					echo ""; \
					exit 1; \
			esac; \
		fi;

.PHONY: _pipx_install_poetry
_pipx_install_poetry:
	@output="$$(pip install poetry --dry-run)"; \
	if echo "$$output" | grep -q computecanada ; then \
		echo ""; \
		echo -e "\e[1;39;41m-- WARNING --\e[0m Compute Canada (DRAC) environment detected: Installing Poetry < 2.0.0"; \
		echo ""; \
		pipx install 'poetry<2.0.0' ; \
	else \
		pipx install poetry ; \
	fi;


.PHONY: poetry-install
poetry-install: ## Install standalone Poetry using pipx. Will ask where to install pipx.
	@echo "Looking for Poetry version...";\
	poetry --version; \
	if [ $$? != "0" ]; then \
		if [ "$(AUTO_INSTALL)" = "true" ]; then \
			ans="y";\
		else \
			echo "Poetry not found..."; \
			echo "Looking for pipx version...";\
			pipx_found=0; \
			pipx --version; \
				if [ $$? != "0" ]; then \
					pipx_found=1; \
					echo "pipx not found..."; \
					echo""; \
					echo -n "Would you like to install pipx and Poetry? [y/N]: "; \
				else \
					echo""; \
					echo -n "Would you like to install Poetry using pipx? [y/N]: "; \
				fi; \
			read ans; \
		fi; \
		case $$ans in \
			[Yy]*) \
				if [ $$pipx_found == "1" ]; then \
					echo""; \
					echo -e "\e[1;39;41m-- WARNING --\e[0m The following pip has been found and will be used to install pipx: "; \
					echo "    -> "$$(which pip); \
					echo""; \
					echo "If you do not have write permission to that environment, using it to install pipx will fail."; \
					echo "If this is the case, you should install pipx using a virtual one."; \
					echo""; \
					echo "See documentation for more information."; \
					echo""; \
					echo -n "Would you like to use the local available pip above, or create virtual environment to install pipx? [local/virtual]: "; \
					read ans_how; \
					case $$ans_how in \
						"LOCAL" | "Local" |"local") \
							make -s poetry-install-local; \
							;; \
						"VIRTUAL" | "Virtual" | "virtual") \
							make -s poetry-install-venv; \
							;; \
						*) \
							echo ""; \
							echo -e "\e[1;39;41m-- WARNING --\e[0m Option $$ans_how not found, exiting process."; \
							echo ""; \
							exit 1; \
					esac; \
				else \
					echo "Installing Poetry"; \
					make -s _pipx_install_poetry; \
				fi; \
				;; \
			*) \
				echo "Skipping installation."; \
				echo " "; \
				;; \
		esac; \
	fi;

PIPX_VENV_PATH := $$HOME/.pipx_venv
.PHONY: poetry-install-venv
poetry-install-venv: ## Install standalone Poetry. Will install pipx in $HOME/.pipx_venv
	@pipx --version; \
	if [ $$? != "0" ]; then \
		echo "Creating virtual environment using venv here : [$(PIPX_VENV_PATH)]"; \
		python3 -m venv $(PIPX_VENV_PATH); \
		echo "Activating virtual environment [$(PIPX_VENV_PATH)]"; \
		source $(PIPX_VENV_PATH)/bin/activate; \
		pip3 install pipx; \
		pipx ensurepath; \
		source $(PIPX_VENV_PATH)/bin/activate && make -s _pipx_install_poetry ; \
	else \
		make -s _pipx_install_poetry ; \
	fi;

.PHONY: poetry-install-local
poetry-install-local: ## Install standalone Poetry. Will install pipx with locally available pip.
	@pipx --version; \
	if [ $$? != "0" ]; then \
		echo "pipx not found; installing pipx"; \
		pip3 install pipx; \
		pipx ensurepath; \
	fi;
	@echo "Installing Poetry"
	@make -s _pipx_install_poetry


.PHONY: poetry-env-info
poetry-env-info: ## Information about the currently active environment used by Poetry
	@poetry env info

.PHONY: poetry-create-env
poetry-create-env: ## Create a Poetry managed environment for the project (Outside of Conda environment).
	@echo "Creating Poetry environment that will use Python $(PYTHON_VERSION)"; \
	poetry env use $(PYTHON_VERSION); \
	poetry env info
	@echo""
	@echo "This environment can be accessed either by using the <poetry run YOUR COMMAND>"
	@echo "command, or activated with the <poetry shell> command."
	@echo""
	@echo "Use <poetry --help> and <poetry list> for more information"
	@echo""

.PHONY: poetry-activate
poetry-activate: ## Print the shell command to activate the project's poetry env.
	poetry env activate

.PHONY: poetry-remove-env
poetry-remove-env: ## Remove current project's Poetry managed environment.
	@if [ "$(AUTO_INSTALL)" = "true" ]; then \
		ans_env="y";\
		env_path=$$(poetry env info -p); \
		env_name=$$(basename $$env_path); \
	else \
		echo""; \
		echo "Looking for poetry environments..."; \
		env_path=$$(poetry env info -p); \
		if [[ "$$env_path" != "" ]]; then \
			echo "The following environment has been found for this project: "; \
			env_name=$$(basename $$env_path); \
			echo""; \
			echo "Env name : $$env_name"; \
			echo "PATH     : $$env_path"; \
			echo""; \
			echo "If the active environment listed above is a Conda environment,"; \
			echo "Choosing to delete it will have no effect; use the target <make conda-clean-env>"; \
			echo""; \
			echo""; \
			echo "If the active environment listed above is a venv environment,"; \
			echo "Choosing to delete it will have no effect; use the bash command $ rm -rf <PATH_TO_VENV>"; \
			echo""; \
			echo -n "Would you like delete the environment listed above? [y/N]: "; \
			read ans_env; \
		else \
			env_name="None"; \
			env_path="None"; \
  		fi; \
	fi; \
	if [[ $$env_name != "None" ]]; then \
		case $$ans_env in \
			[Yy]*) \
				poetry env remove $$env_name || echo "No environment was removed"; \
				;; \
			*) \
				echo "No environment was found/provided - skipping environment deletion"; \
				;;\
		esac; \
	else \
		echo "No environments were found... skipping environment deletion"; \
	fi; \

.PHONY: poetry-uninstall
poetry-uninstall: poetry-remove-env ## Uninstall pipx-installed Poetry and the created environment
	@if [ "$(AUTO_INSTALL)" = "true" ]; then \
		ans="y";\
	else \
		echo""; \
		echo -n "Would you like to uninstall pipx-installed Poetry? [y/N]: "; \
		read ans; \
	fi; \
	case $$ans in \
		[Yy]*) \
			pipx --version ; \
			if [ $$? != "0" ]; then \
				echo "" ; \
				echo "Pipx not found globally, trying with $(PIPX_VENV_PATH) env" ;\
				echo "" ; \
				source $(PIPX_VENV_PATH)/bin/activate && pipx uninstall poetry ; \
			else \
				pipx uninstall poetry ; \
				fi; \
			;; \
		*) \
			echo "Skipping uninstallation."; \
			echo " "; \
			;; \
	esac; \

.PHONY: poetry-uninstall-pipx
poetry-uninstall-pipx: poetry-remove-env ## Uninstall pipx-installed Poetry, the created Poetry environment and pipx
	@if [ "$(AUTO_INSTALL)" = "true" ]; then \
		ans="y";\
	else \
		echo""; \
		echo -n "Would you like to uninstall pipx-installed Poetry and pipx? [y/N]: "; \
		read ans; \
	fi; \
	case $$ans in \
		[Yy]*) \
			pipx uninstall poetry; \
			pip uninstall -y pipx; \
			;; \
		*) \
			echo "Skipping uninstallation."; \
			echo " "; \
			;; \
	esac; \

.PHONY: poetry-uninstall-venv
poetry-uninstall-venv: poetry-remove-env ## Uninstall pipx-installed Poetry, the created Poetry environment, pipx and $HOME/.pipx_venv
	@if [ "$(AUTO_INSTALL)" = "true" ]; then \
		ans="y";\
	else \
		echo""; \
		echo -n "Would you like to uninstall pipx-installed Poetry and pipx? [y/N]: "; \
		read ans; \
	fi; \
	case $$ans in \
		[Yy]*) \
			(source $(PIPX_VENV_PATH)/bin/activate && pipx uninstall poetry); \
			(source $(PIPX_VENV_PATH)/bin/activate && pip uninstall -y pipx); \
			;; \
		*) \
			echo "Skipping uninstallation."; \
			echo " "; \
			;; \
	esac; \

	@if [ "$(AUTO_INSTALL)" = "true" ]; then \
		ans="y";\
	else \
		echo""; \
		echo -n "Would you like to remove the virtual environment located here : [$(PIPX_VENV_PATH)] ? [y/N]: "; \
		read ans; \
	fi; \
	case $$ans in \
		[Yy]*) \
			rm -r $(PIPX_VENV_PATH); \
			;; \
		*) \
			echo "Skipping [$(PIPX_VENV_PATH)] virtual environment removal."; \
			echo ""; \
			;; \
	esac; \

## -- Conda targets ------------------------------------------------------------------------------------------------- ##

.PHONY: conda-install
conda-install: ## Install Conda on your local machine
	@echo "Looking for [$(CONDA_TOOL)]..."; \
	$(CONDA_TOOL) --version; \
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
    			bash ~/miniconda.sh -b -p $${HOME}/.conda; \
    			export PATH=$${HOME}/.conda/bin:$$PATH; \
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

.PHONY: conda-create-env
conda-create-env: conda-install ## Create a local Conda environment based on 'environment.yml' file
	@$(CONDA_TOOL) env create $(CONDA_YES_OPTION) -f environment.yml

.PHONY: conda-env-info
conda-env-info: ## Print information about active Conda environment using <CONDA_TOOL>
	@$(CONDA_TOOL) info

.PHONY: conda-activate
conda-activate: ## Print the shell command to activate the project's Conda env.
	@echo "$(CONDA_TOOL) activate $(CONDA_ENVIRONMENT)"

.PHONY: _conda-poetry-install
_conda-poetry-install:
	@$(CONDA_TOOL) run -n $(CONDA_ENVIRONMENT) python --version;  \
	if [ $$? != "0" ]; then \
		echo "Target environment doesn't seem to exist..."; \
		if [ "$(AUTO_INSTALL)" = "true" ]; then \
				ans="y";\
		else \
			echo ""; \
			echo -n "Do you want to create it? [y/N] "; \
			read ans; \
		fi; \
		case $$ans in \
			[Yy]*) \
				echo "Creating conda environment : [$(CONDA_ENVIRONMENT)]"; \
				make -s conda-create-env; \
				;; \
			*) \
				echo "Exiting..."; \
				exit 1;\
				;; \
		esac;\
	fi;
	$(CONDA_TOOL) run -n $(CONDA_ENVIRONMENT) $(CONDA_TOOL) install $(CONDA_YES_OPTION) -c conda-forge poetry; \
	CURRENT_VERSION=$$($(CONDA_TOOL) run -n $(CONDA_ENVIRONMENT) poetry --version | awk '{print $$NF}' | tr -d ')'); \
	REQUIRED_VERSION="1.6.0"; \
	if [ "$$(printf '%s\n' "$$REQUIRED_VERSION" "$$CURRENT_VERSION" | sort -V | head -n1)" != "$$REQUIRED_VERSION" ]; then \
		echo "Poetry installed version $$CURRENT_VERSION is less than minimal version $$REQUIRED_VERSION, fixing urllib3 version to prevent problems"; \
		$(CONDA_TOOL) run -n $(CONDA_ENVIRONMENT)  poetry add "urllib3<2.0.0"; \
	fi;

.PHONY:conda-poetry-install
conda-poetry-install: ## Install Poetry in the project's Conda environment. Will fail if Conda is not found
	@poetry --version; \
    	if [ $$? != "0" ]; then \
			echo "Poetry not found, proceeding to install Poetry..."; \
			echo "Looking for [$(CONDA_TOOL)]...";\
			$(CONDA_TOOL) --version; \
			if [ $$? != "0" ]; then \
				echo "$(CONDA_TOOL) not found; Poetry will not be installed"; \
			else \
				echo "Installing Poetry with Conda in [$(CONDA_ENVIRONMENT)] environment"; \
				make -s _conda-poetry-install; \
			fi; \
		else \
			echo ""; \
			echo "Poetry has been found on this system :"; \
			echo "    Install location: $$(which poetry)"; \
			echo ""; \
			if [ "$(AUTO_INSTALL)" = "true" ]; then \
				ans="y";\
			else \
				echo -n "Would you like to install poetry in the project's conda environment anyway ? [y/N]: "; \
				read ans; \
			fi; \
			case $$ans in \
				[Yy]*) \
					echo "Installing Poetry with Conda in [$(CONDA_ENVIRONMENT)] environment"; \
					make -s _conda-poetry-install; \
					;; \
				*) \
					echo "Skipping installation."; \
					echo " "; \
					;; \
			esac; \
		fi;

.PHONY: conda-poetry-uninstall
conda-poetry-uninstall: ## Uninstall Poetry located in currently active Conda environment
	$(CONDA_TOOL) run -n $(CONDA_ENVIRONMENT) $(CONDA_TOOL) remove $(CONDA_YES_OPTION) poetry

.PHONY: conda-clean-env
conda-clean-env: ## Completely removes local project's Conda environment
	$(CONDA_TOOL) env remove $(CONDA_YES_OPTION) -n $(CONDA_ENVIRONMENT)

## -- Install targets (All install targets will install Poetry if not found using 'make poetry-install-auto')-------- ##

POETRY_COMMAND := poetry

ifeq ($(DEFAULT_INSTALL_ENV),venv)
POETRY_COMMAND := source $(VENV_ACTIVATE) && poetry
else ifeq ($(DEFAULT_INSTALL_ENV),poetry)
POETRY_COMMAND := poetry
else ifeq ($(DEFAULT_INSTALL_ENV),conda)
POETRY_COMMAND := $(CONDA_TOOL) run -n $(CONDA_ENVIRONMENT) poetry
endif

.PHONY: _check-env
_check-env:
	@if ! [ $(DEFAULT_INSTALL_ENV) ]; then \
		echo -e "\e[1;39;41m-- WARNING --\e[0m No installation environment have been defined." ; \
		echo "" ; \
		echo "Defaulting to Poetry managed environment - Poetry will either use activated environment, or '.venv'," ; \
		echo "if found, or create and manage it's own environment if not." ; \
	elif [ $(DEFAULT_INSTALL_ENV) = "venv" ]; then \
		if [ ! -f $(VENV_ACTIVATE) ]; then \
			make -s venv-create ;\
		fi; \
	elif [ $(DEFAULT_INSTALL_ENV) = "conda" ]; then \
		if ! $(CONDA_TOOL) env list | grep -q $(CONDA_ENVIRONMENT) ; then \
			make -s conda-create-env ; \
		fi; \
	fi;

.PHONY: _remind-env-activate
_remind-env-activate:
	@echo ""
	@echo "Activate your environment using the following command:"
	@echo ""
	@if ! [ $(DEFAULT_INSTALL_ENV) ] || [ $(DEFAULT_INSTALL_ENV) = "poetry" ]; then \
		make -s poetry-env-activate ; \
		echo "" ; \
		echo "You can also use the eval bash command : eval \$$(make poetry-activate)"; \
		echo "" ; \
		echo "The environment can also be used through the 'poetry run <command>' command."; \
		echo "" ; \
		echo "    Ex: poetry run python <path_to_script>"; \
	elif [ $(DEFAULT_INSTALL_ENV) = "venv" ]; then \
		make -s venv-activate ; \
		echo "" ; \
		echo "You can also use the eval bash command : eval \$$(make venv-activate)"; \
	elif [ $(DEFAULT_INSTALL_ENV) = "conda" ]; then \
		make -s conda-activate ; \
		echo "" ; \
		echo "You can also use the eval bash command : eval \$$(make conda-activate)"; \
	fi;
	@echo ""

test-echo:
	@echo "use the eval bash command : eval \$$(make poetry-activate)"

.PHONY: install
install: install-precommit ## Install the application package, developer dependencies and pre-commit hook

.PHONY: install-precommit
install-precommit: install-dev ## Install the pre-commit hooks (also installs developer dependencies)
	@if [ -f .git/hooks/pre-commit ]; then \
		echo "Pre-commit hook found"; \
	else \
	  	echo "Pre-commit hook not found, proceeding to configure it"; \
		$(POETRY_COMMAND) run pre-commit install; \
	fi;

.PHONY: install-dev
install-dev: poetry-install-auto _check-env ## Install the application along with developer dependencies
	@$(POETRY_COMMAND) install --with dev
	@make -s _remind-env-activate

.PHONY: install-with-lab
install-with-lab: poetry-install-auto _check-env ## Install the application and it's dev dependencies, including Jupyter Lab
	@$(POETRY_COMMAND) install --with dev --with lab
	@make -s _remind-env-activate


.PHONY: install-package
install-package: poetry-install-auto _check-env ## Install the application package only
	@$(POETRY_COMMAND) install
	@make -s _remind-env-activate

## -- Versioning targets -------------------------------------------------------------------------------------------- ##

# Use the "dry" target for a dry-run version bump, ex.
# make bump-major dry
BUMP_ARGS ?= --verbose
ifeq ($(filter dry, $(MAKECMDGOALS)), dry)
	BUMP_ARGS := $(BUMP_ARGS) --dry-run --allow-dirty
endif

.PHONY: dry
dry: ## Add the dry target for a preview of changes; ex. 'make bump-major dry'
	@-echo > /dev/null

.PHONY: bump-major
bump-major: ## Bump application major version  <X.0.0>
	$(BUMP_TOOL) $(BUMP_ARGS) bump major

.PHONY: bump-minor
bump-minor: ## Bump application minor version  <0.X.0>
	$(BUMP_TOOL) $(BUMP_ARGS) bump minor

.PHONY: bump-patch
bump-patch: ## Bump application patch version  <0.0.X>
	$(BUMP_TOOL) $(BUMP_ARGS) bump patch

## -- Docker targets ------------------------------------------------------------------------------------------------ ##

## -- Apptainer/Singularity targets --------------------------------------------------------------------------------- ##

## -- Linting targets ----------------------------------------------------------------------------------------------- ##

.PHONY: check-lint
check-lint: ## Check code linting (black, isort, flake8, docformatter and pylint)
	poetry run nox -s check

.PHONY: check-pylint
check-pylint: ## Check code with pylint
	poetry run nox -s pylint

.PHONY: check-complexity
check-complexity: ## Check code cyclomatic complexity with Flake8-McCabe
	poetry run nox -s complexity

.PHONY: fix-lint
fix-lint: ## Fix code linting (black, isort, flynt, docformatter)
	poetry run nox -s fix

.PHONY: precommit
precommit: ## Run Pre-commit on all files manually
	poetry run nox -s precommit

## -- Tests targets ------------------------------------------------------------------------------------------------- ##

.PHONY: test
test: ## Run all tests
	poetry run nox -s test

TEST_ARGS ?=
MARKER_TEST_ARGS = -m "$(TEST_ARGS)"
SPECIFIC_TEST_ARGS = -k "$(TEST_ARGS)"
CUSTOM_TEST_ARGS = "$(TEST_ARGS)"

.PHONY: test-marker
test-marker: ## Run tests using pytest markers. Ex. make test-marker TEST_ARGS="<marker>"
	@if [ -n "$(TEST_ARGS)" ]; then \
		poetry run nox -s test_custom -- -- $(MARKER_TEST_ARGS); \
	else \
		echo "" ; \
    	echo 'ERROR : Variable TEST_ARGS has not been set, please rerun the command like so :' ; \
	  	echo "" ; \
    	echo '            make test-marker TEST_ARGS="<marker>"' ; \
	  	echo "" ; \
    fi
.PHONY: test-specific
test-specific: ## Run specific tests using the -k option. Ex. make test-specific TEST_ARGS="<name-of-test>"
	@if [ -n "$(TEST_ARGS)" ]; then \
  		poetry run nox -s test_custom -- -- $(SPECIFIC_TEST_ARGS); \
	else \
		echo "" ; \
    	echo 'ERROR : Variable TEST_ARGS has not been set, please rerun the command like so :' ; \
	  	echo "" ; \
    	echo '            make test-specific TEST_ARGS="<name-of-the-test>"' ; \
	  	echo "" ; \
    fi

.PHONY: test-custom
test-custom: ## Run tests with custom args. Ex. make test-custom TEST_ARGS="-m 'not offline'"
	@if [ -n "$(TEST_ARGS)" ]; then \
  		poetry run nox -s test_custom -- -- $(CUSTOM_TEST_ARGS); \
	else \
	  	echo "" ; \
    	echo 'ERROR : Variable TEST_ARGS has not been set, please rerun the command like so :' ; \
	  	echo "" ; \
    	echo '            make test-custom TEST_ARGS="<custom-args>"' ; \
	  	echo "" ; \
    fi