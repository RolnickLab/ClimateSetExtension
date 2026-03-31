# Project and Private variables and targets import to override variables for local
# This is to make sure, sometimes the Makefile includes don't work.
-include Makefile.variables
-include Makefile.private

POETRY_COMMAND_WITH_PROJECT_ENV := $(shell command -v poetry 2> /dev/null)
LOCAL_POETRY_PATH := $(shell echo $$HOME/.local/bin/poetry)

ifeq ($(POETRY_COMMAND_WITH_PROJECT_ENV),)
	POETRY_COMMAND_WITH_PROJECT_ENV := $(LOCAL_POETRY_PATH)
endif

ifeq ($(DEFAULT_INSTALL_ENV),venv)
POETRY_COMMAND_WITH_PROJECT_ENV := source $(VENV_ACTIVATE) && $(POETRY_COMMAND_WITH_PROJECT_ENV)
else ifeq ($(DEFAULT_INSTALL_ENV),poetry)
POETRY_COMMAND_WITH_PROJECT_ENV := $(POETRY_COMMAND_WITH_PROJECT_ENV)
endif

# Do not rename these unless you also rename across all other make files in .make/
ENV_COMMAND_TOOL := $(POETRY_COMMAND_WITH_PROJECT_ENV) run
ENV_INSTALL_TOOL := $(POETRY_COMMAND_WITH_PROJECT_ENV) install

ifeq ($(DEFAULT_INSTALL_ENV),conda)
ENV_COMMAND_TOOL := $(CONDA_ENV_TOOL) run -n $(CONDA_ENVIRONMENT)
ENV_INSTALL_TOOL := $(ENV_COMMAND_TOOL) $(POETRY_COMMAND_WITH_PROJECT_ENV) install
endif


## -- Poetry targets ------------------------------------------------------------------------------------------------ ##

.PHONY: poetry-install-auto
poetry-install-auto: ## Install Poetry automatically via pipx
	@echo "Looking for Poetry version..."; \
	$(POETRY_COMMAND_WITH_PROJECT_ENV) --version; \
	if [ $$? != "0" ]; then \
		echo "Poetry not found, proceeding to install Poetry..."; \
		make AUTO_INSTALL=true -s poetry-install-venv; \
	else \
		echo "$$output"; \
		version=$$(echo "$$output" | sed -n 's/.*version \([0-9.]*\).*/\1/p'); \
		if [ -n "$$version" ]; then \
			is_lower=$$(echo "$$version" | awk -F. '{ if ($$1 < 2 || ($$1 == 2 && $$2 < 2) || ($$1 == 2 && $$2 == 2 && $$3 < 1)) print "yes" }'); \
			if [ "$$is_lower" = "yes" ]; then \
				echo ""; \
				echo -e "$(WARNING)  Poetry version $$version is lower than 2.2.1. Some features might not work as expected."; \
				echo ""; \
			fi; \
		fi; \
	fi;

.PHONY: _pipx_install_poetry
_pipx_install_poetry:
	@output="$$(pip install --disable-pip-version-check poetry --dry-run)"; \
	if echo "$$output" | grep -q computecanada ; then \
		echo ""; \
		echo -e "$(WARNING)Compute Canada (DRAC) environment detected: Installing Poetry < 2.0.0"; \
		echo "Some features will not be available - like 'poetry python install' which allows poetry"; \
		echo "to manage python versions automatically. Consider loading the appropriate python module"; \
		echo ""; \
		echo "This will also require the 'pyproject.toml' file to use the classic poetry format."; \
		echo ""; \
		echo "Consider loading the appropriate python module before installing this package with 'make install'"; \
		echo "or switching to 'uv'."; \
		echo ""; \
		pipx install 'poetry<2.0.0' ; \
	else \
		pipx install poetry ; \
	fi;


.PHONY: poetry-install
poetry-install: ## Install Poetry interactively.
	@echo "Looking for Poetry version...";\
	$(POETRY_COMMAND_WITH_PROJECT_ENV) --version; \
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
					echo ""; \
					echo -n "Would you like to install pipx and Poetry? [y/N]: "; \
				else \
					echo ""; \
					echo -n "Would you like to install Poetry using pipx? [y/N]: "; \
				fi; \
			read ans; \
		fi; \
		case $$ans in \
			[Yy]*) \
				if [ $$pipx_found == "1" ]; then \
					echo ""; \
					echo -e "$(WARNING)The following pip has been found and will be used to install pipx: "; \
					echo "    -> "$$(which pip); \
					echo ""; \
					echo "If you do not have write permission to that environment, using it to install pipx will fail."; \
					echo "If this is the case, you should install pipx using a virtual one."; \
					echo ""; \
					echo "See documentation for more information."; \
					echo ""; \
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
							echo -e "$(WARNING)Option $$ans_how not found, exiting process."; \
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


.PHONY: poetry-install-venv
poetry-install-venv: ## Install standalone Poetry. Will install pipx in $HOME/.pipx_venv
	@$(POETRY_COMMAND_WITH_PROJECT_ENV) --version; \
	if [ $$? != "0" ]; then \
		echo "Looking for pipx version..."; \
		pipx --version; \
		if [ $$? != "0" ]; then \
			echo "Looking for previously installed pipx environment..."; \
			source $(PIPX_VENV_PATH)/bin/activate && pipx --version; \
				if [ $$? != "0" ]; then \
					echo "Creating virtual environment using venv here : [$(PIPX_VENV_PATH)]"; \
					python3 -m venv $(PIPX_VENV_PATH); \
					echo "Activating virtual environment [$(PIPX_VENV_PATH)]"; \
					source $(PIPX_VENV_PATH)/bin/activate; \
					pip3 install pipx; \
					pipx ensurepath; \
					source $(PIPX_VENV_PATH)/bin/activate && make -s _pipx_install_poetry ; \
				else \
					echo "Pipx found!"; \
					source $(PIPX_VENV_PATH)/bin/activate && make -s _pipx_install_poetry ; \
				fi; \
		else \
			make -s _pipx_install_poetry ; \
		fi;\
	else \
		echo "Poetry is already installed - skipping"; \
	fi;

.PHONY: poetry-install-local
poetry-install-local: ## Install standalone Poetry. Will install pipx with locally available pip.
	@$(POETRY_COMMAND_WITH_PROJECT_ENV) --version; \
	if [ $$? != "0" ]; then \
		echo "Looking for pipx version..."; \
		pipx --version; \
		if [ $$? != "0" ]; then \
			echo "pipx not found; installing pipx"; \
			pip3 install pipx; \
			pipx ensurepath; \
		fi;\
		echo "Installing Poetry"; \
		make -s _pipx_install_poetry; \
	else \
		echo "Poetry is already installed - skipping"; \
	fi;

.PHONY: poetry-env-info
poetry-env-info: ## Information about the currently active environment used by Poetry
	@$(POETRY_COMMAND_WITH_PROJECT_ENV) env info

.PHONY: poetry-env-set-local
poetry-env-set-local: ## Configure poetry to create env locally for this project.
	@$(POETRY_COMMAND_WITH_PROJECT_ENV) config virtualenvs.in-project true --local

.PHONY: poetry-create-env
poetry-create-env: ## Create a Poetry managed environment for the project (Outside of Conda environment).
	@echo "Searching for python version $(PYTHON_VERSION) ..."
	@available_python=$$($(POETRY_COMMAND_WITH_PROJECT_ENV) python list); \
	if ! echo "$$available_python" | grep -qF "$(PYTHON_VERSION)"; then \
		echo "Python version $(PYTHON_VERSION) not found ..."; \
		$(POETRY_COMMAND_WITH_PROJECT_ENV) python install "$(PYTHON_VERSION)"; \
	fi;
	@if [ $(DEFAULT_INSTALL_ENV) != "conda" ]; then \
		echo "Creating Poetry environment using Python $(PYTHON_VERSION)"; \
		$(POETRY_COMMAND_WITH_PROJECT_ENV) env use $(PYTHON_VERSION); \
		$(POETRY_COMMAND_WITH_PROJECT_ENV) env info; \
		echo ""; \
		echo "This environment can be accessed either by using the <poetry run YOUR COMMAND>"; \
		echo "command, or activated with the <poetry env activate> command."; \
		echo ""; \
		echo "Use <poetry --help> and <poetry list> for more information"; \
		echo ""; \
	else \
		echo "Using project's conda environment"; \
	fi;

.PHONY: poetry-activate
poetry-activate: ## Print the shell command to activate the project's poetry env.
	@$(POETRY_COMMAND_WITH_PROJECT_ENV) env activate

.PHONY: poetry-remove-env
poetry-remove-env: ## Remove current project's Poetry managed environment.
	@if [ "$(AUTO_INSTALL)" = "true" ]; then \
		ans_env="y";\
		env_path=$$($(POETRY_COMMAND_WITH_PROJECT_ENV) env info -p); \
		env_name=$$(basename $$env_path); \
	else \
		echo ""; \
		echo "Looking for poetry environments..."; \
		env_path=$$($(POETRY_COMMAND_WITH_PROJECT_ENV) env info -p); \
		if [[ "$$env_path" != "" ]]; then \
			echo "The following environment has been found for this project: "; \
			env_name=$$(basename $$env_path); \
			echo ""; \
			echo "Env name : $$env_name"; \
			echo "PATH     : $$env_path"; \
			echo ""; \
			echo "If the active environment listed above is a Conda environment,"; \
			echo "Choosing to delete it will have no effect; use the target <make conda-clean-env>"; \
			echo ""; \
			echo ""; \
			echo "If the active environment listed above is a venv environment,"; \
			echo "Choosing to delete it will have no effect; use the bash command $ rm -rf <PATH_TO_VENV>"; \
			echo "or 'make venv-remove'"; \
			echo ""; \
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
				$(POETRY_COMMAND_WITH_PROJECT_ENV) env remove $$env_name || echo "No environment was removed"; \
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
		echo ""; \
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
		echo ""; \
		echo -n "Would you like to uninstall pipx-installed Poetry and pipx? [y/N]: "; \
		read ans; \
	fi; \
	case $$ans in \
		[Yy]*) \
			pipx --version ; \
			if [ $$? != "0" ]; then \
				echo "" ; \
				echo "Pipx not found globally, trying with $(PIPX_VENV_PATH) env" ;\
				echo "" ; \
				source $(PIPX_VENV_PATH)/bin/activate && pipx uninstall poetry && pip uninstall -y pipx; \
			else \
				pipx uninstall poetry ; \
				pip uninstall -y pipx ;\
				fi; \
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
		echo ""; \
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
		echo ""; \
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

## -- Specific install targets (All install targets will install Poetry if the tool is not found) --------------------##

.PHONY: _check-env
_check-env:
	@if ! [ $(DEFAULT_INSTALL_ENV) ]; then \
		echo -e "$(WARNING)No installation environment have been defined." ; \
		echo "" ; \
		echo "[DEFAULT_INSTALL_ENV] is not defined - Poetry will use the currently activated environment." ; \
		echo "If there is no currently active environment (ie. conda or venv)," ; \
		echo "Poetry will create and manage it's own environment." ; \
	elif [ $(DEFAULT_INSTALL_ENV) = "venv" ]; then \
		if [ ! -f $(VENV_ACTIVATE) ]; then \
			make -s venv-create ;\
		fi; \
	elif [ $(DEFAULT_INSTALL_ENV) = "conda" ]; then \
		if ! $(CONDA_ENV_TOOL) env list | grep -q $(CONDA_ENVIRONMENT) ; then \
			make -s conda-create-env ; \
		fi; \
	fi;

.PHONY: _remind-env-activate
_remind-env-activate:
	@echo ""
	@echo "Activate your environment using the following command:"
	@echo ""
	@if ! [ $(DEFAULT_INSTALL_ENV) ] || [ $(DEFAULT_INSTALL_ENV) = "poetry" ]; then \
		make -s poetry-activate ; \
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

.PHONY: install-dev
install-dev: poetry-install-auto _check-env poetry-create-env ## Install the application along with developer dependencies
	@$(ENV_INSTALL_TOOL) --with dev
	@make -s _remind-env-activate

.PHONY: install-all
install-all: poetry-install-auto _check-env poetry-create-env ## Install the application and all it's dependency groups
	@$(ENV_INSTALL_TOOL) --with dev --all-extras
	@make -s _remind-env-activate

.PHONY: install-jupyter
install-jupyter: poetry-install-auto _check-env ## Install the application and it's dev dependencies, including Jupyter Lab
	@$(ENV_INSTALL_TOOL) --with dev --extras lab
	@make -s _remind-env-activate


.PHONY: install-docs
install-docs: poetry-install-auto _check-env poetry-create-env ## Install docs related dependencies
	@$(ENV_INSTALL_TOOL) --with dev --extras docs
	@make -s _remind-env-activate

.PHONY: install-package
install-package: poetry-install-auto _check-env poetry-create-env ## Install the application package only
	@$(ENV_INSTALL_TOOL) --only-root
	@make -s _remind-env-activate
