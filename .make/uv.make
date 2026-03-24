# Project and Private variables and targets import to override variables for local
# This is to make sure, sometimes the Makefile includes don't work.

## -- UV targets ------------------------------------------------------------------------------------------------ ##
ENV_TOOL := $(shell command -v uv 2> /dev/null)
LOCAL_UV_PATH := $(shell echo $$HOME/.local/bin/uv)
ifeq ($(ENV_TOOL),)
	ENV_TOOL := $(LOCAL_UV_PATH)
endif
# Do not rename these unless you also rename across all other make files in .make/
ENV_COMMAND_TOOL := $(ENV_TOOL) run
ENV_INSTALL_TOOL := $(ENV_TOOL) sync -p $(PYTHON_VERSION)

.PHONY: uv-install-auto
uv-install-auto:
	@$(ENV_TOOL) --version; \
	if [ $$? != "0" ]; then \
		make -s uv-install-venv; \
	fi;

.PHONY: uv-install
uv-install: ## Install uv interactively.
	@echo "Looking for uv version...";\
	$(ENV_TOOL) --version; \
	if [ $$? != "0" ]; then \
		if [ "$(AUTO_INSTALL)" = "true" ]; then \
			ans="y";\
		else \
			echo "uv not found..."; \
			echo "Looking for pipx version...";\
			pipx_found=0; \
			pipx --version; \
				if [ $$? != "0" ]; then \
					pipx_found=1; \
					echo "pipx not found..."; \
					echo ""; \
					echo -n "Would you like to install pipx and uv? [y/N]: "; \
				else \
					echo ""; \
					echo -n "Would you like to install uv using pipx? [y/N]: "; \
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
							make -s uv-install-local; \
							;; \
						"VIRTUAL" | "Virtual" | "virtual") \
							make -s uv-install-venv; \
							;; \
						*) \
							echo ""; \
							echo -e "$(WARNING)Option $$ans_how not found, exiting process."; \
							echo ""; \
							exit 1; \
					esac; \
				else \
					echo "Installing uv"; \
					make -s _pipx_install_uv; \
				fi; \
				;; \
			*) \
				echo "Skipping installation."; \
				echo " "; \
				;; \
		esac; \
	fi;

.PHONY: _pipx_install_uv
_pipx_install_uv:
	@pipx install uv



.PHONY: uv-install-venv
uv-install-venv: ## Install standalone uv. Will install pipx in $HOME/.pipx_venv
	@pipx --version; \
	if [ $$? != "0" ]; then \
		echo "Creating virtual environment using venv here : [$(PIPX_VENV_PATH)]"; \
		virtualenv $(PIPX_VENV_PATH); \
		echo "Activating virtual environment [$(PIPX_VENV_PATH)]"; \
		source $(PIPX_VENV_PATH)/bin/activate; \
		pip3 install pipx; \
		pipx ensurepath; \
		source $(PIPX_VENV_PATH)/bin/activate && make -s _pipx_install_uv; \
	else \
		make -s _pipx_install_uv; \
	fi;

.PHONY: uv-install-local
uv-install-local: ## Install standalone uv. Will install pipx with locally available pip.
	@pipx --version; \
	if [ $$? != "0" ]; then \
		echo "pipx not found; installing pipx"; \
		pip3 install pipx; \
		pipx ensurepath; \
	fi;
	@echo "Installing UV"
	@make -s _pipx_install_uv

.PHONY: uv-create-env
uv-create-env: ## Create a virtual environment for uv, using the project's python version.
	@$(ENV_TOOL) python install $(PYTHON_VERSION)
	@$(ENV_TOOL) venv --python $(PYTHON_VERSION)

.PHONY: uv-activate
uv-activate: ## Print out the shell command to activate the project's uv environment.
	@make -s venv-activate

.PHONY: uv-remove-env
uv-remove-env: ## Remove current project's uv managed environment.
	@make -s venv-remove

.PHONY: uv-uninstall
uv-uninstall: uv-remove-env ## Uninstall pipx-installed uv and the created environment
	@if [ "$(AUTO_INSTALL)" = "true" ]; then \
		ans="y";\
	else \
		echo ""; \
		echo -n "Would you like to uninstall pipx-installed uv? [y/N]: "; \
		read ans; \
	fi; \
	case $$ans in \
		[Yy]*) \
			pipx --version ; \
			if [ $$? != "0" ]; then \
				echo "" ; \
				echo "Pipx not found globally, trying with $(PIPX_VENV_PATH) env" ;\
				echo "" ; \
				source $(PIPX_VENV_PATH)/bin/activate && pipx uninstall uv ; \
			else \
				pipx uninstall uv ; \
				fi; \
			;; \
		*) \
			echo "Skipping uninstallation."; \
			echo " "; \
			;; \
	esac; \

.PHONY: uv-uninstall-pipx
uv-uninstall-pipx: uv-remove-env ## Uninstall pipx-installed uv, the created uv environment and pipx
	@if [ "$(AUTO_INSTALL)" = "true" ]; then \
		ans="y";\
	else \
		echo ""; \
		echo -n "Would you like to uninstall pipx-installed uv and pipx? [y/N]: "; \
		read ans; \
	fi; \
	case $$ans in \
		[Yy]*) \
			pipx --version ; \
			if [ $$? != "0" ]; then \
				echo "" ; \
				echo "Pipx not found globally, trying with $(PIPX_VENV_PATH) env" ;\
				echo "" ; \
				source $(PIPX_VENV_PATH)/bin/activate && pipx uninstall uv && pip uninstall -y pipx; \
			else \
				pipx uninstall uv ; \
				pip uninstall -y pipx ;\
				fi; \
			;; \
		*) \
			echo "Skipping uninstallation."; \
			echo " "; \
			;; \
	esac; \

.PHONY: uv-uninstall-venv
uv-uninstall-venv: uv-remove-env ## Uninstall pipx-installed uv, the created uv environment, pipx and $HOME/.pipx_venv
	@if [ "$(AUTO_INSTALL)" = "true" ]; then \
		ans="y";\
	else \
		echo ""; \
		echo -n "Would you like to uninstall pipx-installed uv and pipx? [y/N]: "; \
		read ans; \
	fi; \
	case $$ans in \
		[Yy]*) \
			(source $(PIPX_VENV_PATH)/bin/activate && pipx uninstall uv); \
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

## -- Specific install targets (All install targets will install uv if the tool is not found using) ----------------- ##
.PHONY: _remind-env-activate
_remind-env-activate:
	@echo ""
	@echo "Activate your environment using the following command:"
	@echo ""
	@make -s uv-activate
	@echo ""
	@echo "or use the eval bash command : eval \$$(make uv-activate)"
	@echo ""
	@echo "You can also use the following command line to interact with the environment directly"
	@echo ""
	@echo "  $$ uv run <command>"
	@echo ""

.PHONY: install-dev
install-dev: uv-install-auto ## Install the application along with developer dependencies
	@$(ENV_INSTALL_TOOL) --group dev --all-packages
	@make -s _remind-env-activate

.PHONY: install-all
install-all: uv-install-auto ## Install the application and all it's dependency groups
	@$(ENV_INSTALL_TOOL) --all-extras --all-packages
	@make -s _remind-env-activate

.PHONY: install-jupyterlab
install-jupyterlab: uv-install-auto ## Install Jupyter Lab dependencies
	@$(ENV_INSTALL_TOOL) --extra lab --all-packages
	@make -s _remind-env-activate

.PHONY: install-docs
install-docs: uv-install-auto ## Install docs related dependencies (mkdocs)
	@$(ENV_INSTALL_TOOL) --extra docs --all-packages
	@make -s _remind-env-activate

.PHONY: install-package
install-package: uv-install-auto ## Install the application package only
	@$(ENV_INSTALL_TOOL) --no-dev --all-packages
	@make -s _remind-env-activate
