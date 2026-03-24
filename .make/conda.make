-include Makefile.variables
-include Makefile.private

## -- Conda targets ------------------------------------------------------------------------------------------------- ##

CONDA_ENVIRONMENT_FILE := environment.yml
SLURM_ENV_VAR_PRESENT := env | grep -q "SLURM"

CONDA_ENV_TOOL := $(shell command -v $(CONDA_TOOL) 2> /dev/null)
LOCAL_CONDA_TOOL_PATH := $(shell echo $$HOME/.local/bin/$(CONDA_TOOL))
ifeq ($(CONDA_ENV_TOOL),)
	CONDA_ENV_TOOL := $(LOCAL_CONDA_TOOL_PATH)
endif

.PHONY: conda-install
conda-install: ## General target to install conda like tool - Uses 'CONDA_TOOL' makefile variable
	@echo "### Checking if [$(CONDA_ENV_TOOL)] is installed ..."; \
	$(CONDA_ENV_TOOL) --version; \
	if [ $$? != "0" ]; then \
		echo " "; \
		echo "Your defined Conda tool [$(CONDA_TOOL)] has not been found."; \
		echo " "; \
		echo "If [$(CONDA_TOOL)] or some other Conda tool installed is supposed to already be installed"; \
		echo "on your system, check your [CONDA_TOOL] variable in the Makefile.private for typos."; \
		echo ""; \
		echo "If you already have some version of conda installed, it might not have"; \
		echo "been properly activated (which can also be on purpose when on a compute cluster)."; \
		echo "Consider reloading your shell before you try again."; \
		echo ""; \
		echo "If in doubt, don't install Conda and manually - create and activate"; \
		echo "your own Python environment some other way."; \
		echo ""; \
		echo "This script provides 2 option:"; \
		echo ""; \
		echo "    * You can install 'conda' and 'mamba' through Miniforge3 (https://github.com/conda-forge/miniforge)"; \
		echo "      which is recommended if you need both 'conda' and 'mamba' - Also strongly recommended if you use"; \
		echo "      Pycharm, as 'micromamba' sometimes is not properly recognized by Pycharm."; \
		echo ""; \
		echo "    * Or you can install 'micromamba' (which will also be available as 'mamba' in your shell)"; \
		echo "      from the Micromamba project (https://mamba.readthedocs.io/en/latest/index.html)."; \
		echo "      This is the ultra-lightweight solution."; \
		echo ""; \
		echo -n "Would you like to install one of the tools mentioned above? [y/N]: "; \
		read ans; \
		case $$ans in \
			[Yy]*) \
			  	echo ""; \
				echo -n "Would you like to install miniforge or micromamba? [miniforge/micromamba/none]: "; \
				read conda_provider; \
				case $$conda_provider in \
					"miniforge" | "mini" | "forge" ) \
						echo ""; \
						make -s _miniforge-install;\
					  	;; \
					"micromamba" | "micro" | "mamba" ) \
						echo ""; \
						make -s mamba-install; \
					  	;; \
					"None" | "none" | "no" | "n" ) \
						echo ""; \
						echo "Exiting process"; \
					  	;; \
					*) \
					echo ""; \
					echo "Input is not conform, process is stopping - please try again"; \
				esac; \
				;; \
			*) \
				echo ""; \
				echo "Terminating installation process."; \
				echo ""; \
				;; \
		esac; \
	else \
		echo "Conda tool [$(CONDA_TOOL)] has been found, skipping installation"; \
	fi;

.PHONY: _is_local_bin_on_path
_is_local_bin_on_path:
	@echo ""
	@echo "### Verifying if [$$HOME/.local/bin] is in PATH"
	@if echo $$PATH | tr ':' '\n' | grep -Fxq "$$HOME/.local/bin"; then \
  		echo "";\
  		echo "[$$HOME/.local/bin] found in PATH variable - skipping";\
  		echo "";\
	else \
  	  echo "";\
	  echo "[$$HOME/.local/bin] NOT found in PATH variable"; \
	  echo "";\
	  echo "Adding 'export PATH="\$$HOME/.local/bin:\$$PATH"' to your .bashrc file"; \
	  echo 'export PATH="$$HOME/.local/bin:$$PATH"' >> $$HOME/.bashrc; \
	  echo "";\
	  echo -e "$(WARNING) Consider reloading your shell after this.";\
	  echo "";\
	fi;

.PHONY: _installer-miniforge
_installer-miniforge:
	@make -s _is_local_bin_on_path
	@wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$$(uname)-$$(uname -m).sh"
	@bash Miniforge3-$$(uname)-$$(uname -m).sh -b -p $$HOME/.miniforge3
	@echo ""
	@echo "Adding [conda] and [mamba] to your '\$$HOME/.local/bin' directory"
	@mkdir -p $$HOME/.local/bin
	@ln -s $$HOME/.miniforge3/condabin/conda $$HOME/.local/bin/conda
	@ln -s $$HOME/.miniforge3/condabin/mamba $$HOME/.local/bin/mamba
	@/usr/bin/rm Miniforge3-$$(uname)-$$(uname -m).sh
	@echo ""
	@echo "Please configure the [CONDA_TOOL] variable in your 'Makefile.private file to"
	@echo "either 'CONDA_TOOL := conda' or 'CONDA_TOOL := mamba', depending on which one "
	@echo "you prefer to use"
	@echo ""
	@echo "Consider reloading your shell after this so you can have access to the tool"
	@echo ""

.PHONY: _installer-mamba
_installer-mamba:
	@make -s _is_local_bin_on_path
	@wget -qO- https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
	@echo "Adding [mamba] to your '\$$HOME/.local/bin' directory"
	@mkdir -p $$HOME/.local/bin
	@mv bin/micromamba $$HOME/.local/bin/mamba
	@rm -rf bin/
	@echo ""
	@echo "Please configure the [CONDA_TOOL] variable in your 'Makefile.private file to"
	@echo "'CONDA_TOOL := mamba' for micromamba to be used by the makefile."
	@echo ""
	@echo "Consider reloading your shell after this so you can have access to the tool"
	@echo ""

.PHONY: _slurm-warming
_slurm-warming:
	@echo ""
	@echo "#"
	@echo "#"
	@echo "#"
	@echo ""
	@echo -e "$(WARNING)SLURM Environment variables have been found!!!"
	@echo ""
	@echo "#"
	@echo "#"
	@echo "#"
	@echo ""
	@echo "This indicates you might be on a compute cluster"
	@echo ""
	@echo "It is NOT advisable to execute this command if you are on a Compute Cluster (ie. Mila/DRAC),"
	@echo "as they either have modules available (Mila), or even prohibit the installation "
	@echo "and use of Conda (DRAC) based environments."
	@echo ""
	@echo "Please do not install Conda or similar tools on one the clusters of the"
	@echo "Digital Research Alliance of Canada"
	@echo ""
	@echo "Only proceed if you know what you are doing!!!"
	@echo ""

.PHONY: _miniforge-install
_miniforge-install:
	@echo "### Verifying for SLURM environment variable ..."
	@if $(SLURM_ENV_VAR_PRESENT) ; then \
  		make -s _slurm-warming; \
		echo -n "Are you sure you want to install Miniforge ? [y/N]: "; \
		read ans_slurm; \
		case $$ans_slurm in \
			[Yy]*) \
				echo ""; \
				echo "Installing Miniforge3 - without initializing it."; \
				make -s _installer-miniforge; \
				;; \
			*) \
				echo ""; \
				echo "Terminating installation process."; \
				echo ""; \
				echo "Please activate the required cluster anaconda module or use another way "; \
				echo "to manage your environment before continuing."; \
				echo ""; \
				;; \
		esac; \
	else \
	  echo ""; \
	  echo "Installing and initializing Miniforge"; \
	  echo ""; \
	  make -s _installer-miniforge; \
	  $$HOME/.local/bin/mamba shell init --shell bash --root-prefix=$$HOME/.mamba; \
	  echo ""; \
	  echo "Installation process is complete."; \
	  echo ""; \
	  echo "Only 'mamba' has been initialized. Having 'conda' initialized forces activation the the base environment"; \
	  echo "and can cause conflicts and/or unpredictable behavior with some environments - especially compute clusters."; \
	  echo ""; \
	  echo "To manually initialized, conda use the following command:"; \
	  echo ""; \
	  echo "$HOME/.local/bin/conda init"; \
	  echo ""; \
	fi; \

.PHONY: _mamba-install
_mamba-install:
	@echo "### Verifying for SLURM environment variable ..."
	@if $(SLURM_ENV_VAR_PRESENT) ; then \
		make -s _slurm-warming; \
		echo -n "Are you sure you want to install Micromamba ? [y/N]: "; \
		read ans_slurm; \
		case $$ans_slurm in \
			[Yy]*) \
				echo ""; \
				echo "Installing Micromamba - without initializing it."; \
				make -s _installer-mamba; \
				;; \
			*) \
				echo ""; \
				echo "Terminating installation process."; \
				echo ""; \
				echo "Please activate the required cluster anaconda module or use another way "; \
				echo "to manage your environment before continuing."; \
				echo " "; \
				;; \
		esac; \
	else \
		echo ""; \
		echo "Installing and initializing Micromamba "; \
		echo ""; \
		make -s _installer-mamba; \
		$$HOME/.local/bin/mamba shell init -s bash $$HOME/.micromamba; \
	fi; \

.PHONY: miniforge-install
miniforge-install: ## Install conda and mamba from Miniforge3. (Full functionality for local development)
	@echo "#"; \
	echo "# Miniforge Install process"; \
	echo "#"; \
	echo ""; \
	echo "### Verifying that [conda] is not already installed ..."; \
	conda --version; \
	if [ $$? != "0" ]; then \
		echo ""; \
		echo ""; \
		echo "[conda] has not been found."; \
		echo " "; \
		echo "If [conda] is already supposed to be installed on this system, it might not have"; \
		echo "been properly configured, and/or has not been initialized (which can also be on purpose when "; \
		echo "on a compute cluster). Consider reloading your shell before you try again."; \
		echo ""; \
		echo "If in doubt, don't install Miniforge. Instead, manually create and activate"; \
		echo "your Python environment some other way."; \
		echo ""; \
		echo "### Verifying that [mamba] is not already installed ..."; \
		mamba --version; \
		if [ $$? = "0" ]; then \
			echo -e "$(WARNING)[mamba] has been found - Installing miniforge is probably redundant"; \
		fi; \
		echo ""; \
		echo "### Verifying that [micromamba] is not already installed ..."; \
		micromamba --version; \
		if [ $$? = "0" ]; then \
			echo -e "$(WARNING)[micromamba] has been found - Installing miniforge is probably redundant"; \
		fi; \
		echo ""; \
		echo -n "Would you like to install and initialize Miniforge ? [y/N]: "; \
		read ans; \
		case $$ans in \
			[Yy]*) \
				echo ""; \
			  	make -s _miniforge-install; \
				;; \
			*) \
				echo ""; \
				echo "Terminating installation process."; \
				echo ""; \
				;; \
		esac; \
	else \
		echo "Conda tool [conda] has been found, skipping installation"; \
	fi;

.PHONY: mamba-install
mamba-install: ## Install Micromamba as 'mamba'. (Minimalistic install for env management)
	@echo "#"; \
	echo "# Micromamba Install process"; \
	echo "#"; \
	echo ""; \
	echo "### Verifying that [mamba] is not installed ..."; \
	mamba --version; \
	if [ $$? != "0" ]; then \
		echo ""; \
		echo "### Verifying that [micromamba] is not already installed ..."; \
		micromamba --version; \
		if [ $$? != "0" ]; then \
			echo ""; \
			echo "[mamba] and [micromamba] have not been found."; \
			echo ""; \
			echo "If [mamba] and/or [micromamba] are already supposed to be installed on this system,"; \
			echo "they might not have been properly configured, and/or have not been initialized"; \
			echo "(which can also be on purpose when on a compute cluster). Consider reloading your"; \
			echo "shell before you try again."; \
			echo ""; \
			echo "If in doubt, don't install Micromamba. Instead, manually create and activate"; \
			echo "your Python environment some other way."; \
			echo ""; \
			echo "### Verifying that [conda] is not already installed ..."; \
			conda --version; \
			if [ $$? = "0" ]; then \
				echo -e "$(WARNING)[conda] has been found - Only install micromamba if you really need it"; \
			fi; \
			echo ""; \
			echo -n "Would you like to install and initialize [micromamba] ? [y/N]: "; \
			read ans; \
			case $$ans in \
				[Yy]*) \
					echo ""; \
					make -s _mamba-install; \
					;; \
				*) \
					echo ""; \
					echo "Terminating installation process."; \
					echo ""; \
					;; \
			esac; \
		else \
			echo "[micromamba] has been found, skipping installation"; \
		fi; \
	else \
		echo "[mamba] has been found, skipping installation"; \
	fi; \


.PHONY: conda-create-env
conda-create-env: conda-install ## Create a local Conda environment based on 'environment.yml' file
	@if [ ! -f $(CONDA_ENVIRONMENT_FILE) ]; then \
		$(CONDA_ENV_TOOL) create $(CONDA_YES_OPTION) python=$(PYTHON_VERSION) -c conda-forge -n $(CONDA_ENVIRONMENT); \
		echo "Generating '$(ENV_FILE)' file..."; \
		if [ -f $(CONDA_ENVIRONMENT_FILE) ]; then \
			echo "Warning: $(CONDA_ENVIRONMENT_FILE) already exists. Overwriting..."; \
		fi; \
		( \
			echo "name: $(CONDA_ENVIRONMENT)"; \
			echo "channels:"; \
			echo "  - conda-forge"; \
			echo "dependencies:"; \
			echo -n "  - python=$(PYTHON_VERSION)"; \
		) > $(CONDA_ENVIRONMENT_FILE); \
		echo ""; \
		echo "#"; \
		echo "Done. File content:"; \
		cat $(CONDA_ENVIRONMENT_FILE); \
		echo ""; \
		echo "#"; \
		echo ""; \
	else \
		$(CONDA_ENV_TOOL) env create $(CONDA_YES_OPTION) -f $(CONDA_ENVIRONMENT_FILE); \
	fi;

.PHONY: conda-env-info
conda-env-info: ## Print information about active Conda environment using <CONDA_TOOL>
	@$(CONDA_ENV_TOOL) info

.PHONY: conda-activate
conda-activate: ## Print the shell command to activate the project's Conda env.
	@echo "$(CONDA_ENV_TOOL) activate $(CONDA_ENVIRONMENT)"

.PHONY: conda-clean-env
conda-clean-env: ## Completely removes local project's Conda environment
	$(CONDA_ENV_TOOL) env remove $(CONDA_YES_OPTION) -n $(CONDA_ENVIRONMENT)
