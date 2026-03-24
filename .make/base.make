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

.DEFAULT_GOAL := help

# Basic variables
PROJECT_PATH := $(dir $(abspath $(firstword $(MAKEFILE_LIST))))
MAKEFILE_NAME := $(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))
SHELL := /usr/bin/env bash
BUMP_TOOL := bump-my-version
MAKEFILE_VERSION := 1.3.0
DOCKER_COMPOSE ?= docker compose
AUTO_INSTALL ?=

# Conda variables
# CONDA_TOOL can be overridden in Makefile.private file
CONDA_TOOL := conda
CONDA_ENVIRONMENT ?=
CONDA_YES_OPTION ?=

# Global Variables
PIPX_VENV_PATH := $$HOME/.pipx_venv

# Default variables (if Makefile.variables is missing)
APP_VERSION := 0.0.0
APPLICATION_NAME := core
PYTHON_VERSION := 3.12
DEFAULT_INSTALL_ENV := uv
DEFAULT_BUILD_TOOL := uv
TARGET_GROUPS := lint,test
CONDA_ENVIRONMENT := core

# Targets Colors
_ESC := $(shell printf '\033')
_SECTION := $(_ESC)[1m\033[34m
_TARGET  := $(_ESC)[1m\033[36m
_CYAN := $(_ESC)[36m
_NORMAL  := $(_ESC)[0m
_WARNING := $(_ESC)[1;39;41m

WARNING := $(_WARNING) -- WARNING -- $(_NORMAL)

# Project and Private variables and targets import to override variables for local
# This is to make sure, sometimes the Makefile includes don't work.
-include Makefile.variables
-include Makefile.private

comma := ,
contains = $(if $(filter $(1),$(subst $(comma), ,$(2))),true)
not_in = $(if $(filter $(1),$(subst $(comma), ,$(2))),,true)

INSTALL_ENV_IS_VENV := $(call contains,venv,$(DEFAULT_INSTALL_ENV))
INSTALL_ENV_IS_UV := $(call contains,uv,$(DEFAULT_INSTALL_ENV))
INSTALL_ENV_IS_POETRY := $(call contains,poetry,$(DEFAULT_INSTALL_ENV))
INSTALL_ENV_IS_CONDA := $(call contains,conda,$(DEFAULT_INSTALL_ENV))

BUILD_TOOL_IS_UV := $(call contains,uv,$(DEFAULT_BUILD_TOOL))
BUILD_TOOL_IS_POETRY := $(call contains,poetry,$(DEFAULT_BUILD_TOOL))

CONDA_CONFLICT := $(and $(INSTALL_ENV_IS_CONDA),$(BUILD_TOOL_IS_UV))
UV_CONFLICT := $(and $(INSTALL_ENV_IS_POETRY),$(BUILD_TOOL_IS_UV))
POETRY_CONFLICT := $(and $(INSTALL_ENV_IS_UV),$(BUILD_TOOL_IS_POETRY))
PLEASE_FIX_CONFLICT_MSG := Please fix the conflict in your [Makefile.variables] and/or [Makefile.private] file(s)

IS_MAKEFILE_VARIABLES_MISSING := $(call not_in,Makefile.variables,$(MAKEFILE_LIST))
PLEASE_FIX_MISSING_FILE := Please consider adding a [Makefile.variables] file to your project - See lab-advanced-template for more info

check_configs = $(if $($(1)), \
    $(info ) \
    $(info $(WARNING) $(2)) \
    $(info $(PLEASE_FIX_CONFLICT_MSG)) \
    $(info ) \
)

check_files = $(if $($(1)), \
    $(info ) \
    $(info $(WARNING) $(2)) \
    $(info $(PLEASE_FIX_MISSING_FILE)) \
    $(info ) \
)

# Config Checks
# These run immediately when you type 'make'
$(call check_configs,CONDA_CONFLICT,'conda' environment is enabled while using 'uv')
$(call check_configs,UV_CONFLICT,'poetry' environment is enabled while using 'uv')
$(call check_configs,POETRY_CONFLICT,'uv' environment is enabled while using 'poetry')
$(call check_files,IS_MAKEFILE_VARIABLES_MISSING,The configuration file 'Makefile.variables' is missing - Using default values)


## -- Initialization targets ---------------------------------------------------------------------------------------- ##
.PHONY: project-init
project-init: ## Initialize the project from the template - Only run once!
	@python3 $(PROJECT_PATH).make/scripts/auto_init_script.py

.PHONY: project-init-dry-run
project-init-dry: ## Test run: no changes will be made - Initialize the project from the template
	@python3 $(PROJECT_PATH).make/scripts/auto_init_script.py --dry


## -- Informative targets ------------------------------------------------------------------------------------------- ##

.PHONY: info
info: ## Get project configuration info
	@echo ""
	@echo -e "$(_BLUE)--- Configuration Status ---$(_NORMAL)"
	@echo ""
	@echo -e "$(_CYAN)Application Name$(_NORMAL)         : $(APPLICATION_NAME)"
	@echo -e "$(_CYAN)Application version$(_NORMAL)      : $(APP_VERSION)"
	@echo -e "$(_CYAN)Application Root$(_NORMAL)         : [$(PROJECT_PATH)]"
	@echo -e "$(_CYAN)Application package$(_NORMAL)      : [$(PROJECT_PATH)src/$(APPLICATION_NAME)]"
	@echo -e "$(_CYAN)Environment manager$(_NORMAL)      : $(DEFAULT_INSTALL_ENV)"
	@echo -e "$(_CYAN)Build tool$(_NORMAL)               : $(DEFAULT_BUILD_TOOL)"
	@echo -e "$(_CYAN)Python version$(_NORMAL)           : $(PYTHON_VERSION)"
	@echo -e "$(_CYAN)Active makefile targets$(_NORMAL)  : [$(TARGET_GROUPS)]"
	@echo -e "$(_CYAN)Makefile version$(_NORMAL)         : $(MAKEFILE_VERSION)"


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
	@echo -e "$(_CYAN)Application version$(_NORMAL)  : $(APP_VERSION)"
	@echo -e "$(_CYAN)Makefile version$(_NORMAL)     : $(MAKEFILE_VERSION)"

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

.PHONY: venv-remove
venv-remove: ## Delete the virtualenv '.venv' at the root of the project folder.
	@if [ -d $(VENV_PATH) ]; then \
  	  echo "Current venv folder is [$(VENV_PATH)]"; \
  	  if [ "$(AUTO_INSTALL)" = "true" ]; then \
			ans="y";\
	  else \
	    echo ""; \
		echo -n "Would you like to completely delete this virtual environment? [y/N]: "; \
		read ans; \
	  fi; \
	  case $$ans in \
			[Yy]*) \
				echo ""; \
				echo "Starting deletion process for [$(VENV_PATH)]"; \
				rm -rf $(VENV_PATH); \
				echo ""; \
				echo "-- Deletion complete --"; \
				;; \
			*) \
	    		echo ""; \
				echo "Skipping virtual environment deletion."; \
				echo " "; \
				;; \
		esac; \
  	else \
  	  echo "Venv [$(VENV_PATH)] does not exist, nothing to do"; \
  	fi;

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
	@$(ENV_COMMAND_TOOL) $(BUMP_TOOL) bump $(BUMP_ARGS) major

.PHONY: bump-minor
bump-minor: ## Bump application minor version  <0.X.0>
	@$(ENV_COMMAND_TOOL) $(BUMP_TOOL) bump $(BUMP_ARGS) minor

.PHONY: bump-patch
bump-patch: ## Bump application patch version  <0.0.X>
	@$(ENV_COMMAND_TOOL) $(BUMP_TOOL) bump $(BUMP_ARGS) patch




## -- General Install Targets ---------------------------------------------------------------------------------------- ##

.PHONY: install
install: install-dev install-precommit ## Install the application package, developer dependencies and pre-commit hook

.PHONY: install-precommit
install-precommit: ## Install the pre-commit hook (need to run one of the install targets first)
	@if [ -f .git/hooks/pre-commit ]; then \
		echo "Pre-commit hook found"; \
	else \
	  	echo "Pre-commit hook not found, proceeding to configure it"; \
		$(ENV_COMMAND_TOOL) pre-commit install; \
	fi;

.PHONY: uninstall-precommit
uninstall-precommit: ## Uninstall the pre-commit hook
	@$(ENV_COMMAND_TOOL) pre-commit uninstall
