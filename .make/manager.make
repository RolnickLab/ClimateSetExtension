########################################################################################
#
# MODIFY WITH CARE!!!
# If necessary, override the corresponding variable and/or target, or create new ones
# in one of the following files, depending on the nature of the override :
#
# `Makefile.variables`, `Makefile.targets` or `Makefile.private`,
#
# The only valid reason to modify this file is to fix a bug or to add/remove
# files to include.
#
# REMEMBER!!!
# This is a project level config, any changes here will affect all other users
#
########################################################################################
#
# Necessary make files
#
include .make/base.make
-include Makefile.variables

#
# Optional makefiles targets
#

# Env related
ifneq (,$(INSTALL_ENV_IS_CONDA))
	include .make/conda.make
endif

# Build tool related
ifneq (,$(BUILD_TOOL_IS_UV))
	include .make/uv.make
endif

ifneq (,$(BUILD_TOOL_IS_POETRY))
	include .make/poetry.make
endif

## Linting targets
ifneq (,$(findstring lint,$(TARGET_GROUPS)))
	include .make/lint.make
endif

## Test related targets
ifneq (,$(findstring test,$(TARGET_GROUPS)))
	include .make/test.make
endif

## Docs related targets
ifneq (,$(findstring docs,$(TARGET_GROUPS)))
	include .make/docs.make
endif

#
# Project related makefiles
#
## Custom targets and variables
-include Makefile.targets

## Private variables and targets import to override variables for local
-include Makefile.private

