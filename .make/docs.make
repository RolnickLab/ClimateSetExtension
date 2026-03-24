## -- Docs targets -------------------------------------------------------------------------------------------------- ##
.PHONY: preview-docs
preview-docs: install-docs ## Preview the documentation site locally
	@$(ENV_COMMAND_TOOL) mkdocs serve -a 0.0.0.0:7000


.PHONY: build-docs
build-docs: install-docs ## Build the documentation files locally
	@$(ENV_COMMAND_TOOL) mkdocs build

.PHONY: deploy-docs
deploy-docs: install-docs ## Publish and deploy the documentation to the live Github page
	@echo ""; \
	echo -e "\e[1;39;41m-- WARNING --\e[0m This command will deploy all current changes to the live Github page - Making it publicly available"; \
	echo ""; \
	echo -n "Would you like to deploys the docs? [Y/n]: "; \
	read ans; \
	case $$ans in \
		[Yy]*) \
			echo ""; \
			$(ENV_COMMAND_TOOL) mkdocs gh-deploy; \
			echo ""; \
			;; \
		*) \
			echo ""; \
			echo "Skipping publication to Github Pages."; \
			echo " "; \
			;; \
	esac; \
