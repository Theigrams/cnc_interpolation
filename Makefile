.PHONY: all example test add_path clean help

example: ## Run the example naive agent
	python curve.py

test: ## Run the unit tests
	export PYTHONPATH=$(shell pwd):$$PYTHONPATH && pytest
	
clean: ## Remove all __pycache__ directories
	find . -type d -name __pycache__ -print -exec rm -r {} +

add_path: ## Add the current directory to the PYTHONPATH
	export PYTHONPATH=$(shell pwd):$$PYTHONPATH
	
help: ## Display this help message
	@awk 'BEGIN {FS = ":.*##"; printf "\nUsage:\n  make \033[36m<target>\033[0m\n"} /^[a-zA-Z_-]+:.*?##/ { printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2 } /^##@/ { printf "\n\033[1m%s\033[0m\n", substr($$0, 5) } ' $(MAKEFILE_LIST)