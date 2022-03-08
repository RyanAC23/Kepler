# Modelled after
# https://github.com/simoninireland/introduction-to-epidemics/blob/master/Makefile
SHELL = /bin/bash

RESOURCES = None


ACTIVATE ?= eval "$$(conda shell.bash hook)" && conda activate
ANACONDA_PROJECT ?= CONDA_EXE=$(CONDA_EXE) anaconda-project

ENV ?= kepler
ENV_PATH ?= $(abspath envs/$(ENV))
ACTIVATE_PROJECT ?= $(ACTIVATE) $(ENV_PATH)

# ------- Top-level targets  -------

# Default prints a help message
help:
	@make usage


usage:
	@echo "$$HELP_MESSAGE"


init: _ext/Resources anaconda-project.yaml
	@make _init

_init:
	$(ANACONDA_PROJECT) prepare
	$(ANACONDA_PROJECT) run init  # Custom command: see anaconda-project.yaml

.PHONY: _init

_ext/Resources:
	-git clone $(RESOURCES) $@
	@if [ ! -d "$@" ]; then \
	  echo "$$RESOURCES_ERROR_MESSAGE"; \
	fi


Docs/environment.yaml: anaconda-project.yaml Makefile
	$(ANACONDA_PROJECT) run export 1> $@


sync:
	find . -name ".ipynb_checkpoints" -prune -o \
	       -name "_ext" -prune -o \
	       -name "envs" -prune -o \
	       -name "*.ipynb" -o -name "*.md" \
	       -exec jupytext --sync {} + 2> >(grep -v "is not a paired notebook" 1>&2)
# See https://stackoverflow.com/a/15936384/1088938 for details

clean:
	find . -name "__pycache__" -exec $(RM) -r {} +
	$(RM) -r _htmlcov .coverage .pytest_cache
	$(ACTIVATE) root && conda clean --all -y


reallyclean:
	$(ANACONDA_PROJECT) run clean || true  # Custom command: see anaconda-project.yaml
	$(ANACONDA_PROJECT) clean || true
	$(RM) -r envs


test:
	$(ANACONDA_PROJECT) run test

doc-server:
	sphinx-autobuild --ignore Docs/_build Docs Docs/_build/html


.PHONY: clean realclean init cocalc-init sync doc-server help test


# ----- Usage -----

define HELP_MESSAGE

This Makefile provides several tools to help initialize the project.  It is primarly designed
to help get a CoCalc project up an runnning, but should work on other platforms.

Variables:
   ACTIVATE: (= "$(ACTIVATE)")
                     Command to activate a conda environment as `$$(ACTIVATE) <env name>`
                     Defaults to `conda activate`.
   ANACONDA_PROJECT: (= "$(ANACONDA_PROJECT)")
                     Command to run the `anaconda-project` command.  If you need to first
                     activate an environment, then this should do that.
                     Defaults to `anaconda-project`.
   ENV: (= "$(ENV)")
                     Name of the conda environment user by the project.
                     (Customizations have not been tested.)
                     Defaults to `phys-581-2021`.
   ENV_PATH: (= "$(ENV_PATH)")
                     Path to the conda environment user by the project.
                     (Customizations have not been tested.)
                     Defaults to `envs/$$(ENV)`.
   ACTIVATE_PROJECT: (= "$(ACTIVATE_PROJECT)")
                     Command to activate the project environment in the shell.
                     Defaults to `$$(ACTIVATE)  $$(ENV)`.

Initialization:
   make init         Initialize the environment and kernel.  On CoCalc we do specific things
                     like install mmf-setup, and activate the environment in ~/.bash_aliases.
                     This is done by `make init` if ANACONDA2020 is defined.

Testing:
   make test         Runs the general tests.

Maintenance:
   make clean        Call conda clean --all: saves disk space.
   make reallyclean  delete the environments and kernel as well.

Documentation:
   make doc-server   Build the html documentation server on http://localhost:8000
                     Uses Sphinx autobuild
endef
export HELP_MESSAGE


define RESOURCES_ERROR_MESSAGE

*************************************************************
WARNING: The `_ext/Resources` folder could not be cloned from

  $(RESOURCES)

Likely this is because this repository is private and requires registration in the class.
If you believe that you should have access, please contact your instructor, and provide
your GitLab username.

These resources are not crucial for the project, but are important for the course.
*************************************************************

endef
export RESOURCES_ERROR_MESSAGE
