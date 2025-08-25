ifndef ENGINE
ENGINE := podman
endif

ifndef IMAGE
IMAGE := cfa-uw-modeling
endif

lib_path := /home/RUU7/r_libs

help:
	@echo "Makefile for the $(IMAGE) project."
	@echo ""
	@echo "Usage: make [target]"
	@echo ""
	@echo "Available targets:"
	@echo "  help             : Displays this message."
	@echo "  container_build  : Builds the container using the Dockerfile."
	@echo "  container_run    : Runs the container mounting the current folder. Primarily used for testing container prior to running sims on Azure Batch."
	@echo "  epi_build   			: Builds epimodelcfa package and adds to project repo, should be run prior to container_build if package recently updated."
	@echo "  epi_install 			: Installs the epimodelcfa package into the R library, primarily used for development or local project activities."
	@echo "  epi_clean   			: Cleans up the built package tarball from the project folder."
	@echo "  epi_bic     			: Builds, installs, and cleans the epimodelcfa package."
	@echo ""

epi_build:
	R CMD build /home/RUU7/repos/cfa-epimodel

epi_install:
	R CMD INSTALL --preclean --clean --library=$(lib_path) epimodelcfa_*.tar.gz --dependencies=TRUE

epi_clean:
	rm -f epimodelcfa_*tar.gz

epi_bic: epi_build epi_install epi_clean

container_build: epi_build
	$(ENGINE) build -t $(IMAGE) -f Dockerfile .

container_run:
	$(ENGINE) run --mount type=bind,source=$(PWD),target=/$(IMAGE) \
                -it --rm -w /$(IMAGE) $(IMAGE) bash

.PHONY: container_build container_run epi_install epi_clean help
