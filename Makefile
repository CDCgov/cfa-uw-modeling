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
	@echo "  help              : Displays this message."
	@echo "  container_build   : Builds the container using the Dockerfile."
	@echo "  container_run     : Runs the container mounting the current folder. Most often used for testing container prior to running sims on Azure Batch."
	@echo "  epimodelcfa       : Builds and installs the epimodelcfa package from project repo, should be run prior to container_build."
	@echo ""

container_build:
	$(ENGINE) build -t $(IMAGE) -f Dockerfile .

container_run:
	$(ENGINE) run --mount type=bind,source=$(PWD),target=/$(IMAGE) \
                -it --rm -w /$(IMAGE) $(IMAGE) bash

epimodelcfa:
	R CMD build /home/RUU7/repos/cfa-epimodel
	R CMD INSTALL --preclean --clean --library=$(lib_path) epimodelcfa_*.tar.gz --dependencies=TRUE
	rm -f epimodelcfa_*tar.gz

.PHONY: container_build container_run
