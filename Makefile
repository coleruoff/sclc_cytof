IMAGE_NAME := sclc_project
WORKDIR := /workspace

.PHONY: build reproduce shell clean

build:
	docker build -t $(IMAGE_NAME) .

reproduce:
	docker run --rm -v $(PWD):$(WORKDIR) $(IMAGE_NAME) \
		Rscript run_all_scripts.R

shell:
	docker run -it --rm -v $(PWD):$(WORKDIR) $(IMAGE_NAME) bash

clean:
	docker system prune -f
