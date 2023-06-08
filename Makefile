MAKEFLAGS 				+= --warn-undefined-variables
SHELL 					:= bash
.DEFAULT_GOAL			:= help
URIBASE					:=	http://purl.obolibrary.org/obo
TEST_DATA				:=	testdata
TMP_DATA				:=	data/tmp
NAME					:= $(shell python -c 'import tomli; print(tomli.load(open("pyproject.toml", "rb"))["tool"]["poetry"]["name"])')
VERSION					:= $(shell python -c 'import tomli; print(tomli.load(open("pyproject.toml", "rb"))["tool"]["poetry"]["version"])')
H2_JAR					:= /home/vinicius/.local/share/DBeaverData/drivers/maven/maven-central/com.h2database/h2-1.4.199.jar
help: info
	@echo ""
	@echo "help"
	@echo "make pheval -- this runs the entire pipeline including corpus preparation and pheval run"
	@echo "make semsim -- generate all configured similarity profiles"
	@echo "make semsim-shuffle -- generate new ontology terms to the semsim process"
	@echo "make semsim-scramble -- scramble semsim profile"
	@echo "make semsim-convert -- convert all semsim profiles into exomiser SQL format"
	@echo "make semsim-ingest -- takes all the configured semsim profiles and loads them into the exomiser databases"

	@echo "make clean -- removes corpora and pheval results"
	@echo "make help -- show this help"
	@echo ""

info:
	@echo "Project: $(NAME)"
	@echo "Version: $(VERSION)"

.PHONY: prepare-inputs
prepare-inputs: configurations/phen2gene-1.2.3-default/config.yaml

configurations/phen2gene-1.2.3-default/config.yaml:
	mkdir -p $(shell dirname pwd)/$(shell dirname $@)
	ln -s /home/vinicius/Documents/softwares/Phen2Gene/* $(shell dirname pwd)/$(shell dirname $@)



prepare-inputs: configurations/exomiser-13.2.0-default/config.yaml

configurations/exomiser-13.2.0-default/config.yaml:
	mkdir -p $(shell dirname pwd)/$(shell dirname $@)
	ln -s /home/data/exomiser-data//* $(shell dirname pwd)/$(shell dirname $@)

configurations/exomiser-13.2.0/default/2302_phenotype/2302_phenotype.h2.db: $(TMP_DATA)/semsim1.sql
	test -d configurations/exomiser-13.2.0/default/ || mkdir -p configurations/exomiser-13.2.0/default/
	test -d configurations/exomiser-13.2.0/default/2302_phenotype/ || cp -rf /home/data/exomiser-data/2302_phenotype configurations/exomiser-13.2.0/default/2302_phenotype/
	java -cp $(H2_JAR) org.h2.tools.RunScript -user sa -url jdbc:h2:$(shell pwd)/configurations/exomiser-13.2.0/default/2302_phenotype/2302_phenotype  -script $<

semsim-ingest: configurations/exomiser-13.2.0/default/2302_phenotype/2302_phenotype.h2.db


prepare-inputs: configurations/exomiser-13.2.0-semsim1/config.yaml

configurations/exomiser-13.2.0-semsim1/config.yaml:
	mkdir -p $(shell dirname pwd)/$(shell dirname $@)
	ln -s /home/data/exomiser-data//* $(shell dirname pwd)/$(shell dirname $@)

configurations/exomiser-13.2.0/semsim1/2302_phenotype/2302_phenotype.h2.db: $(TMP_DATA)/semsim1.sql
	test -d configurations/exomiser-13.2.0/semsim1/ || mkdir -p configurations/exomiser-13.2.0/semsim1/
	test -d configurations/exomiser-13.2.0/semsim1/2302_phenotype/ || cp -rf /home/data/exomiser-data/2302_phenotype configurations/exomiser-13.2.0/semsim1/2302_phenotype/
	java -cp $(H2_JAR) org.h2.tools.RunScript -user sa -url jdbc:h2:$(shell pwd)/configurations/exomiser-13.2.0/semsim1/2302_phenotype/2302_phenotype  -script $<

semsim-ingest: configurations/exomiser-13.2.0/semsim1/2302_phenotype/2302_phenotype.h2.db


.PHONY: prepare-corpora


results/exomiser-13.2.0-default/corpus1-scrambled1/results.yml: configurations/exomiser-13.2.0-default/config.yaml
	mkdir -p $(shell pwd)/$(shell dirname $@)
	pheval run \
	 --input-dir $(shell pwd)/configurations/exomiser-13.2.0-default \
	 --testdata-dir $(shell pwd)/corpora/exomiser/corpus1/scrambled1 \
	 --runner exomiserphevalrunner \
	 --tmp-dir data/tmp/ \
	 --version 13.2.0 \
	 --output-dir $(shell pwd)/$(shell dirname $@)

	touch $@

corpora/exomiser/corpus1/scrambled1/corpus.yml: $(TEST_DATA)/template_vcf/template_exome_hg19.vcf.gz
	test -d $(shell dirname $@)/vcf || mkdir -p $(shell dirname $@)/vcf
	test -L $(shell pwd)/corpora/exomiser/corpus1/scrambled1/template_exome_hg19.vcf.gz || ln -s $(shell pwd)/$< $(shell dirname $@)/vcf/
	pheval-utils create-spiked-vcfs \
	 --template-vcf-path $(shell pwd)/$(shell dirname $@)/vcf/template_exome_hg19.vcf.gz \
	 --phenopacket-dir=$(shell pwd)/$(TEST_DATA)/phenopackets/single \
	 --output-dir $(shell pwd)/$(shell dirname $@)/vcf

	test -d $(shell dirname $@)/phenopackets || mkdir -p $(shell dirname $@)/phenopackets
	pheval-utils scramble-phenopackets \
	 --scramble-factor 1 \
	 --output-dir $(shell pwd)/$(shell dirname $@)/phenopackets \
	 --phenopacket-dir=$(shell pwd)/$(TEST_DATA)/phenopackets/single

	touch $@

prepare-corpora: corpora/exomiser/corpus1/scrambled1/corpus.yml

run-corpus1-exomiser:
	$(MAKE) results/exomiser-13.2.0-default/corpus1-scrambled1/results.yml

pheval-run: run-corpus1-exomiser


.PHONY: pheval
pheval:
	$(MAKE) prepare-inputs
	$(MAKE) prepare-corpora
	$(MAKE) pheval-run

include ./resources/custom.Makefile