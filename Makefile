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
	ln -s /home/data/exomiser-data//* $(shell dirname pwd)/$(shell dirname $@)configurations/exomiser-13.2.0/default/2302_phenotype/2302_phenotype.h2.db: $(TMP_DATA)/semsim1.sql
	test -d configurations/exomiser-13.2.0/default/ || mkdir -p configurations/exomiser-13.2.0/default/
	test -d configurations/exomiser-13.2.0/default/2302_phenotype/ || cp -rf /home/data/exomiser-data/2302_phenotype configurations/exomiser-13.2.0/default/2302_phenotype/
	java -cp $(H2_JAR) org.h2.tools.RunScript -user sa -url jdbc:h2:$(shell pwd)/configurations/exomiser-13.2.0/default/2302_phenotype/2302_phenotype  -script $<

semsim-ingest: configurations/exomiser-13.2.0/default/2302_phenotype/2302_phenotype.h2.db


prepare-inputs: configurations/exomiser-13.2.0-semsim1/config.yaml

configurations/exomiser-13.2.0-semsim1/config.yaml:
	mkdir -p $(shell dirname pwd)/$(shell dirname $@)
	ln -s /home/data/exomiser-data//* $(shell dirname pwd)/$(shell dirname $@)configurations/exomiser-13.2.0/semsim1/2302_phenotype/2302_phenotype.h2.db: $(TMP_DATA)/semsim1.sql
	test -d configurations/exomiser-13.2.0/semsim1/ || mkdir -p configurations/exomiser-13.2.0/semsim1/
	test -d configurations/exomiser-13.2.0/semsim1/2302_phenotype/ || cp -rf /home/data/exomiser-data/2302_phenotype configurations/exomiser-13.2.0/semsim1/2302_phenotype/
	java -cp $(H2_JAR) org.h2.tools.RunScript -user sa -url jdbc:h2:$(shell pwd)/configurations/exomiser-13.2.0/semsim1/2302_phenotype/2302_phenotype  -script $<

semsim-ingest: configurations/exomiser-13.2.0/semsim1/2302_phenotype/2302_phenotype.h2.db


.PHONY: prepare-corpora


results/phen2gene-1.2.3-default/corpus1-scrambled0.5/results.yml: configurations/phen2gene-1.2.3-default/config.yaml
	pheval run \
	 --input-dir configurations/phen2gene-1.2.3-default \
	 --testdata-dir $(shell pwd)/corpora/phen2gene/corpus1/scrambled0.5 \
	 --runner phen2genephevalrunner \
	 --tmp-dir data/tmp/ \
	 --output-dir $(shell dirname pwd)/$(shell dirname $@)

	touch $@

corpora/phen2gene/corpus1/scrambled0.5/corpus.yml: $(TEST_DATA)/template_vcf/template_exome_hg19.vcf.gz
	test -d $(shell dirname $@)/phenopackets || mkdir -p $(shell dirname $@)/phenopackets
	pheval-utils scramble-phenopackets \
	 --scramble-factor 0.5 \
	 --output-dir $(shell pwd)/$(shell dirname $@)/phenopackets \
	 --phenopacket-dir=$(shell pwd)/$(TEST_DATA)/phenopackets/single

	touch $@

prepare-corpora: corpora/phen2gene/corpus1/scrambled0.5/corpus.yml

run-corpus1-phen2gene:
	$(MAKE) results/phen2gene-1.2.3-default/corpus1-scrambled0.5/results.yml

pheval-run: run-corpus1-phen2gene

results/exomiser-13.2.0-default/corpus1-scrambled1/results.yml: configurations/exomiser-13.2.0-default/config.yaml
	pheval run \
	 --input-dir configurations/exomiser-13.2.0-default \
	 --testdata-dir $(shell pwd)/corpora/exomiser/corpus1/scrambled1 \
	 --runner exomiserphevalrunner \
	 --tmp-dir data/tmp/ \
	 --output-dir $(shell dirname pwd)/$(shell dirname $@)

	touch $@

corpora/exomiser/corpus1/scrambled1/corpus.yml: $(TEST_DATA)/template_vcf/template_exome_hg19.vcf.gztest -d $(shell dirname $@)/vcf || mkdir -p $(shell dirname $@)/vcf
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

results/exomiser-13.2.0-default/corpus2-scrambled0.3/results.yml: configurations/exomiser-13.2.0-default/config.yaml
	pheval run \
	 --input-dir configurations/exomiser-13.2.0-default \
	 --testdata-dir $(shell pwd)/corpora/exomiser/corpus2/scrambled0.3 \
	 --runner exomiserphevalrunner \
	 --tmp-dir data/tmp/ \
	 --output-dir $(shell dirname pwd)/$(shell dirname $@)

	touch $@

corpora/exomiser/corpus2/scrambled0.3/corpus.yml: $(TEST_DATA)/template_vcf/template_exome_hg19.vcf.gztest -d $(shell dirname $@)/vcf || mkdir -p $(shell dirname $@)/vcf
	test -L $(shell pwd)/corpora/exomiser/corpus2/scrambled0.3/template_exome_hg19.vcf.gz || ln -s $(shell pwd)/$< $(shell dirname $@)/vcf/
	pheval-utils create-spiked-vcfs \
	 --template-vcf-path $(shell pwd)/$(shell dirname $@)/vcf/template_exome_hg19.vcf.gz \
	 --phenopacket-dir=$(shell pwd)/$(TEST_DATA)/phenopackets/single \
	 --output-dir $(shell pwd)/$(shell dirname $@)/vcf

	test -d $(shell dirname $@)/phenopackets || mkdir -p $(shell dirname $@)/phenopackets
	pheval-utils scramble-phenopackets \
	 --scramble-factor 0.3 \
	 --output-dir $(shell pwd)/$(shell dirname $@)/phenopackets \
	 --phenopacket-dir=$(shell pwd)/$(TEST_DATA)/phenopackets/single

	touch $@

prepare-corpora: corpora/exomiser/corpus2/scrambled0.3/corpus.yml

run-corpus2-exomiser:
	$(MAKE) results/exomiser-13.2.0-default/corpus2-scrambled0.3/results.yml

pheval-run: run-corpus2-exomiser


.PHONY: pheval
pheval:
	$(MAKE) prepare-inputs
	$(MAKE) prepare-corpora
	$(MAKE) pheval-run

include ./resources/custom.Makefile