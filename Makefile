MAKEFLAGS 				+= --warn-undefined-variables
SHELL 					:= bash
.DEFAULT_GOAL			:= help
URIBASE					:=	http://purl.obolibrary.org/obo
TEST_DATA				:=	testdata
TMP_DATA				:=	data/tmp
NAME					:= $(shell python -c 'import tomli; print(tomli.load(open("pyproject.toml", "rb"))["tool"]["poetry"]["name"])')
VERSION					:= $(shell python -c 'import tomli; print(tomli.load(open("pyproject.toml", "rb"))["tool"]["poetry"]["version"])')
H2_JAR					:= /home/vinicius/.local/share/DBeaverData/drivers/maven/maven-central/com.h2database/h2-1.4.199.jar
help: status
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

status:
	@echo "Project: $(NAME)"
	@echo "Version: $(VERSION)"

.SECONDEXPANSION:

.PHONY: pheval
pheval: prepare-inputs $$(prepare-corpus) $$(pheval-run)




configurations/exomiser-13.2.0-default/config.yaml:
	mkdir -p $(shell dirname pwd)/$(shell dirname $@)
	ln -s /home/data/exomiser-data//* $(shell dirname pwd)/$(shell dirname $@)



configurations/exomiser-13.2.0/default/2209_phenotype/2209_phenotype.h2.db: configurations/exomiser-13.2.0-default/config.yaml configurations/default/hp-mp.semsim.default.sql configurations/default/hp-mp2.semsim.default.sql
	test -d default/hp-mp-semsim.default || cp -rf configurations/exomiser-13.2.0/default/2209_phenotype/hp-mp-semsim.default
	java -cp $(H2_JAR) org.h2.tools.RunScript -user sa -url jdbc:h2:/configurations/exomiser-13.2.0/default/hp-mp-semsim.default/2209_phenotype -script configurations/default/hp-mp.semsim.default.sql

.PHONY: semsim-ingest
semsim-ingest: configurations/exomiser-13.2.0/default/2209_phenotype/2209_phenotype.h2.db



.PHONY: prepare-inputs
prepare-inputs: configurations/exomiser-13.2.0-default/config.yaml
	echo prepare-inputs $*




results/exomiser-13.2.0/default-corpus1-scrambled1/results.yml: configurations/exomiser-13.2.0-default/config.yaml
	
	pheval run \
	 --input-dir configurations/exomiser-13.2.0-default \
	 --testdata-dir $(shell pwd)/corpora/corpus1/scrambled1 \
	 --runner exomiserphevalrunner \
	 --tmp-dir data/tmp/ \
	 --output-dir $(shell dirname pwd)/$(shell dirname $@) \
	 --config $<
	
	touch $@


corpora/corpus1/scrambled1/corpus.yml: $(TEST_DATA)/template_vcf/template_exome_hg19.vcf.gz
	test -d $(shell dirname $@)/vcf || mkdir -p $(shell dirname $@)/vcf
	test -d $(shell dirname $@)/phenopackets || mkdir -p $(shell dirname $@)/phenopackets
	test -L $(shell pwd)/corpora/corpus1/scrambled1/template_exome_hg19.vcf.gz || ln -s $(shell pwd)/$< $(shell dirname $@)/vcf/
	

.PHONY: prepare-corpus
prepare-corpus +=  corpora/corpus1/scrambled1/corpus.yml



run-corpus1-exomiser:
	$(MAKE) SCRAMBLE_FACTOR=1 results/exomiser-13.2.0/default-corpus1-scrambled1/results.yml

pheval-run += run-corpus1-exomiser




.PHONY: clean
clean:
	rm -rf data/* corpora/* inputs/* configurations/* results/*

include ./resources/custom.Makefile