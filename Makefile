MAKEFLAGS 				+= --warn-undefined-variables
SHELL 					:= bash
.DEFAULT_GOAL			:= help
URIBASE					:=	http://purl.obolibrary.org/obo
TEST_DATA				:=	testdata
TMP_DATA				:=	data/tmp
NAME					:= $(shell python -c 'import tomli; print(tomli.load(open("pyproject.toml", "rb"))["tool"]["poetry"]["name"])')
VERSION					:= $(shell python -c 'import tomli; print(tomli.load(open("pyproject.toml", "rb"))["tool"]["poetry"]["version"])')
H2_JAR					:= /home/vinicius/.local/share/DBeaverData/drivers/maven/maven-central/com.h2database/h2-1.4.199.jar
EXOMISER_VERSION		:= 13.2.0

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

configurations/semsim1/%.owl:
	test -d configurations/semsim1 || mkdir -p configurations/semsim1
	wget $(URIBASE)/$*.owl -O $@

configurations/semsim1/mp.owl:
	test -d configurations/semsim1 || mkdir -p configurations/semsim1
	wget $(URIBASE)/mp/mp-base.owl -O $@

configurations/semsim1/hp-mp-merged.owl: configurations/semsim1/hp.owl configurations/semsim1/mp.owl
	test -d configurations/semsim1 || mkdir -p configurations/semsim1
	robot merge \
	 --input configurations/semsim1/hp.owl \
	 --input configurations/semsim1/mp.owl reason reduce \
	 --output $@

configurations/semsim1/hp-mp-merged2.owl: configurations/semsim1/hp-mp-merged.owl
	test -d configurations/semsim1 || mkdir -p configurations/semsim1
	$(eval TERM1=$(shell bash -c "cat configurations/semsim1/random-hp-terms.txt | cut -d ' ' -f1 "))
	$(eval TERM2=$(shell bash -c "cat configurations/semsim1/random-mp-terms.txt | cut -d ' ' -f2 "))

	robot merge \
	 --input $< remove \
	 --term $(TERM1) \
	 --term $(TERM2) \
	 --axioms logical \
	 --output $@

configurations/semsim1/%.db: configurations/semsim1/%.owl
	docker run -e ROBOT_JAVA_ARGS='-Xmx15G' -e JAVA_OPTS='-Xmx15G' -v $(shell pwd)/:/work -w /work \
	 --rm -ti obolibrary/odkfull semsql make $@

configurations/semsim1/random-hp-terms.txt: configurations/semsim1/hp-mp-merged.owl
	$(eval TERM=$(shell bash -c "cat configurations/semsim1/hp-mp-merged.owl | grep '<!-- http://purl.obolibrary.org/obo/HP_'  | grep '_' |  shuf -n 100  | head -n 10 | rev | cut -d '/' -f 1 | cut -d ' ' -f 2 | rev | tr '_' ':'"))
	echo $(TERM) > $@

configurations/semsim1/random-mp-terms.txt: configurations/semsim1/hp-mp-merged.owl
	$(eval TERM=$(shell bash -c "cat configurations/semsim1/hp-mp-merged.owl | grep '<!-- http://purl.obolibrary.org/obo/MP_'  | grep '_' |  shuf -n 100  | head -n 10 | rev | cut -d '/' -f 1 | cut -d ' ' -f 2 | rev | tr '_' ':'"))
	echo $(TERM) > $@

configurations/semsim1/hp-mp.semsim.tsv: configurations/semsim1/hp-mp-merged.db configurations/semsim1/random-hp-terms.txt configurations/semsim1/random-mp-terms.txt
	$(eval TERM=$(shell bash -c "cat configurations/semsim1/random-hp-terms.txt"))
	$(eval TERM2=$(shell bash -c "cat configurations/semsim1/random-mp-terms.txt"))
	runoak -i $< similarity -p i,p $(TERM) @ $(TERM2) \
	 --autolabel -O csv -o $@

configurations/semsim1/hp-mp2.semsim.tsv: configurations/semsim1/hp-mp-merged2.db configurations/semsim1/random-hp-terms.txt configurations/semsim1/random-mp-terms.txt
	$(eval TERM=$(shell bash -c "cat configurations/semsim1/random-hp-terms.txt"))
	$(eval TERM2=$(shell bash -c "cat configurations/semsim1/random-mp-terms.txt"))
	runoak -i $< similarity -p i,p $(TERM) @ $(TERM2) \
	 --autolabel -O csv -o $@

configurations/semsim1/hp-mp.semsim.scrambled1.tsv: configurations/semsim1/hp-mp.semsim.tsv
	pheval-utils semsim-scramble \
	 --input $< \
	 --output $@ \
	 --scramble-factor 0.5

configurations/semsim1/hp-mp2.semsim.scrambled1.tsv: configurations/semsim1/hp-mp.semsim.tsv
	pheval-utils semsim-scramble \
	 --input $< \
	 --output $@ \
	 --scramble-factor 0.5

#CONVERT SAMPLE
configurations/semsim1/hp-mp.semsim.scrambled1.sql: configurations/semsim1/hp-mp.semsim.scrambled1.tsv
	pheval-utils semsim-convert \
	 --input $< \
	 --output $@ \
	 --subject-prefix HP \
	 --object-prefix MP

configurations/semsim1/hp-mp2.semsim.scrambled1.sql: configurations/semsim1/hp-mp2.semsim.scrambled1.tsv
	pheval-utils semsim-convert \
	 --input $< \
	 --output $@ \
	 --subject-prefix HP \
	 --object-prefix MP


results/phen2gene/default-corpus1-scrambled1/results.yml: configurations/phen2gene-0.0.0/default/config.yml
	pheval run \
	 --input-dir /home/vinicius/Documents/softwares/Phen2Gene/lib \
	 --testdata-dir $(shell pwd)/corpora/corpus1/scrambled1 \
	 --runner phen2genephevalrunner \
	 --tmp-dir data/tmp/ \
	 --output-dir $(shell dirname pwd)/$(shell dirname $@) \
	 --config $<
	
	touch $@



results/exomiser/default-corpus1-scrambled1/results.yml: configurations/exomiser-13.2.0/default/config.yml
	pheval run \
	 --input-dir ./configurations/exomiser-13.2.0/default \
	 --testdata-dir $(shell pwd)/corpora/corpus1/scrambled1 \
	 --runner exomiserphevalrunner \
	 --tmp-dir data/tmp/ \
	 --output-dir $(shell dirname pwd)/$(shell dirname $@) \
	 --config $<
	
	touch $@


configurations/exomiser-13.2.0/default/2209_phenotype/2209_phenotype.h2.db: configurations/exomiser-13.2.0/default/config.yml configurations/semsim1/hp-mp.semsim.scrambled1.sql configurations/semsim1/hp-mp2.semsim.scrambled1.sql
	test -d /home/data/exomiser-data//hp-mp-semsim.scrambled1 || cp -rf configurations/exomiser-13.2.0/default/2209_phenotype/ /home/data/exomiser-data//hp-mp-semsim.scrambled1
	java -cp $(H2_JAR) org.h2.tools.RunScript -user sa -url jdbc:h2://home/data/exomiser-data//hp-mp-semsim.scrambled1/2209_phenotype -script configurations/semsim1/hp-mp.semsim.scrambled1.sql





corpora/corpus1/scrambled1/corpus.yml: $(TEST_DATA)/template_vcf/template_exome_hg19.vcf.gz
	test -d $(shell dirname $@)/vcf || mkdir -p $(shell dirname $@)/vcf
	test -d $(shell dirname $@)/phenopackets || mkdir -p $(shell dirname $@)/phenopackets
	test -L $(shell pwd)/corpora/corpus1/scrambled1/template_exome_hg19.vcf.gz || ln -s $(shell pwd)/$< $(shell dirname $@)/vcf/
	
	pheval-utils create-spiked-vcfs \
	 --template-vcf-path $(shell dirname $@)/vcf/template_exome_hg19.vcf.gz \
	 --phenopacket-dir=$(TEST_DATA)/phenopackets/single \
	 --output-dir $(shell dirname $@)/vcf
	
	pheval-utils scramble-phenopackets \
	 --scramble-factor 1 \
	 --output-dir $(shell dirname $@)/phenopackets \
	 --phenopacket-dir=$(TEST_DATA)/phenopackets/single
	touch $@

configurations/exomiser-$(EXOMISER_VERSION)/default/config.yml: $(EXOMISER_DATA_FOLDER)/config.yml
	test -L $@ || mkdir -p configurations/exomiser-$(EXOMISER_VERSION)/default/ ; ln -s $(EXOMISER_DATA_FOLDER)* $(shell dirname $@)/

.PHONY: semsim
semsim: configurations/semsim1/hp-mp.semsim.tsv configurations/semsim1/hp-mp2.semsim.tsv

.PHONY: semsim-shuffle
semsim-shuffle: configurations/semsim1/random-hp-terms.txt configurations/semsim1/random-mp-terms.txt

.PHONY: semsim-scramble
semsim-scramble: configurations/semsim1/hp-mp.semsim.scrambled1.tsv configurations/semsim1/hp-mp.semsim.scrambled2.tsv configurations/semsim1/hp-mp.semsim.scrambled3.tsv configurations/semsim1/hp-mp2.semsim.scrambled1.tsv configurations/semsim1/hp-mp2.semsim.scrambled2.tsv configurations/semsim1/hp-mp2.semsim.scrambled3.tsv

.PHONY: semsim-convert
semsim-convert: configurations/semsim1/hp-mp.semsim.scrambled1.sql configurations/semsim1/hp-mp.semsim.scrambled2.sql configurations/semsim1/hp-mp.semsim.scrambled3.sql configurations/semsim1/hp-mp2.semsim.scrambled1.sql configurations/semsim1/hp-mp2.semsim.scrambled2.sql configurations/semsim1/hp-mp2.semsim.scrambled3.sql

.PHONY: semsim-ingest
semsim-ingest: configurations/exomiser-$(EXOMISER_VERSION)/default/2209_phenotype/2209_phenotype.h2.db

.PHONY: prepare-inputs1
prepare-inputs1: configurations/exomiser-$(EXOMISER_VERSION)/default/config.yml
	echo prepare-inputs $*

.PHONY: prepare-corpus1
prepare-corpus1: corpora/corpus1/scrambled1/corpus.yml

.PHONY: run-corpus1
run-corpus1:
	$(MAKE) SCRAMBLE_FACTOR=1 results/exomiser-$(EXOMISER_VERSION)/default-corpus1-scrambled1/results.yml
	# $(MAKE) SCRAMBLE_FACTOR=1 results/phen2gene/default-corpus1-scrambled1/results.yml

.PHONY: pheval
pheval: prepare-inputs1 prepare-corpus1 run-corpus1

.PHONY: clean
clean:
	rm -rf data/* corpora/* inputs/* configurations/* results/*

include custom.Makefile