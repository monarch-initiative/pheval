MAKEFLAGS 				+= --warn-undefined-variables
SHELL 					:= bash
.DEFAULT_GOAL			:= help
URIBASE					:=	http://purl.obolibrary.org/obo
RAW_DATA_FOLDER			:=	data/raw
TEST_DATA				:=	testdata
TMP_DATA				:=	data/tmp
SCRAMBLE_FACTOR			:=	0.5
NAME					:= $(shell python -c 'import tomli; print(tomli.load(open("pyproject.toml", "rb"))["tool"]["poetry"]["name"])')
VERSION					:= $(shell python -c 'import tomli; print(tomli.load(open("pyproject.toml", "rb"))["tool"]["poetry"]["version"])')


help: status
	@echo ""
	@echo "make semsim -- setup data required for semsim"
	@echo "make shuffle-semsim -- generate new ontology terms to the semsim process"
	@echo "make pheval -- pheval.exomiser run"
	@echo "make clean -- removes corpora and pheval exomiser results"
	@echo "make help -- show this help"
	@echo ""

status:
	@echo "Project: $(NAME)"
	@echo "Version: $(VERSION)"

$(RAW_DATA_FOLDER)/%.owl:
	test -d $(RAW_DATA_FOLDER) || mkdir -p $(RAW_DATA_FOLDER)
	wget $(URIBASE)/$*.owl -O $@

$(RAW_DATA_FOLDER)/mp.owl:
	test -d $(RAW_DATA_FOLDER) || mkdir -p $(RAW_DATA_FOLDER)
	wget $(URIBASE)/mp/mp-base.owl -O $@

$(RAW_DATA_FOLDER)/hp-mp-merged.owl: $(RAW_DATA_FOLDER)/hp.owl $(RAW_DATA_FOLDER)/mp.owl
	test -d $(RAW_DATA_FOLDER) || mkdir -p $(RAW_DATA_FOLDER)
	robot merge --input $(RAW_DATA_FOLDER)/hp.owl --input $(RAW_DATA_FOLDER)/mp.owl reason reduce --output $@

$(RAW_DATA_FOLDER)/hp-mp-merged2.owl: $(RAW_DATA_FOLDER)/hp-mp-merged.owl
	test -d $(RAW_DATA_FOLDER) || mkdir -p $(RAW_DATA_FOLDER)
	$(eval TERM1=$(shell bash -c "cat $(RAW_DATA_FOLDER)/random-hp-terms.txt | cut -d ' ' -f1 "))
	$(eval TERM2=$(shell bash -c "cat $(RAW_DATA_FOLDER)/random-mp-terms.txt | cut -d ' ' -f2 "))

	robot merge --input $< \
	remove --term $(TERM1) --term $(TERM2) --axioms logical --output $@

$(RAW_DATA_FOLDER)/%.db: $(RAW_DATA_FOLDER)/%.owl
	semsql make $@



$(RAW_DATA_FOLDER)/random-hp-terms.txt: $(RAW_DATA_FOLDER)/hp-mp-merged.owl
	$(eval TERM=$(shell bash -c "cat $(RAW_DATA_FOLDER)/hp-mp-merged.owl | grep '<!-- http://purl.obolibrary.org/obo/HP_'  | grep '_' |  shuf -n 100  | head -n 100 | rev | cut -d '/' -f 1 | cut -d ' ' -f 2 | rev | tr '_' ':'"))
	echo $(TERM) > $@

$(RAW_DATA_FOLDER)/random-mp-terms.txt: $(RAW_DATA_FOLDER)/hp-mp-merged.owl
	#TODO: GET BACK TO HP
	$(eval TERM=$(shell bash -c "cat $(RAW_DATA_FOLDER)/hp-mp-merged.owl | grep '<!-- http://purl.obolibrary.org/obo/HP_'  | grep '_' |  shuf -n 100  | head -n 100 | rev | cut -d '/' -f 1 | cut -d ' ' -f 2 | rev | tr '_' ':'"))
	echo $(TERM) > $@



$(RAW_DATA_FOLDER)/hp-mp.semsim.tsv: $(RAW_DATA_FOLDER)/hp-mp-merged.db $(RAW_DATA_FOLDER)/random-hp-terms.txt $(RAW_DATA_FOLDER)/random-mp-terms.txt
	$(eval TERM=$(shell bash -c "cat $(RAW_DATA_FOLDER)/random-hp-terms.txt"))
	$(eval TERM2=$(shell bash -c "cat $(RAW_DATA_FOLDER)/random-mp-terms.txt"))
	runoak -i $< similarity -p i,p $(TERM) @ $(TERM2) -O csv -o $@

$(RAW_DATA_FOLDER)/hp-mp2.semsim.tsv: $(RAW_DATA_FOLDER)/hp-mp-merged2.db $(RAW_DATA_FOLDER)/random-hp-terms.txt $(RAW_DATA_FOLDER)/random-mp-terms.txt
	$(eval TERM=$(shell bash -c "cat $(RAW_DATA_FOLDER)/random-hp-terms.txt"))
	$(eval TERM2=$(shell bash -c "cat $(RAW_DATA_FOLDER)/random-mp-terms.txt"))
	runoak -i $< similarity -p i,p $(TERM) @ $(TERM2) -O csv -o $@


.PHONY: shuffle-semsim

shuffle-semsim: $(RAW_DATA_FOLDER)/random-hp-terms.txt $(RAW_DATA_FOLDER)/random-mp-terms.txt

semsim: $(RAW_DATA_FOLDER)/hp-mp.semsim.tsv $(RAW_DATA_FOLDER)/hp-mp2.semsim.tsv


.PHONY: pheval

# pheval: results/exomiser13/corpus-1-r%/config1/results.json
# pheval: prepare-corpus-corpus-1
pheval: prepare-corpus-1 run-exomiser-corpus-1

inputs/exomiser13/config1:
	test -d inputs/exomiser13/ || mkdir -p inputs/exomiser13/
	ln -s /home/data/exomiser-data $@
	#This actually not correct, it is just a proxy. Add one: altered hpo database


results/exomiser_13_1_0_config1/corpus-1-r%_results: inputs/exomiser13/config1
	pheval run --input-dir $(shell pwd)/$< --testdata-dir $(shell pwd)/corpora/corpus-1-r$* --runner exomiserphevalrunner --tmp-dir $(shell pwd)/$(TMP_DATA)/ --output-dir $(shell pwd)/results --config $(shell pwd)/$</config.yml


corpora/corpus-1-r%/corpus.yml: #corpora/corpus-1/corpus.yml
	test -d corpora/corpus-1-r$*/phenopackets || mkdir -p corpora/corpus-1-r$*/phenopackets
	test -d corpora/corpus-1-r$*/vcf || mkdir -p corpora/corpus-1-r$*/vcf
	test -f $(shell dirname $@)/vcf/template_exome_hg19.vcf.gz || ln -s $(shell pwd)/corpora/vcf/* $(shell dirname $@)/vcf

	pheval-utils scramble-phenopackets --scramble-factor $* --output-dir corpora/corpus-1-r$*/phenopackets --phenopacket-dir=$(TEST_DATA)/phenopackets/single --output-file-suffix=test
	test -d corpora/corpus-1-r$*/vcf || mkdir -p corpora/corpus-1-r$*/vcf
	touch $@


corpora/vcf/template_exome_hg19.vcf.gz:
	mkdir -p corpora/vcf
	test -f $@ || ln -s $(shell pwd)/$(TEST_DATA)/template_vcf/template_exome_hg19.vcf.gz $@
	pheval-utils create-spiked-vcfs --template-vcf-path $@ --phenopacket-dir=$(TEST_DATA)/phenopackets/single --output-dir $(shell dirname $@)

prepare-corpus-%: corpora/vcf/template_exome_hg19.vcf.gz
	$(MAKE) corpora/corpus-$*-r0.9/corpus.yml
	$(MAKE) corpora/corpus-$*-r0.7/corpus.yml
	$(MAKE) corpora/corpus-$*-r0.5/corpus.yml

run-exomiser-corpus-%: corpora/vcf/template_exome_hg19.vcf.gz
	$(MAKE) results/exomiser_13_1_0_config1/corpus-1-r0.9_results
	$(MAKE) results/exomiser_13_1_0_config1/corpus-1-r0.7_results
	$(MAKE) results/exomiser_13_1_0_config1/corpus-1-r0.5_results

.PHONY: clean
clean:
	rm -rf corpora/* inputs/* results/*


#illustrates how exomiser should work
#randomiser for the phenopacket should work
