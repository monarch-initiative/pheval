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
	@echo "make shuffle-semsim -- generate new ontology terms to the semsim process"
	@echo "make semsim-rust -- Calculates semsim using the new rust method"
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



$(RAW_DATA_FOLDER)/random-hp-terms.txt:
	$(eval TERM=$(shell bash -c "cat $(RAW_DATA_FOLDER)/hp-mp-merged.owl | grep '<!-- http://purl.obolibrary.org/obo/HP_'  | grep '_' |  shuf -n 100  | head -n 100 | rev | cut -d '/' -f 1 | cut -d ' ' -f 2 | rev | tr '_' ':'"))
	echo $(TERM) > $@

$(RAW_DATA_FOLDER)/random-mp-terms.txt:
	$(eval TERM=$(shell bash -c "cat $(RAW_DATA_FOLDER)/hp-mp-merged.owl | grep '<!-- http://purl.obolibrary.org/obo/MP_'  | grep '_' |  shuf -n 100  | head -n 100 | rev | cut -d '/' -f 1 | cut -d ' ' -f 2 | rev | tr '_' ':'"))
	echo $(TERM) > $@



$(RAW_DATA_FOLDER)/hp-mp.semsim.tsv: $(RAW_DATA_FOLDER)/hp-mp-merged.db $(RAW_DATA_FOLDER)/random-hp-terms.txt $(RAW_DATA_FOLDER)/random-mp-terms.txt
	$(eval TERM=$(shell bash -c "cat $(RAW_DATA_FOLDER)/random-hp-terms.txt"))
	$(eval TERM2=$(shell bash -c "cat $(RAW_DATA_FOLDER)/random-mp-terms.txt"))
	runoak -i $< similarity -p i,p $(TERM) @ $(TERM2) -O csv -o $@

$(RAW_DATA_FOLDER)/hp-mp2.semsim.tsv: $(RAW_DATA_FOLDER)/hp-mp-merged2.db $(RAW_DATA_FOLDER)/random-hp-terms.txt $(RAW_DATA_FOLDER)/random-mp-terms.txt
	$(eval TERM=$(shell bash -c "cat $(RAW_DATA_FOLDER)/random-hp-terms.txt"))
	$(eval TERM2=$(shell bash -c "cat $(RAW_DATA_FOLDER)/random-mp-terms.txt"))
	runoak -i $< similarity -p i,p $(TERM) @ $(TERM2) -O csv -o $@


results/pheno_rustsim1.semsim1.tsv:
	test -d results/ || mkdir -p results
	runoak -i rustsim:sqlite:obo:phenio similarity -p i,p hand @ finger -O csv -o $@

.PHONY: semsim-rust
semsim-rust: results/pheno_rustsim1.semsim1.tsv


.PHONY: shuffle-semsim

shuffle-semsim: $(RAW_DATA_FOLDER)/random-hp-terms.txt $(RAW_DATA_FOLDER)/random-mp-terms.txt

semsim: $(RAW_DATA_FOLDER)/hp-mp.semsim.tsv $(RAW_DATA_FOLDER)/hp-mp2.semsim.tsv
