MAKEFLAGS 				+= --warn-undefined-variables
SHELL 					:= bash
.DEFAULT_GOAL			:= help
URIBASE					:=	http://purl.obolibrary.org/obo
TMP_DATA				:=	data/tmp
ROOT_DIR				:=	$(shell pwd)
PHENOTYPE_DIR			:= $(ROOT_DIR)/data/phenotype
RUNNERS_DIR				:= $(ROOT_DIR)/runners
NAME					:= $(shell python -c 'import tomli; print(tomli.load(open("pyproject.toml", "rb"))["tool"]["poetry"]["name"])')
VERSION					:= $(shell python -c 'import tomli; print(tomli.load(open("pyproject.toml", "rb"))["tool"]["poetry"]["version"])')
SEMSIM_BASE_URL			:= https://storage.googleapis.com/data-public-monarchinitiative/semantic-similarity/latest



help: info
	@echo ""
	@echo "help"
	@echo "make setup -- this runs the download and extraction of genomic, phenotypic and runners data"
	@echo "make pheval -- this runs pheval pipeline corpus preparation and run"
	@echo "make all -- this runs the entire pipeline including setup, corpus preparation and pheval run"
	@echo "make help -- show this help"
	@echo ""

info:
	@echo "Project: $(NAME)"
	@echo "Version: $(VERSION)"

.PHONY: prepare-inputs






configurations/template-1.0.0/config.yaml:
	mkdir -p $(ROOT_DIR)/$(shell dirname $@)/
	cp $(RUNNERS_DIR)/configurations/template-1.0.0.config.yaml $(ROOT_DIR)/$(shell dirname $@)/config.yaml







.PHONY: prepare-corpora


$(TMP_DATA)/semsim/%.tsv:
	mkdir -p $(TMP_DATA)/semsim/
	wget $(SEMSIM_BASE_URL)/$*.tsv -O $@


$(TMP_DATA)/semsim/%.sql:
	mkdir -p $(TMP_DATA)/semsim/
	wget $(SEMSIM_BASE_URL)/$*.sql -O $@


results/run_data.txt:
	touch $@

results/gene_rank_stats.svg: results/run_data.txt
	pheval-utils benchmark-comparison -r $< -o $(ROOT_DIR)/$(shell dirname $@)/results --gene-analysis -y bar_cumulative
	mv $(ROOT_DIR)/gene_rank_stats.svg $@

.PHONY: pheval-report
pheval-report: results/gene_rank_stats.svg


results/template-1.0.0/results.yml: configurations/template-1.0.0/config.yaml corpora/lirical/default/corpus.yml



	rm -rf $(ROOT_DIR)/$(shell dirname $@)
	mkdir -p $(ROOT_DIR)/$(shell dirname $@)




	pheval run \
	 --input-dir $(ROOT_DIR)/configurations/template-1.0.0 \
	 --testdata-dir $(ROOT_DIR)/corpora/lirical/default \
	 --runner templatephevalrunner \
	 --tmp-dir data/tmp/ \
	 --version 1.0.0 \
	 --output-dir $(ROOT_DIR)/$(shell dirname $@)

	touch $@
	echo -e "$(ROOT_DIR)/corpora/lirical/default/phenopackets\t$(ROOT_DIR)/$(shell dirname $@)" >> results/run_data.txt

.PHONY: pheval-run
pheval-run: results/template-1.0.0/results.yml
corpora/lirical/default/corpus.yml:
	test -d $(ROOT_DIR)/corpora/lirical/default/ || mkdir -p $(ROOT_DIR)/corpora/lirical/default/

	test -L $(ROOT_DIR)/corpora/lirical/default/template_exome_hg19.vcf.gz || ln -s $(ROOT_DIR)/testdata/template_vcf/template_exome_hg19.vcf.gz $(ROOT_DIR)/corpora/lirical/default/template_exome_hg19.vcf.gz
	pheval-utils create-spiked-vcfs \
	 --template-vcf-path $(ROOT_DIR)/corpora/lirical/default/template_exome_hg19.vcf.gz  \
	 --phenopacket-dir=$(ROOT_DIR)/corpora/lirical/default/phenopackets \
	 --output-dir $(ROOT_DIR)/corpora/lirical/default/vcf
	touch $@



.PHONY: pheval
pheval:
	$(MAKE) prepare-inputs
	$(MAKE) prepare-corpora
	$(MAKE) pheval-run
	$(MAKE) pheval-report

.PHONY: all
all:
	$(MAKE) setup
	$(MAKE) pheval

include ./resources/custom.Makefile