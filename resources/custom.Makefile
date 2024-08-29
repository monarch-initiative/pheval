# This includes all the preprocessing not directly related to the PhEval Running

PHENOTYPE_VERSIONS		:=	2309 2402
EXOMISER_VERSIONS		:=	13.3.0 14.0.0

.PHONY: setup
setup: download-phenotype download-exomiser download-phen2gen download-gado

.PHONY: download-phenotype
download-phenotype: $(addprefix $(PHENOTYPE_DIR)/,$(addsuffix _hg19.sha256,$(PHENOTYPE_VERSIONS))) $(addprefix $(PHENOTYPE_DIR)/,$(addsuffix _hg38.sha256,$(PHENOTYPE_VERSIONS))) $(addprefix $(PHENOTYPE_DIR)/,$(addsuffix _phenotype.sha256,$(PHENOTYPE_VERSIONS)))

$(PHENOTYPE_DIR)/%.sha256: $(TMP_DATA)/phenogeno_%.zip
	unzip $< -d $(PHENOTYPE_DIR)

$(TMP_DATA)/phenogeno_%.zip:
	mkdir -p $(PHENOTYPE_DIR)
	mkdir -p $(TMP_DATA)
	wget https://data.monarchinitiative.org/exomiser/data/$*.zip -O $@

.PHONY: download-exomiser
download-exomiser: $(addprefix $(RUNNERS_DIR)/exomiser-,$(EXOMISER_VERSIONS))

$(TMP_DATA)/exomiser-cli-%-distribution.zip:
	mkdir -p $(TMP_DATA)
	wget https://github.com/exomiser/Exomiser/releases/download/$*/exomiser-cli-$*-distribution.zip -O $@

$(RUNNERS_DIR)/exomiser-%: $(TMP_DATA)/exomiser-cli-%-distribution.zip
	mkdir -p $(RUNNERS_DIR)
	unzip $< -d $(RUNNERS_DIR)
	mv $(RUNNERS_DIR)/exomiser-cli-$* $(RUNNERS_DIR)/exomiser-$*
	cp $(RUNNERS_DIR)/configurations/preset-exome-analysis.yml $(RUNNERS_DIR)/exomiser-$*/


.PHONY: download-phen2gen
download-phen2gen: $(RUNNERS_DIR)/Phen2Gene
$(RUNNERS_DIR)/Phen2Gene:
	mkdir -p $(RUNNERS_DIR)
	cd $(RUNNERS_DIR)
	git clone https://github.com/WGLab/Phen2Gene.git $@
	yes "$@" | bash $@/setup.sh


.PHONY: download-gado
download-gado: $(RUNNERS_DIR)/gado
$(RUNNERS_DIR)/gado:
	mkdir -p $@
	wget https://github.com/molgenis/systemsgenetics/releases/download/v1.0.4/GadoCommandline-1.0.1-dist.zip -O $@/gado.1.0.1.zip
	unzip $@/gado.1.0.1.zip -d $@
	wget https://molgenis26.gcc.rug.nl/downloads/genenetwork/v2.1/genenetwork_bonf_spiked.zip -O $@/genenetwork_bonf_spiked.zip
	unzip $@/genenetwork_bonf_spiked.zip -d $@
	wget https://molgenis26.gcc.rug.nl/downloads/genenetwork/v2.1/predictions_auc_bonf.txt -O $@/predictions_auc_bonf.txt
	wget https://molgenis26.gcc.rug.nl/downloads/genenetwork/v2.1/hpo_prediction_genes.txt -O $@/hpo_prediction_genes.txt
	wget https://molgenis26.gcc.rug.nl/downloads/genenetwork/v2.1/hp.obo -O $@/hp.obo


$(TMP_DATA)/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz:
	wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz -O $@
	gunzip $@

$(TMP_DATA)/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz:
	wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -O $@
	gunzip $@

$(TMP_DATA)/all_phenopackets/all_phenopackets.zip:
	mkdir -p $(TMP_DATA)/all_phenopackets/
	wget https://github.com/monarch-initiative/phenopacket-store/releases/download/0.1.12/all_phenopackets.zip -O $@
	unzip $@ -d $(ROOT_DIR)/$(shell dirname $@)/
	mkdir -p $(TMP_DATA)/all_phenopackets/unpacked_phenopackets
	mv $(TMP_DATA)/all_phenopackets/*/* $(TMP_DATA)/all_phenopackets/unpacked_phenopackets

.PHONY: clean
clean:
	rm -rf $(RUNNERS_DIR) $(TMP_DATA)
