RAW_DIR=raw
DATA_DIR=$(RAW_DIR)/data

EXOMISER=13.1.0
HG=hg19
PHENOTYPE=2209
SCRAMBLE-FACTOR=0.5

EXOMISER_FILE=exomiser-cli-$(EXOMISER)

#PATH VARIABLES
EXOMISER_DOWNLOAD=$(RAW_DIR)/$(EXOMISER_FILE)-distribution.zip
EXOMISER_UNZIP=$(RAW_DIR)/$(EXOMISER_FILE)
#HG VARIABLES
HG_DOWNLOAD=$(DATA_DIR)/$(PHENOTYPE)_$(HG).zip
HG_UNZIP=$(DATA_DIR)/$(PHENOTYPE)_$(HG)

#PHENOTYPE VARIABLES
PHENOTYPE_DOWNLOAD=$(DATA_DIR)/$(PHENOTYPE)_phenotype.zip
PHENOTYPE_UNZIP=$(DATA_DIR)/$(PHENOTYPE)_phenotype

#MAPPING VARIABLES
MAPPING_TABLES=HP_HP_MAPPINGS HP_MP_MAPPINGS HP_ZP_MAPPINGS
#SCRAMBLE VARIABLES 
SCRAMBLE_DIR=inputs/data/run_exomiser$(EXOMISER)_scrambled_$(SCRAMBLE-FACTOR)
EXOMISER_DIR=inputs/data/run_exomiser$(EXOMISER)

.PHONY: help
help :
	@echo "build		:Run download unzip and build targets"
	@echo "download	:Download exomiser, human genome and phenotypic data"
	@echo "unzip		:Unzip all download files from download target."
	@echo "dump		:Dump tables from phenotypic H2 database."
	@echo "clean		:Remove all downloaded unzipped and dumped data."
	@echo "flags		"
	@echo "		 EXOMISER: sets exomiser version (default: 13.1.0)"
	@echo "		 HG: sets human genome version (default: hg19)"
	@echo "		 PHENOTYPE: sets phenotype version (default: 2209)"
	@echo "		 SCRAMBLE-FACTOR: sets scramble factor (default: 0.5)"
	@echo "																											"
	@echo "		 Example: make build EXOMISER="13.0.1" HG="hg19" PHENOTYPE="2209""
	@echo "		 Simplified example version using default values: make build"


.PHONY: download
download: $(EXOMISER_DOWNLOAD) $(HG_DOWNLOAD) $(PHENOTYPE_DOWNLOAD)

.PHONY: unzip
unzip: $(EXOMISER_UNZIP) $(HG_UNZIP) $(PHENOTYPE_UNZIP)

.DEFAULT_GOAL=help

$(EXOMISER_DOWNLOAD):
	wget "https://data.monarchinitiative.org/exomiser/latest/$(EXOMISER_FILE)-distribution.zip" -P $(RAW_DIR)

$(EXOMISER_UNZIP): $(EXOMISER_DOWNLOAD)
	mkdir -p $(RAW_DIR)
	unzip -o $(RAW_DIR)/$(EXOMISER_FILE)-distribution.zip -d $(RAW_DIR)

$(HG_DOWNLOAD):
	mkdir -p $(DATA_DIR)
	wget https://data.monarchinitiative.org/exomiser/latest/$(PHENOTYPE)_$(HG).zip -P $(DATA_DIR)

$(HG_UNZIP): $(HG_DOWNLOAD)
	unzip -o $< -d $@

$(PHENOTYPE_DOWNLOAD):
	wget https://data.monarchinitiative.org/exomiser/latest/$(PHENOTYPE)_phenotype.zip -P $(DATA_DIR)

$(PHENOTYPE_UNZIP): $(PHENOTYPE_DOWNLOAD)
	unzip -o $< -d $@

.PHONY: dump
dump:$(EXOMISER_DIR)/*.tsv

$(EXOMISER_DIR)/%.tsv:
	for prefix in $(MAPPING_TABLES); do \
		echo "CALL CSVWRITE('$(EXOMISER_DIR)/$${prefix}.tsv',  'SELECT * FROM EXOMISER.$${prefix}',  'charset=UTF-8 fieldSeparator=\t')" > dump.sql; \
		java -Xms1024m -Xmx8192m -Dh2.bindAddress=127.0.0.1 -cp "./$(RAW_DIR)/lib/h2.jar" org.h2.tools.RunScript -url jdbc:h2:file:./$(DATA_DIR)/$(PHENOTYPE)_phenotype/$(PHENOTYPE)_phenotype/$(PHENOTYPE)_phenotype -script dump.sql -user sa; \
		rm dump.sql ; \
	done


$(SCRAMBLE_DIR)/%.tsv:
	for prefix in $(MAPPING_TABLES); do \
		pheval scramble -T $${prefix} -S $(SCRAMBLE-FACTOR) -D $(EXOMISER_DIR)/ ; \
		echo "CALL CSVWRITE('$(SCRAMBLE_DIR)/$$prefix.tsv',  'SELECT * FROM EXOMISER.$${prefix}_SCRAMBLE',  'charset=UTF-8 fieldSeparator=\t')" > dump.sql; \
		java -Xms1024m -Xmx8192m -Dh2.bindAddress=127.0.0.1 -cp "./$(RAW_DIR)/lib/h2.jar" org.h2.tools.RunScript -url jdbc:h2:file:./$(DATA_DIR)/$(PHENOTYPE)_phenotype/$(PHENOTYPE)_phenotype/$(PHENOTYPE)_phenotype -script dump.sql -user sa; \
		rm dump.sql ; \
	done


.PHONY: scramble
scramble:dump $(SCRAMBLE_DIR)/*.tsv

.PHONY: cleandownload
cleandownload:
	rm -fv $(HG_DOWNLOAD) $(PHENOTYPE_DOWNLOAD)

.PHONY: clean
clean:
	rm -rf $(EXOMISER_DOWNLOAD) $(EXOMISER_UNZIP)

.PHONY: build
build: download unzip dump scramble