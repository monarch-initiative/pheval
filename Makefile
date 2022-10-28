DATA_DIR=data
DATABASES=HP_HP_MAPPINGS HP_MP_MAPPINGS HP_ZP_MAPPINGS
LIB_DIR=lib

EXOMISER=13.1.0
HG=hg19
PHENOTYPE=2209
SCRAMBLE-FACTOR=0.5

EXOMISER_FILE=exomiser-cli-$(EXOMISER)

#PATH VARIABLES
EXOMISER_DOWNLOAD=$(LIB_DIR)/$(EXOMISER_FILE)-distribution.zip
EXOMISER_UNZIP=$(LIB_DIR)/$(EXOMISER_FILE)

#HG VARIABLES
HG_DOWNLOAD=$(DATA_DIR)/$(PHENOTYPE)_$(HG).zip
HG_UNZIP=$(DATA_DIR)/$(PHENOTYPE)_$(HG)

#PHENOTYPE VARIABLES
PHENOTYPE_DOWNLOAD=$(DATA_DIR)/$(PHENOTYPE)_phenotype.zip
PHENOTYPE_UNZIP=$(DATA_DIR)/$(PHENOTYPE)_phenotype

MAPPING_TABLES=HP_HP_MAPPINGS HP_MP_MAPPINGS HP_ZP_MAPPINGS

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
	@echo "		 PHENOTYPE: sets phenotype version (default: 2209)"
	@echo "																											"
	@echo "		 Example: make build EXOMISER="13.0.1" HG="hg19" PHENOTYPE="2209""
	@echo "		 Simplified example version using default values: make build"

.PHONY: build
build: download unzip dump_db

.PHONY: download
download: $(EXOMISER_DOWNLOAD) $(HG_DOWNLOAD) $(PHENOTYPE_DOWNLOAD)

.PHONY: unzip
unzip: $(EXOMISER_UNZIP) $(HG_UNZIP) $(PHENOTYPE_UNZIP)

.PHONY: dump
dump: $(DATA_DIR)/mapping_db $(DATA_DIR)/mapping_db/$(MAPPING_TABLES).csv


.DEFAULT_GOAL = help

$(EXOMISER_DOWNLOAD):
	wget "https://data.monarchinitiative.org/exomiser/latest/$(EXOMISER_FILE)-distribution.zip" -P $(LIB_DIR)

$(EXOMISER_UNZIP): $(EXOMISER_DOWNLOAD)
	mkdir -p $(LIB_DIR)
	unzip -o $(LIB_DIR)/$(EXOMISER_FILE)-distribution.zip -d $(LIB_DIR)

$(HG_DOWNLOAD):
	mkdir -p $(DATA_DIR)
	wget https://data.monarchinitiative.org/exomiser/latest/$(PHENOTYPE)_$(HG).zip -P $(DATA_DIR)

$(HG_UNZIP): $(HG_DOWNLOAD)
	unzip -o $< -d $@

$(PHENOTYPE_DOWNLOAD):
	wget https://data.monarchinitiative.org/exomiser/latest/$(PHENOTYPE)_phenotype.zip -P $(DATA_DIR)

$(PHENOTYPE_UNZIP): $(PHENOTYPE_DOWNLOAD)
	unzip -o $< -d $@


$(DATA_DIR)/mapping_db:
	mkdir -p $(DATA_DIR)/mapping_db

$(DATA_DIR)/mapping_db/%.csv: $(DATA_DIR)/mapping_db
	echo "CALL CSVWRITE('$@',  'SELECT * FROM EXOMISER.$*',  'charset=UTF-8 fieldSeparator=;')" > dump.sql; \
	java -Xms256m -Xmx2048m -Dh2.bindAddress=127.0.0.1 -cp "./$(LIB_DIR)/h2.jar" org.h2.tools.RunScript -url jdbc:h2:file:./$(DATA_DIR)/$(PHENOTYPE)_phenotype/$(PHENOTYPE)_phenotype/$(PHENOTYPE)_phenotype -script dump.sql -user sa; \
	rm dump.sql;

.PHONY: cleandownload
cleandownload:
	rm -fv $(HG_DOWNLOAD) $(PHENOTYPE_DOWNLOAD)


.PHONY: clean
clean:
	rm -rf $(EXOMISER_DOWNLOAD) $(EXOMISER_UNZIP)