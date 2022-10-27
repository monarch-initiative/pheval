DATA_DIR=data
DATABASES=HP_HP_MAPPINGS HP_MP_MAPPINGS HP_ZP_MAPPINGS
LIB_DIR=lib

ifndef EXOMISER
override EXOMISER=13.1.0
endif

ifndef HG
override HG=hg19
endif

ifndef PHENOTYPE
override PHENOTYPE=2209
endif

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


build: download unzip dump_db

download: $(EXOMISER_DOWNLOAD) $(HG_DOWNLOAD) $(PHENOTYPE_DOWNLOAD)

unzip: $(EXOMISER_UNZIP) $(HG_UNZIP) $(PHENOTYPE_UNZIP)

dump: dump_db

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

.DEFAULT_GOAL = help

$(EXOMISER_DOWNLOAD):
	wget "https://data.monarchinitiative.org/exomiser/latest/$(EXOMISER_FILE)-distribution.zip" -P $(LIB_DIR)

$(EXOMISER_UNZIP):
	mkdir -p $(LIB_DIR)
	unzip -o $(LIB_DIR)/$(EXOMISER_FILE)-distribution.zip -d $(LIB_DIR)

$(HG_DOWNLOAD):
	mkdir -p $(DATA_DIR)
	wget https://data.monarchinitiative.org/exomiser/latest/$(PHENOTYPE)_$(HG).zip -P $(DATA_DIR)

$(HG_UNZIP):
	unzip -o $(HG_DOWNLOAD) -d $(HG_UNZIP)

$(PHENOTYPE_DOWNLOAD):
	wget https://data.monarchinitiative.org/exomiser/latest/$(PHENOTYPE)_phenotype.zip -P $(DATA_DIR)

$(PHENOTYPE_UNZIP):
	unzip -o $(PHENOTYPE_DOWNLOAD) -d $(PHENOTYPE_UNZIP)

dump_db:
	mkdir -p $(DATA_DIR)/mapping_db
	for i in $(MAPPING_TABLES); do \
		echo "CALL CSVWRITE('./$(DATA_DIR)/mapping_db/$$i.csv',  'SELECT * FROM EXOMISER.$$i',  'charset=UTF-8 fieldSeparator=;')" > dump.sql; \
		java -Xms256m -Xmx2048m -Dh2.bindAddress=127.0.0.1 -cp "./$(LIB_DIR)/h2.jar" org.h2.tools.RunScript -url jdbc:h2:file:./$(DATA_DIR)/$(PHENOTYPE)_phenotype/$(PHENOTYPE)_phenotype/$(PHENOTYPE)_phenotype -script dump.sql -user sa; \
		rm dump.sql;	\
	done

clean:
	rm -rf $(EXOMISER_DOWNLOAD) $(EXOMISER_UNZIP)