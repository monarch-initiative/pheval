# This includes all the preprocessing not directly related to the PhEval Running

.PHONY: prepare-inputs1
prepare-inputs1: configurations/exomiser-$(EXOMISER_VERSION)/default/config.yml
	echo prepare-inputs $*

.PHONY: prepare-corpus1
prepare-corpus1: corpora/corpus1/scrambled1/corpus.yml

.PHONY: run-corpus1
run-corpus1:
	$(MAKE) SCRAMBLE_FACTOR=1 results/exomiser-$(EXOMISER_VERSION)/default-corpus1-scrambled1/results.yml
	$(MAKE) SCRAMBLE_FACTOR=1 results/phen2gene/default-corpus1-scrambled1/results.yml
