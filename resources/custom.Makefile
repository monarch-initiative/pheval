# This includes all the preprocessing not directly related to the PhEval Running

configurations/exomiser-13.2.0-semsim1/%.owl:
	test -d configurations/semsim1 || mkdir -p configurations/semsim1
	wget $(URIBASE)/$*.owl -O $@


configurations/semsim1/hp.owl:
	test -d configurations/semsim1 || mkdir -p configurations/semsim1
	wget $(URIBASE)/hp/hp-base.owl -O $@

configurations/semsim1/mp.owl:
	test -d configurations/semsim1 || mkdir -p configurations/semsim1
	wget $(URIBASE)/mp/mp-base.owl -O $@

configurations/semsim1/hp-mp-merged.owl: configurations/semsim1/hp.owl configurations/semsim1/mp.owl
	test -d configurations/semsim1 || mkdir -p configurations/semsim1
	robot merge \
	 --input $< \
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
$(TMP_DATA)/semsim1.sql: configurations/semsim1/hp-mp.semsim.scrambled1.tsv
	test -d $(TMP_DATA) || mkdir -p $(TMP_DATA)
	pheval-utils semsim-convert \
	 --input $< \
	 --output $@ \
	 --subject-prefix HP \
	 --object-prefix MP

$(TMP_DATA)/semsim2.sql: configurations/semsim1/hp-mp2.semsim.scrambled1.tsv
	test -d $(TMP_DATA) || mkdir -p $(TMP_DATA)
	pheval-utils semsim-convert \
	 --input $< \
	 --output $@ \
	 --subject-prefix HP \
	 --object-prefix MP

.PHONY: semsim
semsim: configurations/semsim1/hp-mp.semsim.tsv configurations/semsim1/hp-mp2.semsim.tsv

.PHONY: semsim-shuffle
semsim-shuffle: configurations/semsim1/random-hp-terms.txt configurations/semsim1/random-mp-terms.txt

.PHONY: semsim-scramble
semsim-scramble: configurations/semsim1/hp-mp.semsim.scrambled1.tsv configurations/semsim1/hp-mp.semsim.scrambled2.tsv configurations/semsim1/hp-mp.semsim.scrambled3.tsv configurations/semsim1/hp-mp2.semsim.scrambled1.tsv configurations/semsim1/hp-mp2.semsim.scrambled2.tsv configurations/semsim1/hp-mp2.semsim.scrambled3.tsv

.PHONY: semsim-convert
semsim-convert: configurations/semsim1/hp-mp.semsim.scrambled1.sql configurations/semsim1/hp-mp.semsim.scrambled2.sql configurations/semsim1/hp-mp.semsim.scrambled3.sql configurations/semsim1/hp-mp2.semsim.scrambled1.sql configurations/semsim1/hp-mp2.semsim.scrambled2.sql configurations/semsim1/hp-mp2.semsim.scrambled3.sql