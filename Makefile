include config.mk

## all         : Make everything.
.PHONY : all
all :  simulate learn base constant test best join writeup


## simulate    : Simulate clinical trail data.
.PHONY : simulate
simulate : $(DATA)

define sim_template
$$(DAT_DIR)/data-$(1).rds : R/simulate.R R/sim-functs.R
	mkdir -p $$(DAT_DIR)
	Rscript $$< --dependencies $$(lastword $$^) --output $$@ --scenario $(1) --n_samples 100
endef

$(foreach scenario, $(SCENARIOS), $(eval $(call sim_template,$(scenario))))


## learn       : Apply Q-learning to simulated data.
.PHONY : learn
learn : $(MODELS)

define learn_template
$$(RES_DIR)/q-$(1)-%.rds : R/learn.R $$(DAT_DIR)/data-%.rds
	mkdir -p $$(RES_DIR)
	Rscript $$(wordlist 1, 2, $$^) --output $$@ --method $(1)
endef

$(foreach mod, $(MODEL_TYPES), $(eval $(call learn_template,$(mod))))


## base        : Simulate baseline patient condtions for testing treatment regimes.
.PHONY : base
base : $(DATA_BASE)

define sim_test_template
$$(DAT_DIR)/data-base-$(1).rds : R/simulate.R R/sim-functs.R
	mkdir -p $$(DAT_DIR)
	Rscript $$< --dependencies $$(lastword $$^) --output $$@ --scenario $(1) --seed 20180927 --baseline-only --n_subjects 2000
endef

$(foreach scenario, $(SCENARIOS), $(eval $(call sim_test_template,$(scenario))))


## best        : Simulate best possible treatment sequences for baseline data.
.PHONY : best
best : $(DATA_BEST)

$(RES_DIR)/data-best-%.rds : R/sim-best.R $(DAT_DIR)/data-base-%.rds R/sim-functs.R
	mkdir -p $(RES_DIR)
	Rscript $(wordlist 1, 2, $^) --dependencies $(lastword $^) --output $@


## constant    : Simulate constant dose sequences for baseline data.
.PHONY : constant
constant : $(DATA_CONSTANT)

$(RES_DIR)/data-constant-%.rds : R/sim-constant.R $(DAT_DIR)/data-base-%.rds R/sim-functs.R
	mkdir -p $(RES_DIR)
	Rscript $(wordlist 1, 2, $^) --dependencies $(lastword $^) --output $@


## test        : Simulate dose sequences from models for baseline data.
.PHONY : test
test : $(DATA_TEST)

define test_template
$$(RES_DIR)/data-$(1)-%.rds : R/sim-q.R $$(DAT_DIR)/data-base-%.rds $$(RES_DIR)/q-$(1)-%.rds R/sim-functs.R
	Rscript $$(wordlist 1, 3, $$^) --dependencies $$(lastword $$^) --output $$@
endef

$(foreach mod, $(MODEL_TYPES), $(eval $(call test_template,$(mod))))


## join        : Join all test data.
.PHONY : join
join : $(RES_DIR)/data-all.rds

$(RES_DIR)/data-all.rds : R/join.R $(DATA_BEST) $(DATA_CONSTANT) $(DATA_TEST) 
	Rscript $^ --output $@


## importances  : Get variable importances for each model.
.PHONY : importances
importances : $(RES_DIR)/data-importance.rds

$(RES_DIR)/data-importance.rds : R/importances.R $(MODELS)
	Rscript $^ --output $@


## writeup     : Generate writeup.
.PHONY : writeup
writeup : $(DOC_DIR)/writeup.pdf join

$(DOC_DIR)/writeup.tex : $(DOC_DIR)/writeup.Rnw $(RES_DIR)/data-all.rds $(RES_DIR)/data-importance.rds 
	Rscript -e "pacman::p_load(knitr); knit(input = '$<', output = '$@')"

$(DOC_DIR)/writeup.pdf : $(DOC_DIR)/writeup.tex
	latexmk -pdf -jobname=$(basename $@) -pdflatex="pdflatex -interaction=nonstopmode" -use-make $^


## clean-cache : Remove knitr cache
.PHONY : clean-cache
clean-cache : 
	rm -fR $(CACHE_DIR)
	rm -fR $(FIG_DIR)
	rm -fR $(DOC_DIR)/writeup.pdf 
	rm -fR $(DOC_DIR)/writeup.tex 


## clean       : Remove auto-generated files.
.PHONY : clean
clean : clean-cache
	rm -fR $(DAT_DIR)
	rm -fR $(RES_DIR)


## help        : Show arguments to make and what they do.
.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<
