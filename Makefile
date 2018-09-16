include config.mk

## all         : Make everything.
.PHONY : all
all : simulate learn


## simulate    : Simulate clinical trail data.
.PHONY : simulate
simulate : $(DATA)

define sim_template
$$(DATA_DIR)/data-$(1).rds : R/simulate.R R/sim-functs.R 
	mkdir -p $$(DATA_DIR)
	Rscript $$< --dependencies $$(lastword $$^) --output $$@ --scenario $(1)
endef

$(foreach scenario, $(SCENARIOS), $(eval $(call sim_template,$(scenario))))


## learn       : Apply Q-learning to simulated data.
.PHONY : learn
learn : $(MODELS)

define learn_template
$$(RESULTS_DIR)/q$(1)-%.rds : R/learn-$(1).R $$(DATA_DIR)/data-%.rds R/q-functs.R
	mkdir -p $$(RESULTS_DIR)
	Rscript $$(wordlist 1, 2, $$^) --dependencies $$(lastword $$^) --output $$@ 
endef

$(foreach mod, $(MODEL_TYPES), $(eval $(call learn_template,$(mod))))


## clean      : Remove auto-generated files.
.PHONY : clean
clean :
	rm -fR $(DATA_DIR)
	rm -fR $(RESULTS_DIR)


## help        : Show arguments to make and what they do.
.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<
