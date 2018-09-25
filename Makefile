include config.mk

## all         : Make everything.
.PHONY : all
all : simulate learn base best


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
$$(RESULTS_DIR)/q-$(1)-%.rds : R/learn-$(1).R $$(DATA_DIR)/data-%.rds R/q-functs.R
	mkdir -p $$(RESULTS_DIR)
	Rscript $$(wordlist 1, 2, $$^) --dependencies $$(lastword $$^) --output $$@
endef

$(foreach mod, $(MODEL_TYPES), $(eval $(call learn_template,$(mod))))


## base        : Simulate baseline patient condtions for testing treatment regimes.
.PHONY : base
base : $(DATA_BASE)

define sim_test_template
$$(DATA_DIR)/data-base-$(1).rds : R/simulate.R R/sim-functs.R
	mkdir -p $$(DATA_DIR)
	Rscript $$< --dependencies $$(lastword $$^) --output $$@ --scenario $(1) --baseline-only
endef

$(foreach scenario, $(SCENARIOS), $(eval $(call sim_test_template,$(scenario))))


## best        : Simulate best possible treatment sequences for baseline data.
.PHONY : best
best : $(DATA_BEST)

$(RESULTS_DIR)/data-best-%.rds : R/sim-best.R $(DATA_DIR)/data-base-%.rds R/sim-functs.R
	mkdir -p $(RESULTS_DIR)
	Rscript $(wordlist 1, 2, $^) --dependencies $(lastword $^) --output $@


## constant    : Simulate constant dose sequences for baseline data.
.PHONY : constant
constant : $(DATA_CONSTANT)

$(RESULTS_DIR)/data-constant-%.rds : R/sim-constant.R $(DATA_DIR)/data-base-%.rds R/sim-functs.R
	mkdir -p $(RESULTS_DIR)
	Rscript $(wordlist 1, 2, $^) --dependencies $(lastword $^) --output $@


## test        : Simulate constant dose sequences for baseline data.
.PHONY : test
test : $(DATA_TEST)

define test_template
$$(RESULTS_DIR)/data-$(1)-%.rds : R/sim-q.R $$(DATA_DIR)/data-base-%.rds $$(RESULTS_DIR)/q-$(1)-%.rds R/sim-functs.R
	Rscript $$(wordlist 1, 3, $$^) --dependencies $$(lastword $$^) --output $$@
endef

$(foreach mod, $(MODEL_TYPES), $(eval $(call test_template,$(mod))))


## clean       : Remove auto-generated files.
.PHONY : clean
clean :
	rm -fR $(DATA_DIR)
	rm -fR $(RESULTS_DIR)


## help        : Show arguments to make and what they do.
.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<
