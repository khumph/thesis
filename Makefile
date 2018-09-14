include config.mk

## all         : Make everything.
.PHONY : all
all : simulate learn


## simulate    : Simulate clinical trail data.
.PHONY : simulate
simulate : $(DATA)

$(DUMMIES) :
	touch $(DUMMIES)

$(DATA_DIR)/data-%.rds : R/simulate.R $(DATA_DIR)/dummy-% R/sim-functs.R 
	mkdir -p $(DATA_DIR)
	Rscript $< --dependencies $(lastword $^) --output $@ --scenario $*


## learn       : Apply Q-learning to simulated data.
.PHONY : learn
learn : $(CART_MODELS)

$(RESULTS_DIR)/QCART-%.RData : R/learn-CART.R $(DATA_DIR)/data-%.rds R/q-functs.R 
	mkdir -p $(RESULTS_DIR)
	Rscript $(wordlist 1, 2, $^) --dependencies $(lastword $^) --output $@ 


## remove      : Remove auto-generated files.
.PHONY : remove
remove :
	rm -fR $(DATA_DIR)
	rm -fR $(RESULTS_DIR)


## help        : Show arguments to make and what they do.
.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<
