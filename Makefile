include config.mk

## simulate    : Simulate clinical trail data.
.PHONY : simulate
simulate : $(DATA)

$(DUMMIES) :
	touch $(DUMMIES)

$(DATA_DIR)/data-%.rds : R/simulate.R $(DATA_DIR)/dummy-% R/sim.R 
	mkdir -p $(DATA_DIR)
	Rscript $< --dependencies $(lastword $^) --output $@ \
	--scenario $*

## remove      : Remove auto-generated files.
.PHONY : remove
remove :
	rm -fR $(DATA_DIR)
	rm -fR $(RESULTS_DIR)

## help        : Show arguments to make and what they do.
.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<
