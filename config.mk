SRC_DIR=R
WRITEUP_DIR=docs
RESULTS_DIR=results
DATA_DIR=data

SCENARIOS=simple int noise noise_pred
MODEL_TYPES=rpart mars

DATA=$(addsuffix .rds, $(addprefix $(DATA_DIR)/data-, $(SCENARIOS)))

CART_MODELS=$(patsubst $(DATA_DIR)/data-%.rds, $(RESULTS_DIR)/qrpart-%.RData, $(DATA))
MARS_MODELS=$(patsubst $(DATA_DIR)/data-%.rds, $(RESULTS_DIR)/qmars-%.RData, $(DATA))
MODELS=$(CART_MODELS) $(MARS_MODELS)