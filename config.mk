SRC_DIR = R
WRITEUP_DIR = docs
RESULTS_DIR = results
DATA_DIR = data

SCENARIOS = simple int noise noise_pred
MODEL_TYPES = rpart mars rf

DATA = $(addsuffix .rds, $(addprefix $(DATA_DIR)/data-, $(SCENARIOS)))

MODELS = $(foreach mod, $(MODEL_TYPES), \
  $(patsubst $(DATA_DIR)/data-%.rds, $(RESULTS_DIR)/q$(mod)-%.rds, $(DATA)))