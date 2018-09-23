SRC_DIR = R
WRITEUP_DIR = docs
RESULTS_DIR = results
DATA_DIR = data

SCENARIOS = simple int noise noise-pred
MODEL_TYPES = rpart mars rf

DATA = $(addsuffix .rds, $(addprefix $(DATA_DIR)/data-, $(SCENARIOS)))
DATA_BASE = $(addsuffix .rds, $(addprefix $(DATA_DIR)/data-base-, $(SCENARIOS)))
DATA_BEST = $(patsubst $(DATA_DIR)/data-base-%.rds, $(RESULTS_DIR)/data-best-%.rds, $(DATA_BASE))
DATA_CONSTANT = $(patsubst $(DATA_DIR)/data-base-%.rds, $(RESULTS_DIR)/data-constant-%.rds, $(DATA_BASE))

MODELS = $(foreach mod, $(MODEL_TYPES), \
  $(patsubst $(DATA_DIR)/data-%.rds, $(RESULTS_DIR)/q-$(mod)-%.rds, $(DATA)))