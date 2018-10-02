SRC_DIR = R
DAT_DIR = data
RES_DIR = results
DOC_DIR = docs
CACHE_DIR = cache
FIG_DIR = figure
PDF_DIR = docs-typeset

SCENARIOS = simple noise noise-pred simple-int noise-int noise-pred-int
MODEL_TYPES = rpart mars rf

DATA = $(addsuffix .rds, $(addprefix $(DAT_DIR)/data-, $(SCENARIOS)))
DATA_BASE = $(addsuffix .rds, $(addprefix $(DAT_DIR)/data-base-, $(SCENARIOS)))
DATA_BEST = $(patsubst $(DAT_DIR)/data-base-%.rds, $(RES_DIR)/data-best-%.rds, $(DATA_BASE))
DATA_CONSTANT = $(patsubst $(DAT_DIR)/data-base-%.rds, $(RES_DIR)/data-constant-%.rds, $(DATA_BASE))

MODELS = $(foreach mod, $(MODEL_TYPES), \
  $(patsubst $(DAT_DIR)/data-%.rds, $(RES_DIR)/q-$(mod)-%.rds, $(DATA)))

DATA_TEST = $(foreach mod, $(MODEL_TYPES), \
  $(patsubst $(DAT_DIR)/data-base-%.rds, $(RES_DIR)/data-$(mod)-%.rds, $(DATA_BASE)))

N_SAMPLES = 100 # Number of repeated training samples to take
N_TEST_OBS = 2000 # Number of test observations to be simulated for each regime
