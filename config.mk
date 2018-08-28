SRC_DIR=R
WRITEUP_DIR=docs
RESULTS_DIR=results
DATA_DIR=data

SCENARIOS=simple int noise noise_pred
DUMMIES=$(addprefix $(DATA_DIR)/dummy-, $(SCENARIOS))
DATA=$(addsuffix .rds, $(addprefix $(DATA_DIR)/data-, $(SCENARIOS)))
