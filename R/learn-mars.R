"Perform Q-learning using MARS to fit Q-functions

Usage:
  learn-mars.R <input> (--dependencies <files> --output <file>)
  learn-mars.R -h | --help

Arguments:
  -h --help               Show this screen
  <input>                 .rds data file to use for learning
  --dependencies <files>  Files this program depends on
  --output <file>         Path to .rds file to save fitted models in
" -> doc

pacman::p_load(tidyverse, caret, tictoc, earth)
opts <- docopt::docopt(doc)

main <- function(data_file, output_file, dependencies) {

  lapply(dependencies, source)

  dat <- readRDS(data_file)

  n_stages <- length(unique(dat$month))
  n_samples <- length(unique(dat$samp))

  form <- makeFormula(dat)

  # predictors in PREDICTORS are allowed to interact with predictors in PARENTS
  # but no other interactions are allowed

  PREDICTORS <- form$predictor_names
  PARENTS <-  "dose"
  allowFunct <- function(degree, pred, parents, namesx) {
    if (degree <= 1)
      return(TRUE)
    predictor <- namesx[pred]
    parents <- namesx[parents != 0]
    if ((any(predictor %in% PREDICTORS) &&
         any(parents %in% PARENTS)) ||
        (any(predictor %in% PARENTS) &&
         any(parents %in% PREDICTORS))) {
      return(TRUE)
    }
    FALSE
  }

  set.seed(20170128)
  q <- Qlearn(
    formula = form$formula,
    data = dat,
    n_stages = n_stages,
    method = 'gcvEarth',
    trControl = trainControl(method = "none"),
    nk = 200,
    allowed = allowFunct,
    tuneGrid = expand.grid(degree = 2)
  )

  saveRDS(q, file = output_file, compress = F)
}


tic()
main(opts$input, opts$output, opts$dependencies)
toc()
