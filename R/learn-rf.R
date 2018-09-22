"Perform Q-learning using random forest to fit Q-functions

Usage:
  learn-rf.R <input> (--dependencies <files> --output <file>)
  learn-rf.R -h | --help

Arguments:
  -h --help               Show this screen
  <input>                 .rds data file to use for learning
  --dependencies <files>  Files this program depends on
  --output <file>         Path to .rds file to save fitted models in
" -> doc

pacman::p_load(tidyverse, caret, tictoc, ranger, e1071)
opts <- docopt::docopt(doc)

main <- function(data_file, output_file, dependencies) {

  lapply(dependencies, source)

  dat <- readRDS(data_file)

  n_stages <- length(unique(dat$month))
  n_samples <- length(unique(dat$samp))

  form <- makeFormula(dat)

  tic()
  set.seed(20170128)
  q <- Qlearn(
    formula = form$formula,
    data = dat,
    n_stages = n_stages,
    method = 'ranger',
    tuneGrid = expand.grid(mtry = length(form$predictor_names),
                           min.node.size = 5, splitrule = "variance"),
    trControl = trainControl(method = "none"),
    importance = "impurity",
    always.split.variables = "dose",
    num.trees = 250
  )
  toc()

  saveRDS(q, file = output_file, compress = F)
}


main(opts$input, opts$output, opts$dependencies)
