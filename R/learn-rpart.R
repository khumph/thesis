"Perform Q-learning using CART/rpart to fit Q-functions

Usage:
  learn-rpart.R <input> (--dependencies <files> --output <file>)
  learn-rpart.R -h | --help

Arguments:
  -h --help               Show this screen
  <input>                 .rds data file to use for learning
  --dependencies <files>  Files this program depends on
  --output <file>         Path to .rds file to save fitted models in
" -> doc

pacman::p_load(tidyverse, caret, rpart, tictoc)
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
    method = 'rpart',
    control = rpart.control(
      cp = 1e-5,
      maxcompete = 0,
      maxsurrogate = 0,
      minbucket = 5)
  )
  toc()

  saveRDS(q, file = output_file, compress = F)
}


main(opts$input, opts$output, opts$dependencies)
