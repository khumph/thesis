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

  dat <- read_rds(data_file)

  n_stages <- length(unique(dat$month))
  n_samples <- length(unique(dat$samp))

  form <- makeFormula(dat)

  tic("Fitting all samples")
  set.seed(20170128)
  Q_list <- lapply(
    1:n_samples,
    function(x) {
      tic(paste("Fitting sample", x, "of", n_samples))
      q <- Qlearn(
        formula = form$formula,
        data = dat[dat$samp == x, ],
        n_stages = n_stages,
        method = 'rpart',
        control = rpart.control(
          cp = 1e-5,
          maxcompete = 0,
          maxsurrogate = 0,
          minbucket = 5)
      )
      toc()
      return(q)
    }) %>% set_names(paste0("rpart", 1:n_samples))
  toc()

  write_rds(Q_list, path = output_file)
}


main(opts$input, opts$output, opts$dependencies)
