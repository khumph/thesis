"Perform Q-learning using CART to fit Q-functions

Usage:
  learn-CART.R <input> (--dependencies <files> --output <file>)
  learn-CART.R -h | --help

Arguments:
  -h --help               Show this screen
  <input>                 .rds data file to use for learning
  --dependencies <files>  Files this program depends on
  --output <file>         Path to .RData file to save fitted models in
" -> doc

pacman::p_load(tidyverse, caret, rpart)
opts <- docopt::docopt(doc)

main <- function(data_file, output_file, dependencies) {

  lapply(dependencies, source)

  dat <- read_rds(data_file)

  n_stages = length(unique(dat$month))
  n_samples = length(unique(dat$samp))

  form <-
    paste("Qhat ~ ", paste(
      c("dose", "tumor_mass", "toxicity", grep("\\X|\\Z|\\V", names(dat), value = T)),
      collapse = " + "
    )) %>% as.formula()

  ptm <- proc.time()
  Q_rpart <- lapply(
    1:n_samples,
    function(x) {
      Qlearn(
        formula = form,
        data = dat[dat$samp == x, ],
        n_stages = n_stages,
        method = 'rpart',
        control = rpart.control(
          cp = 1e-5,
          maxcompete = 0,
          maxsurrogate = 0,
          minbucket = 5)
      )
    }) %>% set_names(paste0("CART", 1:n_samples))
  print(proc.time() - ptm)

  save(Q_rpart, file = output_file)
}


main(opts$input, opts$output, opts$dependencies)
