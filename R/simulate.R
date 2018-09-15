"Simulate adaptive clinical trial data

Usage:
  simulate.R (--dependencies <files> --output <file>) [--n_subjects <num> --n_samples <num> --n_stages <num> --seed <num> --scenario <scenario>]
  simulate.R -h | --help

Arguments:
  -h --help               Show this screen
  --dependencies <files>  Files this program depends on
  --output <file>         Path to output .rds file
  --n_subjects <num>      Total number of participants to randomize [default: 1000]
  --n_samples <num>       The number of repeated samples of size n to take [default: 1]
  --n_stages <num>        How many stages of treatment to simulate [default: 3]
  --seed <num>            Number to use as seed for random number generation [default: 20170411]
  --scenario <scenario>   Which scenario to simulate [default: simple]

Possible scenarios are:
  simple      Only the useful variables are simualted
  int         Only useful variables are simulated, with interactions between them and treatment
  noise       Same as 'simple' with 100 noise variables simulated
  noise_pred  Same as 'noise' with 10 of the noise variables correlated with the outcome
" -> doc

pacman::p_load(tidyverse)
opts <- docopt::docopt(doc)

main <- function(seed, n_subjects, n_samples, n_stages, scenario, output_file,
                 dependencies) {

  if (!(scenario %in% c("simple", "int", "noise", "noise_pred"))) {
    stop(paste0("'", scenario, "'", " is not an alias of any scenario."))
  } else if (scenario == "simple") {
    int <- F
    noise <- F
    noise_pred <- F
  } else if (scenario == "int") {
    int <- T
    noise <- F
    noise_pred <- F
  } else if (scenario == "noise") {
    int <- F
    noise <- T
    noise_pred <- F
  } else if (scenario == "noise_pred") {
    int <- F
    noise <- T
    noise_pred <- T
  }

  lapply(dependencies, source)

  dat <- map_df(
    1:n_samples,
    ~ sim(
      n_subjects = n_subjects,
      int = int,
      noise = noise,
      noise_pred = noise_pred,
      n_stages = n_stages,
      seed = seed + .x
    ) %>% mutate(samp = .x)
  )

  write_rds(dat, path = output_file)
}


main(as.numeric(opts$seed), as.numeric(opts$n_subjects), as.numeric(opts$n_samples),
     as.numeric(opts$n_stages), opts$scenario, opts$output, opts$dependencies)
