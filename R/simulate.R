"Simulate adaptive clinical trial data

Usage:
  simulate.R (--dependencies <files> --output <file>) [--seed <num> --n_subjects <num> --n_stages <num> --scenario <scenario> --baseline-only]
  simulate.R -h | --help

Arguments:
  -h --help               Show this screen
  --dependencies <files>  Files this program depends on
  --output <file>         Path to output .rds file
  --seed <num>            Seed for random number generation [default: 20170411]
  --n_subjects <num>      Participants to randomize [default: 1000]
  --n_stages <num>        Stages of treatment to simulate [default: 3]
  --scenario <scenario>   Scenario to simulate [default: simple]
  --baseline-only         Simulate only baseline condtions (for testing treatment regimes)

Possible scenarios are:
  simple          Only useful variables are simualted
  noise           Same as 'simple' with 100 noise variables simulated
  noise-pred      'noise' with 10 noise variables correlated with the outcome
  simple-int      Same as 'simple' with two variables that interact with treatment
  noise-int       Same as 'noise' with two variables that interact with treatment
  noise-pred-int  Same as 'noise-pred' with two variables that interact with treatment
" -> doc

pacman::p_load(tictoc)
opts <- docopt::docopt(doc)

main <- function(seed, n_subjects, n_stages, scenario, output_file,
                 dependencies, baseline_only) {

  scenario_split <- unlist(strsplit(scenario, '-'))
  int <- tail(scenario_split, 1) == 'int'
  if (int) {
    scenario <- paste(head(scenario_split, -1), collapse = '-')
  }

  if (!(scenario %in% c("simple", "noise", "noise-pred"))) {
    stop(paste0("'", scenario, "'", " is not an alias of any scenario."))
  } else if (scenario == "simple") {
    noise <- F
    noise_pred <- F
  } else if (scenario == "noise") {
    noise <- T
    noise_pred <- F
  } else if (scenario == "noise-pred") {
    noise <- T
    noise_pred <- T
  }

  lapply(dependencies, source)

  set.seed(seed)

  dat <- data.frame(
    ID = seq_len(n_subjects),
    tumor_mass = runif(n_subjects, min = 0, max = 2),
    toxicity = runif(n_subjects, min = 0, max = 2),
    cW = rep(1, n_subjects),
    cM = rep(1, n_subjects)
  )

  if (int) {
    dat$X1 <- runif(n_subjects, min = 0, max = 1)
    dat$X2 <- runif(n_subjects, min = 0, max = 1)
    dat$cW <- replace(dat$cW, dat$X1 > 0.5, 1.5)
    dat$cM <- replace(dat$cM, dat$X2 > 0.5, 1.5)
  }

  if (noise | noise_pred) {
    Z1 <- lapply(1:5, function(x) rnorm(n_subjects, mean = 1))
    Z2 <- lapply(1:5, function(x) rnorm(n_subjects, mean = -1))
    V <- lapply(1:90, function(x) rnorm(n_subjects))

    names(Z1) <- paste0("Z", 1:5)
    names(Z2) <- paste0("Z", 6:10)
    names(V) <- paste0("V", 1:90)

    dat <- cbind(dat, Z1, Z2, V)
  }

  if (noise_pred) {
    dat$noise_chng <- dat$Z1 + dat$Z2 + dat$Z3 + dat$Z4 + dat$Z5 +
      dat$Z6 + dat$Z7 + dat$Z8 + dat$Z9 + dat$Z10
  }
  if (!noise_pred) {
    dat$noise_chng <- rep(0, n_subjects)
  }

  if (baseline_only) {
    saveRDS(dat, file = output_file, compress = F)
    return(invisible())
  }

  out <- data.frame()
  for (i in 0:(n_stages - 1)) {
    dat$month <- rep(i, nrow(dat))
    dat$dose <- runif(nrow(dat), min = 0, max = 1)
    dat <- Mnext(dat)
    dat <- Wnext(dat)
    dat$beta <- 1 / lambda(dat$M_next, dat$W_next, Z = dat$noise_chng)
    dat$dead <- replace(dat$dead, dat$beta < 1, T)
    if (i < (n_stages - 1)) {
      dat$reward <- replace(log(i + dat$beta), !dat$dead, 0)
    } else {
      dat$reward <- log(i + dat$beta)
    }
    out <- rbind(out, dat)
    if (i < (n_stages - 1)) {
      dat$tumor_mass <- dat$M_next
      dat$toxicity <- dat$W_next
    }
  }

  out$Qhat <- out$reward

  saveRDS(out, file = output_file, compress = F)
}

tic()
main(
  seed = as.numeric(opts$seed),
  n_subjects = as.numeric(opts$n_subjects),
  n_stages = as.numeric(opts$n_stages),
  scenario = opts$scenario,
  output_file = opts$output,
  dependencies = opts$dependencies,
  baseline_only = opts[['baseline-only']]
)
toc()
