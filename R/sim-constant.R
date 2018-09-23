"Simulate clinical trial data for constant doses, given baseline conditions

Usage:
  sim-const.R <input> (--dependencies <files> --output <file>) [--n_stages <num>]
  sim-const.R -h | --help

Arguments:
  -h --help               Show this screen
  <input>                 .rds file with baseline patient condtions
  --dependencies <files>  Files this program depends on
  --output <file>         Path to output .rds file
  --n_stages <num>        How many stages of treatment to simulate [default: 3]
" -> doc

pacman::p_load(tictoc)
opts <- docopt::docopt(doc)

main <- function(baseline_file, n_stages, output_file, dependencies) {

  lapply(dependencies, source)

  dat <- readRDS(baseline_file)

  n_subjects <- nrow(dat)
  doses <- seq(from = 0.1, to = 1, by = 0.1)

  dat <- dat[rep(dat$ID, each = length(doses)), ]
  dat$dose <- rep(doses, n_subjects)

  out <- data.frame()
  for (i in 0:(n_stages - 1)) {
    dat$month <- rep(i, nrow(dat))
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
    dat$tumor_mass <- dat$M_next
    dat$toxicity <- dat$W_next
  }

  saveRDS(out, file = output_file, compress = F)
}

tic()
main(baseline_file = opts$input,
     n_stages = as.numeric(opts$n_stages),
     output_file = opts$output,
     dependencies = opts$dependencies)
toc()
