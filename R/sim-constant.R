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
  dat$mod <- as.character(dat$dose)

  out <- data.frame()
  for (i in 0:(n_stages - 1)) {
    if (i > 0) {
      dat$tumor_mass <- dat$M_next
      dat$toxicity <- dat$W_next
    }
    dat$month <- rep(i, nrow(dat))
    dat$W_next <- updateW(dat$tumor_mass, dat$toxicity, dat$dose, c = dat$cW)
    dat$W_next <- replace(dat$W_next, dat$dead, NA_real_)
    dat$M_next <- updateM(dat$tumor_mass, dat$toxicity, dat$dose, c = dat$cM)
    dat$M_next <- replace(dat$M_next, dat$dead, NA_real_)
    dat$beta <- 1 / lambda(dat$M_next, dat$W_next, Z = dat$noise_chng)
    dat$dead <- dat$beta < 1
    dat$reward <- log(i + dat$beta)
    if (i < (n_stages - 1)) {
      dat$reward <- replace(dat$reward, !dat$dead, 0)
    }
    out <- rbind(out, dat)
  }

  saveRDS(out, file = output_file, compress = F)
}

tic()
main(baseline_file = opts$input,
     n_stages = as.numeric(opts$n_stages),
     output_file = opts$output,
     dependencies = opts$dependencies)
toc()
