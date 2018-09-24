"Simulate clinical trial data for the best predicted dose sequence, given baseline conditions and fitted models

Usage:
  sim-q.R <base-file> <q-file> (--dependencies <files> --output <file>) [--n_stages <num>]
  sim-q.R -h | --help

Arguments:
  -h --help               Show this screen
  <base-file>             .rds file with baseline patient condtions
  <q-file>                .rds file with fitted models
  --dependencies <files>  Files this program depends on
  --output <file>         Path to output .rds file
  --n_stages <num>        How many stages of treatment to simulate [default: 3]
" -> doc

pacman::p_load(data.table, tictoc)
opts <- docopt::docopt(doc)

main <- function(baseline_file, q_file, n_stages, output_file, dependencies) {

  lapply(dependencies, source)

  dat <- readRDS(baseline_file)
  dat <- as.data.table(dat)

  q <- readRDS(q_file)

  out <- data.frame()
  for (i in 0:(n_stages - 1)) {
    n <- nrow(dat)
    dat$month <- rep(i, n)

    doses <- seq(from = 0, to = 1, by = 0.01)
    dat_expand <- dat[rep(seq_len(n), each = length(doses)), ]
    dat_expand$dose <- rep(doses, n)
    dat_expand$pred <- predict(q[[paste0('month', i)]], dat_expand)

    # filter to the first minimum, also the also the min with the
    # smallest dose, because does are in ascending order:
    dat <- dat_expand[dat_expand[, .I[which.max(pred)], by = ID]$V1]

    dat <- Mnext(dat)
    dat <- Wnext(dat)
    dat$beta <- 1 / lambda(dat$M_next, dat$W_next, Z = dat$noise_chng)
    dat$dead <- replace(dat$dead, dat$beta < 1, T)
    if (i < (n_stages - 1)) {
      dat$reward <- replace(log(i + dat$beta), !dat$dead, 0)
    }
    if (i == n_stages - 1) {
      dat$reward <- log(i + dat$beta)
    }
    out <- rbind(out, dat)
    if (i < (n_stages - 1)) {
      dat <- dat[!(dat$dead), ]
      dat$tumor_mass <- dat$M_next
      dat$toxicity <- dat$W_next
    }
  }

  saveRDS(out, file = output_file, compress = F)
}

tic()
main(baseline_file = opts[['base-file']],
     q_file = opts[['q-file']],
     n_stages = as.numeric(opts$n_stages),
     output_file = opts$output,
     dependencies = opts$dependencies)
toc()
