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

  q_list <- readRDS(q_file)
  n_samples <- length(q_list)

  dat_list <- lapply(seq_len(n_samples), function(s) {

    tic(paste('Sample', s, 'of', n_samples))
    dat$samp <- rep(s, nrow(dat))
    out <- data.frame()
    for (i in 0:(n_stages - 1)) {

      if (i > 0) {
        dat <- dat[!(dat$dead), ]
        dat$tumor_mass <- dat$M_next
        dat$toxicity <- dat$W_next
      }

      n <- nrow(dat)
      dat$month <- rep(i, n)
      doses <- seq(from = 0, to = 1, by = 0.01)
      dat_expand <- dat[rep(seq_len(n), each = length(doses)), ]
      dat_expand$dose <- rep(doses, n)
      dat_expand$pred <- predict(q_list[[s]][[paste0('month', i)]], dat_expand)

      # subset to smallest dose which produces the longest expected survival
      dat <- dat_expand[dat_expand[, .I[which.max(pred)], by = ID]$V1]

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
    toc()
    return(out)
  })

  out <- rbindlist(dat_list)

  saveRDS(out, file = output_file, compress = F)
}

tic()
main(baseline_file = opts[['base-file']],
     q_file = opts[['q-file']],
     n_stages = as.numeric(opts$n_stages),
     output_file = opts$output,
     dependencies = opts$dependencies)
toc()
