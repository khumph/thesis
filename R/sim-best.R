"Simulate clinical trial data for the best possible three stage treatment, given baseline conditions

Usage:
  sim-best.R <input> (--dependencies <files> --output <file>) [--n_stages <num>]
  sim-best.R -h | --help

Arguments:
  -h --help               Show this screen
  <input>                 .rds file with baseline patient condtions
  --dependencies <files>  Files this program depends on
  --output <file>         Path to output .rds file
  --n_stages <num>        How many stages of treatment to simulate [default: 3]

Note: more than three stages may cause you to run out of memory, and/or take a very long time
" -> doc

pacman::p_load(data.table, tictoc)
opts <- docopt::docopt(doc)

filter_best <- function(x) {
  # select rows/dose sequences with longest expected survival for each person
  x <- x[x[, .I[beta == max(beta)], by = ID]$V1]
  # of those that give longest survival, select row with smallest total dose
  x[ , totDose := rowSums(.SD, na.rm = T), .SDcols = grep("dose", names(x), value = T)]
  x[x[, .I[totDose == min(totDose)], by = ID]$V1]
}

main <- function(baseline_file, n_stages, output_file, dependencies) {

  lapply(dependencies, source)

  dat <- readRDS(baseline_file)
  dat <- as.data.table(dat)

  dat[, ('M0') := dat$tumor_mass]
  dat[, ('W0') := dat$toxicity]
  dat <- dat[, c('ID', 'W0', 'M0', 'cW', 'cM', 'noise_chng')]

  # Too many subjects at once will exhaust memory, split into groups
  dat[ , grp := (dat$ID - 1) %/% 20]

  possible_doses <- seq(0, 1, 0.01)

  x <- lapply(
    unique(dat$grp),
    function(y) {
      tic(paste("Group", y + 1, "of", length(unique(dat$grp))))
      x <- dat[grp == y]
      dead_already <- data.table()
      for (i in 0:(n_stages - 1)) {
        n <- nrow(x)
        x <- x[rep(seq_len(n), each = length(possible_doses)), ]

        dose_nm <- paste0('dose', i)
        M_now <- paste0('M', i)
        W_now <- paste0('W', i)
        M_next <- paste0('M', i + 1)
        W_next <- paste0('W', i + 1)
        x[ , (dose_nm) := rep(possible_doses, n)]
        x[ , (M_next) := updateM(x[[M_now]], x[[W_now]], x[[dose_nm]], c = cM)]
        x[ , (W_next) := updateW(x[[M_now]], x[[W_now]], x[[dose_nm]], c = cW)]
        x[ , beta := 1 / lambda(x[[M_next]], x[[W_next]], Z = x$noise_chng)]
        x[ , dead := x$beta <= 1]

        if (any(x$dead) & i < (n_stages - 1)) {
          # select people who will die regardless of dose
          deads <- x[x[, .I[all(dead)], by = ID]$V1]
          deads[ , (paste0('reward', i)) := log(i + deads$beta)]
          if (nrow(deads) > 0) {
            # store the best dose sequences for people who will die
            deads <- filter_best(deads)
            dead_already <- rbind(dead_already, deads, fill = T)
          }
          # next dose can't bring people back from dead, so remove dead
          x <- x[!(dead)]
          x[ , (paste0('reward', i)) := rep(0, nrow(x))]
        }

        if (i == (n_stages - 1)) {
          x[ , (paste0('reward', i)) := log(i + x$beta)]
        }

        if (all(x$dead)) {
          dead_already
        }
      }
      x <- filter_best(x)
      x <- rbind(x, dead_already, fill = T)
      toc()
      return(x)
    }
  )

  x <- rbindlist(x, fill = T)

  x[, c('M3', 'W3') := NULL]

  x <- reshape(x, varying = grep("\\d$", names(x), value = T),
          timevar = "month", idvar = "ID",
          direction = "long", sep = '')

  saveRDS(x, file = output_file, compress = F)
}

tic('Total time')
main(baseline_file = opts$input,
     n_stages = as.numeric(opts$n_stages),
     output_file = opts$output,
     dependencies = opts$dependencies)
toc()
