"Join clinical trial data created using different treatment regimes

Usage:
  join.R <inputs>... (--output <file>)
  join.R -h | --help

Arguments:
  -h --help        Show this screen
  <inputs>         .rds data files to join
  --output <file>  Path to output .rds file
" -> doc

pacman::p_load(data.table, tictoc)
opts <- docopt::docopt(doc)

main <- function(inputs, output) {

  x <- lapply(inputs, readRDS)
  dat_names <- gsub('^results/data-', '', inputs)
  dat_names <- gsub(paste0('.rds$'), '', dat_names)
  x <- mapply(function(X, Y) {
    X$mod <- Y[1]
    X$scenario <- paste(tail(Y, -1), collapse = '_')
    if (Y[1] == "best") {
      names(X)[names(X) == 'W'] <- 'toxicity'
      names(X)[names(X) == 'M'] <- 'tumor_mass'
    }
    if (Y[1] == "constant") {
      X$mod <- as.character(X$dose)
    }
    return(X)
  }, X = x, Y = strsplit(dat_names, split = "-"))
  x <- rbindlist(x, fill = T)

  x <- x[ , tot_reward := sum(reward, na.rm = T), by = .(ID, mod, samp, scenario)]

  saveRDS(x, file = output, compress = F)
}

tic()
main(opts$inputs, opts$output)
toc()

