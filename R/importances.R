"Create a data frame of variable importances for each model used in Q-learning

Usage:
  importances.R <inputs>... (--output <file>)
  importances.R -h | --help

Arguments:
  -h --help        Show this screen
  <inputs>         .rds files of models
  --output <file>  Path to output .rds file
" -> doc

pacman::p_load(data.table, caret, earth, rpart, tictoc)
opts <- docopt::docopt(doc)

main <- function(inputs, output_file) {

  q_lists <- lapply(inputs, readRDS)

  inputs_split <- strsplit(inputs, split = 'rf|rpart|mars')
  scenarios <- lapply(inputs_split,
                      function(x) {gsub(paste0('^-|.rds$'), '', tail(x, 1))})

  imps_list <- mapply(function(q_list, scenario) {
    samp_list <- lapply(seq_along(q_list), function(samp) {
      month_list <- lapply(2:0, function(month) {
        mod <- q_list[[samp]][[paste0('month', month)]]
        if (class(mod$finalModel) == "earth") {
          dat <- caret::varImp(mod$finalModel, useModel = T, value = "rss")
          dat$importance <- dat$Overall
          dat$Overall <- NULL
          dat$var_nm <- rownames(dat)
          dat$month <- month
          dat$samp <- samp
          dat$scenario <- scenario
          dat$mod <- 'mars'
        } else if (class(mod$finalModel) == "ranger") {
          imps <- ranger::importance(mod$finalModel)
          dat <- data.table(var_nm = names(imps),
                            importance = imps / max(imps) * 100,
                            month = month,
                            samp = samp,
                            scenario = scenario,
                            mod = "rf")
        } else if (class(mod) == "rpart") {
          imps <- mod$variable.importance
          imps <- (imps / max(imps)) * 100
          var_nms <- names(imps)
          var_nms <- union(var_nms, c("tumor_mass", "toxicity", "dose"))
          if (grepl('int', scenario)) {
            var_nms <- union(var_nms, c("X1", "X2"))
          }
          if (grepl('noise', scenario)) {
            var_nms <- union(var_nms, c(paste0("V", 1:90), paste0("Z", 1:10)))
          }
          imps <- c(imps, rep(0, length(var_nms) - length(imps)))
          dat <- data.table(var_nm = var_nms,
                            importance = imps,
                            month = month,
                            samp = samp,
                            scenario = scenario,
                            mod = "rpart")
        }
        return(dat)
      })
      return(data.table::rbindlist(month_list))
    })
    return(data.table::rbindlist(samp_list))
  },
  q_lists, scenarios, SIMPLIFY = F
  )

  imps <- data.table::rbindlist(imps_list, fill = T)

  saveRDS(imps, output_file, compress = F)
}

tic()
main(opts$inputs, opts$output)
toc()
