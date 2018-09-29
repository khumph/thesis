"Perform Q-learning using tree-based methods to fit Q-functions

Usage:
  learn.R <input> (--output <file> --method <method>)
  learn.R -h | --help

Arguments:
  -h --help               Show this screen
  <input>                 .rds data file to use for learning
  --output <file>         Path to .rds file to save fitted models in
  --method <method>       String indicating method to fit Q-functions

Possible methods are:
  rpart  Fit Q-functions using CART regression trees
  mars   Fit Q-functions using MARS
  rf     Fit Q-funcitons using random forests
" -> doc

pacman::p_load(data.table, tictoc)
opts <- docopt::docopt(doc)
if (opts$method == 'mars') {
  pacman::p_load(caret, earth)
}
if (opts$method == 'rpart') {
  pacman::p_load(rpart)
}
if (opts$method == 'rf') {
  pacman::p_load(ranger, e1071)
}

main <- function(data_file, output_file, method) {

  dat_all <- readRDS(data_file)

  n_stages <- length(unique(dat_all$month))
  n_samples <- length(unique(dat_all$samp))

  predictor_names <- c("tumor_mass", "toxicity",
                       grep("\\X|\\Z|\\V", names(dat_all), value = T))
  treatment_name <- "dose"
  form <- as.formula(paste("Qhat ~ ", paste(c(treatment_name, predictor_names),
                                    collapse = " + ")))

  if (method == 'mars') {
    # predictors in PREDICTORS are allowed to interact with predictors in PARENTS
    # but no other interactions are allowed

    PREDICTORS <- predictor_names
    PARENTS <-  treatment_name
    allowFunct <- function(degree, pred, parents, namesx) {
      if (degree <= 1)
        return(TRUE)
      predictor <- namesx[pred]
      parents <- namesx[parents != 0]
      if ((any(predictor %in% PREDICTORS) &&
           any(parents %in% PARENTS)) ||
          (any(predictor %in% PARENTS) &&
           any(parents %in% PREDICTORS))) {
        return(TRUE)
      }
      FALSE
    }
  }

  q_list <- lapply(seq_len(n_samples), function(s) {

    tic(paste('Sample', s, 'of', n_samples))

    dat <- dat_all[dat_all$samp == s, ]

    set.seed(20170128)
    mod_list <- list()
    for (i in (n_stages - 1):0) {

      x <- dat[dat$month == i, ]
      x <- data.table::as.data.table(x)

      # remove people who are already dead at this stage, can't optimize dose if dead
      x <- x[!is.na(x$M_next), ]

      if (method == "rpart") {
        mod <- rpart::rpart(form, x,
                            control = rpart::rpart.control(cp = 1e-5, maxcompete = 0,
                                                           maxsurrogate = 0, minbucket = 5))
        mod <- rpart::prune(mod, cp = mod$cptable[which.min(mod$cptable[, "xerror"]), "CP"])
      }
      if (method == 'mars') {
        mod <- caret::train(form = form, data = x, method = 'gcvEarth',
                            trControl = caret::trainControl(method = "none"),
                            nk = 200,
                            allowed = allowFunct,
                            tuneGrid = expand.grid(degree = 2),
                            na.action = na.omit)
      }
      if (method == 'rf') {
        mod <- caret::train(form = form, data = x, method = 'ranger',
                            tuneGrid = expand.grid(mtry = length(predictor_names),
                                                   min.node.size = 5,
                                                   splitrule = "variance"),
                            trControl = caret::trainControl(method = "none"),
                            importance = "impurity",
                            always.split.variables = "dose",
                            num.trees = 250,
                            na.action = na.omit)
      }
      mod_list[[paste0('month', i)]] <- mod

      n <- nrow(x)
      doses <- seq(from = 0, to = 1, by = 0.01)
      x_expand <- x[rep(seq_len(n), each = length(doses)), ]
      x_expand[ , dose := rep(doses, n)]
      x_expand[ , pred := predict(mod, x_expand)]
      max_preds <- x_expand[ , max(pred), by = ID]

      if (i > 0) {
        ind <- (dat$month == i - 1) & (dat$ID %in% max_preds$ID)
        dat[ind, ]$Qhat <- dat[ind, ]$reward + max_preds$V1
      }
    }
    toc()
    return(mod_list)
  })

  saveRDS(q_list, file = output_file, compress = F)
}

tic('All samples')
main(opts$input, opts$output, opts$method)
toc()
