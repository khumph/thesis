renderDocument <- function(sim_path) {
  rmarkdown::render("./code/testbed.Rmd", params = list(sim = sim_path))
}

renderDocument("sim.R")
renderDocument("sim-noise.R")
renderDocument("sim-int.R")
