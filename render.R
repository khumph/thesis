library(pacman)
p_load(stringr)

renderDocument <- function(int, noise) {
  rmarkdown::render("testbed.Rmd", params = list(int = int, noise = noise),
                    output_file = str_c("testbed-",
                                        ifelse(noise, "noise",
                                               ifelse(int, "int",
                                                      "orig")), ".html"))
}

renderDocument(int = F, noise = F)
renderDocument(int = T, noise = F)
renderDocument(int = F, noise = T)

