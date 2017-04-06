getImps <- function(nums) {
  if {classQ$mod_list[[nums]]$finalmodel == "earth"} {
    
  }
  varImp(Q$mod_list[[nums]])$importance %>%
    rownames_to_column("X") %>% 
    mutate(month = nums - 1) %>% 
    rename(importance = Overall) %>% 
    filter(importance > 0)
}

Q <- Q_list$MARS


Q$mod_list[[1]]$finalModel %>% plotmo()
mod <- Q$mod_list[[2]]$finalModel
mod %>% summary()
evimp(mod)

imps <- map_df(1:6, ~ getImps(.))

ggplot(imps, aes(x = as.factor(month), y = importance, fill = X)) +
  geom_bar(stat = "identity", position = "dodge")

varImp(Q$mod_list[[5]]$finalModel, useModel = T)

varImp(Q$mod_list[[1]])$importance
