# path definition
path = "C:/Users/monar/Google Drive/Arbeit/2021_Drube_et_al/210929_other statistics"

setwd(path)

# packages
wants <- c("openxlsx",
           "rstatix", 
           "multcomp", 
           "stringr", 
           "readxl", 
           "dplyr")
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])

lapply(wants, require, character.only = TRUE)

##### import data #####
file_names = list.files(pattern="*.xlsx", full.names = TRUE)
file_names <- c(file_names[which(str_detect(file_names, "Losartan_input"))],
                file_names[which(str_detect(file_names, "Tolvaptan.xlsx"))])
experiments <- substr(file_names, 3, nchar(file_names)-5)

data <- vector("list")
i = 1
for (file in substr(file_names, 3, nchar(file_names))){
  sheet = 3
  data[[experiments[i]]] <- read.xlsx(
    file, sheet= sheet)
  i = i +1
}

#### functions ####
stars_annotation <- function(p_values){
  stars <- c()
  for (p in p_values){
    if (p > 0.05){stars <- c(stars, " ")}
    if (between(p, 0.01, 0.05)){stars <- c(stars, "*")}
    if (between(p, 0.001, 0.01)){stars <- c(stars, "**")}
    if (p < 0.001){stars <- c(stars, "***")}
  }
  return(stars)
}
name_contrasts <- function(contrast_matrix, levels){
  contrast_names <- c()
  for (line in 1:length(contrast_matrix[,1])){
    pos = 1
    temp <- c()
    for (t in contrast_matrix[line,]){
      if (!t == 0){
        temp <- c(temp, levels[pos])
      }
      pos = pos+1
    }
    contrast_names <- c(contrast_names, paste(temp[1], temp[2], sep = " - "))
  }
  rownames(contrast_matrix) <- contrast_names
  return(contrast_matrix)
}
baseStim_dunnett <- function(data_from_list,
                             contrast_matrix, 
                             base_or_stim, 
                             side,
                             bonferroni){
  data_temp <- as.data.frame(data_from_list)
  level_names <-levels(as.factor(data_from_list$condition))
  contrasts <- name_contrasts(contrast_matrix, level_names)
  
  data_temp <- subset(data_temp, state == base_or_stim)
  data_temp[,1] <- as.factor(data_temp[,1])
  
  aov_temp <- aov(BRET ~ condition, data_temp)
  if (bonferroni == TRUE){
    mp_temp <- summary(glht(aov_temp, linfct = mcp(condition = contrasts),
                            alternative = side,
                            test = adjusted("bonferroni")))
    print("bonf")
  }else{
    mp_temp <- summary(glht(aov_temp, linfct = mcp(condition = contrasts),
                            alternative = side))
    print("no bonf")
  }
  
  
  #details
  contrast_temp <- names(mp_temp$test$coefficients)
  data_results <- as.data.frame(contrast_temp)
  data_results$effect_size <- mp_temp$test$coefficients
  
  if (bonferroni == TRUE){
    conf_int <- confint(glht(aov_temp, 
                             linfct = mcp(condition = contrasts),
                             alternative = side,
                             test = adjusted("bonferroni")))
  }else{
    conf_int <- confint(glht(aov_temp, 
                             linfct = mcp(condition = contrasts),
                             alternative = side))
  }
  data_results$lower <- conf_int$confin[,2]
  data_results$upper <- conf_int$confin[,3]
  data_results$tstat <- mp_temp$test$tstat
  data_results$df <- mp_temp$df
  data_results$pvalue <- mp_temp$test$pvalues
  data_results$signi <- stars_annotation(mp_temp$test$pvalues)
  
  return(data_results)
}
save_results <- function(data_list, 
                         contrasts, 
                         base_or_stim, 
                         side, 
                         bonferroni){
  results <- vector("list")
  n = 0
  for (set in data_list){
    n = n +1
    results[[names(data)[n]]] <- baseStim_dunnett(set, 
                                                  contrasts, 
                                                  base_or_stim, 
                                                  side,
                                                  bonferroni)
  }
  return(results)
}
results_export <- function(results, filename){
  n = 0
  for (set in results){
    n = n +1
    temp_df <- as.data.frame(set)
    temp_df <- cbind(a = rep(names(results)[n], length(set[,1])), temp_df)
    if (n == 1){
      results_export <- as.data.frame(temp_df)
    }else{
      results_export <- rbind(results_export, temp_df)
    }
  }
  write.xlsx(results_export, filename)
  return(results_export)
}


##### Losartan analysis #####
levels(as.factor(data[[1]]$condition))
# for base only: test against respective dQ baseline
contrasts_base_los <-  rbind(c(-1, 0, 1, 0, 0, 0),
                             c(-1, 0, 0, 0, 1, 0),
                             c(0, -1, 0, 1, 0, 0),
                             c(0, -1, 0, 0, 0, 1))
Los_data <-list(data[[1]])
base_los <- save_results(Los_data, contrasts_base_los, "base", "greater", FALSE)

# stim comparison --> dQ vs. dQ+Los; GRK2 vs. GRK2+Los; GRK6 vs. GRK6+Los
contrasts_stim_los <- rbind(c(1, -1, 0, 0, 0, 0),
                            c(0, 0, 1, -1, 0, 0),
                            c(0, 0, 0, 0, 1, -1))
stim_los <- save_results(Los_data, contrasts_stim_los, "stim", "two.sided", TRUE)



##### Tolvaptan analysis #####
levels(as.factor(data[[2]]$condition))
# for base only: test against respective dQ baseline
contrasts_base_tol <-  rbind(c(1, 0, -1, 0, 0, 0, 0, 0),
                             c(0, 0, -1, 0, 1, 0, 0, 0),
                             c(0, 0, -1, 0, 0, 0, 1, 0),
                             c(0, 1, 0, -1, 0, 0, 0, 0),
                             c(0, 0, 0, -1, 0, 1, 0, 0),
                             c(0, 0, 0, -1, 0, 0, 0, 1))
Tol_data <-list(data[[2]])
base_tol <- save_results(Tol_data, contrasts_base_tol, "base", "greater", FALSE)
# stim comparison 
contrasts_stim_tol <- rbind(c(-1, 1, 0, 0, 0, 0, 0, 0),
                            c(0, 0, -1, 1, 0, 0, 0, 0),
                            c(0, 0, 0, 0, -1, 1, 0, 0),
                            c(0, 0, 0, 0, 0, 0, -1, 1))
stim_tol <- save_results(Tol_data, contrasts_stim_tol, "stim", "two.sided", TRUE)


##### Export #####
Export <- list(base_los[[1]], stim_los[[1]], base_tol[[1]], stim_tol[[1]])
names(Export) <- c("Losartan_base", "Losartan_stim", "Tolvaptan_base", "Tolvaptan_stim")
results_export(Export, "Losartan_Tolvaptan_results.xlsx")
