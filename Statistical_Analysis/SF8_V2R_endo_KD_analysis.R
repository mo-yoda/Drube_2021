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
file_names <- file_names[which(str_detect(file_names, "SF7"))] 
experiments <- substr(file_names, 3, nchar(file_names)-5)

data <- vector("list")
i = 1
for (file in substr(file_names, 3, nchar(file_names))){
  sheet = 3
  data[[experiments[i]]] <- read.xlsx(
    file, sheet= sheet)
  i = i +1
}

##### functions ####
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

base_greater_dunnett <- function(data_from_list, contrast_matrix){
  data_temp <- as.data.frame(data_from_list)
  level_names <-levels(as.factor(data_from_list$condition))
  contrasts <- name_contrasts(contrast_matrix, level_names)
  
  data_temp <- subset(data_temp, state == "base")
  data_temp[,1] <- as.factor(data_temp[,1])
  
  aov_temp <- aov(BRET ~ condition, data_temp)
  dunn_temp <- summary(glht(aov_temp, linfct = mcp(condition = contrasts), alternative = "greater"))
  
  #details
  contrast_temp <- names(dunn_temp$test$coefficients)
  data_results <- as.data.frame(contrast_temp)
  data_results$effect_size <- dunn_temp$test$coefficients
  conf_int <- confint(glht(aov_temp, 
                           linfct = mcp(condition = contrasts),
                           alternative = "greater"))
  data_results$lower <- conf_int$confin[,2]
  data_results$upper <- conf_int$confin[,3]
  data_results$tstat <- dunn_temp$test$tstat
  data_results$df <- dunn_temp$df
  data_results$pvalue <- dunn_temp$test$pvalues
  data_results$signi <- stars_annotation(dunn_temp$test$pvalues)
  
  return(data_results)
}
save_results <- function(data_list, contrasts){
  results <- vector("list")
  n = 0
  for (set in data_list){
    n = n +1
    results[[names(data)[n]]] <- base_greater_dunnett(set, contrasts)
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


##### analysis ####

# order of factor names for AT1R:  dGRK2/3/6+EV, dQ+EV, dQ+GRK5, dQ+GRK5-KD
contrast_AT1R<- rbind(c(0, -1, 1, 0),
                       c(0, -1, 0, 1),
                       c(1, -1, 0, 0))
AT1R_data <- list(data[["SF7_GRK5_endo_KD"]], data[["SF7_GRK6_endo_KD"]])
names(AT1R_data) <- names(data)[1:2]
results_AT1R <- save_results(AT1R_data, contrast_AT1R)
results_export(results_AT1R, "Results_SF9_AT1R_KD_endo.xlsx")

# order of factor names for V2R: "Control+EV" "dQ+EV" "dQ+GRK2-KD" "dQ+GRK6-KD"
contrast_V2R<- rbind(c(0, -1, 1, 0),
                     c(0, -1, 0, 1),
                     c(1, -1, 0, 0))

V2R_data <- list(data[["SF7_V2R_KD"]])
names(V2R_data) <- names(data)[3]
results_V2R <- save_results(V2R_data, contrast_V2R)
results_export(results_V2R, "Results_SF8_V2R_KD.xlsx")



