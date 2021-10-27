# path definition
path = "C:/path/to/folder/"

setwd(path)

# packages
wants <- c("openxlsx",
           "ggpubr", 
           "rstatix", 
           "multcomp", 
           "stringr", 
           "readxl", 
           "dplyr")
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])

lapply(wants, require, character.only = TRUE)

# import data -> xlsx files with SF1 in filename
file_names = list.files(pattern="*.xlsx", full.names = TRUE)
file_names <- file_names[which(str_detect(file_names, "SF1"))] 
experiments <- substr(file_names, 3, nchar(file_names)-5)

data <- vector("list")
i = 1
for (file in substr(file_names, 3, nchar(file_names))){
  sheet = which(str_detect(excel_sheets(file), "_format"))
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
aov_dunnett <- function(data_from_list){
  data_temp <- as.data.frame(data_from_list)
  data_temp[,length(data_temp)-1] <- as.factor(data_temp[,length(data_temp)-1])
  data_temp$result <- data_temp[,length(data_temp)]

  data_model <- aov(result ~ condition, data = data_temp)
  
  ggqqplot(residuals(data_model))
  # levene_test(data_temp[,2] ~ data_temp[,1], data = data_temp)
  
  data_Dunn <- summary(glht(data_model, linfct = mcp(condition = "Dunnett")))
  
  # details
  contrast<-names(data_Dunn$test$coefficients)
  data_results <- as.data.frame(contrast)
  data_results$effect_size <- data_Dunn$test$coefficients

  conf_int <- confint(glht(data_model, 
                           linfct = mcp(condition = "Dunnett")))
  data_results$lower <- conf_int$confin[,2]
  data_results$upper <- conf_int$confin[,3]
  
  data_results$tstat <- data_Dunn$test$tstat
  data_results$df <- data_Dunn$df
  data_results$pvalue <- data_Dunn$test$pvalues
  data_results$signi <- stars_annotation(data_Dunn$test$pvalues)

  return(data_results)
}
save_results <- function(data_list){
  results <- vector("list")
  n = 0
  for (set in data_list){
    n = n +1
    results[[experiments[n]]] <- aov_dunnett(set)
  }
  return(results)
}
results_export <- function(results){
  n = 0
  for (set in results){
    n = n +1
    temp_df <- as.data.frame(set)
    print(temp_df)
    temp_df <- cbind(a = rep(experiments[n], length(set[,1])), temp_df)
    if (n == 1){
      results_export <- as.data.frame(temp_df)
    }else{
      results_export <- rbind(results_export, temp_df)
    }
  }
  write.xlsx(results_export, "Results_SupplF1.xlsx")
  return(results_export)
}

# perfom analysis
results <- save_results(data)
export <- results_export(results)

