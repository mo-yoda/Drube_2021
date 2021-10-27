# path definition
path = "C:/Users/monar/Google Drive/Arbeit/2021_Drube_et_al/210909_Drube et al._vehicle vs. stim"
setwd(path)

# format two.sided


# packages
wants <- c("openxlsx", "readxl", "stringr", "lsr", "dplyr")
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
lapply(wants, require, character.only = TRUE)


# import data (_format and _mean)
sheets <- excel_sheets("Vehicle_vs_Stim_Master_format.xlsx")
sheets_format <- which(str_detect(sheets, "_format"))

data_format <- vector("list")

for (sheet in sheets_format){
  data_format[[str_remove(sheets[sheet], "_format")]] <- read.xlsx(
    "Vehicle_vs_Stim_Master_format.xlsx", sheet= sheet)
}


#### preprocess data ######

# transform excluded values in _format to NA
r = 0
for (experiment in data_format){
  r = r+1
  data_format[[names(data_format)[r]]][1] <- as.factor(experiment[[1]])
  for (i in 2:length(experiment)){
    if (class(experiment[[i]]) == "character"){
      t = 1
      # print("column detected")
      # print(i)
      for (value in experiment[[i]]){
        if (is.na(as.numeric(value))){
          # print("value detected")
          # print(value)
          # print(data_format[[names(data_format)[r]]][t,i])
          data_format[[names(data_format)[r]]][t,i] <- NA
          # print(data_format[[names(data_format)[r]]][t,i])
        }
        t = t + 1
      }
      data_format[[names(data_format)[r]]][i] <- as.numeric(experiment[[i]])
    }
  }
}


#### analysis #####

analysis_format <- function(data_list, alternative){
  exp <- c()
  GPCR <- c()
  bArr <- c()
  condition <- c()
  p_value <- c()
  effect_size <- c()
  t_stat <- c()
  df <- c()
  lower_CI <- c()
  upper_CI <- c()
  signi <- c()
  
  for (experiment in names(data_list)){
    temp_data <- data_list[[experiment]]

    for (i in 2:length(temp_data)){
      exp <- c(exp, experiment)
      
      condition <- c(condition, names(temp_data)[i])
      
      temp_test <- t.test(temp_data[,i] ~ temp_data[,1], alternative = alternative)
      
      p_value <- c(p_value, temp_test$p.value) # p value
      t_stat <- c(t_stat, temp_test$statistic) # t statistic
      df <- c(df, temp_test$parameter) # degrees of freedom
      lower_CI  <- c(lower_CI, temp_test$conf.int[1]) # lower CI
      upper_CI  <- c(upper_CI, temp_test$conf.int[2]) # upper CI
      # Cohens's D
      # for effect size, calculate pooled sd (difference of means/pooled sd = effect size)
      pooledSD <- sqrt(
        ((sd(subset(temp_data, well == "stim")[,i], na.rm = TRUE))^2 + 
            (sd(subset(temp_data, well == "vehicle")[,i], na.rm = TRUE))^2)
        /2)
      
      effect_size <- c(effect_size,
                       (as.numeric(temp_test$estimate[1])
                        - as.numeric(temp_test$estimate[2]))
                       /pooledSD) # effect size

      if (temp_test$p.value > 0.05){
        signi <- c(signi, " ")}
      
      if (between(temp_test$p.value, 0.01, 0.05)){
        signi <- c(signi, "*")}
      
      if (between(temp_test$p.value, 0.001, 0.01)){
        signi <- c(signi, "**")}
      
      if (temp_test$p.value < 0.001){
        signi <- c(signi, "***")}
    }
  }
  GPCR <- exp %>%
    str_remove("-bArr1")  %>%
    str_remove("-bArr2")
  
  bArr <- exp %>%
    str_sub(-5)
  
  data_results <- as.data.frame(exp)
  data_results$GPCR <- GPCR
  data_results$bArr <- bArr
  data_results$condition <- condition
  data_results$effect_size <- effect_size
  data_results$lower_CI <- lower_CI
  data_results$upper_CI <- upper_CI
  data_results$t_statistic <- t_stat
  data_results$df <- df
  data_results$p_value <- p_value
  data_results$signi <- signi
  
  return(data_results)
  
}

data_format_results <- analysis_format(data_format, "two.sided")



#### export #####

remove(results_wb)
results_wb<-createWorkbook()
addWorksheet(results_wb, "results")
writeData(results_wb, "results", data_format_results, 
          startRow = 1, 
          startCol = 1)
saveWorkbook(results_wb, 
             file = "data_format_results.xlsx", 
             overwrite = TRUE)

