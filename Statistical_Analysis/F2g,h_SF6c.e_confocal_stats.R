# path definition
# home PC
path = "C:/Users/monar/Google Drive/Arbeit/2021_Drube_et_al/211011_confocal statistics"
setwd(path)

# packages
wants <- c("openxlsx", 
           "readxl",
           "stringr",
           "rstatix")
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])

lapply(wants, require, character.only = TRUE)


# import data
sheets <- excel_sheets("microscopy_data.xlsx")
sheets_format <- which(str_detect(sheets, "_format"))

data <- vector("list")

for (sheet in sheets_format){
  data[[str_remove(sheets[sheet], "_format")]] <- read.xlsx(
    "microscopy_data.xlsx", sheet= sheet)
}

# functions

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

paired_test <- function(data_from_list){
  data_temp <- as.data.frame(data_from_list)
  data_temp$condition <- as.factor(data_temp$condition)
  data_temp$state <- as.factor(data_temp$state)
  
  levels <- levels(data_temp$condition)
  base_stim <- as.data.frame(levels)
  row = 0
  for (i in levels){
    temp = subset(data_temp, condition == as.character(i))
    row = row +1
    temp_test <-t.test(temp$colocalisation~temp$state, paired = TRUE)
    
    base_stim[row, 2] <- cohens_d(temp, colocalisation ~ state, paired = TRUE)$effsize
    base_stim[row, 3] <- temp_test$conf.int[1]
    base_stim[row, 4] <- temp_test$conf.int[2]
    base_stim[row, 5] <- temp_test$statistic
    base_stim[row, 6] <- temp_test$parameter
    base_stim[row, 7] <- temp_test$p.value
    base_stim[row, 8] <- stars_annotation(temp_test$p.value)
    
    names(base_stim) <- c("condition", "effect_size", 
                          "lower_95%_CI", "upper_95%_CI",
                          "t_statistic", "degrees_of_freedom",
                          "p_value", "significance")
    results_test <- base_stim
  }
  return(base_stim)
}

exclude_NA_rows <- function(data_from_list){
  data_temp <- as.data.frame(data_from_list)
  if (any(is.na(data_temp)) == TRUE){
    na_rows <- which(is.na(data_temp$colocalisation))
    exclude_rows <- c(na_rows)
    for (row in na_rows){
      exclude_rows <- c(exclude_rows, which(data_temp$ID==data_temp[row,]$ID))
    }
    exclude_rows <- unique(exclude_rows)
    return(data_temp[-exclude_rows,])
    
  }else{
    return(data_from_list)
  }
}

save_results <- function(data_list){
  results <- vector("list")
  n = 0
  for (set in data_list){
    n = n +1
    results[[names(data)[n]]] <- paired_test(exclude_NA_rows(set))
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



results <- save_results(data)
export <- results_export(results, "Results_confocal.xlsx") 





