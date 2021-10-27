path = "C:/path/to/folder/"
setwd(path)

wants <- c("multcomp", "ez", "ggpubr", "rstatix", "openxlsx", "stringr", "dplyr")
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
lapply(wants, require, character.only = TRUE)

file_names = list.files(pattern="*.csv", full.names = TRUE)
all_files = lapply(file_names, read.csv)

#each csv as own object
#temp = list.files(pattern="*.csv")
#for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))

# matrix positions CON, dQ, GRK2, GRK3, GRK5, GRK6 
#(factors are sorted alphabetically)
contrast_matrix <- rbind("dQ -GRK2" = c(0, -1, 1, 0, 0, 0),
                         "dQ -GRK3" = c(0, -1, 0, 1, 0, 0),
                         "dQ -GRK5" = c(0, -1, 0, 0, 1, 0),
                         "dQ -GRK6" = c(0, -1, 0, 0, 0, 1),
                         "dQ -CON" = c(1, -1, 0, 0, 0, 0))



##### analysis ####
index = 0
results_paired <- vector("list", length(all_files))
results_baseline <- vector("list")
results_base_detail <- vector("list")

for (experiment in all_files){
  index = index +1
  print(file_names[index])
  
  # formatting
  colnames(experiment)[1] = "id"
  experiment$id <- as.factor(experiment$id)
  experiment$state <- as.factor(experiment$state)
  experiment$condition <- as.factor(experiment$condition)
  experiment$BRET <- as.numeric(experiment$BRET)
  #print(summary(experiment))

  # check for outliers
  temp_outlier <- experiment %>%
    group_by(state,condition) %>%
    identify_outliers(BRET)
  
  if (!(length(temp_outlier$is.outlier) == 0)){
    print(paste("OUTLIER identified in ", file_names[index]))
    print(temp_outlier)
  }
  
  # shapiro test for each combination of factor levels 
  if (all(summary(experiment$condition)>5)){
    temp_shapiro<-  experiment %>%
      group_by(state, condition) %>%
      shapiro_test(BRET)
  }
    
  # mixed model ANOVA
  temp_mixed_aov<-ezANOVA(data = experiment, 
                          dv = .(BRET), 
                          wid = .(id), 
                          within = .(state), 
                          between = .(condition), 
                          return_aov = TRUE)
  temp_p_value<-temp_mixed_aov$ANOVA[3,5]
  # print(c(temp_mixed_aov$ANOVA[1], temp_mixed_aov$ANOVA[5]))
  
  #give warning if p is not significant for interaction effect
  if (temp_p_value > 0.05){
    print(paste("Note that interaction of condition and state was not significant for ",
                file_names[index],
                " with a p value of ",
                temp_p_value))
    if (temp_mixed_aov$ANOVA[1,5] < 0.05){
      print(paste("the main effect of condition is significant with a p value of",
                  temp_mixed_aov$ANOVA[1,5]))
    }
    if (temp_mixed_aov$ANOVA[2,5] < 0.05){
      print(paste("the main effect of state is significant with a p value of",
                  temp_mixed_aov$ANOVA[2,5]))
    }
  }
  
  # look into simple effects 
  # A)split on the larger variable (condition) -> compare baseline vs. stimulated (within subjects, PAIRED!)
  ### would be paired pairwise.t.tests if more than two levels would exist (here only baseline and stimulated)
  levels <- levels(experiment$condition)
  base_stim <- as.data.frame(levels)
  row = 0
  for (i in levels){
    temp = subset(experiment, condition == as.character(i))
    row = row +1
    base_stim[row, 2] <-t.test(temp$BRET~temp$state, paired = TRUE)$p.value
    names(base_stim)[2] <- "p_value"
    #print(base_stim)
    results_paired[[index]]<- base_stim
  }
  #print(results_paired[index])

  # B)split on within variable to compare different baselines (between subjects are compared)
  # only baselines are interesting for now
  temp = subset(experiment, state == "baseline")
  
  #check outliers
  temp_outlier2 <- temp %>%
    group_by(state,condition) %>%
    identify_outliers(BRET)
  
  baseline_aov<-(aov(BRET ~ condition, temp))
  
  # always test (even if ANOVA is not significant)
  if (summary(baseline_aov)[[1]]$`Pr(>F)`[1] < 2){ 
    # dunnett + greater
    base_only <- summary(
      glht(baseline_aov,
           linfct = mcp(condition = contrast_matrix),
           alternative = "greater"))

    contrast<-names(base_only$test$coefficients)
    post_hoc <- as.data.frame(contrast)
    post_hoc$p_value <- base_only$test$pvalues
    
    # also save details on post hoc tests
    post_hoc_detail <- as.data.frame(contrast)
    
    # effect size (in R here estimate)
    post_hoc_detail$effect_size <- base_only$test$coefficients
    
    # also save conf. int of mean difference 
    conf_int <- confint(glht(baseline_aov,
                             linfct = mcp(condition = contrast_matrix),
                             alternative = "greater"))
    post_hoc_detail$lower <- conf_int$confin[,2]
    post_hoc_detail$upper <- conf_int$confin[,3]
    
    # test statistic
    post_hoc_detail$tstat <- base_only$test$tstat
    # degrees of freedom
    post_hoc_detail$df <- base_only$df
    # p values
    post_hoc_detail$pvalue <-  base_only$test$pvalues
    
    results_base_detail[[index]] <- post_hoc_detail
    
    # if p value <2 save in post hoc --> ALWAYS save post hoc
    if (any(post_hoc$p_value<2) == TRUE){
      results_baseline[[index]]<- post_hoc
    } else {
    print("No difference between baselines was discovered, with a p value of")
    print(summary(baseline_aov)[[1]]$`Pr(>F)`[1])
  }
  }
}
  
names(results_baseline)<-file_names
names(results_paired)<-file_names
names(results_base_detail)<-file_names

# save as xlsx
# base vs. stim
remove(results_wb)
results_wb<-createWorkbook()
for (file in file_names){
  experiment = substr(file, 3, nchar(file)-11)
  addWorksheet(results_wb, experiment)
  writeData(results_wb, experiment, results_paired[[file]], startRow = 1, startCol = 1)
}

saveWorkbook(results_wb, file = "Base_vs_Stim_results.xlsx", overwrite = TRUE)

# baseline comparison
remove(results_wb)
results_wb<-createWorkbook()
for (file in file_names){
  if (is.null(results_baseline[[file]]) == FALSE){
    experiment = substr(file, 3, nchar(file)-11)
    addWorksheet(results_wb, experiment)
    writeData(results_wb, experiment, results_baseline[[file]], startRow = 1, startCol = 1)
  }
}

saveWorkbook(results_wb, file = "Baseline_comparison_results.xlsx", overwrite = TRUE)


# nice-looking export
#### p value overview baseline comparison ####

experiment_names <- file_names %>%
  str_remove(".csv") %>%
  str_remove("./") %>%
  str_remove("_format") %>%
  str_remove("_BL") %>%
  str_remove("_baseline")

gpcr <- experiment_names %>%
  str_remove("_barr1") %>%
  str_remove("_barr2") %>%
  str_remove("_bA1") %>%
  str_remove("_bA2") 

p_value_BL_compare <- data.frame(experiment_names)
p_value_BL_compare$GPCR <- gpcr

#create arrestin column
barr <- str_extract(experiment_names, "barr1")
bA1 <- str_extract(experiment_names, "bA1")
#find out which arrestin in which row
for (row in 1:length(barr)){
  if (is.na(barr[row])){
    if (!is.na(bA1[row])){
      barr[row]<- "barr1"
    }else{
      barr[row]<- "barr2"
    }
  }
}

p_value_BL_compare$barrestin<-barr

#-> collect p values from results_baseline in matrix "p_values_contrast"
p_values_contrast <- matrix(nrow = length(experiment_names), 
                         ncol = length(contrast_matrix[,1])*2)
colnames <- c()
for (con in contrast){
  colnames <- c(colnames, con, paste(con, "_sig"))
}
colnames(p_values_contrast) <- colnames

#[row, column]
row = 0
for (experiment in results_baseline){
  row = row+1
  col = 1
  for (p in experiment$p_value){
    p_values_contrast[row, col] <- p
    if (p > 0.05){sig <- c(" ")}
    if (between(p, 0.01, 0.05)){sig <- c("*")}
    if (between(p, 0.001, 0.01)){sig <- c("**")}
    if (p < 0.001){sig <- c("***")}
    p_values_contrast[row, (col+1)] <- sig
    col = col +2
  }
}

p_values_contrast <- as.data.frame(p_values_contrast)
start = length(p_value_BL_compare)
for (c in 1:length(p_values_contrast)){
  p_value_BL_compare[, (start + c)] <- p_values_contrast[,c]
}
names(p_value_BL_compare)[4:13] <- names(p_values_contrast )

#save p value table
remove(results_wb)
results_wb<-createWorkbook()
addWorksheet(results_wb, "p_values")
writeData(results_wb, "p_values", p_value_BL_compare, 
          startRow = 1, 
          startCol = 1)
saveWorkbook(results_wb, 
             file = "P_values_BL_Compare_overview.xlsx", 
             overwrite = TRUE)


#### p value overview baseline comparison DETAILS ####
details_baseline_overview<- data.frame(experiment_names)
details_baseline_overview$GPCR<-gpcr
details_baseline_overview$barrestin<-barr


#-> collect details from results_base_detail in matrix "base_details"
base_details <- matrix(nrow = length(experiment_names), 
                       ncol = length(contrast_matrix[,1])*7)
colnames <- c()
for (con in contrast){
  colnames <- c(colnames, 
                paste(con, "_df"),
                paste(con, "_effect_size"),
                paste(con, "_lower"),
                paste(con, "_upper"),
                paste(con, "_tstat"),
                paste(con, "_pvalue"),
                paste(con, "_sig"))
}
colnames(base_details) <- colnames

#[row, column]
row = 0
for (experiment in results_base_detail){
  row = row+1
  col = 1
  for (df in experiment$df){
    base_details[row, col] <- df
    col = col +7
  }
  
  col = 2
  for (es in experiment$effect_size){
    base_details[row, col] <- es
    col = col +7
  }
  
  col = 3
  for (lower in experiment$lower){
    base_details[row, col] <- lower
    col = col +7
  }
  col = 4
  for (upper in experiment$upper){
    base_details[row, col] <- upper
    col = col +7
  }
  col = 5
  for (tstat in experiment$tstat){
    base_details[row, col] <- tstat
    col = col +7
  }
  col = col + 3
}

col = 1
n = 1
for (con in contrast){
  base_details[,col+5] <- p_values_contrast[,n]
  base_details[,col+6] <- p_values_contrast[,n +1]
  
  col = col + 7
  n = n +2
}


base_details_df <- as.data.frame(p_value_BL_compare$experiment_names) 
names(base_details_df) <- "experiment_names"

base_details_df$GPCR <- gpcr
base_details_df$barr<- barr

s = length(base_details_df)
base_details <- as.data.frame(base_details)

for (c in 1:length(base_details)){
  base_details_df[, (s + c)] <- base_details[,c]
}
names(base_details_df)[4:length(base_details_df)] <- names(base_details)

# exclude experiments not shown here
exclude <- c("V2R_dFLR_bA1",
             "V2R_dFLR_bA2")


# get row numbers
exclude_rows <- c()
for (exp in 1:length(exclude)){
  temp <- which(base_details_df$experiment_names == exclude[exp])
  exclude_rows <- c(exclude_rows, temp)
}

base_details_df_select<-base_details_df[-exclude_rows,]


#save details table
remove(results_wb)
results_wb<-createWorkbook()
addWorksheet(results_wb, "details")
writeData(results_wb, "details", base_details_df_select, 
          startRow = 1, 
          startCol = 1)
saveWorkbook(results_wb, 
             file = "baseline_detail_statistics.xlsx", 
             overwrite = TRUE)
