path = "C:/Users/monar/Google Drive/Arbeit/Drube et al. statistics/COMPARE data/Formatted compare data/"
setwd(path)

wants <- c("multcomp", "ez", "ggpubr", "rstatix", "openxlsx", "stringr", "dplyr")
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
lapply(wants, require, character.only = TRUE)

file_names = list.files(pattern="*.csv", full.names = TRUE)
all_files = lapply(file_names, read.csv)
names(all_files) <- file_names 

################
# high inequality in number of Ns -> Ns are joined for high number Ns:
# join n for con and dQ
###############


### unequal n ####

# find exp with unequal n
list_join <- c()
i = 1
for (experiment in all_files){
  # print(names(all_files)[i])
  colnames(experiment)[1] = "id"
  experiment$id <- as.factor(experiment$id)
  experiment$state <- as.factor(experiment$state)
  experiment$condition <- as.factor(experiment$condition)
  experiment$BRET <- as.numeric(experiment$BRET)
  
  t <- as.numeric(
    str_remove(
      str_sub(
        summary(experiment)[1:3, 2],
        -4),
      ":"))
  
  if (!t[1] == t[2] || !t[1] == t[3]){
    list_join <- c(list_join, names(all_files)[i])
  }
  
  print(summary(experiment)[1:3, 2])
  i = i+1
}


check.integer <- function(N){
  !grepl("[^[:digit:]]", format(N,  digits = 20, scientific = FALSE))
}

# join unequal n
for (name in list_join){
  print(name)
  temp_join <- all_files[[name]]
  temp_join$condition <- as.factor(temp_join$condition)
  
  rows_na <- which(is.na(temp_join$BRET))
  if (! is.na(rows_na[1])){
    temp_join <- temp_join[-c(rows_na),]
  }
  
  if (name == "./AT1R_barr2_CON_BL_format.csv"){
    BRET_joined <- c(
      mean(c(temp_join$BRET[1:3])),
      mean(c(temp_join$BRET[4:6])),
      mean(c(temp_join$BRET[7:10])),
      mean(c(temp_join$BRET[11:12])),
      mean(c(temp_join$BRET[13:14])),
      temp_join$BRET[15:18],
      mean(c(temp_join$BRET[19:21])),
      mean(c(temp_join$BRET[22:24])),
      mean(c(temp_join$BRET[25:28])),
      mean(c(temp_join$BRET[29:30])),
      mean(c(temp_join$BRET[31:32])),
      temp_join$BRET[33:length(temp_join$BRET)]
    )
    id_joined <- c(seq(1,9), seq(1,9))
    state_joined <- c(rep("baseline", 9), rep("stimulated", 9))
    condition_joined <- 
      c(rep(c(
        rep("barr", 3),
        rep("Goe", 3),
        rep("dFLR", 3))
        , 2))
  }
  
  join_n <- summary(temp_join$condition)/6
  print(join_n)
  
  if (check.integer(join_n)[1] == TRUE){
    BRET_joined <- c()
    id_joined <- c()
    state_joined <- c()
    condition_joined <- c()
    pair <- c(1, 2, 3, 1, 2, 3)
    t = 0
    
    for (cond in names(join_n)){
      print(cond)
      sub_join <- subset(temp_join, condition == cond)
      n = join_n[cond][1]
      id <- pair + t
      t = t+3
      z = 1
      for (i in 1:6){
        id_joined <- c(id_joined, id[i])
        BRET_joined <- c(BRET_joined, mean(sub_join$BRET[z:(z+n-1)]))
        state_joined <- c(state_joined, sub_join$state[z])
        condition_joined <- c(condition_joined, as.character(sub_join$condition[z]))
        
        #print(temp_join$BRET[z])
        #print(temp_join$BRET[z+n-1])
        #print(BRET_joined)
        z = z+n
      }
    }
  }
  result_join <- data.frame(id_joined, 
                            condition_joined, 
                            state_joined, 
                            BRET_joined)
  
  names(result_join) <- names(temp_join)
  print(result_join)
  
  title <- str_remove(name, "_format.csv")
  
  all_files[[paste(title, "JOINED", sep ="_")]] <- result_join
  file_names <- c(file_names, paste(title, "JOINED", sep ="_"))
}


##### analysis ####
index = 0
results_paired <- vector("list", length(all_files))
results_baseline <- vector("list")
results_stimulated <- vector("list")

results_paired_detail <- vector("list")
results_baseline_detail <- vector("list")
results_stimulated_detail <- vector("list")


for (experiment in all_files){
  index = index +1
  print(file_names[index])
  
  # formatting
  colnames(experiment)[1] = "id"
  experiment$id <- as.factor(experiment$id)
  experiment$state <- as.factor(experiment$state)
  experiment$condition <- as.factor(experiment$condition)
  experiment$BRET <- as.numeric(experiment$BRET)
  # print(summary(experiment))
  
  # check for outliers
  temp_outlier <- experiment %>%
    group_by(state,condition) %>%
    identify_outliers(BRET)
  
  if (!(length(temp_outlier$is.outlier) == 0)){
    print(paste("OUTLIER identified in ", file_names[index]))
    print(temp_outlier)
  }
  
  # shapiro test for each combination of factor levels 
  #(can only be done, for n = 3)
  if (all(summary(experiment$condition)>5)){
    temp_shapiro<-  experiment %>%
      group_by(state, condition) %>%
      shapiro_test(BRET)
  }
  print(temp_shapiro)
  
  # remove rows with NA
  rows_w_na <- which(is.na(experiment$BRET))
  if (! is.na(rows_w_na[1])){
      experiment <- experiment[-c(rows_w_na),]
  }

  
  # mixed model ANOVA
  temp_mixed_aov<-ezANOVA(data = experiment, 
                          dv = .(BRET), 
                          wid = .(id), 
                          within = .(state), 
                          between = .(condition), 
                          return_aov = TRUE)
  temp_p_value<-temp_mixed_aov$ANOVA[3,5]
    
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
  
  # if interaction is significant, look into simple effects 
  # A)split on the larger variable (condition) -> compare baseline vs. stimulated (within subjects, PAIRED!)
  ### would be paired pairwise.t.tests if more than two levels would exist (here only baseline and stimulated)
  levels <- levels(experiment$condition)
  base_stim <- as.data.frame(levels)
  base_stim_detail <- as.data.frame(levels)
  row = 0
  for (i in levels){
    temp = subset(experiment, condition == as.character(i))
    row = row +1
    temp_test <- t.test(temp$BRET~temp$state, paired = TRUE)
    
    base_stim[row, 2] <-temp_test$p.value
    names(base_stim)[2] <- "p_value"
    
        results_paired[[index]]<- base_stim
    names(results_paired)[index] <- names(all_files)[index]
    
    # also save details
    base_stim_detail[row, 2] <-temp_test$parameter
    names(base_stim_detail)[2] <- "df"
    
    # for effect size, calculate pooled sd (difference of means/pooled sd = effect size)
    temp_posd <- (sd(temp$BRET[1:3]) + sd(temp$BRET[4:6]))/2
    base_stim_detail[row, 3] <- temp_test$estimate/temp_posd
    names(base_stim_detail)[3] <- "effect_size"
    
    base_stim_detail[row, 4] <- temp_test$conf.int[1]
    base_stim_detail[row, 5] <- temp_test$conf.int[2]
    names(base_stim_detail)[4] <- "lower"
    names(base_stim_detail)[5] <- "upper"
    
    base_stim_detail[row, 6] <- temp_test$statistic
    names(base_stim_detail)[6] <- "tstat"
    
    base_stim_detail[row, 7] <- temp_test$p.value
    names(base_stim_detail)[7] <- "p.value"

    results_paired_detail[[index]]<- base_stim_detail
    names(results_paired_detail)[index] <- names(all_files)[index]
    
  }
  
  # B)split on within variable to compare different baselines (between subjects are compared)
  # only baselines are interesting for now
  states <- levels(experiment$state)
  for (i in states){
    temp = subset(experiment, state == as.character(i))
    
    # check outliers
    temp_outlier2 <- temp %>%
      group_by(state,condition) %>%
      identify_outliers(BRET)
    
    if (!(length(temp_outlier$is.outlier) == 0)){
       print(paste("outlier identified in baselines of ", file_names[index]))
       print(temp_outlier2)
       print(summary(temp$condition))
     }
    
    # check distribution
     if (shapiro_test(temp$BRET)$p.value < 0.05){
       print(paste("Note that ",
                   as.character(i),
                   " values for ",
                   file_names[index],
                   " are not normally distributed with a p value of ",
                   shapiro_test(temp$BRET)$p.value))
     }
    
    lev_test <- levene_test(BRET ~ condition, data = temp)
    
    if (lev_test$p < 0.05)
    {print("LEVENE: Variance of baselines is not homogenous among")
      print(file_names[index])
      print(lev_test$p)
    }
    
    temp_aov<-(aov(BRET ~ condition, temp))
    
    temp_post_hoc <- 
      summary(
        glht(temp_aov,
             linfct = mcp(condition = "Tukey")))
    
    contrast <-names(temp_post_hoc$test$coefficients)
    post_hoc <- as.data.frame(contrast)
    post_hoc$p_value <- temp_post_hoc$test$pvalues

    # also save details of post hoc
    post_hoc_detail <- as.data.frame(contrast)

    # effect size (in R here estimate)
    post_hoc_detail$effect_size <- temp_post_hoc$test$coefficients
    
    # also save conf. int of mean difference 
    conf_int <- confint(glht(temp_aov,
                             linfct = mcp(condition = "Tukey")))
    post_hoc_detail$lower <- conf_int$confin[,2]
    post_hoc_detail$upper <- conf_int$confin[,3]
    
    # test statistic
    post_hoc_detail$tstat <- temp_post_hoc$test$tstat
    # degrees of freedom
    post_hoc_detail$df <- temp_post_hoc$df
    
    post_hoc_detail$pvalue <-  temp_post_hoc$test$pvalues
    
    if (i == "baseline"){
      results_baseline[[index]]<- post_hoc
      names(results_baseline)[index] <- names(all_files)[index]
      
      results_baseline_detail[[index]]<- post_hoc_detail
      names(results_baseline_detail)[index] <- names(all_files)[index]
    }
    if (i == "stimulated"){
      results_stimulated[[index]]<- post_hoc
      names(results_stimulated)[index] <- names(all_files)[index]
      
      results_stimulated_detail[[index]]<- post_hoc_detail
      names(results_stimulated_detail)[index] <- names(all_files)[index]
    }
  }
  
}


#### save p values ####

experiment_names <- file_names %>%
  str_remove(".csv") %>%
  str_remove("./") %>%
  str_remove("_format") %>%
  str_remove("_BL") %>%
  str_remove("_baseline")

cond <- experiment_names %>%
  str_remove("_barr1") %>%
  str_remove("_barr2") %>%
  str_remove("AT1R_")

#create arrestin column
barr <- str_extract(experiment_names, "barr1")
#find out which arrestin in which row
for (row in 1:length(barr)){
  if (is.na(barr[row])){
      barr[row]<- "barr2"
   }
}


####p value overview paired comparison####
p_value_paired <- data.frame(experiment_names)
p_value_paired$condition<-cond
p_value_paired$barrestin<-barr

p_value_paired_detail <- p_value_paired

#-> collect p values from results_paired in matrix "p_values_contrast"
p_values_contrast <- matrix(nrow = length(experiment_names), 
                            ncol = length(levels)*2)
colnames <- c()
for (condition in levels){
  colnames <- c(colnames, condition, paste(condition, "_sig"))
}
colnames(p_values_contrast) <- colnames

#[row, column]
row = 0
for (experiment in results_paired){
  row = row+1
  col = 1
  for (p in experiment$p_value){
    p_values_contrast[row, col] <- p
    
    if (p > 0.05){sig <- c(" ")}
    if (between(p, 0.01, 0.05)){sig <- c("*")}
    if (between(p, 0.001, 0.01)){sig <- c("**")}
    if (between(p, 0.0001, 0.001)){sig <- c("***")}
    if (p < 0.0001){sig <- c("****")}
    
    p_values_contrast[row, (col+1)] <- sig
    
    col = col +2
  }
}
p_values_contrast <- as.data.frame(p_values_contrast)
start = length(p_value_paired)
for (c in 1:length(p_values_contrast)){
  p_value_paired[, (start + c)] <- p_values_contrast[,c]
}
names(p_value_paired)[4:9] <- names(p_values_contrast)


# DETAILS
#-> collect p values from results_paired_detail in matrix "paired_details"
paired_details <- matrix(nrow = length(experiment_names), 
                            ncol = length(levels)*7)
colnames <- c()
for (con in levels){
  colnames <- c(colnames, 
                paste(con, "_df"),
                paste(con, "_effect_size"),
                paste(con, "_lower"),
                paste(con, "_upper"),
                paste(con, "_tstat"),
                paste(con, "_pvalue"),
                paste(con, "_sig"))
}
colnames(paired_details) <- colnames

#[row, column]
row = 0
for (experiment in results_paired_detail){
  row = row+1
  col = 1
  for (df in experiment$df){
    paired_details[row, col] <- df
    col = col +7
  }
  
  col = 2
  for (es in experiment$effect_size){
    paired_details[row, col] <- es
    col = col +7
  }
  
  col = 3
  for (lower in experiment$lower){
    paired_details[row, col] <- lower
    col = col +7
  }
  col = 4
  for (upper in experiment$upper){
    paired_details[row, col] <- upper
    col = col +7
  }
  col = 5
  for (tstat in experiment$tstat){
    paired_details[row, col] <- tstat
    col = col +7
  }
  col = col + 3
}

col = 1
n = 4
for (con in contrast){
  paired_details[,col+5] <- p_value_paired[,n]
  paired_details[,col+6] <- p_value_paired[,n +1]
  
  col = col + 7
  n = n +2
}

paired_details_df <- as.data.frame(p_value_paired$experiment_names) 
names(paired_details_df) <- "experiment_names"

paired_details_df$condition <- p_value_paired$condition
paired_details_df$barr<- barr

s = length(paired_details_df)
paired_details <- as.data.frame(paired_details)

for (c in 1:length(paired_details)){
  paired_details_df[, (s + c)] <- paired_details[,c]
}
names(paired_details_df)[4:length(paired_details_df)] <- names(paired_details)



####p value overview only baseline and only stim####

p_value_baseline <- data.frame(p_value_paired[1:3])
p_value_stimulated <- data.frame(p_value_paired[1:3])

p_value_baseline_details <- data.frame(p_value_paired[1:3])
p_value_stimulated_details <- data.frame(p_value_paired[1:3])

#-> collect p values from results_baseline  in matrix "p_values_base" 
p_values_base <- matrix(nrow = length(experiment_names), 
                            ncol = length(levels)*2)
colnames <- c()
for (condition in contrast){
  colnames <- c(colnames, condition, paste(condition, "_sig"))
}
colnames(p_values_base) <- colnames

#[row, column]
row = 0
for (experiment in results_baseline){
  row = row+1
  col = 1
  for (p in experiment$p_value){
    p_values_base[row, col] <- p
    
    if (p > 0.05){sig <- c(" ")}
    if (between(p, 0.01, 0.05)){sig <- c("*")}
    if (between(p, 0.001, 0.01)){sig <- c("**")}
    if (between(p, 0.0001, 0.001)){sig <- c("***")}
    if (p < 0.0001){sig <- c("****")}
    p_values_base[row, (col+1)] <- sig
    
    col = col +2
  }
}
p_values_base <- as.data.frame(p_values_base)
start = length(p_value_baseline)
for (c in 1:length(p_values_base)){
  p_value_baseline[, (start + c)] <- p_values_base[,c]
}
names(p_value_baseline)[4:9] <- names(p_values_base)


#-> collect p values from results_stimulated  in matrix "p_values_stim" 
p_values_stim <- matrix(nrow = length(experiment_names), 
                        ncol = length(levels)*2)
colnames <- c()
for (condition in contrast){
  colnames <- c(colnames, condition, paste(condition, "_sig"))
}
colnames(p_values_stim) <- colnames

#[row, column]
row = 0
for (experiment in results_stimulated){
  row = row+1
  col = 1
  for (p in experiment$p_value){
    p_values_stim[row, col] <- p
    
    if (p > 0.05){sig <- c(" ")}
    if (between(p, 0.01, 0.05)){sig <- c("*")}
    if (between(p, 0.001, 0.01)){sig <- c("**")}
    if (between(p, 0.0001, 0.001)){sig <- c("***")}
    if (p < 0.0001){sig <- c("****")}
    p_values_stim[row, (col+1)] <- sig
    
    col = col +2
  }
}
p_values_stim <- as.data.frame(p_values_stim)
start = length(p_value_stimulated)
for (c in 1:length(p_values_stim)){
  p_value_stimulated[, (start + c)] <- p_values_stim[,c]
}
names(p_value_stimulated)[4:9] <- names(p_values_stim)


#-> collect p values from results_baseline_detail  in matrix "base_details" 

base_details <- matrix(nrow = length(experiment_names), 
                        ncol = length(levels)*7)
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
for (experiment in results_baseline_detail){
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
n = 4
for (con in contrast){
  base_details[,col+5] <- p_value_baseline[,n]
  base_details[,col+6] <- p_value_baseline[,n +1]
  
  col = col + 7
  n = n +2
}

base_details_df <- as.data.frame(p_value_baseline$experiment_names) 
names(base_details_df) <- "experiment_names"

base_details_df$condition <- p_value_baseline$condition
base_details_df$barr<- barr

s = length(base_details_df)
base_details <- as.data.frame(base_details)

for (c in 1:length(paired_details)){
  base_details_df[, (s + c)] <- base_details[,c]
}
names(base_details_df)[4:length(base_details_df)] <- names(base_details)



#-> collect p values from results_stimulated_detail  in matrix "stim_details" 

stim_details <- matrix(nrow = length(experiment_names), 
                       ncol = length(levels)*7)
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
colnames(stim_details) <- colnames

#[row, column]
row = 0
for (experiment in results_stimulated_detail){
  row = row+1
  col = 1
  for (df in experiment$df){
    stim_details[row, col] <- df
    col = col +7
  }
  
  col = 2
  for (es in experiment$effect_size){
    stim_details[row, col] <- es
    col = col +7
  }
  
  col = 3
  for (lower in experiment$lower){
    stim_details[row, col] <- lower
    col = col +7
  }
  col = 4
  for (upper in experiment$upper){
    stim_details[row, col] <- upper
    col = col +7
  }
  col = 5
  for (tstat in experiment$tstat){
    stim_details[row, col] <- tstat
    col = col +7
  }
  col = col + 3
}

col = 1
n = 4
for (con in contrast){
  stim_details[,col+5] <- p_value_stimulated[,n]
  stim_details[,col+6] <- p_value_stimulated[,n +1]
  
  col = col + 7
  n = n +2
}

stim_details_df <- as.data.frame(p_value_stimulated$experiment_names) 
names(stim_details_df) <- "experiment_names"

stim_details_df$condition <- p_value_stimulated$condition
stim_details_df$barr<- barr

s = length(stim_details_df)
stim_details <- as.data.frame(stim_details)

for (c in 1:length(paired_details)){
  stim_details_df[, (s + c)] <- stim_details[,c]
}
names(stim_details_df)[4:length(stim_details_df)] <- names(stim_details)




#### save all results ####

#save p value table
remove(results_wb)
results_wb<-createWorkbook()
addWorksheet(results_wb, "p_values_paired")
writeData(results_wb, "p_values_paired", p_value_paired, 
          startRow = 1, 
          startCol = 1)
saveWorkbook(results_wb, 
             file = "P_values_paired_overview.xlsx", 
             overwrite = TRUE)

remove(results_wb)
results_wb<-createWorkbook()
addWorksheet(results_wb, "p_values_base")
writeData(results_wb, "p_values_base", p_value_baseline, 
          startRow = 1, 
          startCol = 1)
saveWorkbook(results_wb, 
             file = "P_values_baseline_overview.xlsx", 
             overwrite = TRUE)

remove(results_wb)
results_wb<-createWorkbook()
addWorksheet(results_wb, "p_values_stim")
writeData(results_wb, "p_values_stim", p_value_stimulated, 
          startRow = 1, 
          startCol = 1)
saveWorkbook(results_wb, 
             file = "P_values_stimulated_overview.xlsx", 
             overwrite = TRUE)

# details
remove(results_wb)
results_wb<-createWorkbook()
addWorksheet(results_wb, "details")
writeData(results_wb, "details", paired_details_df, 
          startRow = 1, 
          startCol = 1)
saveWorkbook(results_wb, 
             file = "paired_details_overview.xlsx", 
             overwrite = TRUE)

remove(results_wb)
results_wb<-createWorkbook()
addWorksheet(results_wb, "details")
writeData(results_wb, "details", base_details_df, 
          startRow = 1, 
          startCol = 1)
saveWorkbook(results_wb, 
             file = "baseline_details_overview.xlsx", 
             overwrite = TRUE)

remove(results_wb)
results_wb<-createWorkbook()
addWorksheet(results_wb, "details")
writeData(results_wb, "details", stim_details_df, 
          startRow = 1, 
          startCol = 1)
saveWorkbook(results_wb, 
             file = "stimulated_details_overview.xlsx", 
             overwrite = TRUE)

