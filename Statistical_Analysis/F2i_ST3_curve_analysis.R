##### data preparation ####
path = "C:/path/to/folder/"
setwd(path)

wants <- c("multcomp", "rstatix", "ggpubr", 
           "openxlsx", "rlist", "stringr", 
           "dplyr")
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
lapply(wants, require, character.only = TRUE)

# import data
file_names = list.files(pattern="*.csv", full.names = TRUE)
all_files = lapply(file_names, read.csv)

experiment_names <- c()
for (file in file_names){
  experiment_names <- c(experiment_names, substr(file, 3, nchar(file)-11) )
}


# matrix positions CON, dQ, GRK2, GRK3, GRK5, GRK6 
# (factors are sorted alphabetically)
contrast_matrix <- rbind("dQ -GRK2" = c(0, -1, 1, 0, 0, 0),
                         "dQ -GRK3" = c(0, -1, 0, 1, 0, 0),
                         "dQ -GRK5" = c(0, -1, 0, 0, 1, 0),
                         "dQ -GRK6" = c(0, -1, 0, 0, 0, 1),
                         "dQ -CON" = c(1, -1, 0, 0, 0, 0),
                         "CON -GRK2" = c(1, 0, 1, 0, 0, 0),
                         "CON -GRK3" = c(1, 0, 0, 1, 0, 0),
                         "CON -GRK5" = c(1, 0, 0, 0, 1, 0),
                         "CON -GRK6" = c(1, 0, 0, 0, 0, 1))


##### analysis ####

index = 0
results_curve <- vector("list")
results_detail <- vector("list")

for (experiment in all_files)
{experiment$condition <- as.factor(experiment$condition)
  index = index +1
  print(experiment_names[index])
  #print(summary(experiment))
  
    #check for outliers
  temp_outlier <- experiment %>%
    group_by(condition) %>%
    identify_outliers(fold_change)
  
  #give warning if there is an outlier, remove it for n > 3
  if (!(length(temp_outlier$is.outlier) == 0))
  {print(paste("outlier identified in", experiment_names[index]))
    print(temp_outlier)
    n = length(experiment$fold_change)
    if (n>3){
      outlier<-temp_outlier$fold_change
      row <- c()
      for (i in 1:length(outlier)){
        row[i] = which(experiment$fold_change == outlier[i])
      }
      experiment_temp<-experiment
      experiment<-experiment_temp[-row,]
      print("outliers were removed in rows")
      print(row)
    }}
  

  # anova
  temp_aov <- aov(fold_change ~ condition, data = experiment)
  temp_p_value <- summary(temp_aov)[[1]]$`Pr(>F)`[1]
  
  # give warning if p is not significant
  if (temp_p_value > 0.05)
  {print(paste("Note that ANOVA was not significant for", experiment_names[index]))
    print(paste("with a p value of", temp_p_value))
    }
    
  #check distribution
  temp_shapiro<-shapiro_test(residuals(temp_aov))
  if (temp_shapiro$p.value < 0.05)
  {print("Note that normality assumption is violated, since shapiro gave a p value of")
    print(temp_shapiro$p.value)}

  #check homogeneity of variance assumption, if it is violated, perform Welch ANOVA
  levene_test(fold_change ~ condition, data = experiment)
  
  if (levene_test(fold_change ~ condition, data = experiment)[4] < 0.05)
  {print("Variance is not homogenous among")
    print(experiment_names[index])}
  
  mult_comp <- summary(glht(temp_aov, 
                            linfct = mcp(condition = contrast_matrix)), 
                            test = adjusted("bonferroni"))
  
  contrast <- names(mult_comp$test$coefficients)
  post_hoc <- as.data.frame(contrast)
  post_hoc$p_value <- mult_comp$test$pvalues 
  
  
  results_curve[[experiment_names[index]]]<- post_hoc
  
  
  # save also unadjusted p values
  mult_comp_noadj <- summary(glht(temp_aov, 
                            linfct = mcp(condition = contrast_matrix)), 
                       test = adjusted("none"))
  
  contrast<-names(mult_comp_noadj$test$coefficients)
  post_hoc <- as.data.frame(contrast)
  post_hoc$p_value <- mult_comp_noadj$test$pvalues 
  
  name <- paste(experiment_names[index], "noadj", sep = "_")
  
  results_curve[[name]]<- post_hoc
  
  
  #save more details of testing
  post_hoc_detail <- as.data.frame(contrast)
  # effect size
  post_hoc_detail$effect_size <- mult_comp$test$coefficients
  # also save conf. int of mean difference 
  conf_int <- confint(glht(temp_aov, 
                           linfct = mcp(condition = contrast_matrix)))
  post_hoc_detail$lower <- conf_int$confin[,2]
  post_hoc_detail$upper <- conf_int$confin[,3]
  # test statistic
  post_hoc_detail$tstat <- mult_comp$test$tstat
  # degrees of freedom
  post_hoc_detail$df <- mult_comp$df
  post_hoc_detail$pvalue <-  mult_comp_noadj$test$pvalues
  post_hoc_detail$pvalue_adjusted <- mult_comp$test$pvalues
  results_detail[[experiment_names[index]]]<- post_hoc_detail
  }


##### save results ####

# save as xlsx
remove(results_wb)
results_wb<-createWorkbook()
for (experiment in list.names(results_curve)){
  file = paste(experiment, ".xlsx", sep = "")
  addWorksheet(results_wb, file)
  writeData(results_wb, file, results_curve[[experiment]], 
            startRow = 1, 
            startCol = 1)
}

saveWorkbook(results_wb, file = "Curve_results.xlsx", overwrite = TRUE)

# save for visualizations
results_curve_table<- data.frame()
results_curve_table_noadj <- data.frame()
exp_adj<- c()
exp_noadj<- c()

noadjusted_list <- str_extract(list.names(results_curve), "noadj")

for (i in 1:length(results_curve)){
  if (is.na(noadjusted_list[i])){
    results_curve_table<-rbind(results_curve_table, results_curve[[i]])
    exp_adj <- c(exp_adj, rep(list.names(results_curve)[i], 9))
  }else{
    results_curve_table_noadj<-
      rbind(results_curve_table_noadj, results_curve[[i]])
    exp_noadj <- c(exp_noadj, rep(list.names(results_curve)[i], 9))
  }
}

results_curve_table$experiment <- as.factor(exp_adj)
results_curve_table$contrast <- as.factor(results_curve_table$contrast)

results_curve_table_noadj$experiment <- as.factor(exp_noadj)
results_curve_table_noadj$contrast <- as.factor(results_curve_table$contrast)


# save adjusted for visualization
remove(results_wb)
results_wb<-createWorkbook()
addWorksheet(results_wb, "results")
writeData(results_wb, "results", results_curve_table, 
          startRow = 1, 
          startCol = 1)
saveWorkbook(results_wb, 
             file = "Curve_results_overview.xlsx", 
             overwrite = TRUE)

# save non - adjusted for visualization
remove(results_wb)
results_wb<-createWorkbook()
addWorksheet(results_wb, "results")
writeData(results_wb, "results", results_curve_table_noadj, 
          startRow = 1, 
          startCol = 1)
saveWorkbook(results_wb, 
             file = "Curve_results_overview_noadj.xlsx", 
             overwrite = TRUE)




##### check violations ####
results_curve_corrected <- results_curve

# normality violation (without outliers)
# M1 bA1
M1_A1 <- read.csv("M1R_barr1_format.csv")
M1_A1$condition <- as.factor(M1_A1$condition)
aov_M1_A1<- aov(fold_change ~ condition, data = M1_A1)
ggqqplot(residuals(aov_M1_A1))
# -> remove last value for CON

M1_A1_subset <- M1_A1[1:23, ]
aov_M1_A1<- aov(fold_change ~ condition, data = M1_A1_subset)
ggqqplot(residuals(aov_M1_A1))
summary(aov_M1_A1)
mult_comp_M1_correct <- summary(glht(aov_M1_A1, 
                          linfct = mcp(condition = contrast_matrix)), 
                     test = adjusted("bonferroni"))

# save this in results
contrast<-names(mult_comp_M1_correct$test$coefficients)
post_hoc <- as.data.frame(contrast)
post_hoc$p_value <- mult_comp_M1_correct$test$pvalues

results_curve_corrected$M1R_barr1<-post_hoc

#also save unadjusted
mult_comp_M1_correct_noadj <- summary(glht(aov_M1_A1, 
                                     linfct = mcp(condition = contrast_matrix)), 
                                    test = adjusted("none"))
post_hoc$p_value <- mult_comp_M1_correct_noadj$test$pvalues
results_curve_corrected$M1R_barr1_noadj <-post_hoc

# also save detail
post_hoc_detail_M1_barr1 <- post_hoc_detail
post_hoc_detail_M1_barr1$effect_size <- mult_comp_M1_correct_noadj$test$coefficients

M1_barr1_conf <- confint(glht(aov_M1_A1, 
                              linfct = mcp(condition = contrast_matrix)))
post_hoc_detail_M1_barr1$lower <- M1_barr1_conf$confin[,2]
post_hoc_detail_M1_barr1$upper <- M1_barr1_conf$confin[,3]

post_hoc_detail_M1_barr1$tstat <- mult_comp_M1_correct_noadj$test$tstat
post_hoc_detail_M1_barr1$df <- mult_comp_M1_correct_noadj$df

post_hoc_detail_M1_barr1$pvalue <-   mult_comp_M1_correct_noadj$test$pvalues
post_hoc_detail_M1_barr1$pvalue_adjusted <- mult_comp_M1_correct$test$pvalues

results_detail$M1R_barr1 <- post_hoc_detail_M1_barr1


# M3_bA2
M3_A2 <- read.csv("M3R_barr2_format.csv")
M3_A2$condition <- as.factor(M3_A2$condition)
aov_M3_A2<- aov(fold_change ~ condition, data = M3_A2)
ggqqplot(residuals(aov_M3_A2))
boxplot(fold_change ~ condition, M3_A2)
# still ok


# V2R_bA2
V2R_bA2 <- read.csv("V2R_bA2_format.csv")
V2R_bA2$condition <- as.factor(V2R_bA2$condition)
aov_V2R_bA2<- aov(fold_change ~ condition, data = V2R_bA2)
ggqqplot(residuals(aov_V2R_bA2))
# one "outlier".. but only n =3
plot(V2R_bA2)
which.max(V2R_bA2$fold_change)
V2R_bA2$fold_change[16]
boxplot(fold_change ~ condition, V2R_bA2)

# b2AR bA1
b2AR_bA1 <- read.csv("b2AR_curves_bA1_format.csv")
b2AR_bA1$condition <- as.factor(b2AR_bA1$condition)
boxplot(fold_change ~ condition, b2AR_bA1)
boxplot(fold_change ~ condition, b2AR_bA1[-16,])
# exclude first n of CON
b2AR_bA1 <- b2AR_bA1[-16,]
aov_b2AR_bA1<- aov(fold_change ~ condition, data = b2AR_bA1)
ggqqplot(residuals(aov_b2AR_bA1))
summary(aov_b2AR_bA1)
mult_comp_b2AR_bA1_correct <- summary(glht(aov_b2AR_bA1, 
                          linfct = mcp(condition = contrast_matrix)), 
                     test = adjusted("bonferroni"))
# save this in results
contrast<-names(mult_comp_b2AR_bA1_correct$test$coefficients)
post_hoc <- as.data.frame(contrast)
post_hoc$p_value <- mult_comp_b2AR_bA1_correct$test$pvalues

results_curve_corrected$b2AR_curves_bA1 <-post_hoc

#also save unadjusted
mult_b2AR_bA1_correct_noadj <- 
  summary(glht(aov_b2AR_bA1,
               linfct = mcp(condition = contrast_matrix)), 
          test = adjusted("none"))
post_hoc$p_value <- mult_b2AR_bA1_correct_noadj$test$pvalues
results_curve_corrected$b2AR_curves_bA1_noadj <-post_hoc

# also save details
post_hoc_detail_b2AR_bA1 <- post_hoc_detail
post_hoc_detail_b2AR_bA1$effect_size <- mult_b2AR_bA1_correct_noadj$test$coefficients

b2AR_bA1_conf <- confint(glht(aov_b2AR_bA1, 
                              linfct = mcp(condition = contrast_matrix)))
post_hoc_detail_b2AR_bA1$lower <- b2AR_bA1_conf $confin[,2]
post_hoc_detail_b2AR_bA1$upper <- b2AR_bA1_conf $confin[,3]

post_hoc_detail_b2AR_bA1$tstat <- mult_b2AR_bA1_correct_noadj$test$tstat
post_hoc_detail_b2AR_bA1$df <- mult_b2AR_bA1_correct_noadj$df

post_hoc_detail_b2AR_bA1$pvalue <-   mult_b2AR_bA1_correct_noadj$test$pvalues
post_hoc_detail_b2AR_bA1$pvalue_adjusted <- mult_comp_b2AR_bA1_correct$test$pvalues

results_detail$b2AR_curves_bA1 <- post_hoc_detail_b2AR_bA1


##### save corrected results ####

# save as xlsx
remove(results_wb)
results_wb<-createWorkbook()
for (experiment in list.names(results_curve)){
  file = paste(experiment, ".xlsx", sep = "")
  addWorksheet(results_wb, file)
  writeData(results_wb, file, results_curve_corrected[[experiment]], 
            startRow = 1, 
            startCol = 1)
}

saveWorkbook(results_wb, file = "Curve_results_corrected.xlsx", overwrite = TRUE)


# save for visualizations
results_curve_table_cor<- data.frame()
results_curve_table_cor_noadj <- data.frame()
exp_adj<- c()
exp_noadj<- c()

noadjusted_list <- str_extract(list.names(results_curve_corrected), "noadj")

for (i in 1:length(results_curve_corrected)){
  if (is.na(noadjusted_list[i])){
    results_curve_table_cor<-rbind(results_curve_table_cor, results_curve_corrected[[i]])
    exp_adj <- c(exp_adj, rep(list.names(results_curve_corrected)[i], 9))
  }else{
    results_curve_table_cor_noadj<-
      rbind(results_curve_table_cor_noadj, results_curve_corrected[[i]])
    exp_noadj <- c(exp_noadj, rep(list.names(results_curve_corrected)[i], 9))
  }
}

results_curve_table_cor$experiment <- as.factor(exp_adj)
results_curve_table_cor$contrast <- as.factor(results_curve_table_cor$contrast)

results_curve_table_cor_noadj$experiment <- as.factor(exp_noadj)
results_curve_table_cor_noadj$contrast <- as.factor(results_curve_table_cor$contrast)


# save adjusted for visualization
remove(results_wb)
results_wb<-createWorkbook()
addWorksheet(results_wb, "results")
writeData(results_wb, "results", results_curve_table_cor, 
          startRow = 1, 
          startCol = 1)
saveWorkbook(results_wb, 
             file = "Curve_results_overview_corrected.xlsx", 
             overwrite = TRUE)

# save non - adjusted save adjusted for visualization
remove(results_wb)
results_wb<-createWorkbook()
addWorksheet(results_wb, "results")
writeData(results_wb, "results", results_curve_table_cor_noadj, 
          startRow = 1, 
          startCol = 1)
saveWorkbook(results_wb, 
             file = "Curve_results_overview_corrected_noadj.xlsx", 
             overwrite = TRUE)




###p value table of corrected adjusted and non-adjusted####

p_value_table <- data.frame(experiment_names)

# create list of gpcrs from experiment names 
#by removing arrestin specification
gpcr <- experiment_names %>%
  str_remove("_barr1") %>%
  str_remove("_barr2") %>%
  str_remove("_bA1") %>%
  str_remove("_bA2") %>%
  str_remove("_curves") %>%
  str_remove("_2nd")

p_value_table$GPCR<-gpcr

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

p_value_table$barrestin<-barr

#-> collect p values from 
#results_curve_table_cor and results_curve_table_cor_noadj in one table

c = length(names(p_value_table))
for (title in contrast){
  c = c+1
  temp_noadj <- 
    subset(results_curve_table_cor_noadj, contrast == title)$p_value
  temp_adj <- subset(results_curve_table_cor, contrast == title)$p_value
  
  p_value_table[c] <- temp_noadj
  names(p_value_table)[c] <- paste(title, "noadj", sep="_")
  
  
  c = c+1
  p_value_table[c] <- temp_adj
  names(p_value_table)[c] <- paste(title, "adj", sep="_")
  
  temp_adj_stars <- c()
  for (p in temp_adj){
    print(p)
    if (p > 0.05){
      temp_adj_stars <- c(temp_adj_stars, " ")
      print("n.s.")}
    
    if (between(p, 0.01, 0.05)){
      temp_adj_stars <- c(temp_adj_stars, "*")
      print("*")}
    
    if (between(p, 0.001, 0.01)){
      temp_adj_stars <- c(temp_adj_stars, "**")
      print("**")}
    
    if (p < 0.001){
      temp_adj_stars <- c(temp_adj_stars, "***")
      print("***")}

  }
  c = c+1
  p_value_table[c] <- temp_adj_stars
}


#remove experiments which are not presented in paper
exclude <- c("V2R_dFLR_bA1",
             "V2R_dFLR_bA2")
#get row numbers
exclude_rows <- c()
for (exp in 1:length(exclude)){
  temp <- which(p_value_table$experiment_names == exclude[exp])
  exclude_rows <- c(exclude_rows, temp)
}

p_value_table_select<-p_value_table[-exclude_rows,]

#save p value table
remove(results_wb)
results_wb<-createWorkbook()
addWorksheet(results_wb, "p_values")
writeData(results_wb, "p_values", p_value_table_select, 
          startRow = 1, 
          startCol = 1)
saveWorkbook(results_wb, 
             file = "P_values_overview.xlsx", 
             overwrite = TRUE)



#-> collect details from results_detail in one table
c = length(names(details_table))
effect_sizes <- c()
lower <- c()
upper <- c()
tstat <- c()
df <- c()

for (title in names(results_detail)){
  temp <- results_detail[[title]]
  effect_sizes <- c(effect_sizes, temp$effect_size)
  lower <- c(lower, temp$lower)
  upper <- c(upper, temp$upper)
  tstat <- c(tstat, temp$tstat)
  df <- c(df, temp$df)
}

details_table <- data.frame(results_curve_table_cor$experiment)
names(details_table) <- "experiment"
details_table$contrast <- results_curve_table_cor$contrast

details_table$effect_sizes <- effect_sizes 
details_table$lower <- lower
details_table$upper <- upper
details_table$tstat <- tstat
details_table$df <- df


# add details to p value table
details_table_overview <- data.frame(experiment_names)
details_table_overview$GPCR<-gpcr
details_table_overview$barrestin<-barr
# take p values from p_value_table, afterwards exclude experiments not in paper

c = length(names(details_table_overview))
n = 4
for (title in contrast){
  c = c+1
  temp_effect_sizes <- 
    subset(details_table, contrast == title)$effect_sizes
  temp_lower <- 
    subset(details_table, contrast == title)$lower
  temp_upper <- 
    subset(details_table, contrast == title)$upper
  temp_tstat <- 
    subset(details_table, contrast == title)$tstat
  temp_df <- 
    subset(details_table, contrast == title)$df

  details_table_overview[c] <- temp_effect_sizes
  names(details_table_overview)[c] <- paste(title, "effect_sizes", sep="_")
  
  details_table_overview[c+1] <- temp_lower
  names(details_table_overview)[c+1] <- paste(title, "lower", sep="_")
  
  details_table_overview[c+2] <- temp_upper
  names(details_table_overview)[c+2] <- paste(title, "upper", sep="_")
  
  details_table_overview[c+3] <- temp_tstat 
  names(details_table_overview)[c+3] <- paste(title, "tstat", sep="_")
  
  details_table_overview[c+4] <- temp_df
  names(details_table_overview)[c+4] <- paste(title, "df", sep="_")
  
  details_table_overview[c+5] <- p_value_table[n]
  details_table_overview[c+6] <- p_value_table[n+1]
  details_table_overview[c+7] <- p_value_table[n+2]
  n = n + 3
  c = c+7
}

#remove experiments which are not presented in paper
exclude <- c("V2R_dFLR_bA1",
             "V2R_dFLR_bA2")
#get row numbers
exclude_rows <- c()
for (exp in 1:length(exclude)){
  temp <- which(details_table_overview$experiment_names == exclude[exp])
  exclude_rows <- c(exclude_rows, temp)
}

details_table_overview_select<-details_table_overview[-exclude_rows,]

#save p value table
remove(results_wb)
results_wb<-createWorkbook()
addWorksheet(results_wb, "details")
writeData(results_wb, "details", details_table_overview_select, 
          startRow = 1, 
          startCol = 1)
saveWorkbook(results_wb, 
             file = "details_table_overview.xlsx", 
             overwrite = TRUE)






