# path definition
path = "C:/Users/monar/Google Drive/Arbeit/2021_Drube_et_al/210929_other statistics"
setwd(path)

# packages
wants <- c("openxlsx",
           "rstatix", 
           "multcomp", 
           "stringr", 
           "readxl", 
           "dplyr",
           "ggpubr",
           "rstatix")
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])

lapply(wants, require, character.only = TRUE)



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

##### import data Figure 1d #####
data<- read.xlsx("Control_cells_GRK_expression.xlsx", sheet= 2)

##### analysis #####
data$antibody <- as.factor(data$antibody)

aov <- aov(intensity ~ antibody, data = data)
result <- summary(glht(aov, linfct = mcp(antibody = "Tukey")))

contrasts <- names(result$test$coefficients)
data_results <- as.data.frame(contrasts)
data_results$effect_size <- result$test$coefficients
conf_int <- confint(glht(aov, 
                         linfct = mcp(antibody = "Tukey")))
data_results$lower <- conf_int$confin[,2]
data_results$upper <- conf_int$confin[,3]
data_results$tstat <- result$test$tstat
data_results$df <- result$df
data_results$pvalue <- result$test$pvalues
data_results$signi <- stars_annotation(result$test$pvalues)

##### export #####
write.xlsx(data_results, "F1_Control_expression.xlsx")



##### import data - SF2c fold over endogenous #####
data<- read.xlsx("Quantifizierung Fold over endogenous SupplFig2.xlsx", sheet= 2)

##### analysis #####
data$antibody <- as.factor(data$antibody)

aov <- aov(intensity ~ antibody, data = data)
summary(aov)

result <- summary(glht(aov, linfct = mcp(antibody = "Tukey")))

contrasts <- names(result$test$coefficients)
data_results <- as.data.frame(contrasts)
data_results$effect_size <- result$test$coefficients
conf_int <- confint(glht(aov, 
                         linfct = mcp(antibody = "Tukey")))
data_results$lower <- conf_int$confin[,2]
data_results$upper <- conf_int$confin[,3]
data_results$tstat <- result$test$tstat
data_results$df <- result$df
data_results$pvalue <- result$test$pvalues
data_results$signi <- stars_annotation(result$test$pvalues)


##### export #####
write.xlsx(data_results, "SF2_fold_overexpression.xlsx")



##### import data - F1f GRK-YFP intensity #####
data<- read.xlsx("GRK-YFP_Fig1f.xlsx", sheet= 3) # _mean


##### analysis #####
data$GRK <- as.factor(data$GRK)

aov <- aov(intensity ~ GRK, data = data)
summary(aov)

result <- summary(glht(aov, linfct = mcp(GRK = "Tukey")))

contrasts <- names(result$test$coefficients)
data_results <- as.data.frame(contrasts)
data_results$effect_size <- result$test$coefficients
conf_int <- confint(glht(aov, 
                         linfct = mcp(GRK = "Tukey")))
data_results$lower <- conf_int$confin[,2]
data_results$upper <- conf_int$confin[,3]
data_results$tstat <- result$test$tstat
data_results$df <- result$df
data_results$pvalue <- result$test$pvalues
data_results$signi <- stars_annotation(result$test$pvalues)


##### export #####
write.xlsx(data_results, "F1f_GRK-YFP.xlsx")

