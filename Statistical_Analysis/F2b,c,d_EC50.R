# path for results export
path = "C:/path/to/folder/"

setwd(path)

# packages
wants <- c("multcomp",
           "dplyr",
           "openxlsx")
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])

lapply(wants, require, character.only = TRUE)

####### functions #######
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

# generates data from mean + sd -> is needed for aov() input
gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}


# function for analysis + details

aov_dunnett <- function(simulated_df){
  simulated_df$factor <- as.factor(simulated_df[,3])
  simulated_aov <- aov(y ~ factor, data = simulated_df)
  print(summary(simulated_aov))
  dunn <-summary(glht(simulated_aov, linfct = mcp(factor = "Dunnett")))
  
  contrast <- names(dunn$test$coefficients)
  results <- as.data.frame(contrast)
  results$effect_size <- dunn$test$coefficients
  
  conf_int <- confint(glht(simulated_aov, 
                           linfct = mcp(factor = "Dunnett")))
  results$lower <- conf_int$confin[,2]
  results$upper <- conf_int$confin[,3]
  
  results$tstat <- dunn$test$tstat
  results$df <- dunn$df
  results$pvalue <- dunn$test$pvalues
  results$signi <- stars_annotation(dunn$test$pvalues)
  return(results)
}
save_results <- function(data_list){
  results <- vector("list")
  n = 0
  for (set in data_list){
    n = n +1
    results[[names(data_list)[n]]] <- aov_dunnett(set)
  }
  return(results)
}
results_export <- function(results){
  n = 0
  for (set in results){
    n = n +1
    temp_df <- as.data.frame(set)
    print(temp_df)
    temp_df <- cbind(a = rep(names(results)[n], length(set[,1])), temp_df)
    if (n == 1){
      results_export <- as.data.frame(temp_df)
    }else{
      results_export <- rbind(results_export, temp_df)
    }
  }
  write.xlsx(results_export, "Results_EC50.xlsx")
  return(results_export)
}

##### EC50 data #####
# mean + SEM of EC50 as calculated in prism
condition <- c(rep("con", 3), rep("GRK2",3), rep("GRK6", 3))

# isoprenaline
mean <- c(0.09473, 0.03679, 0.02605)
SEM <- c(0.047510154, 0.020005187, 0.015003024)
SD <- SEM/sqrt(3)
n <- c(3, 3, 3)
iso <- as.data.frame(mean)
iso$SD <- SD
iso$n <- n
iso_simulated <- gen_data(iso$mean, iso$SD, iso$n)
iso_simulated$condition <- as.factor(condition)

# epinephrine
mean <- c(0.3485, 0.3234, 0.09021)
SEM <- c(0.232614423, 0.2249934, 0.060275368)
SD <- SEM/sqrt(3)
n <- c(3, 3, 3)
epi <- as.data.frame(mean)
epi$SD <- SD
epi$n <- n
epi_simulated <- gen_data(epi$mean, epi$SD, epi$n)
epi_simulated$condition <- as.factor(condition)

# norepinephrine
mean <- c(7.109, 5.928, 1.648)
SEM <- c(4.95539736, 3.555900308, 1.057243813)
SD <- SEM/sqrt(3)
n <- c(3, 3, 3)
norepi <- as.data.frame(mean)
norepi$SD <- SD
norepi$n <- n
norepi_simulated <- gen_data(norepi$mean, norepi$SD, norepi$n)
norepi_simulated$condition <- as.factor(condition)

# list from input data
data_list <- vector("list")
data_list[["iso"]] <- iso_simulated
data_list[["epi"]] <- epi_simulated
data_list[["norepi"]] <- norepi_simulated

##### analysis #####
analysis <- save_results(data_list)
results_export(analysis)











