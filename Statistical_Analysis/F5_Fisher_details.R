# path definition
path = "C:/path/to/folder/"

setwd(path)

wants <- c("openxlsx", 
           "pheatmap", 
           "rstatix", 
           "stringr",
           "xlsx")
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
lapply(wants, require, character.only = TRUE)


#removing duplicates position data (created by the GPCR-barr pairs)
remove_duplicates <- function(P_pos_in_table, P_pos_table){
  rows_to_remove <- c()
  for (u in 1:length(P_pos_in_table)){
    temp <- P_pos_in_table[u]
    # which rows have same position
    temp_list <- which(P_pos_in_table == temp)
    # remove same row fron list which is currently looked at
    temp_list <- temp_list[-which(temp_list == u)]
    # remove already tested rows from list
    if (length(which(temp_list < (u + 0.5))) > 0){
      temp_list <- temp_list[-which(temp_list < (u + 0.5))]
    }
    decide <- c()
    for (index in temp_list){
      # do they have same GPCR? and structure status? and GRK group?
      decide <- c(P_pos_table$GPCR[index] == P_pos_table$GPCR[u],
                  P_pos_table$IL3_or_Cterm[index] == P_pos_table$IL3_or_Cterm[u],
                  P_pos_table$GRK[index] == P_pos_table$GRK[u])
      if (all(decide) == TRUE){
        rows_to_remove <- c(rows_to_remove, index)
      }
    }
  }
  output <- P_pos_table[-rows_to_remove,]
  return(output)
}

# input: position of PXPP
PXPP_pos <- read.csv("PXPP_pos.csv")
PXPP_pos$P_type <- rep("PXPP", length(PXPP_pos[1]))
PXPP_pos<- PXPP_pos[-which(PXPP_pos$GPCR == "b2V2R"),]
PXPP_pos_select <- remove_duplicates(PXPP_pos$PXPP_position, PXPP_pos)

# defining end and middle
PXPP_end <- PXPP_pos_select[c(which((PXPP_pos_select$PXPP_position/PXPP_pos_select$length) < 0.25),
                              which((PXPP_pos_select$PXPP_position/PXPP_pos_select$length) > 0.75)),]

PXPP_middle <- PXPP_pos_select[-c(which((PXPP_pos_select$PXPP_position/PXPP_pos_select$length) < 0.25),
                                  which((PXPP_pos_select$PXPP_position/PXPP_pos_select$length) > 0.75)),]

#### analysis #### 
# amount of PXPP ~ GRK group
length(subset(PXPP_end, GRK=="GRK23")[,1])
length(subset(PXPP_middle, GRK=="GRK23")[,1])

length(subset(PXPP_end, GRK=="GRK2356")[,1])
length(subset(PXPP_middle, GRK=="GRK2356")[,1])

# creating contingency table ~ GRK group
GRKgroup_compare <- data.frame(GRK23 = c(2,15), GRK2356 = c(10,3), row.names = c("peripheral", "central"))
GRK_fisher <-fisher.test(GRKgroup_compare) #0.0005367


# amount of PXPP ~ barr class
length(subset(PXPP_end, arr_class=="A")[,1])
length(subset(PXPP_middle, arr_class=="A")[,1])

length(subset(PXPP_end, arr_class=="B")[,1])
length(subset(PXPP_middle, arr_class=="B")[,1])

# creating contingency table ~ barr class
barrClass_compare <- data.frame(class_A = c(14,18), class_B = c(4,8), row.names = c("peripheral", "central"))
barr_fisher <- fisher.test(barrClass_compare) #0.7328


# fisher results
fisher_results <- data.frame(group = c("GRK","barr"))
fisher_results$odds_ratio <- c(GRK_fisher$estimate, barr_fisher$estimate)
fisher_results$lower_CI <- c(GRK_fisher$conf.int[1], barr_fisher$conf.int[1])
fisher_results$upper_CI <- c(GRK_fisher$conf.int[2], barr_fisher$conf.int[2])
fisher_results$p_value <- c(GRK_fisher$p.value, barr_fisher$p.value)


#save p value table
remove(results_wb)
results_wb<-createWorkbook()
addWorksheet(results_wb, "GRK_group_input")
writeData(results_wb, "GRK_group_input", GRKgroup_compare, 
          startRow = 1, 
          startCol = 1,
          rowNames = TRUE)
addWorksheet(results_wb, "barr_group_input")
writeData(results_wb, "barr_group_input", barrClass_compare, 
          startRow = 1, 
          startCol = 1,
          rowNames = TRUE)
addWorksheet(results_wb, "Fisher_test_results")
writeData(results_wb, "Fisher_test_results", fisher_results, 
          startRow = 1, 
          startCol = 1)
saveWorkbook(results_wb, 
             file = "Figure5_results.xlsx", 
             overwrite = TRUE)

