path = "C:/path/to/folder/"
setwd(path)

wants <- c("openxlsx", "RColorBrewer", "pheatmap", "rlist", "stringr")
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
lapply(wants, require, character.only = TRUE)

# import results from curve analysis
file_names = list.files(pattern="*.csv", full.names = TRUE)
experiment_names_OG = c()
for (file in file_names){
  experiment_names_OG = c(experiment_names_OG, substr(file, 3, nchar(file)-11))
}

results_files = list.files(pattern="results_overview", full.names = TRUE)
results = lapply(results_files, read.xlsx)
names(results) <- results_files

contrast_matrix <- rbind("dQ -GRK2" = c(0, -1, 1, 0, 0, 0),
                         "dQ -GRK3" = c(0, -1, 0, 1, 0, 0),
                         "dQ -GRK5" = c(0, -1, 0, 0, 1, 0),
                         "dQ -GRK6" = c(0, -1, 0, 0, 0, 1),
                         "dQ -CON" = c(1, -1, 0, 0, 0, 0),
                         "CON -GRK2" = c(1, 0, 1, 0, 0, 0),
                         "CON -GRK3" = c(1, 0, 0, 1, 0, 0),
                         "CON -GRK5" = c(1, 0, 0, 0, 1, 0),
                         "CON -GRK6" = c(1, 0, 0, 0, 0, 1))
number_contr <- length(rownames(contrast_matrix))

# without PD mutants + dFLR; a subset of experiment_names

exclude_exp_OG <- c("V2R_dFLR_bA1", 
                    "V2R_dFLR_bA2",
                    "AT1R_barr1",
                    "AT1R_barr2",
                    "V2R_bA2",
                    "V2R_curves_bA1")

list_noadj<-str_extract(list.names(results), "noadj")

all_data_heatmap<-list()
i = 0
for (result in results){
  i = i +1
  # matrix for heatmap
  # each row one receptor + arr, each column is one contrast
  data_for_heatmap <- matrix(1:number_contr, ncol = number_contr)
  
  if (is.na(list_noadj[i])){
    experiment_names <- experiment_names_OG
    exclude_exp <-exclude_exp_OG
  }else{
    experiment_names <- experiment_names_OG
    exclude_exp <-exclude_exp_OG
    experiment_names <- paste(experiment_names_OG, "noadj", sep = "_")
    exclude_exp <- paste(exclude_exp_OG, "noadj", sep = "_")
  }
  
  for (name in experiment_names){
  data_for_heatmap <- 
    rbind(data_for_heatmap, 
          subset(result, experiment == name)$p_value)
  }
  data_for_heatmap <- data_for_heatmap[-1,]
  rownames(data_for_heatmap) <- experiment_names
  colnames(data_for_heatmap) <- rownames(contrast_matrix)
  
  # exclude experiments from heatmap
  data_for_heatmap_select <- data_for_heatmap
  for (experiment in exclude_exp){
    row = which(rownames(data_for_heatmap_select) == experiment)
    data_for_heatmap_select <- data_for_heatmap_select[-row,]
  }
  name <- substr(names(results)[i], 3, nchar(names(results)[i])-5)
  all_data_heatmap[[name]] <- data_for_heatmap_select
    
}
  
# calculate log10 of all pvalues
all_data_heatmap_log <- list()
n = 0
for (data in all_data_heatmap){
  n = n+1
  all_data_heatmap_log[[list.names(all_data_heatmap)[n]]] <- log10(data)
}


#### heatmap ####
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

color<-colorRampPalette(brewer.pal(8, "RdBu"))(25)

rowlabels <- c("b2ADR   bA2", "b2ADR   bA1",
               "b2V2   bA2", "b2V2   bA1",
               "C5aR   bA1", "C5aR   bA2",
               "M1R   bA1", "M1R   bA2",
               "M2R   bA1", "M2R   bA2",
               "M3R   bA1", "M3R   bA2",
               "M4R   bA1", "M4R   bA2",
               "M5R   bA1", "M5R   bA2",
               "mOR   bA1", "mOR   bA2",
               "PTH1R   bA1", "PTH1R   bA2")

breaks <- c(-6, -3.5, -3, -2.5, -2.25, -2, -1.75, -1.5, -1.25, -1, -0.75, -0.5, -0.25, -0.15, -0.01)
color<- c(color[1], color[3], color[4.5],
                color[5.5], color[7.5], color[9.5], color[12], color[14],
                color[16], color[18], color[21],
                color[22], color[24],
                color[25])

pheatmap(all_data_heatmap_log$Curve_results_overview_corrected_noadj, 
         border_color = NA, 
         cluster_cols = FALSE, 
         cluster_rows = TRUE,
         col = color, 
         breaks = breaks,
         cellwidth = 60,
         cellheight = 25,
         #cutree_rows = 4,
         clustering_distance_rows = "canberra",
         treeheight_row = 35,
         fontsize_row = 12,
         fontsize_col = 12,
         #display_numbers = TRUE,
         main = "Curve_results_overview_corrected_noadj",
         #legend_breaks = c(-6, -3.5, -3, -2.5, -2, -1.5, -1, -0.5),
         legend_breaks = c(-6, -3.5, -3, -2.5,-2.25, -2, -1.5, -1, -0.5),
         legend_labels = c("-6", "-3.5","-3", "-2.5", "----adjusted_p = 0.05", "-2", "-1.5", "-1", "-0.5"),
         #legend_breaks = c(-2.255272, -1.954243),
         #legend_labels = c("--adjusted_p = 0.05", "--adjusted_p = 0.1"),
         labels_row= rowlabels
)

