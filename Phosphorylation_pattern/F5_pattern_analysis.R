path = "C:/Users/monar/Google Drive/Arbeit/Drube et al. phosphorylation pattern/"
setwd(path)

wants <- c("openxlsx", "pheatmap", "rstatix", "stringr", "ggplot2", "hrbrthemes", "RColorBrewer")
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
lapply(wants, require, character.only = TRUE)

count <- read.csv("moving_frame_count.csv")
count$GPCR <- as.factor(count$GPCR)
count$barr <- as.factor(count$barr)
count$GRK <- as.factor(count$GRK)
count$arr_class <- as.factor(count$arr_class)
count$IL3_or_Cterm <- as.factor(count$IL3_or_Cterm)

count$sum_count <- count$p_count + count$n_count
## remove b2V2 from analysis
count <- count[-which(count$GPCR == "b2V2R"),]

P_pos <- read.csv("P_pos.csv")
P_pos<- P_pos[-which(P_pos$GPCR == "b2V2R"),]

PPP_pos <- read.csv("PPP_pos.csv")
PPP_pos$P_type <- rep("PPP", length(PPP_pos[1]))
PPP_pos<- PPP_pos[-which(PPP_pos$GPCR == "b2V2R"),]

PXPP_pos <- read.csv("PXPP_pos.csv")
PXPP_pos$P_type <- rep("PXPP", length(PXPP_pos[1]))
PXPP_pos<- PXPP_pos[-which(PXPP_pos$GPCR == "b2V2R"),]

PXPXXP_pos <- read.csv("PXPXXP_pos.csv")
PXPXXP_pos$P_type <- rep("PXPXXP", length(PXPXXP_pos[1]))
PXPXXP_pos<- PXPXXP_pos[-which(PXPXXP_pos$GPCR == "b2V2R"),]

PXXPXXP_pos <- read.csv("PXXPXXP_pos.csv")
PXXPXXP_pos$P_type <- rep("PXXPXXP", length(PXXPXXP_pos[1]))
PXXPXXP_pos<- PXXPXXP_pos[-which(PXXPXXP_pos$GPCR == "b2V2R"),]


#### removing duplicate position data ####

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

P_pos_select <- remove_duplicates(P_pos$PPP_position, P_pos)
PPP_pos_select <- remove_duplicates(PPP_pos$PPP_position, PPP_pos)
PXPP_pos_select <- remove_duplicates(PXPP_pos$PXPP_position, PXPP_pos)
PXPXXP_pos_select <- remove_duplicates(PXPXXP_pos$PXPXXP_position, PXPXXP_pos)
PXXPXXP_pos_select <- remove_duplicates(PXXPXXP_pos$PXXPXXP_position, PXXPXXP_pos)

rows_to_remove <- c()

#### test removing duplikates from count data ####
rows_to_remove <- c()
for (u in 1:length(count$GPCR)){
  temp <- count$GPCR[u]
  # which rows have same GPCR
  temp_list <- which(count$GPCR == temp)
  # remove same row fron list which is currently looked at
  temp_list <- temp_list[-which(temp_list == u)]
  # remove already tested rows from list
  if (length(which(temp_list < (u + 0.5))) > 0){
    temp_list <- temp_list[-which(temp_list < (u + 0.5))]
  }
  decide <- c()
  for (index in temp_list){
    # do they have  same structure? and same GRK group?
    decide <- c(count$IL3_or_Cterm[index] == count$IL3_or_Cterm[u],
                count$GRK[index] == count$GRK[u])
    if (all(decide) == TRUE){
      rows_to_remove <- c(rows_to_remove, index)
    }
  }
}

count_select <- count[-rows_to_remove,]


# create table containing all position information
all_pos <- data.frame(c(PPP_pos[,1], PXPP_pos[,1], PXPXXP_pos[,1], PXXPXXP_pos[,1]))
for (col in 2:length(PPP_pos)){
  all_pos[col] <- c(PPP_pos[,col], PXPP_pos[,col], PXPXXP_pos[,col], PXXPXXP_pos[,col])
}
names(all_pos) <- c(names(PPP_pos)[1:6], "position", names(PPP_pos)[8])


summary(subset(count, GRK == "GRK23")$arr_class)
summary(subset(count, GRK == "GRK2356")$arr_class)
summary(count$arr_class)

###### make pie chart for Class A and Class B receptors ######

pie_classA <- c(
  as.numeric(summary(subset(count, GRK == "GRK23")$arr_class)[1]),
  as.numeric(summary(subset(count, GRK == "GRK2356")$arr_class)[1]),
  as.numeric(summary(subset(count, GRK == "ns")$arr_class)[1]),
  as.numeric(summary(subset(count, GRK == "precoupling")$arr_class)[1])
)

pie_classB <- c(
  as.numeric(summary(subset(count, GRK == "GRK23")$arr_class)[2]),
  as.numeric(summary(subset(count, GRK == "GRK2356")$arr_class)[2]),
  as.numeric(summary(subset(count, GRK == "ns")$arr_class)[2]),
  as.numeric(summary(subset(count, GRK == "precoupling")$arr_class)[2])
)

###  pie charts
pie(pie_classA, 
    labels = c("GRK23", "GRK2356", "ns", "precoupling"),
    edges = 1000,
    border = "white",
    col = brewer.pal(4, name = "Set1" ),
    main = "Class_A") 

pie_classB <- c(pie_classB[1], pie_classB[2], pie_classB[4], pie_classB[3])
pie(pie_classB, 
    labels = c("GRK23", "GRK2356", "precoupling", "ns"),
    edges = 1000,
    border = "white",
    col = brewer.pal(4, name = "Set1" ),
    main = "Class_B")

##### version 02 #####

# abundance of PPP clusters ~ GRK group
count_select %>%
  ggplot(aes(x = GRK, y = PPP_count)) +
  geom_dotplot(aes(fill = GRK),
               binaxis = "y",
               binwidth = 0.008,
               dotsize = 60,
               stackdir = "center",
               fill = "darkgrey") +
  ylim(0,15) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5) + 
  theme(panel.background = element_rect(fill="white"), 
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))

# abundance of PXPP clusters ~ GRK group
count_select %>%
  ggplot(aes(x = GRK, y = PXPP_count)) +
  geom_dotplot(aes(fill = GRK),
               binaxis = "y",
               binwidth = 0.008,
               dotsize = 60,
               stackdir = "center",
               fill = "darkgrey") +
  ylim(0,15) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5) + 
  theme(panel.background = element_rect(fill="white"), 
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))

# abundance of PXPXXP pattern ~ GRK group
count_select %>%
  ggplot(aes(x = GRK, y = PXPXXP_count)) +
  geom_dotplot(aes(fill = GRK),
               binaxis = "y",
               binwidth = 0.008,
               dotsize = 60,
               stackdir = "center",
               fill = "darkgrey") +
  ylim(0,15) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5) + 
  theme(panel.background = element_rect(fill="white"), 
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))

# abundance of PXXPXXP pattern ~ GRK group
count_select %>%
  ggplot(aes(x = GRK, y = PXXPXXP_count)) +
  geom_dotplot(aes(fill = GRK),
               binaxis = "y",
               binwidth = 0.008,
               dotsize = 60,
               stackdir = "center",
               fill = "darkgrey") +
  ylim(0,15) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5) + 
  theme(panel.background = element_rect(fill="white"), 
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))


# P positions distribution
P_pos_select %>%
  ggplot(aes(x = GRK, y = PPP_position/length)) +
  geom_dotplot(aes(fill = GRK),
               binaxis = "y",
               binwidth = 0.008,
               dotsize = 4.2,
               stackdir = "center",
               fill = "darkgrey") +
  ylim(0,1) +
  ggtitle("P positions") +
  stat_summary(fun = median, geom = "crossbar", width = 0.5) + 
  theme(panel.background = element_rect(fill="white"), 
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))

# PPP position ~ GRK group
PPP_pos_select %>%
  ggplot(aes(x = GRK, y = PPP_position/length)) +
  geom_dotplot(aes(fill = GRK),
               binaxis = "y",
               binwidth = 0.008,
               dotsize = 4.2,
               stackdir = "center",
               fill = "darkgrey") +
  ylim(0,1) +
  ggtitle("PPP positions") +
  stat_summary(fun = median, geom = "crossbar", width = 0.5) + 
  theme(panel.background = element_rect(fill="white"), 
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))

# PPP position  ~ barr class
PPP_pos_select %>%
  ggplot(aes(x = arr_class, y = PPP_position/length)) +
  geom_dotplot(aes(fill = arr_class),
               binaxis = "y",
               binwidth = 0.008,
               dotsize = 4.2,
               stackdir = "center",
               fill = "darkgrey") +
  ylim(0,1) +
  ggtitle("PPP positions") +
  stat_summary(fun = median, geom = "crossbar", width = 0.5) + 
  theme(panel.background = element_rect(fill="white"), 
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))


# PXPP position ~ GRK group
PXPP_pos_select %>%
  ggplot(aes(x = GRK, y = PXPP_position/length)) +
  geom_dotplot(aes(fill = GRK),
               binaxis = "y",
               binwidth = 0.008,
               dotsize = 4.2,
               stackdir = "center",
               fill = "darkgrey") +
  ylim(0,1) +
  ggtitle("PXPP positions") +
  stat_summary(fun = median, geom = "crossbar", width = 0.5) + 
  theme(panel.background = element_rect(fill="white"), 
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))

# PXPP position ~ barr class
PXPP_pos_select %>%
  ggplot(aes(x = arr_class, y = PXPP_position/length)) +
  geom_dotplot(aes(fill = arr_class),
               binaxis = "y",
               binwidth = 0.008,
               dotsize = 4.2,
               stackdir = "center",
               fill = "darkgrey") +
  ylim(0,1) +
  ggtitle("PXPP positions") +
  stat_summary(fun = median, geom = "crossbar", width = 0.5) + 
  theme(panel.background = element_rect(fill="white"), 
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))


# PXPXXP position ~ GRK group
PXPXXP_pos_select %>%
  ggplot(aes(x = GRK, y = PXPXXP_position/length)) +
  geom_dotplot(aes(fill = GRK),
               binaxis = "y",
               binwidth = 0.008,
               dotsize = 4.2,
               stackdir = "center",
               fill = "darkgrey") +
  ylim(0,1) +
  ggtitle("PXPXXP positions") +
  stat_summary(fun = median, geom = "crossbar", width = 0.5) + 
  theme(panel.background = element_rect(fill="white"), 
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))

# PXPXXP position ~ barr class
PXPXXP_pos_select %>%
  ggplot(aes(x = arr_class, y = PXPXXP_position/length)) +
  geom_dotplot(aes(fill = arr_class),
               binaxis = "y",
               binwidth = 0.008,
               dotsize = 4.2,
               stackdir = "center",
               fill = "darkgrey") +
  ylim(0,1) +
  ggtitle("PXPXXP positions") +
  stat_summary(fun = median, geom = "crossbar", width = 0.5) + 
  theme(panel.background = element_rect(fill="white"), 
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))

# PXXPXXP position ~ GRK group
PXXPXXP_pos_select %>%
  ggplot(aes(x = GRK, y = PXXPXXP_position/length)) +
  geom_dotplot(aes(fill = GRK),
               binaxis = "y",
               binwidth = 0.008,
               dotsize = 4.2,
               stackdir = "center",
               fill = "darkgrey") +
  ylim(0,1) +
  ggtitle("PXXPXXP positions") +
  stat_summary(fun = median, geom = "crossbar", width = 0.5) + 
  theme(panel.background = element_rect(fill="white"), 
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))

# PXXPXXP position ~ barr class
PXXPXXP_pos_select %>%
  ggplot(aes(x = arr_class, y = PXXPXXP_position/length)) +
  geom_dotplot(aes(fill = arr_class),
               binaxis = "y",
               binwidth = 0.008,
               dotsize = 4.2,
               stackdir = "center",
               fill = "darkgrey") +
  ylim(0,1) +
  ggtitle("PXXPXXP positions") +
  stat_summary(fun = median, geom = "crossbar", width = 0.5) + 
  theme(panel.background = element_rect(fill="white"), 
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))

##### export of plottet data ####
setwd("C:/Users/monar/Google Drive/Arbeit/Drube et al. phosphorylation pattern/plottet data_for nature comm sub/")
write.xlsx(count_select, "count_select.xlsx")
write.xlsx(P_pos_select, "P_pos_select.xlsx")
write.xlsx(PPP_pos_select, "PPP_pos_select.xlsx")
write.xlsx(PXPP_pos_select, "PXPP_pos_select.xlsx")
write.xlsx(PXPXXP_pos_select, "PXPXXP_pos_select.xlsx")
write.xlsx(PXXPXXP_pos_select, "PXXPXXP_pos_select.xlsx")

