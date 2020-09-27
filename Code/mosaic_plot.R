##### Project: Ultralow risk STO3 ######################################
# mosaic_plot.R
# Functions for Mosaic plots
# Author: Annelie Johansson & Nancy Yu & Adina Iftimi
# Modified last: January 2020 by Annelie Johansson
#########################################################################

library(ggplot2)
library(ggmosaic)

#### FUNCTIONS ####
source("functions_for_mosaic_plot.R")

#### DATA PREPARATION ####

#### Load data -----

# Load numbers for all characteristics/gene modules
table_mosaic <- read.table("table_mosaic.txt", header = TRUE, stringsAsFactors = FALSE)

# Load p-values from Fisher's test
table_pvals <- read.table("table_fishers_test.txt", header = TRUE)

# Make sure we have equal rownames in table_pvals as names of characteristics in table_mosaic
identical(rownames(table_pvals), unique(table_mosaic$Characteristics)) # TRUE

# Change name on categories, to get the same size of all mosaic plots.
table_mosaic$level[which(table_mosaic$Characteristics == "Ki67")] <- c("negative", "positive")
table_mosaic$level[which(table_mosaic$Characteristics == "Size")] <- c("negative", "positive")
table_mosaic$level[which(table_mosaic$Characteristics == "Grade")] <- c("negative", "middle", "positive")


#### FIGURE 1 ####
# Size, Grade, PR, HER2, Ki-67

# Colors in order: light to dark (light for the less aggresive feature)
col_list <- list(PR = c("#e8b9d3", "#e088b9"),
                 HER2 = c("#d1c0d8","#8f5da3"),
                 Ki67 = c("#b8edb8","#0dc478"),
                 Size = c("#b8dced", "#1abce0"), 
                 Grade = c("#f7d1a0","#f7ab47","#e8890d"))
# Control that colors look good
barplot(rep(5,(length(col_list)*2+1)), col = unlist(col_list),
        names = paste0(table_mosaic$Characteristics[1:11], ":", table_mosaic$level[1:11]), las=2)

g_size <- create_gg_mosaic(char = "Size", t = table_mosaic)
g_grade <- create_gg_mosaic(char = "Grade", t = table_mosaic)
g_pr <- create_gg_mosaic(char = "PR", t = table_mosaic)
g_her2 <- create_gg_mosaic(char = "HER2", t = table_mosaic)
g_ki67 <- create_gg_mosaic(char = "Ki67", t = table_mosaic)


pdf(file = "Figure1.pdf", width = 6, height = 6)
multiplot(g_size, 
          g_grade,
          g_pr, 
          g_her2, 
          g_ki67, 
          cols=1)
dev.off()


##### FIGURE 2 ####
# Significant gene modules: AURKA, AKT-MTOR, ERBB2, IGF1, PIK3CA, PTEN, IMMUNE1, STAT1

# Change name on AKT_MTOR to AKT
table_mosaic$Characteristics[which(table_mosaic$Characteristics == "AKT_MTOR")] <- "AKT"
rownames(table_pvals)[which(rownames(table_pvals) == "AKT_MTOR")] <- "AKT"

# Colors in order: light to dark (light for the less aggresive feature)
col_list <- list(AURKA = c("azure2", "azure4"),
                 AKT = c("lightblue1","royalblue1"),
                 ERBB2 = c("#b8edb8","#0dc478"),
                 IGF1 = c("mistyrose", "lightcoral"),
                 PIK3CA = c("khaki1","gold2"),
                 PTEN = c("darkolivegreen1","darkolivegreen4"),
                 IMMUNE1 = c("peachpuff", "sienna1"),
                 STAT1 =  c("plum2", "mediumorchid")
                 )
# Control that colors looks good
modules <- c("AURKA", "AKT", "ERBB2", "IGF1", "PIK3CA", "PTEN", "IMMUNE1", "STAT1")
table_mosaic_temp <- table_mosaic[which(table_mosaic$Characteristics %in% modules),]
barplot(rep(5,(length(col_list)*2)), col = unlist(col_list),
        names = paste0(table_mosaic_temp$Characteristics, ":", table_mosaic_temp$level), las = 2)

g_aurka <- create_gg_mosaic(char = "AURKA", t = table_mosaic)
g_akt_mtor <- create_gg_mosaic(char = "AKT", t = table_mosaic)
g_erbb2 <- create_gg_mosaic(char = "ERBB2", t = table_mosaic)
g_igf1 <- create_gg_mosaic(char = "IGF1",  t = table_mosaic)
g_pik3ca <- create_gg_mosaic(char = "PIK3CA",  t = table_mosaic)
g_pten <- create_gg_mosaic(char = "PTEN",  t = table_mosaic)
g_immune1 <- create_gg_mosaic(char = "IMMUNE1", t = table_mosaic)
g_stat1 <- create_gg_mosaic(char = "STAT1", t = table_mosaic)

pdf(file = "Figure2.pdf", width = 12, height = 4.8)
multiplot(g_aurka,
          g_erbb2,
          g_pik3ca,
          g_immune1,
          # next column
          g_akt_mtor,
          g_igf1,
          g_pten, 
          g_stat1, 
          cols=2)
dev.off()
