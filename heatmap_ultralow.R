##### Project: Ultralow risk STO3 ###############################################
# heatmap_ultralow.R
# Figure 3
# Author: Annelie Johansson, Nancy Yu & Adina Iftimi
# Modified last: January 2020 by Annelie Johansson
#################################################################################

library(RColorBrewer)
library(gplots)
library(ggplot2)
library(reshape2)
library(colorspace)
library(operators) # for the %!in% operation
library(dplyr)
options(stringsAsFactors = FALSE)

### Load data ----

# Load data file containing: annot1, ord_dat1, ord_pat1
load("20170404patient_array.RData")

# Load Supplementery Table 2
eTable2 <- read.table("eTable2_full.txt", header=TRUE, sep="\t")

# Remove probes that are in no gene group
length(which(is.na(eTable2$Group)) ) # 400
r <- which(is.na(eTable2$Group))
eTable2 <- eTable2[-r,]

# Remove gene groups with few probes in them
t <- table(eTable2$Group)
t <- t[order(t, decreasing = TRUE)]
groups_remove <- names(t[which(t < 5)])
length(groups_remove) # 11
groups_remove
# [1] "WNT_BETA_CATENIN_SIGNALING" "ADIPOGENESIS"               "HYPOXIA"                    "SPERMATOGENESIS"           
# [5] "UNFOLDED_PROTEIN_RESPONSE"  "UV_RESPONSE"                "ANDROGEN_RESPONSE"          "ANGIOGENESIS"              
# [9] "COAGULATION"                "MYOGENESIS"                 "NOTCH_SIGNALING"

remove <- which(eTable2$Group  %in% groups_remove)
eTable2 <- eTable2[-remove,]

nrow(eTable2) # 375

## Prepare patient data ----
ord_pat1 <- select(ord_pat1, AgendiaRUNN, ERstatus, MammaPrint_Result, PAM50, grade, size20, mitosis)

# Select only ER-positive patients
pat_er <- subset(ord_pat1, ERstatus == "Positive") 
nrow(pat_er) # 538

# Prepare variables
pat_er$MammaPrint_Result <- factor(pat_er$MammaPrint_Result, levels=c("High Risk", "Low Risk", "Ultra Low Risk"))
pat_er$PAM50 <- factor(pat_er$PAM50, levels=c("Basal", "Her2", "LumB", "LumA", "Normal"))

# Remove patients with missing info of grade and size20
pat_use <- subset(pat_er, !is.na(size20) & !is.na(grade))

# Order patient data as we want it in the heatmap
pat_use <- pat_use[order(pat_use$MammaPrint_Result, pat_use$PAM50, -pat_use$grade, -pat_use$size20, -pat_use$mitosis),]

nrow(pat_use) # 526

## Prepare probes data ----

# 1: Order groups by t-statistics
groups <- unique(eTable2$Group)
t_mean <- c()
for(i in 1:length(groups)){
  t_mean <- c(t_mean, mean(eTable2$tstat[which(eTable2$Group == groups[i])]) )
}
names(t_mean) <- groups
t_mean <- t_mean[order(t_mean)]

eTable2$Group <- factor(eTable2$Group, levels=names(t_mean))
eTable2 <- eTable2[order(eTable2$Group),]

# 2: Order within each group by t-statistics
new_order <- c()
for(i in 1:length(groups)){
  old_order <- which(eTable2$Group == groups[i])
  new_order <- old_order[order(eTable2$tstat[old_order])]
  eTable2[old_order,] <- eTable2[new_order,]
}
eTable2$Group <- factor(eTable2$Group, levels = unique(eTable2$Group))

## Prepare gene expression data ----

# Match gene expression data with patient data
dat_use <- ord_dat1[, match(as.character(pat_use$AgendiaRUNN), colnames(ord_dat1))]  
identical(colnames(dat_use), as.character(pat_use$AgendiaRUNN)) # TRUE, same order

# Order gene expression data after probes to be used for the heatmap
probes_use <- eTable2$ProbeID 
all(probes_use %in% rownames(dat_use)) # TRUE
dat_use <- dat_use[which(rownames(dat_use) %in% probes_use), ]
dat_use <- dat_use[probes_use, ]
identical(rownames(dat_use), probes_use) # TRUE

dim(dat_use) # 375 x 526

## Convert gene expression data to Z-scores                      
# Values > abs(3) were converted to -3 and +3, respectively.
dat_zscores <- t(apply(dat_use, 1, function(x) { (x - mean(x, na.rm = T)) / sd(x, na.rm = T) }) )
summary(as.vector(t(dat_zscores)))
#     Min.  1st Qu.   Median     Mean  3rd Qu. Max.
# -10.87446  -0.69776  -0.07851   0.00000   0.62685   6.92652  
dat_zscores_cutoff <- replace(dat_zscores, dat_zscores > 3, 3)
dat_zscores_cutoff <- replace(dat_zscores_cutoff, dat_zscores < -3, -3)
summary(as.vector(t(dat_zscores_cutoff)))
#     Min.  1st Qu.   Median     Mean  3rd Qu. Max.
# -3.000000 -0.697758 -0.078506 -0.002741  0.626853  3.000000  


## Colors for the patient data ----
luminal_colors <- c("#d92405", "#6c7197", "#eac124", "#3563eb","#739211") # Basal red, HER2 grey, LumB yellow, LumA blue, Normal green
size_colors <- c("#b8a3ff", "#3e2b7c")
grade_colors <- c("#f9b877", "#e57504","#512704")
mammaprint_colors <- c("#212d28", "#137a44","#5fce9e")
white_color <- c("white")

# This is for the clinical labels for each patient (color legend on the top for the columns) 
patient_colororder <- cbind( mammaprint_colors[factor(pat_use$MammaPrint_Result)],
                             rep(white_color, nrow(pat_use)),
                             luminal_colors[factor(pat_use$PAM50)],
                             grade_colors[factor(pat_use$grade)],
                             size_colors[factor(pat_use$size20)])

## Colors for the gene categories ----
names(table(eTable2$Group)) 
# [1] "HISTONE"                           "MYC_TARGETS"                       "REACTIVE_OXIGEN_SPECIES_PATHWAY"  
# [4] "CELL_CYCLE"                        "PI3K_AKT_MTOR"                     "DNA_REPAIR"                       
# [7] "IMMUNE"                            "APOPTOSIS"                         "METABOLIC"                        
# [10] "PROTEIN_SECRETION"                 "ESTROGEN_RESPONSE"                 "P53_PATHWAY"                      
# [13] "HOMEOBOX"                          "EPITHELIAL_STRUCTURE"              "KRAS_SIGNALING"                   
# [16] "EPITHELIAL_MESENCHYMAL_TRANSITION"

col <- rep("", length(unique(eTable2$Group)))
col[which(unique(eTable2$Group) == "HISTONE")] <- "darkgoldenrod2" 
col[which(unique(eTable2$Group) == "MYC_TARGETS")] <- "darkseagreen4" 
col[which(unique(eTable2$Group) == "REACTIVE_OXIGEN_SPECIES_PATHWAY")] <- "deepskyblue" 
col[which(unique(eTable2$Group) == "CELL_CYCLE")] <- "darkorchid1" 
col[which(unique(eTable2$Group) == "PI3K_AKT_MTOR")] <- "slateblue" 
col[which(unique(eTable2$Group) == "DNA_REPAIR")] <- "darkorange2" 
col[which(unique(eTable2$Group) == "IMMUNE")] <- "turquoise3" 
col[which(unique(eTable2$Group) == "APOPTOSIS")] <- "indianred1" 
col[which(unique(eTable2$Group) == "METABOLIC")] <- "steelblue4" 
col[which(unique(eTable2$Group) == "PROTEIN_SECRETION")] <- "lightsalmon"
col[which(unique(eTable2$Group) == "ESTROGEN_RESPONSE")] <- "deeppink2"
col[which(unique(eTable2$Group) == "P53_PATHWAY")] <- "seagreen3" 
col[which(unique(eTable2$Group) == "HOMEOBOX")] <- "tomato" 
col[which(unique(eTable2$Group) == "EPITHELIAL_STRUCTURE")] <- "skyblue1" 
col[which(unique(eTable2$Group) == "KRAS_SIGNALING")] <- "gold2"
col[which(unique(eTable2$Group) == "EPITHELIAL_MESENCHYMAL_TRANSITION")] <- "yellowgreen"

# Control of colors
barplot(rep(5,length(unique(eTable2$Group))), col = col, names.arg = names(table(eTable2$Group)) , las = 2, cex.names = 0.6)

names(col) <- names(table(eTable2$Group))
annot_colors <- col
annot_colororder <- as.matrix(annot_colors[as.factor(eTable2$Group)])

## HEATMAP  ----
source("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

heatcolors <- colorRampPalette(c("blue", "black", "yellow"))(256)

pdf("Figure2_heatmap.pdf", useDingbats = F, width = 17, height = 10)

heatmap.3(as.matrix(dat_zscores_cutoff), col = heatcolors, scale = "none", main = "",
          margins = c(9,20), ColSideColorsSize = 5,  
          ColSideColors = patient_colororder, 
          Colv = F,
          Rowv = F,
          RowSideColors = t(annot_colororder), RowSideColorsSize = 1,
          dendrogram = "none", key=TRUE, keysize = 0.7, symkey = F, symbreaks = T,
          symm = T, density.info = "none", trace = "none", 
          labRow = FALSE, labCol = FALSE)

legend(x = 0.865, y = 0.75, xpd = TRUE, 
       legend=c("70-gene signature", "High Risk", "Low Risk", "Ultralow Risk", "",
                "PAM50 subtype", "Basal", "Her2",  "LumB", "LumA", "Normal", "",
                "Tumor grade", "Grade1", "Grade2", "Grade3", "",
                "Tumor size", "<=20mm", ">20mm", ""),
       fill=c("white", mammaprint_colors, "white", 
              "white", luminal_colors, "white", 
              "white", grade_colors, "white", 
              "white", size_colors, "white"), 
       border = FALSE, bty = "n", y.intersp = 1, cex = 1)

text(x = 0.1, y = 0.255, xpd = TRUE, pos = 2,
    "Histone\n\n
    MYC signaling
    Reactive oxygen species\n\n\n\n
    Cell cycle\n\n\n\n
    PI3K/Akt/mTOR pathway\n
    DNA repair\n\n\n\n
    Immune response\n\n\n\n
    Apoptosis\n\n\n
    Metabolic processes\n\n\n
    Protein secretion\n
    Estrogen response\n\n
    P53 pathway
    Homeobox\n
    Epithelial structure
    KRAS signaling
    Epithelial-mesenchymal transition", cex = 0.65, font = 2)

dev.off()
