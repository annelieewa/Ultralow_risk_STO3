###################################################################################
# Manuscript: Clinical and Molecular Characteristics of ER-Positive Ultralow Risk
#            Breast Cancer Tumors Identified by the 70-Gene Signature
# Author: Annelie Johansson, Nancy Yu & Adina Iftimi
# Modified last: January 2020 by Annelie Johansson
#################################################################################

setwd("/Volumes/Annelie Encrypted/Projects/Ultralow_risk_genes/Analysis/")

library(dplyr)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(reshape2)
library(colorspace)
library(operators) # for the %!in% operation
options(stringsAsFactors = FALSE)

### Load data ------------------------------------------------------------------------------
load("input/20170404patient_array.RData")
patient_data <- select(ord_pat1, AgendiaRUNN, MammaPrint_Result, PAM50,
                       ERstatus, PRstatus, Ki67status, HER2status, grade, size20)
microarray <- ord_dat1

# Load Supplementery Table 4
stable4 <- read.table("output/stable4_heatmap", header = TRUE, sep = "\t")

# Remove probes that are in no gene group
length(which(is.na(stable4$Group)) ) # 400
r <- which(is.na(stable4$Group))
stable4 <- stable4[-r,]

# Remove gene groups with few probes in them
t <- table(stable4$Group)
t <- t[order(t, decreasing = TRUE)]
groups_remove <- names(t[which(t < 5)])
length(groups_remove) # 11
groups_remove
# [1] "WNT_BETA_CATENIN_SIGNALING"
# [2] "ADIPOGENESIS"              
# [3] "HYPOXIA"                   
# [4] "SPERMATOGENESIS"           
# [5] "UNFOLDED_PROTEIN_RESPONSE" 
# [6] "UV_RESPONSE"               
# [7] "ANDROGEN_RESPONSE"         
# [8] "ANGIOGENESIS"              
# [9] "COAGULATION"               
# [10] "MYOGENESIS"                
# [11] "NOTCH_SIGNALING"  

remove <- which(stable4$Group  %in% groups_remove)
stable4 <- stable4[-remove,]

nrow(stable4) # 375

## Fix patient data ---------------------------------------------------------------------------------------------------------

# Select only ER-positive
pat_er <- subset(patient_data, ERstatus == "Positive") 
nrow(pat_er) # 538

# Remove patients with missing information
pat_use <- subset(pat_er, !is.na(size20) & !is.na(grade) & 
                    PRstatus != "Unknown" & HER2status != "Unknown" &
                    Ki67status != "Unknown")
nrow(pat_use) # 495

# Prepare variables, order from poor to good prognosis markers
pat_use$MammaPrint_Result <- factor(pat_use$MammaPrint_Result, levels = c("High Risk", "Low Risk", "Ultra Low Risk"))
pat_use$PAM50 <- factor(pat_use$PAM50, levels = c("Basal", "Her2", "LumB", "LumA", "Normal"))
pat_use$grade <- factor(pat_use$grade, levels = c("3", "2", "1"))
pat_use$size20 <- factor(pat_use$size20, levels = c("1", "0"))
pat_use$PRstatus <- factor(pat_use$PRstatus, levels = c("Negative", "Positive"), labels = c("0", "1"))
pat_use$HER2status <- factor(pat_use$HER2status, levels = c("Positive", "Negative"), labels = c("1", "0"))
pat_use$Ki67status <- factor(pat_use$Ki67status, levels = c("Positive", "Negative"), labels = c("1", "0"))

# Order patient data as we want it in the heatmap
pat_use <- pat_use[order(pat_use$MammaPrint_Result, pat_use$PAM50, pat_use$grade, pat_use$size20,
                         pat_use$PRstatus, pat_use$HER2status, pat_use$Ki67status),]

nrow(pat_use) # 495

## Fix probes data ---------------------------------------------------------------------------------------------------

# 1: Order groups by t-statistics
groups <- unique(stable4$Group)
t_mean <- c()
for(i in 1:length(groups)){
  t_mean <- c(t_mean, mean(stable4$tstat[which(stable4$Group == groups[i])]) )
}
names(t_mean) <- groups
t_mean <- t_mean[order(t_mean)]

stable4$Group <- factor(stable4$Group, levels=names(t_mean))
stable4 <- stable4[order(stable4$Group),]

# 2: Order within each group by t-statistics
new_order <- c()
for(i in 1:length(groups)){
  old_order <- which(stable4$Group == groups[i])
  new_order <- old_order[order(stable4$tstat[old_order])]
  stable4[old_order,] <- stable4[new_order,]
}

stable4$Group <- factor(stable4$Group, levels=unique(stable4$Group))

table(stable4$Group)
# HISTONE 
# 45 
# MYC_TARGETS 
# 9 
# REACTIVE_OXIGEN_SPECIES_PATHWAY 
# 7 
# CELL_CYCLE 
# 68 
# PI3K_AKT_MTOR 
# 14 
# DNA_REPAIR 
# 17 
# IMMUNE 
# 60 
# APOPTOSIS 
# 20 
# METABOLIC 
# 51 
# PROTEIN_SECRETION 
# 7 
# ESTROGEN_RESPONSE 
# 30 
# P53_PATHWAY 
# 6 
# HOMEOBOX 
# 13 
# EPITHELIAL_STRUCTURE 
# 13 
# KRAS_SIGNALING 
# 7 
# EPITHELIAL_MESENCHYMAL_TRANSITION 
# 8 

## Prepare gene expression data ---------------------------------------------------------------------------------------------------------

# Match gene expression data with patient data
dat_use <- microarray[, match(as.character(pat_use$AgendiaRUNN), colnames(microarray))]  
identical(colnames(dat_use), as.character(pat_use$AgendiaRUNN)) # TRUE, same order

# Order gene expression data after probes to be used for the heatmap
probes_use <- stable4$ProbeID 
all(probes_use %in% rownames(dat_use)) # TRUE
dat_use <- dat_use[which(rownames(dat_use) %in% probes_use), ]
dat_use <- dat_use[probes_use, ]
identical(rownames(dat_use), probes_use) # TRUE

dim(dat_use) # 375 495

## Convert GE data to Z-scores                      
# Method: Heatmap was generated by computing Z scores across the cohort. Values > abs(3) were converted to -3 and +3, respectively.
dat_zscores <- t(apply(dat_use, 1, function(x) { (x - mean(x, na.rm = T)) / sd(x, na.rm = T) }) )
summary(as.vector(t(dat_zscores)))
#     Min.    1st Qu.    Median     Mean    3rd Qu.    Max.
# -10.79224  -0.69811  -0.07942   0.00000   0.62637   6.87763 
dat_zscores_cutoff <- replace(dat_zscores, dat_zscores > 3, 3)
dat_zscores_cutoff <- replace(dat_zscores_cutoff, dat_zscores < -3, -3)
summary(as.vector(t(dat_zscores_cutoff)))
#     Min.   1st Qu.   Median     Mean    3rd Qu.   Max.
# -3.000000 -0.69811 -0.07942 -0.00273   0.62637  3.00000 

## Colors for the patient data -----------------------------------------------------------------------------------------------------------------
mammaprint_colors <- c("#212d28", "#137a44","#5fce9e")
luminal_colors <- c("#d92405", "#6c7197", "#eac124", "#3563eb","#739211")
size_colors <- c("#3e2b7c", "#b8a3ff")
grade_colors <- c("#512704", "#e57504", "#f9b877")
PR_colors <- c("lightpink3", "lightpink") 
HER2_colors <- c("lemonchiffon3", "lemonchiffon")
Ki67_colors <- c("slategray3", "slategray1")
white_color <- c("white")

## This is for the clinical labels for each patient (Color legend on the top for the columns) 
patient_colororder <- cbind( mammaprint_colors[factor(pat_use$MammaPrint_Result)],
                             rep(white_color, nrow(pat_use)),
                             luminal_colors[factor(pat_use$PAM50)],
                             grade_colors[factor(pat_use$grade)],
                             Ki67_colors[factor(pat_use$Ki67status)],
                             PR_colors[factor(pat_use$PRstatus)],
                             size_colors[factor(pat_use$size20)],
                             HER2_colors[factor(pat_use$HER2status)]
                             )

# Control of colors
barplot(rep(5,19),
        col = c(mammaprint_colors, luminal_colors, size_colors, grade_colors, PR_colors, HER2_colors, Ki67_colors),
        names.arg = c(levels(pat_use$MammaPrint_Result), levels(pat_use$PAM50), levels(pat_use$size20), levels(pat_use$grade),
        levels(pat_use$PRstatus), levels(pat_use$HER2status), levels(pat_use$Ki67status)), las = 2, cex.names = 0.6)

## Colors for the gene categories -----------------------------------------------------------------------------------------------------------------

names(table(stable4$Group)) 
# [1] "HISTONE"                          
# [2] "MYC_TARGETS"                      
# [3] "REACTIVE_OXIGEN_SPECIES_PATHWAY"  
# [4] "CELL_CYCLE"                       
# [5] "PI3K_AKT_MTOR"                    
# [6] "DNA_REPAIR"                       
# [7] "IMMUNE"                           
# [8] "APOPTOSIS"                        
# [9] "METABOLIC"                        
# [10] "PROTEIN_SECRETION"                
# [11] "ESTROGEN_RESPONSE"                
# [12] "P53_PATHWAY"                      
# [13] "HOMEOBOX"                         
# [14] "EPITHELIAL_STRUCTURE"             
# [15] "KRAS_SIGNALING"                   
# [16] "EPITHELIAL_MESENCHYMAL_TRANSITION"

col <- rep("", length(unique(stable4$Group)))
col[which(unique(stable4$Group) == "HISTONE")] <- "darkgoldenrod2" 
col[which(unique(stable4$Group) == "MYC_TARGETS")] <- "darkseagreen4" 
col[which(unique(stable4$Group) == "REACTIVE_OXIGEN_SPECIES_PATHWAY")] <- "deepskyblue" 
col[which(unique(stable4$Group) == "CELL_CYCLE")] <- "darkorchid1" 
col[which(unique(stable4$Group) == "PI3K_AKT_MTOR")] <- "slateblue" 
col[which(unique(stable4$Group) == "DNA_REPAIR")] <- "darkorange2" 
col[which(unique(stable4$Group) == "IMMUNE")] <- "turquoise3" 
col[which(unique(stable4$Group) == "APOPTOSIS")] <- "indianred1" 
col[which(unique(stable4$Group) == "METABOLIC")] <- "steelblue4" 
col[which(unique(stable4$Group) == "PROTEIN_SECRETION")] <- "lightsalmon"
col[which(unique(stable4$Group) == "ESTROGEN_RESPONSE")] <- "deeppink2"
col[which(unique(stable4$Group) == "P53_PATHWAY")] <- "seagreen3" 
col[which(unique(stable4$Group) == "HOMEOBOX")] <- "tomato" 
col[which(unique(stable4$Group) == "EPITHELIAL_STRUCTURE")] <- "skyblue1" 
col[which(unique(stable4$Group) == "KRAS_SIGNALING")] <- "gold2"
col[which(unique(stable4$Group) == "EPITHELIAL_MESENCHYMAL_TRANSITION")] <- "yellowgreen"

# Control of colors
barplot(rep(5,length(unique(stable4$Group))), col = col, names.arg = names(table(stable4$Group)) , las = 2, cex.names = 0.6)

names(col) <- names(table(stable4$Group))
annot_colors <- col
annot_colororder <- as.matrix(annot_colors[as.factor(stable4$Group)])

## HEATMAP  ----------------------------------------------------------------------------------------------------------
source("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

heatcolors <- colorRampPalette(c("blue", "black", "yellow"))(256)

pdf("output/Figure3.pdf", useDingbats = F, width = 15.5, height = 13)

heatmap.3(as.matrix(dat_zscores_cutoff), col = heatcolors, scale = "none", main = "",
          margins = c(9,20),
          ColSideColorsSize = 5, ColSideColors=patient_colororder,
          Colv = F, Rowv = F,
          RowSideColors = t(annot_colororder), RowSideColorsSize = 1,
          dendrogram = "none", key = TRUE, keysize = 0.7, symkey = F, symbreaks = T,
          symm = T, density.info = "none", trace = "none", 
          labRow = FALSE, labCol = FALSE)

legend(x = 0.85, y = 0.7, xpd = TRUE, 
       legend = c("HER2 status", "Positive", "Negative", "",
                  "Tumor size", ">20mm", "≤20mm", "",
                  "PR status", "Negative", "Positive", "",
                  "Ki-67 status", "Med/High", "Low", "",
                  "Tumor grade", "Grade3", "Grade2", "Grade1", "",
                  "Molecular subtype", "Basal", "HER2-enriched",  "Luminal B", "Luminal A", "Normal-like", ""),
       fill = c("white", HER2_colors, "white",
                "white", size_colors, "white",
                "white", PR_colors, "white",
                "white", Ki67_colors, "white",
                "white",  grade_colors,"white", 
                "white", luminal_colors, "white"), 
       border = FALSE, bty = "n", y.intersp = 1, cex = 1)

legend(x = 0.85, y = 0.15, xpd = TRUE,
       legend = c("70-gene signature", "High Risk", "Low Risk", "Ultralow Risk", ""),
       fill = c("white", mammaprint_colors,  "white"), 
       border = FALSE, bty = "n", y.intersp = 1, cex = 1)

text(x = 0.1, y = 0.25, xpd = TRUE, pos = 2,
    "Histone\n\n
    MYC signaling
    ROS\n\n\n\n
    Cell cycle\n\n\n\n
    PI3K/Akt/mTOR pathway\n
    DNA repair\n\n\n\n
    Immune response\n\n\n\n
    Apoptosis\n\n\n\n
    Metabolic processes\n\n\n
    Protein secretion\n
    Estrogen response\n\n
    P53 pathway
    Homeobox\n
    Epithelial structure
    KRAS signaling
    EMT", cex = 0.85, font = 2)

dev.off()
