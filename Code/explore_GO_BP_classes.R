###################################################################################
# Manuscript: Clinical and Molecular Characteristics of ER-Positive Ultralow Risk
#            Breast Cancer Tumors Identified by the 70-Gene Signature
# Author: Annelie Johansson 
# Modified last: January 2020 by Annelie Johansson
#################################################################################

setwd("/Volumes/Annelie Encrypted/Projects/Ultralow_risk_genes/Analysis/")

#### In this script we try to find GO BP classes with similar functions as the Hallmarks ####
library(sigPathway)

## Hallmark genes and GO biological processes genes from here:
# http://software.broadinstitute.org/gsea/msigdb/index.jsp

# Read in Hallmarks gene sets
hallmarks <- gmtToG("input/h.all.v6.2.symbols.gmt.txt") 
hnames <- sapply(hallmarks, function(x) {x$src})
hprobes <- sapply(hallmarks, function(x) {x$probes})
length(hprobes) # 50
length(hnames) # 50 

# We put together ESTROGEN_RESPONSE_EARLY and ESTROGEN_RESPONSE_LATE
hnames[grep("ESTROGEN", hnames)] <- "HALLMARK_ESTROGEN_RESPONSE"

# We put together KRAS_SIGNALING_DN and KRAS_SIGNALING_UP
hnames[grep("KRAS", hnames)] <- "HALLMARK_KRAS_SIGNALING"

# We put together MYC_TARGETS_V1 and MYC_TARGETS_V2
hnames[grep("MYC", hnames)] <- "HALLMARK_MYC_TARGETS"

# We put together UV_RESPONSE_DN and UV_RESPONSE_UP
hnames[grep("UV_RESPONSE", hnames)] <- "HALLMARK_UV_RESPONSE"

length(unique(hnames)) # 46
hnames <- hnames[-which(duplicated(hnames) == TRUE)]

hnames <- sapply(hnames, function(x) { strsplit(x[1], "HALLMARK_")[[1]][2] }, USE.NAMES = FALSE)
hnames <- hnames[order(hnames)]

#### This script aims to explore/identify GO BP gene groups that correspons to the Hallmarks ####

## Read in GO biological process
c5.proc <- gmtToG("input/c5.bp.v6.2.symbols.gmt.txt") 
c5names <- sapply(c5.proc, function(x) { x$src })
c5probes <- sapply(c5.proc, function(x) { x$probes })
length(c5names) # 4436 gene sets

# Fill up these:
GO_use <- c()
GO_use_save <- c()
c5names_use <- c()
c5probes_use <- NULL

## Function, to check if genes for what we are searching for are in the 'big group(s)'
# Useful when searching for matching GO BP groups
# Returns a list with data-frames a and b:
# - a shows the groups that do NOT have all genes in the big group(s)
# - b shows the groups that have all genes in the big group(s)
fun <- function(big_group, search){
  genes_in_chosen <- c()
  for(i in 1:length(big_group)){
    genes_in_chosen <- c(genes_in_chosen, unlist(c5probes[which(c5names == big_group[i])]) )
  }
  c5_others <- c5names[grep(search, c5names)]
  res <- matrix(nrow=length(grep(search, c5names)), ncol=3)
  for(i in 1:length(grep(search, c5names)) ){
    c5_now_genes <- unlist(c5probes[which(c5names == c5_others[i])])
    a <- length(c5_now_genes)
    b <- length(which(c5_now_genes %in% genes_in_chosen))
    res[i,] <- c(c5_others[i], a, b)
  }
  return(list(a=res[which(res[,2] != res[,3]),], b=res[which(res[,2] == res[,3]),]))
}

## * hnames[1] ADIPOGENESIS ----------------------------------------------------------------------------------------
# Genes up-regulated during adipocyte differentiation (adipogenesis).
length( c5names[grep("ADIPO", c5names)] ) # 1
c5names[grep("ADIPO", c5names)] # GO_ADIPOSE_TISSUE_DEVELOPMENT

GO_use <- c("GO_ADIPOSE_TISSUE_DEVELOPMENT")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "ADIPOGENESIS")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[2] ALLOGRAFT_REJECTION --------------------------------------------------------------------------------
# Genes up-regulated during transplant rejection.
# Graft rejection is an immunologic destruction of transplanted tissues or organs between two members or strains of
# a species differing at the major histocompatibility complex for that species (i.e., HLA in man and H-2 in the mouse).
length( c5names[grep("REJECT", c5names)] ) # 0
length( c5names[grep("TRANSPLANT", c5names)] ) # 0
length( c5names[grep("GRAFT", c5names)] ) # 0

GO_use <- c("NONE")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "ALLOGRAFT_REJECTION")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[3] ANDROGEN_RESPONSE --------------------------------------------------------------------------------------
# Genes defining response to androgens.
length( c5names[grep("ANDROGEN", c5names)] ) # 5

f <- fun(big_group = c("GO_ANDROGEN_RECEPTOR_SIGNALING_PATHWAY", "GO_REGULATION_OF_ANDROGEN_RECEPTOR_SIGNALING_PATHWAY"), search = "ANDROGEN")
nrow(f$a) # 2 --> This is OK
nrow(f$b) # 3

GO_use <- c("GO_ANDROGEN_RECEPTOR_SIGNALING_PATHWAY", "GO_REGULATION_OF_ANDROGEN_RECEPTOR_SIGNALING_PATHWAY")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "ANDROGEN_RESPONSE")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[4] ANGIOGENESIS --------------------------------------------------------------------------------------------
# Genes up-regulated during formation of blood vessels (angiogenesis).
length( c5names[grep("ANGIOGENESIS", c5names)] ) # 9

f <- fun(big_group = c("GO_ANGIOGENESIS"), search = "ANGIOGENESIS")
nrow(f$a) # 1 --> OK
nrow(f$b) # 8

GO_use <- c("GO_ANGIOGENESIS")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "ANGIOGENESIS")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[5] APICAL_JUNCTION --------------------------------------------------------------------------------------
# Genes encoding components of apical junction complex.
# Cell junction
length( c5names[grep("APICAL", c5names)] ) # 3
length( c5names[grep("JUNCTION", c5names)] ) # 17
length( c5names[grep("APICAL_JUNCTION", c5names)] ) # 1
length( c5names[grep("CELL_JUNCTION", c5names)] ) # 6

f <- fun(big_group = c("GO_APICAL_JUNCTION_ASSEMBLY"), search = "APICAL")
nrow(f$a) # 2 --> OK
nrow(f$b) # 1

f <- fun(big_group = c("GO_CELL_JUNCTION_ASSEMBLY", "GO_REGULATION_OF_CELL_JUNCTION_ASSEMBLY"), search = "JUNCTION")
nrow(f$a) # 6 --> OK
nrow(f$b) # 11 -->  OBS - APICAL_JUNCTION_ASSEMBLY covered in here! So it is enough to use these.

GO_use <- c("GO_CELL_JUNCTION_ASSEMBLY", "GO_REGULATION_OF_CELL_JUNCTION_ASSEMBLY")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "APICAL_JUNCTION")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[6] APICAL_SURFACE -------------------------------------------------------------------------------------
# Surface of an epithelial cell that faces the body surface
# Genes encoding proteins over-represented on the apical surface of epithelial cells,
# e.g., important for cell polarity (apical area).
length( c5names[grep("SURFACE", c5names)] ) # 7 --> None really good
length( intersect(c5names[grep("EPITHELIAL", c5names)], c5names[grep("SURFACE", c5names)]) ) # 0
length( c5names[grep("CELL_POLA", c5names)] ) # 8

f <- fun(big_group = c("GO_ESTABLISHMENT_OR_MAINTENANCE_OF_CELL_POLARITY", "GO_REGULATION_OF_ESTABLISHMENT_OR_MAINTENANCE_OF_CELL_POLARITY"), search = "CELL_POLA")
nrow(f$a) # 1 --> OK
nrow(f$b) # 7

GO_use <- c("GO_ESTABLISHMENT_OR_MAINTENANCE_OF_CELL_POLARITY", "GO_REGULATION_OF_ESTABLISHMENT_OR_MAINTENANCE_OF_CELL_POLARITY")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "APICAL_SURFACE")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[7] APOPTOSIS -------------------------------------------------------------------------------------------
# Genes mediating programmed cell death (apoptosis) by activation of caspases.
length(grep("APOPTOTIC_SIGNALING", c5names)) # 37
length(grep("CELL_DEATH", c5names)) # 12

f <- fun(big_group = c("GO_APOPTOTIC_SIGNALING_PATHWAY", "GO_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY"), search = "APOPTOTIC_SIGNALING")
nrow(f$a) # 0 
nrow(f$b) # 37

f <- fun(big_group = c("GO_CELL_DEATH", "GO_REGULATION_OF_CELL_DEATH"), search = "CELL_DEATH")
nrow(f$a) # 0 
nrow(f$b) # 12

# GO_CELL_DEATH covers all genes in APOPTOTIC_SIGNALING!
f <- fun(big_group = c("GO_CELL_DEATH", "GO_REGULATION_OF_CELL_DEATH"), search = "APOPTOTIC_SIGNALING")
nrow(f$a) # 0 = all!
nrow(f$b) # 37
# So we only need to use GO_CELL_DEATH and GO_REGULATION_OF_CELL_DEATH

GO_use <- c("GO_CELL_DEATH", "GO_REGULATION_OF_CELL_DEATH")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "APOPTOSIS")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[8] BILE_ACID_METABOLISM ----------------------------------------------------------------------------------
# Genes involve in metabolism of bile acids and salts
length( c5names[grep("BILE_ACID_MET", c5names)] ) # 1

GO_use <- c("GO_BILE_ACID_METABOLIC_PROCESS")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "BILE_ACID_METABOLISM")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[9] CHOLESTEROL_HOMEOSTASIS -------------------------------------------------------------------------------
# Genes involved in cholesterol homeostasis.
# membrane maintenance, transmembrane signaling pathways, hormones, bile acids, vit D. in liver. 
length( c5names[grep("CHOLESTEROL_HOMEOSTASIS", c5names)] ) # 1

GO_use <- c("GO_REGULATION_OF_CHOLESTEROL_HOMEOSTASIS")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "CHOLESTEROL_HOMEOSTASIS")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[10] COAGULATION --------------------------------------------------------------------------------------------
# Genes encoding components of blood coagulation system; also up-regulated in platelets.
length( c5names[grep("COAG", c5names)] ) # 5

f <- fun(big_group = c("GO_REGULATION_OF_COAGULATION"), search = "COAGULATION")
nrow(f$a) # 2 --> OK
nrow(f$b) # 3

GO_use <- c("GO_REGULATION_OF_COAGULATION")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "COAGULATION")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[11] COMPLEMENT ---------------------------------------------------------------------------------------------
# Genes encoding components of the complement system,
# which is part of the innate immune system.
# 'complement cascade', enhancs antibodies, phagocytes etc. 
length( c5names[grep("COMPLEMENT", c5names)] ) # 2
length( intersect(c5names[grep("COMPLEMENT", c5names)], c5names[grep("COMPLEX", c5names)]) ) # 0

f <- fun(big_group = c("GO_COMPLEMENT_ACTIVATION"), search = "COMPLEMENT")
nrow(f$a) # 0
nrow(f$b) # 2

GO_use <- c("GO_COMPLEMENT_ACTIVATION")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "COMPLEMENT")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[12] DNA_REPAIR -----------------------------------------------------------------------------------------------
# Genes involved in DNA repair.
length( c5names[grep("DNA_REPAIR", c5names)] ) # 6

f <- fun(big_group = c("GO_DNA_REPAIR", "GO_REGULATION_OF_DNA_REPAIR"), search = "DNA_REPAIR")
nrow(f$a) # 0
nrow(f$b) # 6

GO_use <- c("GO_DNA_REPAIR", "GO_REGULATION_OF_DNA_REPAIR")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "DNA_REPAIR")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[13] E2F_TARGETS ----------------------------------------------------------------------------------------------
# Genes encoding cell cycle related targets of E2F transcription factors.
# G1/S transition, cell cycle regulation
length( c5names[grep("E2", c5names)] ) # 0
length( c5names[grep("G1_S", c5names)] ) # 6

f <- fun(big_group = c("GO_CELL_CYCLE_G1_S_PHASE_TRANSITION", "GO_REGULATION_OF_CELL_CYCLE_G1_S_PHASE_TRANSITION"), search = "G1_S")
nrow(f$a) # 0
nrow(f$b) # 6

GO_use <- c("GO_CELL_CYCLE_G1_S_PHASE_TRANSITION", "GO_REGULATION_OF_CELL_CYCLE_G1_S_PHASE_TRANSITION")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "E2F_TARGETS")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[14] EPITHELIAL_MESENCHYMAL_TRANSITION ------------------------------------------------------------------------
# Genes defining epithelial-mesenchymal transition, as in wound healing,
# fibrosis and metastasis.
length( c5names[grep("EPITHELIAL_TO_MESENCHYMAL_TRANSITION", c5names)] ) # 6

f <- fun(big_group = c("GO_EPITHELIAL_TO_MESENCHYMAL_TRANSITION", "GO_REGULATION_OF_EPITHELIAL_TO_MESENCHYMAL_TRANSITION"), search = "EPITHELIAL_TO_MESENCHYMAL_TRANSITION")
nrow(f$a) # 0
nrow(f$b) # 6

GO_use <- c("GO_EPITHELIAL_TO_MESENCHYMAL_TRANSITION", "GO_REGULATION_OF_EPITHELIAL_TO_MESENCHYMAL_TRANSITION")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "EPITHELIAL_MESENCHYMAL_TRANSITION")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[15] ESTROGEN_RESPONSE ------------------------------------------------------------------------------------------
# Genes defining early response to estrogen.
# Genes defining late response to estrogen.
length( c5names[grep("ESTROGEN", c5names)] ) # 7

f <- fun(big_group = c("GO_RESPONSE_TO_ESTROGEN", "GO_INTRACELLULAR_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY", "GO_REGULATION_OF_INTRACELLULAR_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY"), search = "ESTROGEN")
nrow(f$a) # 2 --> OK
nrow(f$b) # 5

GO_use <- c("GO_RESPONSE_TO_ESTROGEN",
            "GO_INTRACELLULAR_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY",
            "GO_REGULATION_OF_INTRACELLULAR_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "ESTROGEN_RESPONSE")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[16] FATTY_ACID_METABOLISM -------------------------------------------------------------------------------------
# Genes encoding proteins involved in metabolism of fatty acids.
# creates energy + important molecules, ATP
length( c5names[grep("FATTY_ACID_MET", c5names)] ) # 8

f <- fun(big_group = c("GO_FATTY_ACID_METABOLIC_PROCESS", "GO_REGULATION_OF_FATTY_ACID_METABOLIC_PROCESS"), search = "FATTY_ACID_MET")
nrow(f$a) # 0
nrow(f$b) # 8

GO_use <- c("GO_FATTY_ACID_METABOLIC_PROCESS", "GO_REGULATION_OF_FATTY_ACID_METABOLIC_PROCESS")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "FATTY_ACID_METABOLISM")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[17] G2M_CHECKPOINT --------------------------------------------------------------------------------------------
# Genes involved in the G2/M checkpoint, as in progression through the cell division cycle.
length( c5names[grep("G2_M", c5names)] ) # 5

f <- fun(big_group = c("GO_CELL_CYCLE_G2_M_PHASE_TRANSITION", "GO_REGULATION_OF_CELL_CYCLE_G2_M_PHASE_TRANSITION")
         , search = "G2_M")
nrow(f$a) # 0
nrow(f$b) # 5

GO_use <- c("GO_CELL_CYCLE_G2_M_PHASE_TRANSITION", "GO_REGULATION_OF_CELL_CYCLE_G2_M_PHASE_TRANSITION")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "G2M_CHECKPOINT")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[18] GLYCOLYSIS -------------------------------------------------------------------------------------------
# Genes encoding proteins involved in glycolysis and gluconeogenesis.
# Glycolysis is the process of breaking down glucose
# Gluconeogenesis is a metabolic pathway that results in the generation of glucose
length( c5names[grep("GLYCOLYSIS", c5names)] ) # 0
length( c5names[grep("GLUCONEOGENESIS", c5names)] ) # 3
length( c5names[grep("GLUCOSE", c5names)] ) # 17

f <- fun(big_group = c("GO_REGULATION_OF_GLUCONEOGENESIS"), search = "GLUCONEOGENESIS")
nrow(f$a) # 0
nrow(f$b) # 3

f <- fun(big_group = c("GO_GLUCOSE_METABOLIC_PROCESS", "GO_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS"), search = "GLUCOSE")
nrow(f$a) # 13 --> OK
nrow(f$b) # 4

GO_use <- c("GO_GLUCOSE_METABOLIC_PROCESS", "GO_REGULATION_OF_GLUCOSE_METABOLIC_PROCESS", "GO_REGULATION_OF_GLUCONEOGENESIS")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "GLYCOLYSIS")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

# hnames[19] HEDGEHOG_SIGNALING ---------------------------------------------------------------------------------------
# Genes up-regulated by activation of hedgehog signaling.
length( c5names[grep("HEDGE",c5names)] ) # 0
length( c5names[grep("SHH",c5names)] ) # 0
length( c5names[grep("SONIC",c5names)] ) # 0

GO_use <- c("NONE")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "HEDGEHOG_SIGNALING")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[20] HEME_METABOLISM -------------------------------------------------------------------------------------
# Genes involved in metabolism of heme (a cofactor consisting of iron and porphyrin)
# and erythroblast differentiation.
# binding O2, electron transfer, oxidation reactions
length( c5names[grep("HEME_MET", c5names)] ) # 1
length( c5names[grep("ERYTHROBLAST", c5names)] ) # 0

GO_use <- c("GO_HEME_METABOLIC_PROCESS")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "HEME_METABOLISM")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[21] HYPOXIA --------------------------------------------------------------------------------------------
# Genes up-regulated in response to low oxygen levels (hypoxia).
length( c5names[grep("HYPOXIA", c5names)] ) # 2

GO_use <- c("GO_REGULATION_OF_CELLULAR_RESPONSE_TO_HYPOXIA")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "HYPOXIA")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[22] IL2_STAT5_SIGNALING -------------------------------------------------------------------------------------
# Genes up-regulated by STAT5 in response to IL2 stimulation.
# STAT = Signaling transducer and Activators of Transcription
# JAK: tyrosine phosphorylation of STAT, activates them, and STAT activates genes. 

# 22 a) Interleukin 2 ----
length( c5names[grep("INTERLEUKIN_2", c5names)] ) # 5
c5names[grep("INTERLEUKIN_2", c5names)]
# [1] "GO_REGULATION_OF_INTERLEUKIN_2_BIOSYNTHETIC_PROCESS"         
# [2] "GO_POSITIVE_REGULATION_OF_INTERLEUKIN_2_PRODUCTION"          
# [3] "GO_REGULATION_OF_INTERLEUKIN_2_PRODUCTION"                   
# [4] "GO_POSITIVE_REGULATION_OF_INTERLEUKIN_2_BIOSYNTHETIC_PROCESS"
# [5] "GO_NEGATIVE_REGULATION_OF_INTERLEUKIN_2_PRODUCTION"
## None really good?

# 22 b) STAT5 (signal transducer and activator of transcription 5) ----
length( c5names[grep("STAT5", c5names)] ) # 2

f <- fun(big_group = c("GO_REGULATION_OF_TYROSINE_PHOSPHORYLATION_OF_STAT5_PROTEIN"), search = "STAT5")
nrow(f$a) # 0
nrow(f$b) # 2 --> OK, and includes the IL2 gene

GO_use <- c("GO_REGULATION_OF_TYROSINE_PHOSPHORYLATION_OF_STAT5_PROTEIN")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "IL2_STAT5_SIGNALING")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[23] IL6_JAK_STAT3_SIGNALING --------------------------------------------------------------------------------
# Genes up-regulated by IL6 via STAT3,  e.g., during acute phase response.
# STAT = Signaling transducer and Activators of Transcription
# JAK --> tyrosine phosphorylation of STAT, activates them, and STAT activates genes. 

# 23 a) Interleukin 6 ----
length( c5names[grep("INTERLEUKIN_6", c5names)] ) # 7

f <- fun(big_group = c("GO_RESPONSE_TO_INTERLEUKIN_6"), search = "INTERLEUKIN_6")
nrow(f$a) # 5
nrow(f$b) # 2

# 23 b) JAK (Janus kinases) ----
length( c5names[grep("JAK", c5names)] ) # 1 --> Not good
length( c5names[grep("JANUS", c5names)] ) # 0

# 23 c) STAT3 (signal transducer and activator of transcription 3) ----
length( c5names[grep("STAT3", c5names)] ) # 2

f <- fun(big_group = c("GO_REGULATION_OF_TYROSINE_PHOSPHORYLATION_OF_STAT3_PROTEIN"), search = "STAT3")
nrow(f$a) # 0
nrow(f$b) # 2

GO_use <- c("GO_RESPONSE_TO_INTERLEUKIN_6", "GO_REGULATION_OF_TYROSINE_PHOSPHORYLATION_OF_STAT3_PROTEIN")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "IL6_JAK_STAT3_SIGNALING")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[24] INFLAMMATORY_RESPONSE ----------------------------------------------------------------------------------
# Genes defining inflammatory response
length( c5names[grep("INFLAMMATORY_RESPONSE", c5names)] ) # 15

f <- fun(big_group = c("GO_INFLAMMATORY_RESPONSE", "GO_REGULATION_OF_INFLAMMATORY_RESPONSE"), search = "INFLAMMATORY_RESPONSE")
nrow(f$a) # 0
nrow(f$b) # 15

GO_use <- c("GO_INFLAMMATORY_RESPONSE", "GO_REGULATION_OF_INFLAMMATORY_RESPONSE")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "INFLAMMATORY_RESPONSE")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[25] INTERFERON_ALPHA_RESPONSE -------------------------------------------------------------------------------
# Genes up-regulated in response to alpha interferon proteins.
# regulate the activity of the immune system, viral infections
length( c5names[grep("INTERFERON_ALPHA", c5names)] ) # 3

f <- fun(big_group = c("GO_RESPONSE_TO_INTERFERON_ALPHA"), search = "INTERFERON_ALPHA")
nrow(f$a) # 2 --> This is OK!
nrow(f$b) # 1

GO_use <- c("GO_RESPONSE_TO_INTERFERON_ALPHA")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "INTERFERON_ALPHA_RESPONSE")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[26] INTERFERON_GAMMA_RESPONSE ---------------------------------------------------------------------------------
# Genes up-regulated in response to IFNG
# interferons = proteins release during viral infections
# innate and adaptive immunity against infections
length( c5names[grep("INTERFERON_GAMMA", c5names)] ) # 11

f <- fun(big_group = c("GO_RESPONSE_TO_INTERFERON_GAMMA", "GO_REGULATION_OF_RESPONSE_TO_INTERFERON_GAMMA"), search = "INTERFERON_GAMMA")
nrow(f$a) # 7 --> This is OK
nrow(f$b) # 4

GO_use <- c("GO_RESPONSE_TO_INTERFERON_GAMMA", "GO_REGULATION_OF_RESPONSE_TO_INTERFERON_GAMMA")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "INTERFERON_GAMMA_RESPONSE")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[27] KRAS_SIGNALING --------------------------------------------------------------------------------------------
# Genes down-regulated by KRAS activation.
# The KRAS gene provides instructions for making a protein called K-Ras, part of the RAS/MAPK pathway.
length( c5names[grep("KRAS",c5names)] ) # 0
length( c5names[grep("_RAS",c5names)] ) # 2
length( c5names[grep("MAPK",c5names)] ) # 9

## Suggestions:
# GO_RAS_PROTEIN_SIGNAL_TRANSDUCTION
# GO_REGULATION_OF_RAS_PROTEIN_SIGNAL_TRANSDUCTION
## Not as specific as we want them...

GO_use <- c("NONE")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "KRAS_SIGNALING")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[28] MITOTIC_SPINDLE ------------------------------------------------------------------------------------------
# Genes important for mitotic spindle assembly.
length( c5names[grep("MITOTIC_SPINDLE", c5names)] ) # 5

f <- fun(big_group = c("GO_MITOTIC_SPINDLE_ASSEMBLY", "GO_MITOTIC_SPINDLE_ORGANIZATION"), search = "MITOTIC_SPINDLE")
nrow(f$a) # 3 --> OK 
nrow(f$b) # 2

GO_use <- c("GO_MITOTIC_SPINDLE_ASSEMBLY", "GO_MITOTIC_SPINDLE_ORGANIZATION")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "MITOTIC_SPINDLE")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[29] MTORC1_SIGNALING -----------------------------------------------------------------------------------------
# Genes up-regulated through activation of mTORC1 complex.
length( c5names[grep("MTOR",c5names)] ) # 0
length( c5names[grep("M_TOR",c5names)] ) # 0
length( c5names[grep("RAPAMYCIN",c5names)] ) # 0
length( c5names[grep("_TOR",c5names)] ) # 4

## Suggestions: (for PI3K_MTOR_AKT later on)
# GO_TOR_SIGNALING
# GO_REGULATION_OF_TOR_SIGNALING

GO_use <- c("NONE")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "MTORC1_SIGNALING")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[30] MYC_TARGETS ------------------------------------------------------------------------------------------------
# A subgroup of genes regulated by MYC
length( c5names[grep("MYC",c5names)] ) # 0

GO_use <- c("NONE")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "MYC_TARGETS")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[31] MYOGENESIS -------------------------------------------------------------------------------------------------
# Genes involved in development of skeletal muscle (myogenesis).
length( c5names[grep("MYOGENESIS",c5names)] ) # 0
length( c5names[grep("SKELETAL_MUSCLE",c5names)] ) # 12
length( c5names[grep("SKELETAL_MUSCLE_TISSUE_DEVELOPMENT",c5names)] ) # 3

f <- fun(big_group = c("GO_REGULATION_OF_SKELETAL_MUSCLE_TISSUE_DEVELOPMENT"), search = "SKELETAL_MUSCLE")
nrow(f$a) #  8 --> This is OK
nrow(f$b) # 4

GO_use <- c("GO_REGULATION_OF_SKELETAL_MUSCLE_TISSUE_DEVELOPMENT")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "MYOGENESIS")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[32] NOTCH_SIGNALING -------------------------------------------------------------------------------------------
# Genes up-regulated by activation of Notch signaling.
length( c5names[grep("NOTCH",c5names)] ) # 5

f <- fun(big_group = c("GO_NOTCH_SIGNALING_PATHWAY", "GO_REGULATION_OF_NOTCH_SIGNALING_PATHWAY"), search = "NOTCH")
nrow(f$a) # 0
nrow(f$b) # 5

GO_use <- c("GO_NOTCH_SIGNALING_PATHWAY", "GO_REGULATION_OF_NOTCH_SIGNALING_PATHWAY")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "NOTCH_SIGNALING")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[33] OXIDATIVE_PHOSPHORYLATION --------------------------------------------------------------------------------------
# Genes encoding proteins involved in oxidative phosphorylation.
# in mitochondria, production of ATP
length( c5names[grep("OXIDATIVE_PHOSPHORYLATION", c5names)] ) # 2

GO_use <- c("GO_OXIDATIVE_PHOSPHORYLATION", "GO_REGULATION_OF_OXIDATIVE_PHOSPHORYLATION")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "OXIDATIVE_PHOSPHORYLATION")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[34] P53_PATHWAY ----------------------------------------------------------------------------------------------
# Genes involved in p53 pathways and networks.
# invovled in cell cycle arrest, cellular senescence or apoptosis....
length( c5names[grep("P53",c5names)] ) # 12

f <- fun(big_group = c("GO_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR", "GO_REGULATION_OF_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR"), search = "P53")
nrow(f$a)  # 0
nrow(f$b) # 12

GO_use <- c("GO_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR", "GO_REGULATION_OF_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "P53_PATHWAY")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )
 
## * hnames[35] PANCREAS_BETA_CELLS  ------------------------------------------------------------------------------------
# Genes specifically up-regulated in pancreatic beta cells.
# Beta cells are a type of cell found in pancreatic islets that synthesize and secrete insulin and amylin
length( c5names[grep("BETA_CELL",c5names)] ) # 0
length( c5names[grep("PANCREAS",c5names)] ) # 2 --> None good
length( c5names[grep("INSULIN_SEC",c5names)] ) # 4

f <- fun(big_group = c("GO_INSULIN_SECRETION", "GO_POSITIVE_REGULATION_OF_INSULIN_SECRETION"), search = "INSULIN_SEC")
nrow(f$a)  # 1 --> This is OK
nrow(f$b) # 3

GO_use <- c("GO_INSULIN_SECRETION", "GO_POSITIVE_REGULATION_OF_INSULIN_SECRETION")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "PANCREAS_BETA_CELLS")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[36] PEROXISOME ----------------------------------------------------------------------------------------------
# Genes encoding components of peroxisome.
# membrane-bound organelle, interacts with mitochondria, e.g. beta-ocidatino of fatty acids and metabolism of ROS etc.
# close contact with endoplasmic reticulum
# lipid metabolism etc. 
length( c5names[grep("PEROXISOME", c5names)] ) # 3

f <- fun(big_group = c("GO_PEROXISOME_ORGANIZATION"), search = "PEROXISOME")
nrow(f$a) # 1 This is OK
nrow(f$b) # 2

GO_use <- c("GO_PEROXISOME_ORGANIZATION")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "PEROXISOME")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[37] PI3K_AKT_MTOR_SIGNALING ----------------------------------------------------------------------------------
# Genes up-regulated by activation of the PI3K/AKT/mTOR pathway.

## 37 a) PI3K (Phosphoinositide/Phosphatidylinositol 3-kinases) ---- 
length( c5names[grep("PI3K",c5names)] ) # 0
length( c5names[grep("PHOSPHATIDYLINOSITOL_3_KINASE",c5names)] ) # 4

f <- fun(big_group = c("GO_PHOSPHATIDYLINOSITOL_3_KINASE_SIGNALING", "GO_REGULATION_OF_PHOSPHATIDYLINOSITOL_3_KINASE_SIGNALING"), search = "PHOSPHATIDYLINOSITOL_3")
nrow(f$a)  # 2 This is OK
nrow(f$b) # 3

## 37 b) AKT (Protein Kinase B) ----
length( c5names[grep("_AKT_",c5names)] ) # 0
length( c5names[grep("PROTEIN_KINASE_B",c5names)] ) # 5

f <- fun(big_group = c("GO_PROTEIN_KINASE_B_SIGNALING", "GO_REGULATION_OF_PROTEIN_KINASE_B_SIGNALING"), search = "PROTEIN_KINASE_B")
nrow(f$a)  # 1 This is OK
nrow(f$b) # 4

## 37 c) mTOR ----
# mammalian/mechanistic target of Rapamycin
length( c5names[grep("MTOR",c5names)] ) # 0
length( c5names[grep("M_TOR",c5names)] ) # 0
length( c5names[grep("RAPAMYCIN",c5names)] ) # 0
length( c5names[grep("_TOR",c5names)] ) # 4

## Suggestions:
# GO_TOR_SIGNALING
# GO_REGULATION_OF_TOR_SIGNALING

f <- fun(big_group = c("GO_TOR_SIGNALING", "GO_REGULATION_OF_TOR_SIGNALING"), search = "_TOR")
nrow(f$a)  # 0
nrow(f$b) # 4

GO_use <- c("GO_PHOSPHATIDYLINOSITOL_3_KINASE_SIGNALING", "GO_REGULATION_OF_PHOSPHATIDYLINOSITOL_3_KINASE_SIGNALING",
            "GO_PROTEIN_KINASE_B_SIGNALING", "GO_REGULATION_OF_PROTEIN_KINASE_B_SIGNALING",
            "GO_TOR_SIGNALING", "GO_REGULATION_OF_TOR_SIGNALING")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "PI3K_AKT_MTOR_SIGNALING")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[38] PROTEIN_SECRETION ----------------------------------------------------------------------------------------
# Genes involved in protein secretion pathway.
length( c5names[grep("PROTEIN_SECRETION",c5names)] ) # 4

f <- fun(big_group = c("GO_PROTEIN_SECRETION", "GO_REGULATION_OF_PROTEIN_SECRETION"), search = "PROTEIN_SECRETION")
nrow(f$a)  # 0
nrow(f$b) # 4

GO_use <- c("GO_PROTEIN_SECRETION", "GO_REGULATION_OF_PROTEIN_SECRETION")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "PROTEIN_SECRETION")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )
 
## * hnames[39] REACTIVE_OXYGEN_SPECIES_PATHWAY ---------------------------------------------------------------------------
# Genes up-regulated by reactive oxigen species (ROS).
length( c5names[grep("REACTIVE_OXYGEN_SPECIES",c5names)] ) # 12

f <- fun(big_group = c("GO_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES", "GO_REGULATION_OF_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES"), search = "REACTIVE_OXYGEN_SPECIES")
nrow(f$a)  # 8 This is OK
nrow(f$b) # 4

GO_use <- c("GO_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES", "GO_REGULATION_OF_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "REACTIVE_OXIGEN_SPECIES_PATHWAY")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[40] SPERMATOGENESIS ------------------------------------------------------------------------------------------
# Genes up-regulated during production of male gametes (sperm), as in spermatogenesis.
length( c5names[grep("SPERM",c5names)] ) # 7
length( c5names[grep("GAMETE",c5names)] ) # 4

## Not really relevant this project, + not as specific as we want?

GO_use <- c("NONE")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "SPERMATOGENESIS")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[41] TGF_BETA_SIGNALING ------------------------------------------------------------------------------------
# Genes up-regulated in response to TGFB1
# involvedin cell growth, differenation, development, apoptosis..... and a part of the immune system.

length( c5names[grep("TRANSFORMING_GROWTH_FACTOR_BETA", c5names)] ) # 7

f <- fun(big_group = c("GO_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA", "GO_REGULATION_OF_CELLULAR_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA_STIMULUS"), search = "TRANSFORMING_GROWTH_FACTOR_BETA")
nrow(f$a) # 2 --> OK
nrow(f$b) # 5

GO_use <- c("GO_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA", "GO_REGULATION_OF_CELLULAR_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA_STIMULUS")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "TGF_BETA_SIGNALING")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[42] TNFA_SIGNALING_VIA_NFKB ------------------------------------------------------------------------------
# Genes regulated by NF-kB in response to TNF
# TNF (tumor necrosis factor) alpha signaling via NF-kappaB (NF-kappaB is a transcription factor)
# in absence of TNFa, NFkB is not active

# 42 a) TNFA ----
length( c5names[grep("TNF", c5names)] ) # 0
length( c5names[grep("TUMOR_NECROSIS_FACTOR", c5names)] ) # 8

f <- fun(big_group = c("GO_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY",
                       "GO_REGULATION_OF_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY"),
         search = "TUMOR_NECROSIS_FACTOR")
nrow(f$a) # 5 --> OK
nrow(f$b) # 3
# These should be OK, see description: http://www.informatics.jax.org/vocab/gene_ontology/GO:0033209

# 42 b) NFkB ----
length( c5names[grep("NF_KAPPAB", c5names)] ) # 13
## Not as specific as we want it to be! 

GO_use <- c("GO_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY",
            "GO_REGULATION_OF_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "TNFA_SIGNALING_VIA_NFKB")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[43] UNFOLDED_PROTEIN_RESPONSE -------------------------------------------------------------------------------
# Genes up-regulated during unfolded protein response, a cellular stress
# response related to the endoplasmic reticulum.
length( c5names[grep("UNFOLDED_PROTEIN_RESPONSE",c5names)] ) # 6

f <- fun(big_group = c("GO_REGULATION_OF_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE"), search = "UNFOLDED_PROTEIN_RESPONSE")
nrow(f$a) # 2 --> This is OK
nrow(f$b) # 4

GO_use <- c("GO_REGULATION_OF_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "UNFOLDED_PROTEIN_RESPONSE")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[44] UV_RESPONSE --------------------------------------------------------------------------------------------
# Genes down-regulated in response to ultraviolet (UV) radiation.
# Genes up-regulated in response to ultraviolet (UV) radiation.
length( c5names[grep("_UV",c5names)] ) # 5

f <- fun(big_group = c("GO_RESPONSE_TO_UV"), search = "_UV")
nrow(f$a) # 0
nrow(f$b) # 5

GO_use <- c("GO_RESPONSE_TO_UV")
c5names_use <- c(c5names_use, "UV_RESPONSE")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[45] WNT_BETA_CATENIN_SIGNALING -----------------------------------------------------------------------------
# Genes up-regulated by activation of WNT signaling through accumulation of beta catenin CTNNB1
length( c5names[grep("WNT",c5names)] ) # 13

f <- fun(big_group = c("GO_WNT_SIGNALING_PATHWAY", "GO_REGULATION_OF_WNT_SIGNALING_PATHWAY"), search = "WNT")
nrow(f$a) # 0
nrow(f$b) # 13

length( c5names[grep("BETA_C",c5names)] ) # 2 --> Not good
length( c5names[grep("CATENIN",c5names)] ) # 4 --> Not good

GO_use <- c("GO_WNT_SIGNALING_PATHWAY", "GO_REGULATION_OF_WNT_SIGNALING_PATHWAY")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "WNT_BETA_CATENIN_SIGNALING")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

## * hnames[46] XENOBIOTIC_METABOLISM -------------------------------------------------------------------------
# Genes encoding proteins involved in processing of drugs and other xenobiotics.
length( c5names[grep("XENO", c5names)] ) # 2

GO_use <- c("GO_RESPONSE_TO_XENOBIOTIC_STIMULUS")
GO_use_save <- c(GO_use_save, GO_use)
c5names_use <- c(c5names_use, "XENOBIOTIC_METABOLISM")
c5probes_use <- c(c5probes_use, list(unlist(c5probes[which(c5names %in% GO_use)])) )

# END #

length(c5names_use) # 46
length(c5probes_use) # 46
# cbind(hnames, c5names_use)
identical(hnames, c5names_use) # TRUE
length(GO_use_save) # 76

save(c5names_use, c5probes_use, file="output/GO_use.RData")




