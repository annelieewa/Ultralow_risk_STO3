###################################################################################
# Manuscript: Clinical and Molecular Characteristics of ER-Positive Ultralow Risk
#            Breast Cancer Tumors Identified by the 70-Gene Signature
# Author: Annelie Johansson
# Date: August 2021
###################################################################################

setwd("/Volumes/Annelie Encrypted/Projects/Ultralow_risk_genes/Analysis/")

library(ggplot2)
library(ggmosaic)
options(stringsAsFactors = FALSE)

# Load data:
tab <- read.table("output/stable3.txt", sep = "\t", header = TRUE, row.names = NULL)
module_names <- tab[seq(1,nrow(tab),3),2]
tab <- tab[-seq(1,nrow(tab),3),]

# save p-values
pvals <- dplyr::select(tab, P1, P2, P3)
  
tab <- tab[,-which(colnames(tab) %in% c("P1", "P2", "P3"))]
rownames(tab)[seq(1,nrow(tab),2)] <- module_names
rownames(tab)[seq(2,nrow(tab),2)] <- paste0(module_names, "_")
colnames(tab)[which(colnames(tab) == "row.names")] <- "Cat"

# extract percentages
fun <- function(x){ strsplit(strsplit(x, "[(]")[[1]][2], "[)]")[[1]][1] }
tab$UL <- sapply(tab$UL, fun)
tab$ERpos <- sapply(tab$ERpos, fun)
tab$LumA <- sapply(tab$LumA, fun)
tab$LumB <- sapply(tab$LumB, fun)


#### Multiple plot function ----
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout_ If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


get_stripes <- function(pos){
  x = c(pos+0.12, pos-0.12, rep(pos-0.35, 10))
  xend = c(rep(pos+0.35, 10), c(pos+0.12, pos-0.12))
  y = c(rep(0,3), seq(10, 90, 10))
  yend = c(seq(10,100,10), c(100,100))
  return(list(x=x, xend=xend, y=y, yend=yend))
}

##### Create stacked barplots
fun_barplot <- function(module){
  
  ind <- grep(module, rownames(tab))
  
  UL <- cbind(tab[ind, c("Cat", "UL")], c("UL", "UL"))
  ERpos <- cbind(tab[ind, c("Cat", "ERpos")], c("ERpos", "ERpos"))
  LumA <- cbind(tab[ind, c("Cat", "LumA")], c("LumA", "LumA"))
  LumB <- cbind(tab[ind, c("Cat", "LumB")], c("LumB", "LumB"))
  colnames(UL) <- colnames(ERpos) <- colnames(LumA) <- colnames(LumB) <- c("Level", "Perc", "Cat")
  
  temp <- rbind(UL, ERpos, LumA, LumB)

  pvals_temp <- pvals[which(rownames(tab) == module),]
  
  temp$Level <- factor(temp$Level, levels = temp$Level[1:2])
  temp$Cat <- factor(temp$Cat, levels = c("UL", "ERpos", "LumA", "LumB"))
  temp$Perc <- as.numeric(temp$Perc)
  
  a <- 100 - (temp$Perc[seq(1, nrow(temp), 2)])/2
  b <- temp$Perc[seq(2, nrow(temp), 2)]/2
  pos <- c()
  for(j in 1:length(a)){
      pos <- c(pos, a[j], b[j])
  }

  pos[which(pos == 100)] <- 200
  pos[which(pos == 0)] <- -100
  pos[which(pos > 0 & pos < 2.5)] <- 4
  pos[which(pos > 97.5 & pos < 100)] <- 96
  temp$pos <- pos
  
  stripes2 <- NULL
  stripes3 <- NULL
  stripes4 <- NULL
  if(pvals_temp$P1 == "<0.001" | pvals_temp$P1 < 0.05){ stripes2 <- get_stripes(2) }
  if(pvals_temp$P2 == "<0.001" | pvals_temp$P2 < 0.05){ stripes3 <- get_stripes(3) }
  if(pvals_temp$P3 == "<0.001" | pvals_temp$P3 < 0.05){ stripes4 <- get_stripes(4) }
  
  pvals_text <- pvals_temp
  pvals_text[which(pvals_text != "<0.001")] <- paste0("P=",pvals_text[which(pvals_text != "<0.001")])
  pvals_text[which(pvals_text == "<0.001")] <- paste0("P",pvals_text[which(pvals_text == "<0.001")])
  
  p <- ggplot(data = temp, aes(x = Cat, y = Perc)) +
    geom_col(aes(fill = Level), width = 0.7) +
    scale_y_continuous(limits = c(0,104),
                       breaks = c(0, 25, 50, 75, 100),
                       labels = paste0(c(0, 25, 50, 75, 100), "%")) +
    annotate("segment", x = stripes2$x, xend = stripes2$xend, y = stripes2$y, yend = stripes2$yend, col="grey", alpha=1) +
    annotate("segment", x = stripes3$x, xend = stripes3$xend, y = stripes3$y, yend = stripes3$yend, col="grey", alpha=1) +
    annotate("segment", x = stripes4$x, xend = stripes4$xend, y = stripes4$y, yend = stripes4$yend, col="grey", alpha=1) +
    labs(y = "", x = "") +
    ggtitle(titles[which(names(titles) == module)]) +
    #ggtitle(module) +
    scale_x_discrete(labels=c("Ultralow", "ER-positive", "Luminal A", "Luminal B")) +
    geom_text(aes(y = pos, label = formatC(Perc, format="f", digits=1), group = Perc), color = "white") +
    scale_fill_manual(values = as.character(unlist(col_list[module]))) +
    theme_bw() + theme(legend.title = element_blank(),
                       panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    annotate("text", x = 1, y = 104, label = "Reference", size = 3) + 
    annotate("text", x = 2, y = 104, label = pvals_text$P1, size = 3) + 
    annotate("text", x = 3, y = 104, label = pvals_text$P2, size = 3) +
    annotate("text", x = 4, y = 104, label = pvals_text$P3, size = 3)
  return(p)
}


# in order: light to dark
col_list <- list(AURKA = c("azure3", "azure4"),
                 AKT_MTOR = c("lightblue2","royalblue1"),
                 ERBB2 = c("lightgreen","#0dc478"),
                 IGF1 = c("lightpink", "lightcoral"),
                 PIK3CA = c("khaki2","gold2"),
                 PTEN = c("olivedrab2","darkolivegreen4"),
                 IMMUNE1 = c("lightsalmon", "sienna1"),
                 IMMUNE2 =  c("plum2", "mediumorchid")
                 )
# Check if colors looks good
barplot(rep(5,(length(col_list)*2)), col=unlist(col_list))

titles <- c("AKT-MTOR (pathway)",
            "AURKA (proliferation)",
            "BETAC (beta-catenin, pathway)",
            "CASP3 (apoptosis)",
            "E2F3 (pathway)",
            "ERBB2 (HER2-signaling)",
            "ESR1 (ER-signaling)",
            "IGF1 (pathway)",
            "IMMUNE1 (immune response)",
            "IMMUNE2 (immune response)",
            "MAPK (pathway)",
            "MYC (pathway)",
            "PIK3CA-mutations",
            "PTEN-loss",
            "RAS (pathway)",
            "SRC (pathway)",
            "STROMA1 (stromal invasion)",
            "STROMA2 (stromal invasion)",
            "VEGF (angiogenesis)")
names(titles) <- module_names

g_aurka <- fun_barplot("AURKA")
g_akt_mtor <- fun_barplot("AKT_MTOR")
g_erbb2 <- fun_barplot("ERBB2")
g_igf1 <- fun_barplot("IGF1")
g_pik3ca <- fun_barplot("PIK3CA")
g_pten <- fun_barplot("PTEN")
g_immune1 <- fun_barplot("IMMUNE1")
g_immune2 <- fun_barplot("IMMUNE2")

pdf(file = "output/Figure2.pdf", width = 9, height = 10)
multiplot(g_aurka,
          g_erbb2,
          g_pik3ca,
          g_immune1,
          # next column
          g_akt_mtor,
          g_igf1,
          g_pten, 
          g_immune2, 
          cols=2)
dev.off()
