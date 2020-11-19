##### Project: Ultralow risk STO3 ######################################
# mosaic_plot.R
# Functions for mosaic plots
# Author: Annelie Johansson & Nancy Yu & Adina Iftimi
# Modified last: January 2020 by Annelie Johansson
#########################################################################

library(ggplot2)
library(ggmosaic)

#### FUNCTIONS ####

#### Multiple plot function ----
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout, if present, 'cols' is ignored.
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
multiplot <- function(..., plotlist = NULL, file, cols = 1 , layout = NULL) {
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
  
  if (numPlots == 1) {
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


#### Function to create each plot, for each characteristic/module ----
# Plots in following order: ER-pos, UL, LumA, LumB
create_gg_mosaic <- function(char, t){
  
  # get numbers from table_mosaic for this characteristic
  t_temp <- t[which(table_mosaic$Characteristics == char), ]
  
  # get p-values (UL vs ERp , UL vs LumA, UL vs LumB)
  pvals <- table_pvals[which(rownames(table_pvals) == char),]
  
  # get colors for the plot
  col <- col_list[which(names(col_list) == char)]
  
  # Process the data for the mosaic plot
  t_temp$level <- factor(t_temp$level, levels = c(t_temp$level))
  characteristic <- rep(t_temp$level, 4)
  type <- c(rep("UL", nrow(t_temp)), rep("LumA", nrow(t_temp)),
            rep("LumB", nrow(t_temp)), rep("ERp", nrow(t_temp)))
  counts <- c(t_temp$UL, t_temp$LumA, t_temp$LumB, t_temp$ERp)
  t_temp_d <- data.frame(type, characteristic, counts)
  
  type_n <- with(t_temp_d, rep(x = type, times = counts))
  characteristic_n <- with(t_temp_d, rep(x = characteristic, times = counts))
  
  t_temp_dn <- data.frame(type_n, characteristic_n)
  
  t_temp_dn$type_n <- factor(t_temp_dn$type_n, levels = c("UL", "ERp", "LumA", "LumB"))
  levels(t_temp_dn$type_n) <- c("Ultralow", "ER positive", "Luminal A", "Luminal B")
  
  # Create figure - add stripes for significant p-values
  if(pvals[1] < 0.05) { stripes_ERp <- get_stripes(0.116, 0.584)
  }else{ stripes_ERp <- list(X = 0, Xend = 0, Y = 0, Yend = 0) }
  if(pvals[2] < 0.05) { stripes_LumA <- get_stripes(0.593, 0.858)
  }else{ stripes_LumA <- list(X = 0, Xend = 0, Y = 0, Yend = 0) }
  if(pvals[3] < 0.05) { stripes_LumB <- get_stripes(0.869, 1)
  }else{ stripes_LumB <- list(X = 0, Xend = 0, Y = 0, Yend = 0) }
  g <- ggplot(data = t_temp_dn) +
    geom_mosaic(aes(x = product(characteristic_n, type_n), fill = characteristic_n)) + 
    scale_fill_manual(values = as.vector(unlist(col)), name = char) +  
    scale_y_continuous(breaks = c(0,1), labels = c("0%","100%")) +
    theme_light() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
    guides(fill = guide_legend(reverse=TRUE)) +
    ylab("") +
    annotate("segment", x=stripes_ERp$X, xend = stripes_ERp$Xend, y = stripes_ERp$Y, yend = stripes_ERp$Yend, col = "grey", alpha = 1) +
    annotate("segment", x=stripes_LumA$X, xend = stripes_LumA$Xend, y = stripes_LumA$Y, yend = stripes_LumA$Yend, col = "grey", alpha = 1) +
    annotate("segment", x=stripes_LumB$X, xend = stripes_LumB$Xend, y = stripes_LumB$Y, yend = stripes_LumB$Yend, col = "grey", alpha = 1)
  
  return(g)
}

#### Function to mark significant differences in the plot (p<0.05) ----
# OBS: not the prettiest code, but it works for now for these numbers....
get_stripes <- function(min, max){
  if(min >= max){
    print("Try again: max value [2] must be larger than min value [1].")
  }else{
    
    a <- seq(from = 0, to = 1, by = 0.025)
    # X1
    X1 <- a[which(a > min & a <= max)]
    # Y1
    Y1 <- rep(0, length(X1))
    # Y2
    if(as.character(min) %in% a){
      Y2 <- seq(0, 0.95, 0.05)
    }else{
      y2 <- rev(a[which(a < min)])[1]
      Y2 <- seq((min - y2)*2, 1, 0.05)
    }
    # X2
    X2 <- rep(min, length(Y2))
    
    X <- c(X1, X2)
    Y <- c(Y1, Y2)
    
    if(as.character(max) %in% a){
      b <- rev(seq(1, 0.05, -0.05))
    }else{
      y3 <- a[which(a > max)][1]
      b <- rev(seq(1 - (y3 - max)*2, 0, -0.05) )
    }
    
    # Yend
    if(max == 1){
      Yend1_a <- c(rev(b[1:(length(Y1)-1)]) , 0)
      Yend1_b <- b[(length(Y1)):(length(b)-1)]
    }else{
      Yend1_a <- rev(b[1:(length(Y1))])
      Yend1_b <- b[(length(Y1)+1):length(b)]
    }
    Yend1 <- c(Yend1_a, Yend1_b)
    
    Yend2 <- rep(1, length(X1))
    Yend <- c(Yend1, Yend2)
    
    # Xend
    Xend1 <- rep(max, length(X2))
    Xend2 <- c(rev(X1)) 
    Xend <- c(Xend1, Xend2)
    
    if( any(Y == Yend) | any(X == Xend) ){
      r <- which(Y == Yend | X == Xend)
      X <- X[-r]
      Y <- Y[-r]
      Xend <- Xend[-r]
      Yend <- Yend[-r]
    }
    
    return(list(X = X, Xend = Xend, Y = Y, Yend = Yend))
  }
}