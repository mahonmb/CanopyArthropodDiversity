##Canopy arthropod diversity analyses
##R script 1/4
##Setup for the scripts; read in R packages, set theme, etc.

##Load in packages
library(tidyverse)
library(PARTITIONR)
library(vegan)
library(raster)
library(iNEXT)
library(cowplot)
library(FD)

##set a ggplot theme to be used throughout the plotting
theme_nmds <- function(){
  theme_set(theme_classic())+
    theme(
      panel.background=element_rect(fill=NA,color="black"),
      legend.background=element_blank(),
      axis.text = element_text(size=16,color="black"),
      axis.title.y=element_text(vjust=0.42),
      panel.grid.major = element_blank(),
      panel.grid.minor =element_blank(),
      legend.title = element_text(size=18),
      legend.text = element_text(size=16),
      legend.key = element_rect(colour = NA, fill = NA),
      strip.text.x = element_text(size=14),
      strip.background = element_rect(color="white", fill="white"),
      axis.title=element_text(size=17)
    )
}


fortify.cca <- function(model, data, axes = 1:6,
                        display = c("sp", "wa", "lc", "bp", "cn"), ...) {
  ## extract scores
  scrs <- scores(model, choices = axes, display = display, ...)
  ## handle case of only 1 set of scores
  if (length(display) == 1L) {
    scrs <- list(scrs)
    nam <- switch(display,
                  sp = "species",
                  species = "species",
                  wa = "sites",
                  sites = "sites",
                  lc = "constraints",
                  bp = "biplot",
                  cn = "centroids",
                  stop("Unknown value for 'display'"))
    names(scrs) <- nam
  }
  miss <- vapply(scrs, function(x ) all(is.na(x)), logical(1L))
  scrs <- scrs[!miss]
  nams <- names(scrs)
  nr <- vapply(scrs, FUN = NROW, FUN.VALUE = integer(1))
  df <- do.call('rbind', scrs)
  rownames(df) <- NULL
  df <- as.data.frame(df)
  df <- cbind(Score = factor(rep(nams, times = nr)),
              Label = unlist(lapply(scrs, rownames), use.names = FALSE),
              df)
  df
}

arrowMul <- function(arrows, data, at = c(0, 0), fill = 0.75) {
  u <- c(range(data[,1], range(data[,2])))
  u <- u - rep(at, each = 2)
  r <- c(range(arrows[, 1], na.rm = TRUE), range(arrows[, 2], na.rm = TRUE))
  rev <- sign(diff(u))[-2]
  if (rev[1] < 0)
    u[1:2] <- u[2:1]
  if (rev[2] < 0)
    u[3:4] <- u[4:3]
  u <- u/r
  u <- u[is.finite(u) & u > 0]
  fill * min(u)
}
