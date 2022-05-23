##Canopy arthropod diversity analyses
##R script 3/4
##dbRDA analysis and plots
##NOTE: best models/model comparisons were conducted using the PRIMER PERMANOVA+
##software from Marti Anderson. That code is available upon request. The purposes
##of this script are for figure generationonly.


##Read in datasets
antsSite <- read.csv("Data/SiteAnts.csv")
antsFSite <- read.csv("Data/FG_CanopyAntsStand.csv")
spidersSite <- read.csv("Data/CanopySpidersPRIMER.csv")
spidersFSite <- read.csv("Data/CanopySpiderFamPRIMER.csv")[-13,]

##stand-level environmental data
latlon <- read.csv("G:/Shared drives/Canopy Study/Site Info/Canopy_Stand_Data.csv")

##generate bray-curtis matrices
responseA <- vegdist(sqrt(antsSite[,-1]), method = "bray")
responseAF <- vegdist(sqrt(antsFSite[,-1]), method = "bray")
responseS <- vegdist(sqrt(spidersSite[,-1]), method = "bray")
responseSF <- vegdist(sqrt(spidersFSite[,-1]), method = "bray")

##best models, per PRIMER PERMANOVA+ software
dbrdaA <- dbrda(responseA ~ Precip_Seas,
                data = latlon)
dbrdaAF <- dbrda(responseAF ~ Precip_Seas + Isothermality + lndscp_GYRATE_AM,
                 data = latlon)
dbrdaS <- dbrda(responseS ~ Isothermality + Temp_Warm + lndscp_CLUMPY,
                data = latlon)
dbrdaSF <- dbrda(responseSF ~ Isothermality + TRich + lndscp_GYRATE_AM +
                   lndscp_CLUMPY,
                 data = latlon)



##Extract biplot vectors for each of the year terms (main plus NLCDs); drop
##first three rows which are the biplots for the main effects of landuse
##Fortify - this is all from Gavin Simpson's autoplot.rda function
obj <- fortify.cca(dbrdaA, axes = c(1:2))

## sort out x, y aesthetics
vars <- colnames(obj)[3:4]

## subset out the layers wanted
obj <- obj[obj[["Score"]] %in% c("sites", "biplot", "centroids"), , drop = FALSE]

## scale biplot vectors to length 
want <- obj[["Score"]] == "biplot"
## use 'mul' for adjusting vector length of the species!!
mul <- arrowMul(obj[want, vars, drop = FALSE],
                obj[!want, vars, drop = FALSE])

obj[want,vars] <- mul * obj[want, vars]

##Generate the biplot vectors (bpvecs) 
bpvecs <- obj[want,]

##Generate the site scores dataset (with identifying info)
sitescores <- obj %>% filter(Score == "sites")
sitescores$Eco <- latlon$ECO
sitescores$Site <- latlon$SITE

sitescoresA <- sitescores
mulA <- mul
bpvecsA <- bpvecs
bpvecsA$Label = c("Precipitation\nSeasonality")
A <- ggplot(sitescoresA)+
  ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = mulA), color = "black")+
  scale_y_continuous(limits = c(-1.2, 1.2))+
  scale_x_continuous(limits = c(-1.2, 1.2))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_point(aes(x = dbRDA1, y = MDS1, color = Eco, shape = Site),
             size = 4)+
  geom_segment(data = bpvecsA,
               aes(x = 0, y = 0,
                   xend = dbRDA1, yend = MDS1),
               arrow = arrow(length = unit(0.3, "cm")),
               size = 1)+
  ggrepel::geom_text_repel(data = bpvecsA,
                           aes(x = dbRDA1, y = MDS1, label = Label),
                           size = 3.5, max.time = 0,
                           box.padding = 0.1, point.padding = 0.1,
                           fontface = "bold",
                           #alpha = 0.7,
                           bg.color = "white",
                           bg.r = .1,
                           nudge_x = -0.2,
                           nudge_y = 0.2)+
  xlab("dbRDA1 (100% of fitted, 16.3% of total variation)")+
  scale_color_manual(values = c("#FAA300","#770C83"), guide = "none")+
  guides(shape = guide_legend(nrow = 2))+
  theme_nmds()+
  theme(legend.position = "bottom",
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13))
A2 <- A + theme(legend.position = "none")

##
obj <- fortify.cca(dbrdaAF, axes = c(1:2))

## sort out x, y aesthetics
vars <- colnames(obj)[3:4]

## subset out the layers wanted
obj <- obj[obj[["Score"]] %in% c("sites", "biplot", "centroids"), , drop = FALSE]

## scale biplot vectors to length 
want <- obj[["Score"]] == "biplot"
## use 'mul' for adjusting vector length of the species!!
mul <- arrowMul(obj[want, vars, drop = FALSE],
                obj[!want, vars, drop = FALSE])

obj[want,vars] <- mul * obj[want, vars]

##Generate the biplot vectors (bpvecs) 
bpvecs <- obj[want,]

##Generate the site scores dataset (with identifying info)
sitescores <- obj %>% filter(Score == "sites")
sitescores$Eco <- latlon$ECO
sitescores$Site <- latlon$SITE

sitescoresAF <- sitescores
mulAF <- mul
bpvecsAF <- bpvecs
bpvecsAF$Label <- c("Precipitation\nSeasonality",
                    "Isothermality",
                    "Patch Connectedness")
bpvecsAF$dbRDA12 <- bpvecsAF$dbRDA1 - c(.3,0,0)
bpvecsAF$dbRDA22 <- bpvecsAF$dbRDA2 - c(0,0.1,0.1)

B <- ggplot(sitescoresAF)+
  ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = mulAF), color = "black")+
  scale_y_continuous(limits = c(-1.2, 1.2))+
  scale_x_continuous(limits = c(-1.2, 1.2))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_point(aes(x = dbRDA1, y = dbRDA2, color = Eco, shape = Site),
             size = 4)+
  geom_segment(data = bpvecsAF,
               aes(x = 0, y = 0,
                   xend = dbRDA1, yend = dbRDA2),
               arrow = arrow(length = unit(0.3, "cm")),
               size = 1)+
  ggrepel::geom_text_repel(data = bpvecsAF,
                           aes(x = dbRDA12, y = dbRDA22, label = Label),
                           size = 3.5, max.time = 0,
                           box.padding = 0.1, point.padding = 0.1,
                           fontface = "bold",
                           #alpha = 0.7,
                           bg.color = "white",
                           bg.r = .1)+
  xlab("dbRDA1 (66.5% of fitted, 43.5% of total variation)")+
  ylab("dbRDA2 (25.9% of fitted, 16.9% of total variation)")+
  scale_color_manual(values = c("#FAA300","#770C83"), guide = "none")+
  theme_nmds()+
  theme(legend.position = "none",
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13))


##
obj <- fortify.cca(dbrdaS, axes = c(1:2))

## sort out x, y aesthetics
vars <- colnames(obj)[3:4]

## subset out the layers wanted
obj <- obj[obj[["Score"]] %in% c("sites", "biplot", "centroids"), , drop = FALSE]

## scale biplot vectors to length 
want <- obj[["Score"]] == "biplot"
## use 'mul' for adjusting vector length of the species!!
mul <- arrowMul(obj[want, vars, drop = FALSE],
                obj[!want, vars, drop = FALSE])

obj[want,vars] <- mul * obj[want, vars]

##Generate the biplot vectors (bpvecs) 
bpvecs <- obj[want,]

##Generate the site scores dataset (with identifying info)
sitescores <- obj %>% filter(Score == "sites")
sitescores$Eco <- latlon$ECO
sitescores$Site <- latlon$SITE

sitescoresS <- sitescores
mulS <- mul
bpvecsS <- bpvecs
bpvecsS$Label <- c("Isothermality",
                   "Max Temperature\nWarmest Month",
                    "Fragmentation\nIndex")
bpvecsS$dbRDA12 <- bpvecsS$dbRDA1 - c(.45,0,0)
bpvecsS$dbRDA22 <- bpvecsS$dbRDA2 - c(-0.05,.15,0.15)
C <- ggplot(sitescoresS)+
  ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = mulS), color = "black")+
  scale_y_continuous(limits = c(-1.3, 1.3))+
  scale_x_continuous(limits = c(-1.3, 1.3))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_point(aes(x = dbRDA1, y = dbRDA2, color = Eco, shape = Site),
             size = 4)+
  geom_segment(data = bpvecsS,
               aes(x = 0, y = 0,
                   xend = dbRDA1, yend = dbRDA2),
               arrow = arrow(length = unit(0.3, "cm")),
               size = 1)+
  ggrepel::geom_text_repel(data = bpvecsS,
                           aes(x = dbRDA12, y = dbRDA22, label = Label),
                           size = 3.5, max.time = 0,
                           box.padding = 0.1, point.padding = 0.1,
                           fontface = "bold",
                           #alpha = 0.7,
                           bg.color = "white",
                           bg.r = .1)+
  xlab("dbRDA1 (43.8% of fitted, 19.1% of total variation)")+
  ylab("dbRDA2 (33.7% of fitted, 14.7% of total variation)")+
  scale_color_manual(values = c("#FAA300","#770C83"), guide = "none")+
  theme_nmds()+
  theme(legend.position = "none",
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13))


##
obj <- fortify.cca(dbrdaSF, axes = c(1:2))

## sort out x, y aesthetics
vars <- colnames(obj)[3:4]

## subset out the layers wanted
obj <- obj[obj[["Score"]] %in% c("sites", "biplot", "centroids"), , drop = FALSE]

## scale biplot vectors to length 
want <- obj[["Score"]] == "biplot"
## use 'mul' for adjusting vector length of the species!!
mul <- arrowMul(obj[want, vars, drop = FALSE],
                obj[!want, vars, drop = FALSE])

obj[want,vars] <- mul * obj[want, vars]

##Generate the biplot vectors (bpvecs) 
bpvecs <- obj[want,]

##Generate the site scores dataset (with identifying info)
sitescores <- obj %>% filter(Score == "sites")
sitescores$Eco <- latlon$ECO
sitescores$Site <- latlon$SITE

sitescoresSF <- sitescores
mulSF <- mul
bpvecsSF <- bpvecs
bpvecsSF$Label <- c("Isothermality",
                    "Tree Richness",
                    "Patch Connectedness",
                    "Fragmentation\nIndex").
bpvecsSF$dbRDA12 <- bpvecsSF$dbRDA1 - c(-.35,0,0,0.2)
bpvecsSF$dbRDA22 <- bpvecsSF$dbRDA2 - c(-.05,.05,-0.1,-0.05)
D <- ggplot(sitescoresSF)+
  ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = mulSF), color = "black")+
  scale_y_continuous(limits = c(-1, 1))+
  scale_x_continuous(limits = c(-1, 1))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_point(aes(x = dbRDA1, y = dbRDA2, color = Eco, shape = Site),
             size = 4)+
  geom_segment(data = bpvecsSF,
               aes(x = 0, y = 0,
                   xend = dbRDA1, yend = dbRDA2),
               arrow = arrow(length = unit(0.3, "cm")),
               size = 1)+
  ggrepel::geom_text_repel(data = bpvecsSF,
                           aes(x = dbRDA12, y = dbRDA22, label = Label),
                           size = 3.5, max.time = 0,
                           box.padding = 0.1, point.padding = 0.1,
                           fontface = "bold",
                           #alpha = 0.7,
                           bg.color = "white",
                           bg.r = .1)+
  xlab("dbRDA1 (46.6% of fitted, 27.8% of total variation)")+
  ylab("dbRDA2 (36.1% of fitted, 21.6% of total variation)")+
  scale_color_manual(values = c("#FAA300","#770C83"), guide = "none")+
  theme_nmds()+
  theme(legend.position = "none",
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13))

f2 <- cowplot::plot_grid(A2,B,C,D, labels = c("A", "B", "C", "D"),
                         label_x = .15,
                         label_size = 20,
                         nrow = 2,
                         rel_heights = c(1,1),
                         rel_widths = c(1,1))

legend_B = cowplot::get_legend(A) 

legends <- cowplot::plot_grid(legend_B, nrow = 1)

f2_legend <- cowplot::plot_grid(f2, legends, nrow = 2, rel_heights = c(3.5,0.5))

##save plot
ggsave("composition_plots.tiff", f2_legend, dpi = 200, 
       height = 10.5, width = 9.5)
