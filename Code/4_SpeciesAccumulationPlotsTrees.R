##Canopy arthropod diversity analyses
##R script 4/4
##Species accumulation plots by tree species

##start with ants
treelvl <- canants[, c(6, 7, 9:31)] %>%
  group_by(TREE.SP, TREE.NO) %>%
  summarize_all(.funs = sum)

##Only enough replicates to make meaningful species accumulation curves for
##Beech, Hickory, Maple, Red Oak, White Oak, and Tulip poplar

##Estimate and observed species richness for specific tree species
specpool(subset(treelvl[ , c(1,3:25)], TREE.SP == "Beech")[,-1])
specpool(subset(treelvl[ , c(1,3:25)], TREE.SP == "Hickory")[,-1])
specpool(subset(treelvl[ , c(1,3:25)], TREE.SP == "Maple")[,-1])
specpool(subset(treelvl[ , c(1,3:25)], TREE.SP == "Red Oak")[,-1])
specpool(subset(treelvl[ , c(1,3:25)], TREE.SP == "White Oak")[,-1])
specpool(subset(treelvl[ , c(1,3:25)], TREE.SP == "Tulip")[,-1])

##Transfer these values to a dataframe to be used in the plot below
antlabs <- data.frame(Trees = c(12,
                                6,
                                22.5,
                                13,
                                20,
                                14.5),
                      Rich = c(13.75,
                               16.75,
                               17,
                               18.75,
                               19.75,
                               11),
                      label = c("O = 14, E = 20",
                                "O = 16, E = 27",
                                "O = 18, E = 19",
                                "O = 18, E = 41",
                                "O = 19, E = 22",
                                "O = 11, E = 13"))

##Generate species accumulation data for specific tree species
beechants <- specaccum(subset(treelvl[ , c(1,3:25)], TREE.SP == "Beech")[,-1])
hickoryants <- specaccum(subset(treelvl[ , c(1,3:25)], TREE.SP == "Hickory")[,-1])
mapleants <- specaccum(subset(treelvl[ , c(1,3:25)], TREE.SP == "Maple")[,-1])
roakants <- specaccum(subset(treelvl[ , c(1,3:25)], TREE.SP == "Red Oak")[,-1])
woakants <- specaccum(subset(treelvl[ , c(1,3:25)], TREE.SP == "White Oak")[,-1])
tulipants <- specaccum(subset(treelvl[ , c(1,3:25)], TREE.SP == "Tulip")[,-1])

##Convert these to dataframes to be used for plotting purposes
beecha <- data.frame(Trees = beechants$sites,
                     Rich = beechants$richness,
                     SD = beechants$sd,
                     Species = rep("Beech",time = length(beechants$sites)))
hickorya <- data.frame(Trees = hickoryants$sites,
                       Rich = hickoryants$richness,
                       SD = hickoryants$sd,
                       Species = rep("Hickory",time = length(hickoryants$sites)))
maplea <- data.frame(Trees = mapleants$sites,
                     Rich = mapleants$richness,
                     SD = mapleants$sd,
                     Species = rep("Maple",time = length(mapleants$sites)))
roaka <- data.frame(Trees = roakants$sites,
                    Rich = roakants$richness,
                    SD = roakants$sd,
                    Species = rep("Red Oak",time = length(roakants$sites)))
woaka <- data.frame(Trees = woakants$sites,
                    Rich = woakants$richness,
                    SD = woakants$sd,
                    Species = rep("White Oak",time = length(woakants$sites)))
tulipa <- data.frame(Trees = tulipants$sites,
                     Rich = tulipants$richness,
                     SD = tulipants$sd,
                     Species = rep("Tulip",time = length(tulipants$sites)))

##bind all of the datasets together
ants.tree.accum <- bind_rows(maplea, woaka, beecha, tulipa,
                                         roaka, hickorya)
##add a origin (0,0) for all tree species
ants.tree.accum <- bind_rows(ants.tree.accum,
                             data.frame(Trees = 0, Rich = 0,
                                        Species = unique(ants.tree.accum$Species)))

##generate the species accumulation plot for ants
antac <- ggplot(ants.tree.accum, aes(x = Trees, y = Rich)) + 
  scale_x_continuous(limits = c(0,26), breaks = seq(0,25,5))+
  geom_line(aes(color = Species), size = 1)+
  geom_point(aes(color = Species, shape = Species), size = 3)+
  geom_text(data = antlabs, aes(label = label))+
  ylab("Species Richness")+
  xlab("Number of Trees Sampled")+
  #geom_segment(x = 7.5, y = 40, xend = 8.5, yend = 38.25)+
  #geom_segment(x = 7.5, y = 25, xend = 8.25, yend = 23)+
  scale_color_manual(values = c("dodgerblue2", "#E31A1C", # red
                                "green4",
                                "#6A3D9A", # purple
                                "#FF7F00", # orange
                                "gold1"))+
  theme_nmds()+
  theme(legend.position = c(0.8, 0.25),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.grid.major.y = element_line(color = "grey"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14))



##spiders
treelvl.s <- canspids[, c(6,7,9:105)] %>%
  group_by(TREE.SP, TREE.NO) %>%
  summarize_all(.funs = sum)

##Estimate and observed species richness for specific tree species
specpool(subset(treelvl.s[ , c(1,3:99)], TREE.SP == "Beech")[,-1])
specpool(subset(treelvl.s[ , c(1,3:99)], TREE.SP == "Hickory")[,-1])
specpool(subset(treelvl.s[ , c(1,3:99)], TREE.SP == "Maple")[,-1])
specpool(subset(treelvl.s[ , c(1,3:99)], TREE.SP == "RedOak")[,-1])
specpool(subset(treelvl.s[ , c(1,3:99)], TREE.SP == "WhiteOak")[,-1])
specpool(subset(treelvl.s[ , c(1,3:99)], TREE.SP == "Tulip")[,-1])

##Transfer these values to a dataframe to be used in the plot below
spidlabs <- data.frame(Trees = c(4, 11.5,
                                 22.75, 19,
                                 13.5, 15),
                       Rich = c(40, 23,
                                55, 35,
                                64, 39),
                       label = c("O = 38, E = 74","O = 26, E = 52",
                                 "O = 65, E = 102","O = 35, E = 50",
                                 "O = 64, E = 159","O = 36, E = 52"))

##Estimate and observed species richness for specific tree species
beechs <- specaccum(subset(treelvl.s[ , c(1,3:99)], TREE.SP == "Beech")[,-1])
hickorys <- specaccum(subset(treelvl.s[ , c(1,3:99)], TREE.SP == "Hickory")[,-1])
maples <- specaccum(subset(treelvl.s[ , c(1,3:99)], TREE.SP == "Maple")[,-1])
roaks <- specaccum(subset(treelvl.s[ , c(1,3:99)], TREE.SP == "RedOak")[,-1])
woaks <- specaccum(subset(treelvl.s[ , c(1,3:99)], TREE.SP == "WhiteOak")[,-1])
tulips <- specaccum(subset(treelvl.s[ , c(1,3:99)], TREE.SP == "Tulip")[,-1])

##Convert these to dataframes to be used for plotting purposes
beechss <- data.frame(Trees = beechs$sites,
                      Rich = beechs$richness,
                      SD = beechs$sd,
                      Species = rep("Beech",time = length(beechs$sites)))
hickoryss <- data.frame(Trees = hickorys$sites,
                        Rich = hickorys$richness,
                        SD = hickorys$sd,
                        Species = rep("Hickory",time = length(hickorys$sites)))
mapless <- data.frame(Trees = maples$sites,
                      Rich = maples$richness,
                      SD = maples$sd,
                      Species = rep("Maple",time = length(maples$sites)))
roakss <- data.frame(Trees = roaks$sites,
                     Rich = roaks$richness,
                     SD = roaks$sd,
                     Species = rep("Red Oak",time = length(roaks$sites)))
woakss <- data.frame(Trees = woaks$sites,
                     Rich = woaks$richness,
                     SD = woaks$sd,
                     Species = rep("White Oak",time = length(woaks$sites)))
tulipss <- data.frame(Trees = tulips$sites,
                      Rich = tulips$richness,
                      SD = tulips$sd,
                      Species = rep("Tulip",time = length(tulips$sites)))

##bind all of the datasets together
spid.tree.accum <- bind_rows(mapless, woakss, beechss, tulipss,
                                         roakss, hickoryss)

##add a origin (0,0) for all tree species
spid.tree.accum <- bind_rows(spid.tree.accum,
                             data.frame(Trees = 0, Rich = 0,
                                        Species = unique(spid.tree.accum$Species)))

##generate the species accumulation plots for spiders
spiderac <- ggplot(spid.tree.accum, aes(x = Trees, y = Rich)) + 
  scale_x_continuous(limits = c(0,26), breaks = seq(0,25,5))+
  geom_line(aes(color = Species), size = 1)+
  geom_point(aes(color = Species, shape = Species), size = 3)+
  geom_text(data = spidlabs, aes(label = label))+
  ylab("Species Richness")+
  xlab("Number of Trees Sampled")+
  #geom_segment(x = 7.5, y = 40, xend = 8.5, yend = 38.25)+
  #geom_segment(x = 7.5, y = 25, xend = 8.25, yend = 23)+
  scale_color_manual(values = c("dodgerblue2", "#E31A1C", # red
                                "green4",
                                "#6A3D9A", # purple
                                "#FF7F00", # orange
                                "gold1"))+
  theme_nmds()+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.grid.major.y = element_line(color = "grey"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))

##generate the paneled figure
tspac <- cowplot::plot_grid(antac+ 
                              theme(axis.text.x = element_blank(),
                                    # axis.ticks.x = element_blank(),
                                    axis.title.x = element_blank(),
                                    plot.margin = margin(b = 0)), 
                            spiderac,
                            nrow = 2,
                            labels = c("A", "B"),
                            label_x = .15,
                            label_size = 20,
                            align = "v",
                            rel_heights = c(1,1.2))

##save
ggsave("TreeAccum_plots.tiff", tspac, dpi = 300,
       height = 7.506, width = 4.25)
