##Canopy arthropod diversity analyses
##R script 2/4
##Hierarchical partitioning of diversity for Ants and Spiders


##Read in datasets
#Read in canopy ant species dataset (taxonomic)
#change Tree number and sample number to factors
canants <- read.csv("CanopyAnts.csv", header = TRUE)

canants$TREE.NO <- as.factor(canants$TREE.NO)
canants$SAMPLE <- as.factor(canants$SAMPLE)

#Read in canopy ant functional group dataset (functional)
#change Tree number and sample number to factors
ant.c.f <- read.csv("CanopyAnts.csv", header = TRUE)

ant.c.f$TREE.NO <- as.factor(ant.c.f$TREE.NO)
ant.c.f$SAMPLE <- as.factor(ant.c.f$SAMPLE)


#Read in canopy spider species dataset (taxonomic)
#change Tree number and sample number to factors
canspids <- read.csv("CanopySpiders.csv", header = TRUE)
canspids$TREE.NO <- as.factor(canspids$TREE.NO)
canspids$SAMPLE <- as.factor(canspids$SAMPLE)

#Read in canopy spider family dataset (functional)
#change Tree number and sample number to factors
canspidsf <- read.csv("CanopySpiderFwTree.csv", header = TRUE)
canspidsf$TREE.NO <- as.factor(canspidsf$TREE.NO)
canspidsf$SAMPLE <- as.factor(canspidsf$SAMPLE)


##Conduct hierarchical partitioning of diversity
##Note that values generated here will vary slightly from those presented in
##the manuscript, due to randomization procedures

##Start with ant species
part.0 <- partition(canants,
                    levels = c("SAMPLE", "TREE.NO", "HAB_ST", 
                               "SITE", "ECO"),
                    low.level = 2,
                    q = 0,
                    method = "sample",
                    perms = 1000
)

summary(part.0, p.value = "two-sided")

#Need to save the output as an RDS, so it gives the same values
##Standard effect sizes for this model
(part.0$Div[2] - mean(part.0$Rand.Alpha)) / sd(part.0$Rand.Alpha)
(part.0$Div[7] - mean(part.0$Rand.Beta.Mult[[1]])) / sd(part.0$Rand.Beta.Mult[[1]])
(part.0$Div[8] - mean(part.0$Rand.Beta.Mult[[2]])) / sd(part.0$Rand.Beta.Mult[[2]])
(part.0$Div[9] - mean(part.0$Rand.Beta.Mult[[3]])) / sd(part.0$Rand.Beta.Mult[[3]])
(part.0$Div[10] - mean(part.0$Rand.Beta.Mult[[4]])) / sd(part.0$Rand.Beta.Mult[[4]])


##ant functional groups
part.0.af <- partition(ant.c.f,
                       levels = c("SAMPLE", "TREE.NO", "HAB_ST", 
                                  "SITE", "ECO"),
                       low.level = 2,
                       q = 0,
                       method = "sample",
                       perms = 1000)

summary(part.0.af, p.value = "two-sided")

#Need to save the output as an RDS, so it gives the same values
##Standrd effect sizes for this model
(part.0.af$Div[2] - mean(part.0.af$Rand.Alpha)) / sd(part.0.af$Rand.Alpha)
(part.0.af$Div[7] - mean(part.0.af$Rand.Beta.Mult[[1]])) / sd(part.0.af$Rand.Beta.Mult[[1]])
(part.0.af$Div[8] - mean(part.0.af$Rand.Beta.Mult[[2]])) / sd(part.0.af$Rand.Beta.Mult[[2]])
(part.0.af$Div[9] - mean(part.0.af$Rand.Beta.Mult[[3]])) / sd(part.0.af$Rand.Beta.Mult[[3]])
(part.0.af$Div[10] - mean(part.0.af$Rand.Beta.Mult[[4]])) / sd(part.0.af$Rand.Beta.Mult[[4]])



##spider species
part.0.s <- partition(canspids,
                      levels = c("SAMPLE", "TREE.NO", "HAB_ST", 
                                 "SITE", "ECO"),
                      low.level = 2,
                      q = 0,
                      method = "sample",
                      perms = 1000)

summary(part.0.s, p.value = "two-sided")

#Need to save the output as an RDS, so it gives the same values
##Standard effect sizes for this model
(part.0.s$Div[2] - mean(part.0.s$Rand.Alpha)) / sd(part.0.s$Rand.Alpha)
(part.0.s$Div[7] - mean(part.0.s$Rand.Beta.Mult[[1]])) / sd(part.0.s$Rand.Beta.Mult[[1]])
(part.0.s$Div[8] - mean(part.0.s$Rand.Beta.Mult[[2]])) / sd(part.0.s$Rand.Beta.Mult[[2]])
(part.0.s$Div[9] - mean(part.0.s$Rand.Beta.Mult[[3]])) / sd(part.0.s$Rand.Beta.Mult[[3]])
(part.0.s$Div[10] - mean(part.0.s$Rand.Beta.Mult[[4]])) / sd(part.0.s$Rand.Beta.Mult[[4]])



##Spider family (functional guild) analysis
part.0.sf <- partition(canspidsf, levels = c("SAMPLE", "TREE.NO", "HAB_ST", 
                                             "SITE", "ECO"),
                       low.level = 2,
                       q = 0,
                       method = "sample",
                       perms = 1000)

summary(part.0.sf, p.value = "two-sided")
(part.0.sf$Div[2] - mean(part.0.sf$Rand.Alpha)) / sd(part.0.sf$Rand.Alpha)
(part.0.sf$Div[7] - mean(part.0.sf$Rand.Beta.Mult[[1]])) / sd(part.0.sf$Rand.Beta.Mult[[1]])
(part.0.sf$Div[8] - mean(part.0.sf$Rand.Beta.Mult[[2]])) / sd(part.0.sf$Rand.Beta.Mult[[2]])
(part.0.sf$Div[9] - mean(part.0.sf$Rand.Beta.Mult[[3]])) / sd(part.0.sf$Rand.Beta.Mult[[3]])
(part.0.sf$Div[10] - mean(part.0.sf$Rand.Beta.Mult[[4]])) / sd(part.0.sf$Rand.Beta.Mult[[4]])


##Code to generate the Standard effect sizes for each taxa and endpoint for
##richness
##Data comes from the PARTITIONR analyses;
##The SES allow us to compare relative strength of deviations across taxa
##(differences are scaled to sd of the null distribution; z-score)
##So that the larger the difference, the greater the deviation from the
##observed value
##This indicates that at local scales (within/among trees), deviations of the
##diversity of ants is greater than the deviations of spiders; this flips at
##the stand level, at which point deviations in spiders are greater than ants
##This indicates that -
##at local scales: ants are less aggregated within trees and more aggregated
##                 among trees than by chance than spiders
##at broad scales: spiders are more aggregated across stands, sites, and regions

##further, after correcting for differences in species pool size, diversities
##were different between ants and spiders, dictated by spatial scale
##Meaning, that the net outcome of species assembly processes is NOT consistent
## in terms of their effects on beta diversity across scales and between taxa
##Suggesting that variation in beta diversity across scales is driven in part by
##differences in mechanisms of community assembly, dispersal ability, or density-
##dependent interactions


##These SES values were from our run of the above randomizations, the SES values
##generated subsequently may be slightly different than those presented here
SESfigdat <- data.frame(SES = c(-0.01, 0.52, 3.94, 24.67, -19.60,
                                0, 1.27, 0.15, 15.23, -13.01,
                                0.92, 2.98, 5.97, 7.01, -6.56,
                                0.19, 0.36, 4.27, 7.12, -6.80),
                        Levels = as.factor(rep(c("Beta 4 (Region)", "Beta 3 (Site)",
                                                 "Beta 2 (Stand)", "Beta 1 (Among Trees)",
                                                 "Alpha (Within Trees)"), times = 4)),
                        Taxa = rep(c("Ants","Spiders"), each = 10),
                        Endpoint = rep(c("Taxonomic","Functional"), each = 5, times = 2))

SESfigdat$Endpoint = factor(SESfigdat$Endpoint, levels=c("Taxonomic","Functional"))

ggplot(SESfigdat, aes(x = Levels, y = SES, shape = Taxa, color = Taxa))+
  facet_wrap(.~Endpoint)+
  # geom_line(aes(group = Taxa))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_segment(aes(xend=Levels, y=0, yend=SES))+
  geom_point(size = 3)+
  ylab("SES (Diversity Deviations)")+
  scale_color_manual(values = c("#FE938C","#427AA1"))+
  coord_flip()+
  theme_nmds()+
  theme(legend.position = c(0.375,0.125),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size=14),
        strip.background = element_rect(color="black", fill="white"),
        axis.text.x = element_text(size = 14),
        legend.title = element_blank())

#Save
ggsave("SES.tiff",
       dpi = 300,
       width = 7.6,
       height = 4.2)














