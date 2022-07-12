library(ggplot2); library(lme4); library(lmerTest)
library(MuMIn); library(car); library(blme)
library(sjPlot)
options(na.action = "na.fail")

##Read in dataset
stnddiv <- read.csv("Data/StandDivDat_3719.csv")
colnames(stnddiv)

##scale all variables, so that standardized regression coefficients can be obtained
stnddiv[,c(6:31)] <- scale(stnddiv[,c(6:31)])

#Ant models
##Taxonomic Diversity
###Multiplicative beta
ant.bm <- lm(AntBm ~ X6o5KM_GYRATE_MN + X6o5KM_GYRATE_AM + X6o5KM_PARA_MN + X6o5KM_CLUMPY + 
               TRich + Isothermality + Temp.DQ + Temp.Warm + PrecipSeas + Precip.WQ, 
             data = stnddiv)

ant.bm.d <- summary(model.avg(dredge(ant.bm,
                                     m.lim = c(0,6),
                                     rank = "AICc")))

ant.bm.d

ant.bm.best <- lm(AntBm ~ X6o5KM_GYRATE_AM + PrecipSeas + TRich, 
                  data = stnddiv)

summary(ant.bm.best)

###Alpha
ant.al <- lm(AntAl ~ X6o5KM_GYRATE_MN + X6o5KM_GYRATE_AM + X6o5KM_PARA_MN + X6o5KM_CLUMPY + 
               TRich + Isothermality + Temp.DQ + Temp.Warm + PrecipSeas + Precip.WQ, 
             data = stnddiv)

ant.al.d <- summary(model.avg(dredge(ant.al,
                                     m.lim = c(0,6),
                                     rank = "AICc")))
ant.al.d

ant.al.best <- lm(AntAl ~ PrecipSeas, data = stnddiv)

summary(ant.al.best)


##Functional diversity
###Multiplicative beta
ant.f.bm <- lm(AFBm ~ X6o5KM_GYRATE_MN + X6o5KM_GYRATE_AM + X6o5KM_PARA_MN + X6o5KM_CLUMPY + 
                 TRich + Isothermality + Temp.DQ + Temp.Warm + PrecipSeas + Precip.WQ, 
               data = stnddiv)

ant.f.bm.d <- summary(model.avg(dredge(ant.f.bm,
                                    m.lim = c(0,6),
                                    rank = "AICc")))

ant.f.bm.d

ant.f.bm.best <- lm(AFBm ~ Temp.Warm, data = stnddiv)
summary(ant.f.bm.best)


###Alpha
ant.f.al <- lm(AFAl ~ X6o5KM_GYRATE_MN + X6o5KM_GYRATE_AM + X6o5KM_PARA_MN + X6o5KM_CLUMPY + 
                 TRich + Isothermality + Temp.DQ + Temp.Warm + PrecipSeas + Precip.WQ, 
               data = stnddiv)

ant.f.al.d <- summary(model.avg(dredge(ant.f.al,
                                      m.lim = c(0,6),
                                      rank = "AICc")))
ant.f.al.d

ant.f.al.best <- (lm(AFAl ~ PrecipSeas + X6o5KM_GYRATE_AM, data = stnddiv))
summary(ant.f.al.best)


#Spider Models
##Taxonomic diversity
###Multiplicative Beta

spider.bm <- lm(SpiBm ~ X6o5KM_GYRATE_MN + X6o5KM_GYRATE_AM + X6o5KM_PARA_MN + X6o5KM_CLUMPY + 
                 TRich + Isothermality + Temp.DQ + Temp.Warm + PrecipSeas + Precip.WQ, 
             data = stnddiv)

spider.bm.d <- summary(model.avg(dredge(spider.bm,
                                  m.lim = c(0,6),
                                  rank = "AICc")))

spider.bm.d

spider.bm.best <- lm(SpiBm ~ Isothermality, data = stnddiv)

summary(spider.bm.best)

###Alpha
spider.al <- lm(SpiAl ~ X6o5KM_GYRATE_MN + X6o5KM_GYRATE_AM + X6o5KM_PARA_MN + X6o5KM_CLUMPY + 
                 TRich + Isothermality + Temp.DQ + Temp.Warm + PrecipSeas + Precip.WQ, 
               data = stnddiv)

spider.al.d <- summary(model.avg(dredge(spider.al,
                                      m.lim = c(0,6),
                                      rank = "AICc")))

spider.al.d

spider.al.best <- lm(SpiAl ~ 1, data = stnddiv)

summary(spider.al.best)


##Functional Diversity
###Multiplicative Beta

spider.f.bm <- lm(SFBm ~ X6o5KM_GYRATE_MN + X6o5KM_GYRATE_AM + X6o5KM_PARA_MN + X6o5KM_CLUMPY + 
                 TRich + Isothermality + Temp.DQ + Temp.Warm + PrecipSeas + Precip.WQ, 
               data = stnddiv)

spider.f.bm.d <- summary(model.avg(dredge(spider.f.bm,
                                    m.lim = c(0,6),
                                    rank = "AICc")))

spider.f.bm.d


spider.f.bm.best <- lm(SFBm ~ PrecipSeas + TRich, data = stnddiv)
summary(spider.f.bm.best)


###Alpha
spider.f.al <- lm(SFAl ~ X6o5KM_GYRATE_MN + X6o5KM_GYRATE_AM + X6o5KM_PARA_MN + X6o5KM_CLUMPY + 
                 TRich + Isothermality + Temp.DQ + Temp.Warm + PrecipSeas + Precip.WQ, 
               data = stnddiv)

spider.f.al.d <- summary(model.avg(dredge(spider.f.al,
                                      m.lim = c(0,6),
                                      rank = "AICc")))
spider.f.al.d

spider.f.al.best <- lm(SFAl ~ PrecipSeas + X6o5KM_PARA_MN, data = stnddiv)
summary(spider.f.al.best)

################################################################################
##Moran's I test for spatial autocorrelation
library(ape)
library(geosphere)
library(ncf)

colnames(stnddiv)
x = stnddiv$Long*10
y = stnddiv$Lat

#### IMPORTANT BELOW

##Ant BetaM - CLEAR
ncf.cor <- correlog(x, y, residuals(ant.bm.best),
                    increment=0.4,
                    resamp=999)

p.adjust(ncf.cor$p, method = "BH", n = length(ncf.cor$p))


##Ant Alpha - CLEAR
ncf.cor <- correlog(x, y, residuals(ant.al.best),
                    increment=0.4,
                    resamp=999)

p.adjust(ncf.cor$p, method = "BH", n = length(ncf.cor$p))


##Spider BetaM - CLEAR
ncf.cor <- correlog(x, y, residuals(spider.bm.best),
                    increment=0.4,
                    resamp=999)

p.adjust(ncf.cor$p, method = "BH", n = length(ncf.cor$p))


##Spider Alpha - CLEAR
ncf.cor <- correlog(x, y, residuals(spider.al.best),
                    increment=0.4,
                    resamp=999)

p.adjust(ncf.cor$p, method = "BH", n = length(ncf.cor$p))


##Ant F BetaM - CLEAR
ncf.cor <- correlog(x, y, residuals(ant.f.bm.best),
                    increment=0.4,
                    resamp=999)

p.adjust(ncf.cor$p, method = "BH", n = length(ncf.cor$p))


##Ant F Alpha - CLEAR
ncf.cor <- correlog(x, y, residuals(ant.f.al.best),
                    increment=0.4,
                    resamp=999)

p.adjust(ncf.cor$p, method = "BH", n = length(ncf.cor$p))


##Spider F BetaM - CLEAR
ncf.cor <- correlog(x, y, residuals(spider.f.bm.best),
                    increment=0.4,
                    resamp=999)

p.adjust(ncf.cor$p, method = "BH", n = length(ncf.cor$p))


##Spider F Alpha - CLEAR
ncf.cor <- correlog(x, y, residuals(spider.f.al.best),
                    increment=0.4,
                    resamp=999)

p.adjust(ncf.cor$p, method = "BH", n = length(ncf.cor$p))

