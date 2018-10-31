# Chapter 2 analysis 10-27-2018 #

rm(list=ls()) 

#### Load Libraries ####
library(Rmisc)
library(ggplot2)
library(plyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(arm)
library(coefplot)
library(reshape2)
library(tidyr)
library(cowplot)
library(car)
library(sjstats)

#### Load: Bootstraps ####
## load previous boots
# # CI boots
# load("Boots-aug18/CI/CI-1.Rdata") # mortality
# load("Boots-aug18/CI/CI-2.Rdata") # seed set
# load("Boots-aug18/CI/CI-3.Rdata") # seed bank carryover
# load("Boots-aug18/CI/CI-4.Rdata") # germination
# load("Boots-aug18/CI/CI-6.Rdata") # plot level lambda
# 
# # Parameter boots (need to redo mortality with set seed sometime)
# load("Boots-aug18/Param-sims/BS-m.Rdata") # mortality
# load("Boots-aug18/Param-sims/BS-F.Rdata") # seed set
# load("Boots-aug18/Param-sims/BS-sb.Rdata") # seed bank carryover
# load("Boots-aug18/Param-sims/BS-g.Rdata") # germination
# 
# # Parameter boots without problem species for Germ and SB
# load("Boots-aug18/Param-sims/BS-sb-AG.Rdata") # seed bank carryover
# load("Boots-aug18/Param-sims/BS-g-CL.Rdata") # germination

#### Prep: Vital Rates ####

###
# Grass
###
grass <- read.csv("Data/Cleaned Data for Analysis/grass.csv")
names(grass)[3] <- "cover.g"

treat <- read.csv("Data/Marina-Treatment-30.csv")

###
# Germ and mortality data. plots 85, 86, 87 have been removed due to burrow damage
###
dem.16 <- read.csv("Data/Cleaned Data for Analysis/dem-data-16.csv")
dem.16 <- filter(dem.16, Plot != 87, Plot != 85, Plot != 86)
dem.17 <- read.csv("Data/Cleaned Data for Analysis/dem-data-17.csv")

dem <- rbind(dem.16, dem.17)
dem <- merge(dem, grass, by = c("Year", "Plot", "Subplot", "Treat.Code"), all.x = T)
dem <- dem[dem$Subplot != "Thatch", -c(8,9,12)] # get rid of thatch for now and extra columns
dem <- filter(dem, !(Year == 2017 & Treat.Code == "D")) # shelters in 2017 did not apply a drought and had the complete opposite effect on plants due to the warmer greenhouse effect
dem$Subplot <- factor(dem$Subplot, levels = c("No Grass", "Grass"))
dem$p.mort <- dem$tot.mort/dem$germ.proj
rm(dem.16, dem.17)

###
# Flowering/seed set data 
###

# 10/28/2018 NOTE: changed code for final-flo-seed.csv to take number of flowers per individual from date of highest number flowering instead of date of highest number setting seed, patterns are the same but with far fewer outliers)

flo.seed <- read.csv("Data/Cleaned Data for Analysis/final-flo-seed.csv")
flo.seed <- flo.seed[flo.seed$Subplot != "Thatch",]
flo.seed <- filter(flo.seed, !(Year == 2017 & Treat.Code == "D")) 
flo.seed$Subplot <- factor(flo.seed$Subplot, levels = c("No Grass", "Grass"))
flo.seed <- filter(flo.seed, n.seed.ind != 0) # only removes 3 scenarios
flo.seed$n.seed.ind <- flo.seed$n.seed.ind * flo.seed$p.viable # adjust for viability
flo.seed <- flo.seed[,-16]

###
# Seed survival
###
sb <- read.csv("Data/Cleaned Data for Analysis/seed-carryover-plot.csv")[,c(1:9)]
sb <- filter(sb, !(Year == 2017 & Treat.Code == "D"))
names(sb)[1] <- "Species"
sb <- sb[,-c(5,6)]

#### Prep: Traits ####
trait.w <- read.csv("Data/Cleaned Data for Analysis/final-traits-w.csv")
trait.w <- trait.w[,-c(5,7,9)]

sla.13 <- read.csv("Data/Cleaned Data for Analysis/final-sla-13c.csv")[,c(1,3,5)]

PCA.G <- prcomp(trait.w[,c(18, 20)], scale = T) # Greenhouse trait
PCA.F <- prcomp(trait.w[, c(7, 8)], scale = T) # Field trait
PCA.s13 <- prcomp(sla.13[, c(2, 3)], scale = T) # Full SLA + D13C field

# biplot(PCA.G)
# biplot(PCA.F)
# biplot(PCA.s13)

# summary(PCA.G)
# summary(PCA.F)
# summary(PCA.s13)

trait.w$PC.G <- PCA.G$x[,1] 
trait.w$PC.F <- PCA.F$x[,1] 
trait.w$PC.s13 <- PCA.s13$x[c(2,3,5,6,7,9),1] 

# merge datasets with traits
dem <- merge(dem, trait.w[,c(1,23:25)], by = "Species")
sb <- merge(sb, trait.w[,c(1,23:25)], by = "Species")
flo.seed <- merge(flo.seed, trait.w[,c(1,23:25)], by = "Species")

rm(PCA.F, PCA.G, PCA.s13, sla.13, grass)

#### Prep: Plot-level Lambda ####

full <- merge(dem, flo.seed[,-c(6:13)], by = c("Year","Plot","Treat.Code", "Subplot","Species", "PC.F", "PC.G", "PC.s13"), all = T)

full <- merge(full, sb, by = c("Year","Plot","Treat.Code","Species","PC.F", "PC.G", "PC.s13"), all.x = T)

full$p.germ <- full$germ.tot/full$viable
full$L.sb <- full$p.surv*(1-full$p.germ)
full$L.sa <- ifelse(full$germ.proj > 0 , ((full$germ.proj - full$tot.mort)/full$germ.proj), 0)
full$L.seeds <- ifelse(full$n.seed.ind >= 0, full$L.sa * full$n.seed.ind, 0)

full$L <- ifelse(full$L.sa == 0,  full$L.sb, full$L.sa*full$L.seeds) #45 missing 


#### M1: Mortality ####
dem$Year <- as.factor(dem$Year)

###
# Overall
###

# Model with random slope for treatment and subplot previously determined to be best model

# Original Model
m1.4 <- glmer(cbind(tot.mort, germ.proj - tot.mort) ~ Treat.Code * Subplot + (1|Year) + (1|Plot) + (1|Species), family = binomial, dem, glmerControl(calc.derivs = F))
plot(fitted(m1.4), resid(m1.4))
qqnorm(resid(m1.4))
qqline(resid(m1.4), col = 2, lwd = 2, lty = 2)
summary(m1.4)
overdisp(m1.4)

# model with what I think are correct random effects but so many of them that it's absorbing most effects
m1.4 <- glmer(cbind(tot.mort, germ.proj - tot.mort) ~ Treat.Code * Subplot + (1|Year/Plot/Subplot/Species), family = binomial, dem, glmerControl(calc.derivs = F))
plot(fitted(m1.4), resid(m1.4))
qqnorm(resid(m1.4))
qqline(resid(m1.4), col = 2, lwd = 2, lty = 2)
summary(m1.4)
overdisp(m1.4)

# model with OLRE
m1.4 <- glmer(cbind(tot.mort, germ.proj - tot.mort) ~ Treat.Code * Subplot + (1|obs), family = binomial, dem, glmerControl(calc.derivs = F))
plot(fitted(m1.4), resid(m1.4))
qqnorm(resid(m1.4))
qqline(resid(m1.4), col = 2, lwd = 2, lty = 2)
summary(m1.4) # gives the exact same answers as the correct random effect model
overdisp(m1.4)


###
# Traits
###

# original model
m1.trait <- glmer(cbind(tot.mort, germ.proj - tot.mort) ~ Treat.Code * Subplot * PC.F + (1|Year) + (1|Plot), family = binomial, dem, glmerControl(calc.derivs = F))
plot(fitted(m1.trait), resid(m1.trait))
qqnorm(resid(m1.trait)) 
qqline(resid(m1.trait), col = 2, lwd = 2, lty = 2) 
summary(m1.trait) 
overdisp(m1.trait)

# model with what I think are correct random effects
m1.trait <- glmer(cbind(tot.mort, germ.proj - tot.mort) ~ Treat.Code * Subplot * PC.F + (1|Year/Plot/Subplot/Species), family = binomial, dem, glmerControl(calc.derivs = F))
plot(fitted(m1.trait), resid(m1.trait))
qqnorm(resid(m1.trait)) 
qqline(resid(m1.trait), col = 2, lwd = 2, lty = 2) 
summary(m1.trait) 
overdisp(m1.trait)

# model with OLRE
m1.trait <- glmer(cbind(tot.mort, germ.proj - tot.mort) ~ Treat.Code * Subplot * PC.F + (1|obs), family = binomial, dem, glmerControl(calc.derivs = F))
plot(fitted(m1.trait), resid(m1.trait))
qqnorm(resid(m1.trait)) 
qqline(resid(m1.trait), col = 2, lwd = 2, lty = 2) 
summary(m1.trait) 

## compare fitted and predicted values for different models
par(mfrow=c(1,2))
plot(fitted(m1.trait) ~ dem$p.mort, col="darkgrey", 
     xlim=c(0, 1), ylim=c(0, 1), 
     xlab="Y (response)", ylab="Fitted Values")
abline(a=0, b=1, col="red")
plot(fitted(m1.trait2) ~ dem$p.mort, col="darkgrey", 
     xlim=c(0, 1), ylim=c(0, 1), 
     xlab="Y (response)", ylab="Fitted Values")
abline(a=0, b=1, col="red")
summary(m1.trait2) 

library(DHARMa)
require(DHARMa)

fittedModel <- m1.trait
simulationOutput <- simulateResiduals(fittedModel = fittedModel)
plot(simulationOutput)

fittedModel <- m1.trait2
simulationOutput <- simulateResiduals(fittedModel = fittedModel)
plot(simulationOutput)


#### M2: Seed Set ####

###
# Overall
###

# Model with random slope for Subplot previously determined to be best model
m2.3 <- lmer(log(n.seed.ind) ~ Treat.Code * Subplot + (1|Year) + (1|Plot) + (1 + Subplot|Species), flo.seed)
plot(fitted(m2.3), resid(m2.3))
qqnorm(resid(m2.3))
qqline(resid(m2.3), col = 2, lwd = 2, lty = 2)
summary(m2.3)

###
# Traits 
###

# Model without 3 way interaction previously determined to be best model
m2.trait <- lmer(log(n.seed.ind) ~ Treat.Code * Subplot * PC.F - Treat.Code : Subplot : PC.F  + (1|Plot) + (1|Year), flo.seed) 
plot(fitted(m2.trait), resid(m2.trait))
qqnorm(resid(m2.trait)) 
qqline(resid(m2.trait), col = 2, lwd = 2, lty = 2) 
summary(m2.trait)

#### M3: Seed Carryover ####

m3.trait <- glmer(cbind(Count, avg.num.per.sub - Count) ~  PC.F + (1|Year) + (1|Plot), family = binomial, data = sb)
plot(fitted(m3.trait), resid(m3.trait))
qqnorm(resid(m3.trait))
qqline(resid(m3.trait), col = 2, lwd = 2, lty = 2)
summary(m3.trait)
hist(resid(m3.trait))

# Model without AGHE previously determined to change sign of relationship; all species tested
m3.trait.AG <- glmer(cbind(Count, avg.num.per.sub - Count) ~  PC.F + (1|Year) + (1|Plot), family = binomial, data = sb[sb$Species != "Agoseris heterophylla",])
plot(fitted(m3.trait.AG), resid(m3.trait.AG))
qqnorm(resid(m3.trait.AG))
qqline(resid(m3.trait.AG), col = 2, lwd = 2, lty = 2)
summary(m3.trait.AG) 

#### M4: Germination ####

m4.trait <- glmer(cbind(germ.tot, viable-germ.tot) ~ PC.F + (1|Year) + (1|Plot/obs), family = binomial, data = dem)
plot(fitted(m4.trait), resid(m4.trait))
qqnorm(resid(m4.trait))
qqline(resid(m4.trait), col = 2, lwd = 2, lty = 2)
summary(m4.trait)

hist(resid(m4.trait))

m4.trait <- glmer(germ.tot ~ PC.F + (1|Year) + (1|Plot), family = poisson, weights = viable, data = dem, glmerControl(calc.derivs = F, optCtrl=list(maxfun=2e5), optimizer="bobyqa"))
plot(fitted(m4.trait), resid(m4.trait))
qqnorm(resid(m4.trait))
qqline(resid(m4.trait), col = 2, lwd = 2, lty = 2)
summary(m4.trait)

hist(resid(m4.trait))
overdisp(m4.trait)

overdisp_fun <- function(model) {
    rdf <- df.residual(model)
    rp <- residuals(model,type="pearson")
    Pearson.chisq <- sum(rp^2)
    prat <- Pearson.chisq/rdf
    pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
    c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(m4.trait)
m4.trait2 <- glmer.nb(germ.tot ~ PC.F + (1|Year) + (1|Plot), data = dem, glmerControl(calc.derivs = F, optCtrl=list(maxfun=2e5), optimizer="bobyqa"))
plot(fitted(m4.trait2), resid(m4.trait2))
qqnorm(resid(m4.trait2))
qqline(resid(m4.trait2), col = 2, lwd = 2, lty = 2)
summary(m4.trait2)

hist(resid(m4.trait2))

overdisp(m4.trait)

dem$obs <- 1:nrow(dem)
m4.trait <- glmer(germ.tot ~ PC.F + (1|Year) + (1|Plot/obs), family = poisson, weights = viable, data = dem, glmerControl(calc.derivs = F, optCtrl=list(maxfun=2e5), optimizer="bobyqa"))
plot(fitted(m4.trait), resid(m4.trait))
qqnorm(resid(m4.trait))
qqline(resid(m4.trait), col = 2, lwd = 2, lty = 2)
summary(m4.trait)

library(glmmTMB)
m4.bb <- glmmTMB(p.germ ~ PC.F + (1|Year) + (1|Plot), family = betabinomial, data = dem, glmmTMBControl(optCtrl=list(maxfun=2e5)))
plot(fitted(m4.bb), resid(m4.bb))
qqnorm(resid(m4.bb))
qqline(resid(m4.bb), col = 2, lwd = 2, lty = 2)
summary(m4.bb)

# ## Influential points - no few outliers driving the pattern in clarkia
# m4.cd <- cooks.distance(m4.trait)
# dem.inf <- merge(dem[,c(1:6,9, 12)], m4.cd, by = "row.names")
# test <- filter(dem[,c(1:6,9, 12)], !(Species == "Clarkia purpurea" & Treat.Code == "C"))
# m4.trait <- glmer(cbind(germ.tot, viable-germ.tot) ~ PC.F + (1|Year) + (1|Plot), family = binomial, data = test)
# plot(fitted(m4.trait), resid(m4.trait))
# qqnorm(resid(m4.trait))
# qqline(resid(m4.trait), col = 2, lwd = 2, lty = 2)
# summary(m4.trait)
# 
# ## dealing with non-normality in residuals
# plot(resid(m4.trait) ~ dem$PC.F)
# dem$p.germ <- dem$germ.tot/dem$viable
# hist(asin(sqrt(dem$p.germ)))
# hist(asin((dem$p.germ)^(1/3)))
# m4.trait <- lmer(asin((dem$p.germ)^(1/3)) ~ PC.F * Treat.Code + (1|Year) + (1|Plot), data = dem)
# plot(fitted(m4.trait), resid(m4.trait))
# qqnorm(resid(m4.trait))
# qqline(resid(m4.trait), col = 2, lwd = 2, lty = 2)
# summary(m4.trait)
# 
# BS.g <- bootMer(m4.trait, function(x) predict(x, newdata = df.bt, type = 'response', allow.new.levels = T, re.form = NA), nsim = 100, ncpus = 2, verbose = T) # 49/1000, 4.9%
# 
# compare model predictions

new <- cbind(df.bt, data.frame(pr.germ = apply(BS.g2$t, 2, function(x) mean(x))), data.frame(bi.germ = apply(BS.g$t, 2, function(x) mean(x))))
new <- ddply(new, .(Treat.Code, PC.F), summarize, pr.germ = mean(pr.germ), bi.germ = mean(bi.germ))
Germ.sp.tr <- ddply(dem, .(Species, Treat.Code), summarize, obs.germ = mean(p.germ))
new <- merge(new, trait.w[c(1,24)], by = "PC.F")
new <- merge(new, Germ.sp.tr, by = c("Species", "Treat.Code"))

new <- melt(new, id.vars = c("Species", "PC.F", "Treat.Code"))
ggplot(new[new$variable != "pr.germ",], aes(x = PC.F, y = value, col = variable, group = variable)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~Treat.Code)

## messing with year
# m4.trait <- glmer(cbind(germ.tot, viable-germ.tot) ~ PC.F, random = ~1|Year, family = binomial, data = dem)
# plot(fitted(m4.trait), resid(m4.trait))
# qqnorm(resid(m4.trait))
# qqline(resid(m4.trait), col = 2, lwd = 2, lty = 2)
# summary(m4.trait)
# 
# lattice::qqmath(m4.trait,id=0.1,idLabels=~.obs)
# 
# m.test <- rlmer(asin(sqrt(p.germ)) ~ PC.F + (1|Year) + (1|Plot), data = dem)
# summary(m.test)
# 
# m4.trait <- glmer(germ.tot ~ PC.F + (1|Year) + (1|Plot), family = poisson, data = dem, weights = viable)

# Model without CLPU previously determined to change sign of relationship; all species tested

m4.trait.CL <- glmer(cbind(germ.tot, viable-germ.tot) ~ PC.F * Treat.Code + (1|Year) + (1|Plot), family = binomial, data = dem[dem$Species != "Clarkia purpurea",])
plot(fitted(m4.trait.CL), resid(m4.trait.CL))
qqnorm(resid(m4.trait.CL))
qqline(resid(m4.trait.CL), col = 2, lwd = 2, lty = 2)
summary(m4.trait.CL) 

#### M6: Plot Lambda ####
hist(log(full$L + .8))
full$L.log <- log(full$L + .8)

###
# Overall
###

# Model with random slope for treatment previously determined to be the best model
m6.2 <- lmer(L.log ~ Treat.Code * Subplot + (1|Year) + (1|Plot) + (1 + Treat.Code|Species), full)
plot(fitted(m6.2), resid(m6.2))
qqnorm(resid(m6.2))
qqline(resid(m6.2), col = 2, lwd = 2, lty = 2)
summary(m6.2)


###
# Trait
###

# Model without 3-way interaction previously determined to be the best model
m6.trait <- lmer(L.log ~ Treat.Code * Subplot * PC.F - Treat.Code : Subplot : PC.F + (1|Year) + (1|Plot), full)
#m6.trait <- lmer(L.log ~ Treat.Code * Subplot * PC.F - Treat.Code : Subplot : PC.F + (1|Year) + (1|Plot), full[full$Species != "Clarkia purpurea",])
plot(fitted(m6.trait), resid(m6.trait))
qqnorm(resid(m6.trait))
qqline(resid(m6.trait), col = 2, lwd = 2, lty = 2)
summary(m6.trait)

#### Prep: Bootstrap CIs ####

df <- expand.grid(Treat.Code = c("D", "C", "W"), Subplot = c("Grass", "No Grass"), PC.F = seq(from = min(trait.w$PC.F), to = max(trait.w$PC.F), length.out = 1000)) 

###
# Mortality
###

boot1 <- bootMer(m1.trait, function(x) predict(x, newdata = df, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, ncpus = 2, verbose = T) # 50 did not converge (5%); according to Ben Bolker, this number isn't that bad (https://stackoverflow.com/questions/30235758/bootstrap-failed-using-mixed-model-in-lme4-package)

PI.boot <- data.frame(
           p.mort = apply(boot1$t, 2, function(x) as.numeric(quantile(x, probs = .5, na.rm=TRUE))),
           p.mort.lo = apply(boot1$t, 2, function(x) as.numeric(quantile(x, probs = .025, na.rm=TRUE))),
           p.mort.hi = apply(boot1$t, 2, function(x) as.numeric(quantile(x, probs = .975, na.rm=TRUE))))

CI.1 <- merge(PI.boot, df, by = "row.names")
save(CI.1, file = "Boots-aug18/CI/CI-1.Rdata")

###
# Seed Set
###

boot2 <- bootMer(m2.trait, function(x) predict(x, newdata = df, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, verbose = T)

PI.boot <- data.frame(
           n.seed.ind = apply(boot2$t, 2, function(x) as.numeric(quantile(x, probs = .5, na.rm=TRUE))),
           seed.lo = apply(boot2$t, 2, function(x) as.numeric(quantile(x, probs = .025, na.rm=TRUE))),
           seed.hi = apply(boot2$t, 2, function(x) as.numeric(quantile(x, probs = .975, na.rm=TRUE))))

CI.2 <- merge(PI.boot, df, by = "row.names")
save(CI.2, file = "Boots-aug18/CI/CI-2.Rdata")

###
# Seed Carryover
###

boot3 <- bootMer(m3.trait, function(x) predict(x, newdata = df, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, verbose = T) # 3/1000 did not converge

PI.boot <- data.frame(p.surv = apply(boot3$t, 2, function(x) mean(x)),
           sb.surv.lo = apply(boot3$t, 2, function(x) as.numeric(quantile(x, probs=.05, na.rm=TRUE))),
           sb.surv.hi = apply(boot3$t, 2, function(x) as.numeric(quantile(x, probs=.95, na.rm=TRUE))))

CI.3 <- merge(PI.boot, df, by = "row.names")
save(CI.3, file = "Boots-aug18/CI/CI-3.Rdata")

###
# Germination
###

boot4 <- bootMer(m4.trait, function(x) predict(x, newdata = df, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, verbose = T)  # more than 50 didn't converge...

PI.boot <- data.frame(germ = apply(boot4$t, 2, function(x) mean(x)),
           germ.lo = apply(boot4$t, 2, function(x) as.numeric(quantile(x, probs=.05, na.rm=TRUE))),
           germ.hi = apply(boot4$t, 2, function(x) as.numeric(quantile(x, probs=.95, na.rm=TRUE))))

CI.4 <- merge(PI.boot, df, by = "row.names")
save(CI.4, file = "Boots-aug18/CI/CI-4.Rdata")

###
# Plot Lambda
###

boot6 <- bootMer(m6.trait, function(x) predict(x, newdata = df, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, ncpus = 2, verbose = T)

PI.boot <- data.frame(
           L = apply(boot6$t, 2, function(x) as.numeric(quantile(x, probs = .5, na.rm=TRUE))),
           lower = apply(boot6$t, 2, function(x) as.numeric(quantile(x, probs = .025, na.rm=TRUE))),
           upper = apply(boot6$t, 2, function(x) as.numeric(quantile(x, probs = .975, na.rm=TRUE))))

CI.6 <- merge(PI.boot, df, by = "row.names")
save(CI.6, file = "Boots-aug18/CI/CI-6.Rdata")

rm(boot1, boot2, boot3, boot4, boot6, PI.boot)

#### Prep: Bootstrap Vital Rates ####

# Dataframe for parameter boots
df.bt <- expand.grid(Treat.Code = c("D", "C", "W"), Subplot = c("Grass", "No Grass"), PC.F = unique(trait.w$PC.F)) 

###
# Mortality
###
BS.m <- bootMer(m1.trait, function(x) predict(x, newdata = df.bt, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, ncpus = 2) #52/1000 = 5.2% convergence warnings

save(BS.m, file = "Boots-aug18/Param-sims/BS-m.Rdata")

###
# Seed set
###
BS.F <- bootMer(m2.trait, function(x) predict(x, newdata = df.bt, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, ncpus = 2, verbose = T)

save(BS.F, file = "Boots-aug18/Param-sims/BS-F.Rdata")

###
# Seed carryover
###

# with AGHE
BS.sb <- bootMer(m3.trait, function(x) predict(x, newdata = df.bt, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, ncpus = 2)

save(BS.sb, file = "Boots-aug18/Param-sims/BS-sb.Rdata")

# without AGHE
BS.sb.AG <- bootMer(m3.AG, function(x) predict(x, newdata = df.bt, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, ncpus = 2)

save(BS.sb.AG, file = "Boots-aug18/Param-sims/BS-sb-AG.Rdata")

###
# Germination
###

# with CLPU
BS.g <- bootMer(m4.trait, function(x) predict(x, newdata = df.bt, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, ncpus = 2) # 49/1000, 4.9%

save(BS.g, file = "Boots-aug18/Param-sims/BS-g.Rdata")

# without CLPU
BS.g.CL <- bootMer(m4.CL, function(x) predict(x, newdata = df.bt, type = 'response', allow.new.levels = T, re.form = NA), nsim = 1000, ncpus = 2, verbose = T) # 47/1000, 4.7%

save(BS.g.CL, file = "Boots-aug18/Param-sims/BS-g-CL.Rdata")

#### Prep: Boots L (all models) ####
# NOTE: boots and sims give essentially the same answer, for some reason sims Lambdas are much higher though, I honestly prefer boots though more because I dont get the weird warnings, I calculate them all at once, and there is not "NA" issue

# Create dataframe of parameter means and sds
BS <- list(BS.m, BS.F, BS.g, BS.sb)
params <- c("p.mort.b", "seed.set.b", "p.germ.b", "p.sb.sur.b", "p.mort.sd",  "seed.set.sd", "p.germ.sd", "p.sb.sur.sd")
pm <- as.data.frame(matrix(NA, 36, length(params)))
names(pm) = params
pm <- cbind(df.bt, pm)

for(i in 1:4){
    pm[, i+3] = apply(BS[[i]]$t, 2, function(x) mean(x))
    pm[, i+7] = apply(BS[[i]]$t, 2, function(x) sd(x))
}

# Revalue levels for creating scenarios and matching later
pm <- merge(pm, trait.w[,c(1,24)], by = "PC.F")
pm$Species <- revalue(pm$Species, c("Agoseris heterophylla" =  "AGHE", "Calycadenia pauciflora" = "CAPA", "Clarkia purpurea" = "CLPU", "Hemizonia congesta" = "HECO", "Lasthenia californica" = "LACA", "Plantago erecta" = "PLER"))
pm$Subplot <- revalue(pm$Subplot, c("No Grass" = "N", "Grass" = "G"))
pm$Scenario <- paste(pm$Treat.Code, pm$Subplot, pm$Species, sep = ".")
Scenario <- pm$Scenario
later <- pm[,c(2,3,1,12,13)]

# Create empty Lambda dataframe to fill
L.sim <- as.data.frame(Scenario)
sims = 1001
L.sim[,c(2:sims)] <- NA

# Other dfs for storage and later merger
g.sim <- L.sim
m.sim <- L.sim
s.sim <- L.sim
F.sim <- L.sim

# Calculate L from random draws of parameter distributions
for(j in 2:sims){
  for(i in 1:length(Scenario)){
    g.sim[i,j] <- rbinom(n = 1, size = 100, prob = pm[i, "p.germ.b"])/100 #germ
    s.sim[i,j] <- rbinom(n = 1, size = 100, prob = pm[i, "p.sb.sur.b"])/100 #seed surv
    m.sim[i,j] <- rnorm(n = 1, mean = pm[i, "p.mort.b"], sd = pm[i, "p.mort.sd"])
    F.sim[i,j] <- exp(rnorm(n = 1, mean = pm[i, "seed.set.b"], sd = pm[i, "seed.set.sd"])) #fecund
    L.sim[i,j] = s.sim[i,j]*(1 - g.sim[i,j]) + g.sim[i,j]*(1 - m.sim[i,j])*F.sim[i,j]
  }
}

# melt all dataframes
sim.list <- list(L.sim, g.sim, s.sim, m.sim, F.sim)
for(i in 1:5){ 
  sim.list[[i]] <- melt(sim.list[[i]])
}

# merge all dataframes
names(sim.list[[1]])[3] <- "L"
names(sim.list[[2]])[3] <- "g"
names(sim.list[[3]])[3] <- "s"
names(sim.list[[4]])[3] <- "m"
names(sim.list[[5]])[3] <- "Fe"

L.sim <- merge(sim.list[[1]], sim.list[[2]], by = c("Scenario", "variable"))
L.sim <- merge(L.sim, sim.list[[3]], by = c("Scenario", "variable"))
L.sim <- merge(L.sim, sim.list[[4]], by = c("Scenario", "variable"))
L.sim <- merge(L.sim, sim.list[[5]], by = c("Scenario", "variable"))
L.sim <- merge(L.sim, later, by = "Scenario")

B.L <- L.sim[,c(8:11,3:7)]

B.L <- ddply(B.L, .(Species, Treat.Code, Subplot, PC.F), summarize, L = mean(L), g = mean(g), sb = mean(s), m = mean(m), Fe = mean(Fe))

rm(g.sim, s.sim, m.sim, F.sim, params, pm, i, j, later, Scenario, sim.list, sims)

#### Prep: Boots L (no CL in germ) ####
# Note: Removing AG from seedbank model has no affect on lambdas, only removing CL from germ model

# Create dataframe of parameter means and sds
BS <- list(BS.m, BS.F, BS.g.CL, BS.sb)
params <- c("p.mort.b", "seed.set.b", "p.germ.b", "p.sb.sur.b", "p.mort.sd",  "seed.set.sd", "p.germ.sd", "p.sb.sur.sd")
pm <- as.data.frame(matrix(NA, 36, length(params)))
names(pm) = params
pm <- cbind(df.bt, pm)

for(i in 1:4){
    pm[, i+3] = apply(BS[[i]]$t, 2, function(x) mean(x))
    pm[, i+7] = apply(BS[[i]]$t, 2, function(x) sd(x))
}

# Revalue levels for creating scenarios and matching later
pm <- merge(pm, trait.w[,c(1,24)], by = "PC.F")
pm$Species <- revalue(pm$Species, c("Agoseris heterophylla" =  "AGHE", "Calycadenia pauciflora" = "CAPA", "Clarkia purpurea" = "CLPU", "Hemizonia congesta" = "HECO", "Lasthenia californica" = "LACA", "Plantago erecta" = "PLER"))
pm$Subplot <- revalue(pm$Subplot, c("No Grass" = "N", "Grass" = "G"))
pm$Scenario <- paste(pm$Treat.Code, pm$Subplot, pm$Species, sep = ".")
Scenario <- pm$Scenario
later <- pm[,c(2,3,1,12,13)]

# Create empty Lambda dataframe to fill
L.sim <- as.data.frame(Scenario)
sims = 1001
L.sim[,c(2:sims)] <- NA

# Other dfs for storage and later merger
g.sim <- L.sim
m.sim <- L.sim
s.sim <- L.sim
F.sim <- L.sim

# Calculate L from random draws of parameter distributions
for(j in 2:sims){
  for(i in 1:length(Scenario)){
    g.sim[i,j] <- rbinom(n = 1, size = 100, prob = pm[i, "p.germ.b"])/100 #germ
    s.sim[i,j] <- rbinom(n = 1, size = 100, prob = pm[i, "p.sb.sur.b"])/100 #seed surv
    m.sim[i,j] <- rnorm(n = 1, mean = pm[i, "p.mort.b"], sd = pm[i, "p.mort.sd"])
    F.sim[i,j] <- exp(rnorm(n = 1, mean = pm[i, "seed.set.b"], sd = pm[i, "seed.set.sd"])) #fecund
    L.sim[i,j] = s.sim[i,j]*(1 - g.sim[i,j]) + g.sim[i,j]*(1 - m.sim[i,j])*F.sim[i,j]
  }
}

# melt all dataframes
sim.list <- list(L.sim, g.sim, s.sim, m.sim, F.sim)
for(i in 1:5){ 
  sim.list[[i]] <- melt(sim.list[[i]])
}

# merge all dataframes
names(sim.list[[1]])[3] <- "L"
names(sim.list[[2]])[3] <- "g"
names(sim.list[[3]])[3] <- "s"
names(sim.list[[4]])[3] <- "m"
names(sim.list[[5]])[3] <- "Fe"

L.sim <- merge(sim.list[[1]], sim.list[[2]], by = c("Scenario", "variable"))
L.sim <- merge(L.sim, sim.list[[3]], by = c("Scenario", "variable"))
L.sim <- merge(L.sim, sim.list[[4]], by = c("Scenario", "variable"))
L.sim <- merge(L.sim, sim.list[[5]], by = c("Scenario", "variable"))
L.sim <- merge(L.sim, later, by = "Scenario")

B.L.CL <- L.sim[,c(8:11,3:7)]
B.L.CL <- ddply(B.L.CL, .(Species, Treat.Code, Subplot, PC.F), summarize, L = mean(L), g = mean(g), sb = mean(s), m = mean(m), Fe = mean(Fe))

rm(g.sim, s.sim, m.sim, F.sim, params, pm, i, j, later, Scenario, sim.list, sims)


### Prep: Simulate L (all models) ####

# create dataframe
sim.dem <- treat[,c(1,3,5)]
sim.dem$Subplot <- revalue(sim.dem$Subplot, c("A" = "Grass", "B" = "No Grass"))
sim.dem <- sim.dem[rep(seq_len(nrow(sim.dem)), each = 6),]
sim.dem$Species <- rep(trait.w$Species, 60)
sim.dem <- rbind(sim.dem, sim.dem)
sim.dem$Year <- rep(c(2016, 2017), each = 360)
sim.dem <- merge(sim.dem, trait.w[,c(1,24)], by = "Species")
sim.dem <- unique(sim.dem[, c(1,3,4,6)])
sim.dem$Year <- 2018
sim.dem$Plot <- 100

###
# Mortality
###

m1.t.sim <- simulate(m1.trait, nsim = 1000, newdata = sim.dem, re.form = NA, allow.new.levels = T) # still dont understand this warning, or rather i understand it but no idea why its giving it to me; contacted Ben Bolker about this issue 10/8/18

sim.dem.m <- as.data.frame(matrix(NA, nrow = 36, ncol = 1000))

for(i in 1:1000) {
  sim.dem.m[,i] <- m1.t.sim[[i]][,1]/(m1.t.sim[[i]][,1] + m1.t.sim[[i]][,2])
}

###
# Seed Set - No NaNs
###

sim.dem.s <- exp(simulate(m2.trait, nsim = 1000, newdata = sim.dem, re.form = NA, allow.new.levels = T))

###
# Seed Carryover - No NaNs
###
m3.t.sim <- simulate(m3.trait, nsim = 1000, newdata = sim.dem, re.form = NA, allow.new.levels = T)

sim.dem.b <- as.data.frame(matrix(NA, nrow = 36, ncol = 1000))

for(i in 1:1000) {
  sim.dem.b[,i] <- m3.t.sim[[i]][,1]/(m3.t.sim[[i]][,1] + m3.t.sim[[i]][,2])
}


###
# Germination - No NaNs
###

m4.t.sim <- simulate(m4.trait, nsim = 1000, newdata = sim.dem, re.form = NA, allow.new.levels = T)

sim.dem.g <- as.data.frame(matrix(NA, nrow = 36, ncol = 1000))

for(i in 1:1000) {
  sim.dem.g[,i] <- m4.t.sim[[i]][,1]/(m4.t.sim[[i]][,1] + m4.t.sim[[i]][,2])
}


###
# Put it together in L
###
# mortality sims: sim.dem.m
# germinatin sims: sim.dem.g
# seed carryover: sim.dem.b
# seed set: sim.dem.s

sims <- list(sim.dem.m, sim.dem.g, sim.dem.b, sim.dem.s)
#sims <- list(sim.dem.m, sim.dem.s)

L <- as.data.frame(matrix(NA, nrow = 36, ncol = 1000))


for(i in 1:1000){
  L[,i] <- sims[[3]][,i]*(1-sims[[2]][,i]) + sims[[2]][,i]*(1-sims[[1]][,i])*sims[[4]][,i]
}


L <- cbind(sim.dem, L)

L <- melt(L, id.vars = c("Species", "Treat.Code", "Subplot", "PC.F", "Plot", "Year"))

names(L)[8] <- "L.sim"

L <- L[,c(1,2,3,4,7,8)]

###
# Merge with other rates
###
# Prep germ sims

sim.dem.g.l <- cbind(sim.dem, sim.dem.g)

sim.dem.g.l <- melt(sim.dem.g.l, id.vars = c("Species", "Treat.Code", "Subplot", "PC.F", "Plot", "Year"))

names(sim.dem.g.l)[8] <- "g.sim"

sim.dem.g.l <- sim.dem.g.l[,c(1,2,3,4,7,8)]

# Prep Mortality sims
sim.dem.m.l <- cbind(sim.dem, sim.dem.m)

sim.dem.m.l <- melt(sim.dem.m.l, id.vars = c("Species", "Treat.Code", "Subplot", "PC.F", "Plot", "Year"))

names(sim.dem.m.l)[8] <- "m.sim"

sim.dem.m.l <- sim.dem.m.l[,c(1,2,3,4,7,8)]

# Prep bank data
sim.dem.b.l <- cbind(sim.dem, sim.dem.b)

sim.dem.b.l <- melt(sim.dem.b.l, id.vars = c("Species", "Treat.Code", "Subplot", "PC.F", "Plot", "Year"))

names(sim.dem.b.l)[8] <- "b.sim"

sim.dem.b.l <- sim.dem.b.l[,c(1,2,3,4,7,8)]

# Prep seed set data
sim.dem.s.l <- cbind(sim.dem, sim.dem.s)

sim.dem.s.l <- melt(sim.dem.s.l, id.vars = c("Species", "Treat.Code", "Subplot", "PC.F", "Plot", "Year"))

names(sim.dem.s.l)[8] <- "s.sim"

sim.dem.s.l <- sim.dem.s.l[,c(1,2,3,4,7,8)]

# Merge with L data
sim.dem.full <- merge(L, sim.dem.b.l, by = c("Species", "Treat.Code", "Subplot", "PC.F", "variable"))

sim.dem.full <- merge(sim.dem.full, sim.dem.g.l, by = c("Species", "Treat.Code", "Subplot", "PC.F", "variable"))

sim.dem.full <- merge(sim.dem.full, sim.dem.m.l, by = c("Species", "Treat.Code", "Subplot", "PC.F", "variable"))

sim.dem.s.l$variable <- gsub("sim_", "V", sim.dem.s.l$variable)
sim.dem.full <- merge(sim.dem.full, sim.dem.s.l, by = c("Species", "Treat.Code", "Subplot", "PC.F", "variable"))

sim.dem.full.sum <- ddply(sim.dem.full, .(Species, Treat.Code, Subplot, PC.F), summarize, L = mean(L.sim, na.rm = T))
sim.dem.full.sum$Subplot <- revalue(sim.dem.full.sum$Subplot, c("No Grass" = "N", "Grass" = "G"))

###  Prep: Simulate L (no CL in germ) #### ####

# create dataframe
sim.dem <- treat[,c(1,3,5)]
sim.dem$Subplot <- revalue(sim.dem$Subplot, c("A" = "Grass", "B" = "No Grass"))
sim.dem <- sim.dem[rep(seq_len(nrow(sim.dem)), each = 6),]
sim.dem$Species <- rep(trait.w$Species, 60)
sim.dem <- rbind(sim.dem, sim.dem)
sim.dem$Year <- rep(c(2016, 2017), each = 360)
sim.dem <- merge(sim.dem, trait.w[,c(1,24)], by = "Species")
sim.dem <- unique(sim.dem[, c(1,3,4,6)])
sim.dem$Year <- 2018
sim.dem$Plot <- 100

###
# Mortality
###

m1.t.sim <- simulate(m1.trait, nsim = 1000, newdata = sim.dem, re.form = NA, allow.new.levels = T) # still dont understand this warning, or rather i understand it but no idea why its giving it to me; contacted Ben Bolker about this issue 10/8/18

sim.dem.m <- as.data.frame(matrix(NA, nrow = 36, ncol = 1000))

for(i in 1:1000) {
  sim.dem.m[,i] <- m1.t.sim[[i]][,1]/(m1.t.sim[[i]][,1] + m1.t.sim[[i]][,2])
}

###
# Seed Set - No NaNs
###

sim.dem.s <- exp(simulate(m2.trait, nsim = 1000, newdata = sim.dem, re.form = NA, allow.new.levels = T))

###
# Seed Carryover - No NaNs
###
m3.t.sim <- simulate(m3.trait, nsim = 1000, newdata = sim.dem, re.form = NA, allow.new.levels = T)

sim.dem.b <- as.data.frame(matrix(NA, nrow = 36, ncol = 1000))

for(i in 1:1000) {
  sim.dem.b[,i] <- m3.t.sim[[i]][,1]/(m3.t.sim[[i]][,1] + m3.t.sim[[i]][,2])
}


###
# Germination - No NaNs
###

m4.t.sim <- simulate(m4.trait.CL, nsim = 1000, newdata = sim.dem, re.form = NA, allow.new.levels = T)

sim.dem.g <- as.data.frame(matrix(NA, nrow = 36, ncol = 1000))

for(i in 1:1000) {
  sim.dem.g[,i] <- m4.t.sim[[i]][,1]/(m4.t.sim[[i]][,1] + m4.t.sim[[i]][,2])
}


###
# Put it together in L
###
# mortality sims: sim.dem.m
# germinatin sims: sim.dem.g
# seed carryover: sim.dem.b
# seed set: sim.dem.s

sims <- list(sim.dem.m, sim.dem.g, sim.dem.b, sim.dem.s)
#sims <- list(sim.dem.m, sim.dem.s)

L <- as.data.frame(matrix(NA, nrow = 36, ncol = 1000))


for(i in 1:1000){
  L[,i] <- sims[[3]][,i]*(1-sims[[2]][,i]) + sims[[2]][,i]*(1-sims[[1]][,i])*sims[[4]][,i]
}


L <- cbind(sim.dem, L)

L <- melt(L, id.vars = c("Species", "Treat.Code", "Subplot", "PC.F", "Plot", "Year"))

names(L)[8] <- "L.sim"

L <- L[,c(1,2,3,4,7,8)]

###
# Merge with other rates
###
# Prep germ sims

sim.dem.g.l <- cbind(sim.dem, sim.dem.g)

sim.dem.g.l <- melt(sim.dem.g.l, id.vars = c("Species", "Treat.Code", "Subplot", "PC.F", "Plot", "Year"))

names(sim.dem.g.l)[8] <- "g.sim"

sim.dem.g.l <- sim.dem.g.l[,c(1,2,3,4,7,8)]

# Prep Mortality sims
sim.dem.m.l <- cbind(sim.dem, sim.dem.m)

sim.dem.m.l <- melt(sim.dem.m.l, id.vars = c("Species", "Treat.Code", "Subplot", "PC.F", "Plot", "Year"))

names(sim.dem.m.l)[8] <- "m.sim"

sim.dem.m.l <- sim.dem.m.l[,c(1,2,3,4,7,8)]

# Prep bank data
sim.dem.b.l <- cbind(sim.dem, sim.dem.b)

sim.dem.b.l <- melt(sim.dem.b.l, id.vars = c("Species", "Treat.Code", "Subplot", "PC.F", "Plot", "Year"))

names(sim.dem.b.l)[8] <- "b.sim"

sim.dem.b.l <- sim.dem.b.l[,c(1,2,3,4,7,8)]

# Prep seed set data
sim.dem.s.l <- cbind(sim.dem, sim.dem.s)

sim.dem.s.l <- melt(sim.dem.s.l, id.vars = c("Species", "Treat.Code", "Subplot", "PC.F", "Plot", "Year"))

names(sim.dem.s.l)[8] <- "s.sim"

sim.dem.s.l <- sim.dem.s.l[,c(1,2,3,4,7,8)]

# Merge with L data
sim.dem.full <- merge(L, sim.dem.b.l, by = c("Species", "Treat.Code", "Subplot", "PC.F", "variable"))

sim.dem.full <- merge(sim.dem.full, sim.dem.g.l, by = c("Species", "Treat.Code", "Subplot", "PC.F", "variable"))

sim.dem.full <- merge(sim.dem.full, sim.dem.m.l, by = c("Species", "Treat.Code", "Subplot", "PC.F", "variable"))

sim.dem.s.l$variable <- gsub("sim_", "V", sim.dem.s.l$variable)
sim.dem.full <- merge(sim.dem.full, sim.dem.s.l, by = c("Species", "Treat.Code", "Subplot", "PC.F", "variable"))

sim.dem.full.sum2 <- ddply(sim.dem.full, .(Species, Treat.Code, Subplot, PC.F), summarize, L = mean(L.sim, na.rm = T))
sim.dem.full.sum2$Subplot <- revalue(sim.dem.full.sum2$Subplot, c("No Grass" = "N", "Grass" = "G"))

#### Prep: Elasticity (all models) ####

# Turn mortality into survival first
B.L$sa <- 1 - B.L$m

# Elasticity with respect to germination
B.L$e.g <- (B.L$g * (B.L$sa * B.L$Fe - B.L$sb))/B.L$L

# Elasticity with respect to seedling survival SAME AS FECUNDITY so really dont know how to compare the two most important things
B.L$e.sa <- (B.L$g * B.L$sa * B.L$Fe)/B.L$L

# Elasticity with respect to seed survival
B.L$e.sb <- (B.L$sb * (1 - B.L$g))/B.L$L

# Elasticity with respect to Fecundity
B.L$e.F <- (B.L$g * B.L$sa * B.L$Fe)/B.L$L

###look at elasticity by species ###
elas.sum <- melt(B.L[,c(1:4,11:14)], id.vars = c("Species", "Treat.Code", "Subplot", "PC.F"))

# Fe and sa elas the same so remove one
elas.sum <- filter(elas.sum, variable != "e.sa")

elas.sum.max <- ddply(elas.sum, .(Species, Treat.Code, Subplot, PC.F), summarize, vital = variable[which.max(value)]) # unsurprisingly, all Fecundity... or mortality who knows

#### Prep: Sensitivity ####

# Sensitivity with respect to germination
B.L$s.g <- B.L$sa * B.L$Fe - B.L$sb

# Sensitivity with respect to seedling survival 
B.L$s.sa <- B.L$g * B.L$Fe

# Sensitivity with respect to seed survival
B.L$s.sb <- 1 - B.L$g

# Sensitivity with respect to Fecundity
B.L$s.F <- B.L$g * B.L$sa

###look at Sensitivity by stage ###
ggplot(B.L, aes(y = s.F, x = PC.F, col = Treat.Code, group = Treat.Code)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c("red1", "orange1", "mediumblue"), labels = c("Shelter", "Control", "Watered")) +
  labs(y = "Sensitivity Fecundity", x = "Drought Tolerance") +
  facet_wrap(~ Subplot) +
  geom_smooth(method = "lm", se = F)

ggplot(B.L, aes(y = s.sa, x = PC.F, col = Treat.Code, group = Treat.Code)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c("red1", "orange1", "mediumblue"), labels = c("Shelter", "Control", "Watered")) +
  labs(y = "Sensitivity to survival", x = "Drought Tolerance") +
  facet_wrap(~ Subplot) +
  geom_smooth(method = "lm", se = F)

ggplot(B.L, aes(y = s.g, x = PC.F)) +
  geom_point() +
  theme_classic() +
  labs(y = "Sensitivity to germination", x = "Drought Tolerance") +
  facet_wrap(~ Subplot) +
  geom_smooth(method = "lm", se = F)

ggplot(B.L, aes(y = s.sb, x = PC.F)) +
  geom_point() +
  theme_classic() +
  labs(y = "Sensitivity to seed survival", x = "Drought Tolerance") +
  facet_wrap(~ Subplot) +
  geom_smooth(method = "lm", se = F)

#### M7: Boot Lambda ####

###
# All models
###

B.L$Treat.Code <- factor(B.L$Treat.Code, levels = c("C", "D", "W"))
B.L$Subplot <- factor(B.L$Subplot, levels = c("N", "G"))
B.L$log.L <- log(B.L$L)
m.B1 <- lm(log.L ~ Treat.Code * Subplot * PC.F, B.L)
plot(fitted(m.B1), resid(m.B1))
qqnorm(resid(m.B1))
qqline(resid(m.B1), col = 2, lwd = 2, lty = 2) 
summary(m.B1)

# Get confidence intervals for lambda
B.L.CI <- predict(m.B1, newdata = df, interval = "confidence")
B.L.CI <- merge(B.L.CI, df, by = "row.names")

###
# Without Clarkia for germ model
###
B.L.CL$Treat.Code <- factor(B.L.CL$Treat.Code, levels = c("C", "D", "W"))
B.L.CL$Subplot <- factor(B.L.CL$Subplot, levels = c("N", "G"))
B.L.CL$log.L <- log(B.L.CL$L)
m.B2 <- lm(log.L ~ Treat.Code * Subplot * PC.F, B.L.CL)
plot(fitted(m.B2), resid(m.B2))
qqnorm(resid(m.B2))
qqline(resid(m.B2), col = 2, lwd = 2, lty = 2) 
summary(m.B2)

# Get confidence intervals for lambda
B.L.CL.CI <- predict(m.B2, newdata = df, interval = "confidence")
B.L.CL.CI <- merge(B.L.CL.CI, df, by = "row.names")

#### M8: Sim Lambda ####

# Get confidence intervals for lambda
L.CI <- predict(m8, newdata = df, interval = "confidence")
L.CI <- merge(L.CI, df, by = "row.names")

#### M9: Elasticity ####


#### M10: Sensitivity ####


#### Graph: Mortality Overall ####
# Summarize Data
dem.avg <- summarySE(dem, measurevar = "p.mort", groupvars = c("Treat.Code", "Subplot"), na.rm = T)
dem.avg$Treat.Code <- factor(dem.avg$Treat.Code, levels = c("D", "C", "W"))
dem.avg$Treat.Code <- revalue(dem.avg$Treat.Code, c("D" = "Drought", "C" = "Control", "W" = "Watered"))

# Plot Data
ggplot(dem.avg, aes(y = p.mort, x = Treat.Code, fill = Subplot)) +
  geom_bar(stat = "identity", position = position_dodge(width=0.9), col = "black") +
  geom_errorbar(aes(ymin = p.mort - se, ymax = p.mort + se), width = 0.2, position = position_dodge(width=0.9)) +
   theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.key.size = unit(2, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_fill_manual(values = c("mediumblue", "orange1"), labels = c("No Grass", "Grass")) +
  scale_y_continuous(labels = scales::percent, limits = c(0,.5)) +
  labs(y = "Forb Mortality") 

#### Graph: Mortality PC ####
dem.sum <- summarySE(dem, measurevar = "p.mort", groupvars = c("Treat.Code", "Subplot", "Species", "PC.F"), na.rm = T)

dem.sum$Treat.Code <- factor(dem.sum$Treat.Code, levels = c("D", "C", "W"))


# control and shelter
ggplot(dem.sum[dem.sum$Treat.Code != "W",], aes(y = p.mort, x = PC.F, col = Treat.Code, group = Treat.Code)) +
  geom_point() +
  geom_errorbar(aes(ymin = p.mort - se, ymax = p.mort + se)) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 10), 
        plot.title = element_text(size=30, face="bold", vjust = 2),
        axis.title = element_text(size = 15), 
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = c(0.12, 0.82),
        legend.key.size = unit(1.5, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("red1", "orange1"), labels = c("Shelter", "Control")) +
  scale_y_continuous(labels = scales::percent, limits = c(0,.85)) +
  labs(y = "Mortality", x = "Drought Tolerance") +
  facet_wrap(~ Subplot) +
  geom_line(data = CI.1[CI.1$Treat.Code != "W",], aes(y = p.mort)) +
  geom_ribbon(data = CI.1[CI.1$Treat.Code != "W",], aes(ymin = p.mort.lo, ymax = p.mort.hi), alpha = 0.2, linetype = 0) 

#ggsave(plot.m2.trait.d, filename = "Plots/mortality-trait-d.pdf", width = 150, height = 90, units = "mm", dpi = 600)

# control and water  
ggplot(dem.sum[dem.sum$Treat.Code != "D",], aes(y = p.mort, x = PC.F, col = Treat.Code, group = Treat.Code)) +
  geom_point() +
  geom_errorbar(aes(ymin = p.mort - se, ymax = p.mort + se)) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 10), 
        plot.title = element_text(size=30, face="bold", vjust = 2),
        axis.title = element_text(size = 15), 
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = c(0.12, 0.82),
        legend.key.size = unit(1.5, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("orange1", "mediumblue"), labels = c("Control", "Watered")) +
  scale_y_continuous(labels = scales::percent, limits = c(0,.8)) +
  labs(y = "Mortality", x = "Drought Tolerance") +
  facet_wrap(~ Subplot) +
  geom_line(data = CI.1[CI.1$Treat.Code != "D",], aes(y = p.mort)) +
  geom_ribbon(data = CI.1[CI.1$Treat.Code != "D",], aes(ymin = p.mort.lo, ymax = p.mort.hi), alpha = 0.2, linetype = 0) 

#### Graph: Seed Set Overall ####
flo.seed.avg <- summarySE(flo.seed, measurevar = "n.seed.ind", groupvars = c("Treat.Code", "Subplot"), na.rm = T)
flo.seed.avg$Treat.Code <- factor(flo.seed.avg$Treat.Code, levels = c("D", "C", "W"))
flo.seed.avg$Treat.Code <- revalue(flo.seed.avg$Treat.Code, c("D" = "Drought", "C" = "Control", "W" = "Watered"))

ggplot(flo.seed.avg, aes(y = n.seed.ind, x = Treat.Code, fill = Subplot)) +
  geom_bar(stat = "identity", position = position_dodge(width=0.9), col = "black") +
  geom_errorbar(aes(ymin = n.seed.ind - se, ymax = n.seed.ind + se), width = 0.2, position = position_dodge(width=0.9)) +
   theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 17),
        legend.text = element_text(size = 12),
        legend.position = "right",
        legend.key.size = unit(1.4, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_fill_manual(values = c("mediumblue", "orange1"), labels = c("No Grass", "Grass")) +
  labs(y = "Forb seed set per individual")

#### Graph: Seed Set PC ####
flo.seed.sum <- summarySE(flo.seed, measurevar = "n.seed.ind", groupvars = c("Treat.Code", "Subplot", "Species", "PC.F"), na.rm = T)
flo.seed.sum$Treat.Code <- factor(flo.seed.sum$Treat.Code, levels = c("D", "C", "W"))

ggplot(flo.seed.sum, aes(y = n.seed.ind, x = PC.F, col = Treat.Code, group = Treat.Code)) +
  geom_point() +
  geom_errorbar(aes(ymin = n.seed.ind - se, ymax = n.seed.ind + se)) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15), 
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 11),
        #legend.position = c(.63, .78),
        legend.key.size = unit(1.4, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("red1", "orange1", "mediumblue"), labels = c("Shelter", "Control", "Watered")) +
  labs(y = "Seed Set", x = "Drought Tolerance") +
  facet_wrap(~ Subplot) +
  geom_line(data = CI.2, aes(y = exp(n.seed.ind))) +
  geom_ribbon(data = CI.2, aes(ymin = exp(seed.lo), ymax = exp(seed.hi)), alpha = 0.2, linetype = 0) 

#### Graph: Germ Overall ####
dem.sum.yr <- ddply(dem, .(Species, Year, PC.F), summarize, p.germ = mean(germ.tot/viable))

ggplot(dem, aes(y = germ.tot/viable, x = factor(PC.F))) +
  geom_boxplot() +
  facet_wrap(~Year)

#### Graph: Germ PC ####
germ.sum <- summarySE(dem, measurevar = "p.germ", groupvars = c("Species", "PC.F", "Treat.Code"), na.rm = T)

CI.4 <- ddply(CI.4, .(PC.F, Treat.Code), summarize, p.surv = mean(p.sb.surv), sb.surv.lo = mean(sb.surv.lo), sb.surv.hi = mean(sb.surv.hi))

ggplot(germ.sum, aes(y = p.germ, x = PC.F, col = Treat.Code, group = Treat.Code)) +
  geom_point() +
  geom_errorbar(aes(ymin = p.germ - se, ymax = p.germ + se)) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 20), 
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.key.size = unit(2, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  labs(y = "Germination", x = "Drought Tolerance") +
  geom_line(data = CI.4, aes(y = germ)) +
  geom_ribbon(data = CI.4, aes(y = NULL, ymin = germ.lo, ymax = germ.hi), alpha = 0.2, linetype = 0) 

#### Graph: Seed Surv Overall ####


#### Graph: Seed Surv PC ####


#### Graph: Plot L Overall ####


#### Graph: Plot L PC ####

full.sum <- summarySE(full, groupvars = c("Treat.Code", "Subplot", "Species", "PC.F"), measurevar = "L", na.rm = T)

ggplot(full.sum, aes(y = L, x = PC.F, col = Treat.Code, group = Treat.Code)) +
  geom_point() +
  geom_errorbar(aes(ymin = L - se, ymax = L + se), width = 0.2) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 10), 
        plot.title = element_text(size=30, face="bold", vjust = 2),
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.key.size = unit(2, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("orange1", "red1", "mediumblue"), labels = c("Control", "Drought", "Watered")) +
  labs(y = "Lambda", x = "Drought Tolerance") +
  facet_wrap(~ Subplot) +
  geom_line(data = CI.6, aes(y = exp(L))) +
  geom_ribbon(data = CI.6, aes(y = NULL, ymin = exp(lower), ymax = exp(upper)), alpha = 0.2, linetype = 0) 

#### Graph: BS L PC ####
ggplot(B.L, aes(y = L, x = PC.F, col = Treat.Code, group = Treat.Code)) +
  geom_point() +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 10), 
        plot.title = element_text(size=30, face="bold", vjust = 2),
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.key.size = unit(2, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("orange1", "red1", "mediumblue"), labels = c("Control", "Drought", "Watered")) +
  labs(y = "B1 Lambda", x = "Drought Tolerance") +
  facet_wrap(~ Subplot) +
  geom_line(data = B.L.CI, aes(y = exp(fit)))

#### Graph: BS L PC (no g.CL) ####
ggplot(B.L.CL, aes(y = L, x = PC.F, col = Treat.Code, group = Treat.Code)) +
  geom_point() +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 10), 
        plot.title = element_text(size=30, face="bold", vjust = 2),
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.key.size = unit(2, 'lines'),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("orange1", "red1", "mediumblue"), labels = c("Control", "Drought", "Watered")) +
  labs(y = "B1 Lambda", x = "Drought Tolerance") +
  facet_wrap(~ Subplot) +
  geom_line(data = B.L.CL.CI, aes(y = exp(fit)))
