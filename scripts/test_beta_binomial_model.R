# Try beta binomial model for overdispersion on mortality data 

library(rstan)
library(lme4)
library(sjstats)

dem <- read.csv("../data/dem-103018.csv")
head(dem)

# Fit Marina's original model 
m1.4 <- glmer(cbind(tot.mort, germ.proj - tot.mort) ~ Treat.Code * Subplot * PC.F+ (1+Subplot|Year/Plot/Species), family = binomial, dem, glmerControl(calc.derivs = F))
summary(m1.4)
overdisp(m1.4)

m1.4 <- glmer(cbind(tot.mort, germ.proj - tot.mort) ~ Treat.Code * Subplot * PC.F+ (1|Year/Plot/Species) + (1+Subplot|Plot), family = binomial, dem, glmerControl(calc.derivs = F))
summary(m1.4)
overdisp(m1.4)

# Fit Marina's saturated model
m1.full <- glmer(cbind(tot.mort, germ.proj - tot.mort) ~ Treat.Code * Subplot + (1|Year/Plot/Subplot/Species), family = binomial, dem, glmerControl(calc.derivs = F))
summary(m1.full)
