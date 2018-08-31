
###
# Dataframes needed
###

# dem = mortality and germination data
# flo.seed = seed set data
# sb = seed carryover in seedbank data
# sim.dem = data frame to store simulations

###
# Models needs
###
# m1.trait = mortality model
# m2.trait = seed set model
# m3.trait = seed carryover model
# m4 = germination model

###
# Mortality (NaNs produced here because of 0/0)
###

m1.t.sim <- simulate(m1.trait, nsim = 1000, seed = 1, newdata = sim.dem, re.form = NULL, allow.new.levels = T) # still dont understand this warning, or rather i understand it but no idea why its giving it to me

sim.dem.m <- as.data.frame(matrix(NA, nrow = 720, ncol = 1000))

for(i in 1:1000) {
  sim.dem.m[,i] <- m1.t.sim[[i]][,1]/(m1.t.sim[[i]][,1] + m1.t.sim[[i]][,2])
}


###
# Seed Set
###

sim.dem.s <- exp(simulate(m2.trait, nsim = 1000, seed=1, newdata = sim.dem, re.form = NULL, allow.new.levels = T))


###
# Seed Carryover 
###
m3.t.sim <- simulate(m3.trait, nsim = 1000, seed = 1, newdata = sim.dem, re.form = NULL, allow.new.levels = T)

sim.dem.b <- as.data.frame(matrix(NA, nrow = 720, ncol = 1000))

for(i in 1:1000) {
  sim.dem.b[,i] <- m3.t.sim[[i]][,1]/(m3.t.sim[[i]][,1] + m3.t.sim[[i]][,2])
}


###
# Germination
###

m4.t.sim <- simulate(m4, nsim = 1000, seed = 1, newdata = sim.dem, re.form = NULL, allow.new.levels = T)

sim.dem.g <- as.data.frame(matrix(NA, nrow = 720, ncol = 1000))

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

L <- as.data.frame(matrix(NA, nrow = 720, ncol = 1000))

for(i in 1:1000){
  L[,i] <- sims[[3]][,i]*(1-sims[[2]][,i]) + sims[[2]][,i]*(1-sims[[1]][,i])*sims[[4]][,i]
}

sim.dem$L.sim <- apply(L, 1, function(x) mean(x, na.rm = T)) #127 NAN because of zeroes, have to go back and figure that out

###
# Model L
###
## Usng simulated parameters ##

sim.dem$Subplot <- factor(sim.dem$Subplot, levels = c("No Grass", "Grass"))
sim.dem$Treat.Code <- factor(sim.dem$Treat.Code, levels = c("C", "D", "W"))

hist(log(sim.dem$L.sim))
m6.1 <- lmer(log(L.sim) ~ Treat.Code * Subplot * PC.F + (1|Year) + (1|Plot), sim.dem)
plot(fitted(m6.1), resid(m6.1))
qqnorm(resid(m6.1)) 
qqline(resid(m6.1), col = 2, lwd = 2, lty = 2) 
summary(m6.1)

###
# Graph it
###
sim.L.sum <- summarySE(sim.dem, measurevar = "L.sim", groupvars = c("Treat.Code", "Subplot", "PC.F"), na.rm = T)

ggplot(sim.L.sum, aes(y = L.sim, x = PC.F, col = Treat.Code, group = Treat.Code)) +
  geom_point() +
  geom_errorbar(aes(ymin = L.sim - se, ymax = L.sim + se), width = 0.2) + 
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
  geom_smooth(method = "lm", se = F)