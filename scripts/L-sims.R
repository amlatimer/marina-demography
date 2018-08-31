library(lme4)

###
# Dataframes needed
###

# dem = mortality and germination data
dem <- read.csv("../data/dem.csv")
# flo.seed = seed set data
flo.seed <- read.csv("../data/flo-seed.csv")
# sb = seed carryover in seedbank data
sb <- read.csv("../data/sb.csv")
# sim.dem = data frame to store simulations
sim.dem <- read.csv("../data/sim-dem.csv")


###
# Models needs
###
# m1.trait = mortality model
load("../data/m1-trait.rds")
# m2.trait = seed set model
load("../data/m2-trait.rds")
# m3.trait = seed carryover model
load("../data/m3-trait.rds")
# m4 = germination model
load("../data/m4.rds")


###
# Mortality (NaNs produced here because of 0/0)
###

summary(m1.trait)

## To get a bunch of simulations of one particular treatment combination, assuming that the uncertainty about this includes uncertainty about what the plot and year random effects are. 

# new data for one species in a particular combination of treatments, here watered, grassy plots in new unknown plot and year
testdata <- data.frame(Species = "Agoseris heterophylla", Plot = 100, Year = 2018, PC.F = -0.876, Subplot = "Grass", Treat.Code = "W")
testsim <- simulate(m1.trait, nsim=100, newdata = testdata, re.form=NA, allow.new.levels=T)
head(m1.testsim)
# Seems to work, although it's more complex because each simulation generates both a number successes (dead plants) and failures (survivors). 

## Simulate data for many combinations of treatments in the same way 

# new data for one species in all combinations of treatments
testdata <- data.frame(Species=rep("Agoseris heterophylla", 6), Plot=rep(100,6), Year = rep(2018, 6), PC.F = rep(-0.876, 6), Subplot = rep(c("Grass", "No Grass"), rep(3, 2)), Treat.Code = rep(c("W", "C", "D"), 2))
testsim <- simulate(m1.trait, nsim=100, newdata = testdata, re.form=NA, allow.new.levels=T)
head(m1.testsim)
# seems to work, but confusing to keep track of which treatments are in which row? 

## Just for fun, compare the mortality rates for Agoseris in Grassy plots for Drought vs Control conditions 
nsims <- 1000
testdata_drought <- data.frame(rep(Species = "Agoseris heterophylla", nsims), Plot = rep(100,  nsims), Year = rep(2018, nsims), PC.F = rep(-0.876, nsims), Subplot = rep("Grass", nsims), Treat.Code = rep("D", nsims))
testdata_control <- data.frame(rep(Species = "Agoseris heterophylla", nsims), Plot = rep(100, nsims), Year = rep(2018, nsims), PC.F = rep(-0.876, nsims), Subplot = rep("Grass", nsims), Treat.Code = rep("C", nsims))
testsim_drought <- simulate(m1.trait, nsim=1, newdata = testdata_drought, re.form = ~(1|Plot) + (1|Year), allow.new.levels=T)
testsim_control <- simulate(m1.trait, nsim=1, newdata = testdata_control, re.form = ~(1|Plot) + (1|Year), allow.new.levels=T)

sim.dem.drought <- testsim_drought[[1]][,1] / (testsim_drought[[1]][,1] + testsim_drought[[1]][,2])
sim.dem.control <- testsim_control[[1]][,1] / (testsim_control[[1]][,1] + testsim_control[[1]][,2])

boxplot(cbind(sim.dem.control, sim.dem.drought))
# Huh seems to work ok. Higher mortality rates in drought, but a lot of noise

# Make sure it isn't different if you specify random effects as NULL or as the model form, assuming you use new levels for all random effects 
testsim1 <- simulate(m1.trait, nsim=1, newdata = testdata_control, re.form = NULL, allow.new.levels=T)
testsim2 <- simulate(m1.trait, nsim=1, newdata = testdata_control, re.form = ~(1|Plot) + (1|Year), allow.new.levels=T)
boxplot(cbind(testsim1[[1]][,1], testsim2[[1]][,1]))
mean(testsim1[[1]][,1]); mean(testsim2[[1]][,1])
var(testsim1[[1]][,1]); var(testsim2[[1]][,1])
# seems to come out about the same, variance is similar and either one can be higher if you repeat the procedure a bunch of times. 

# So, I think the way to go is to make a data frame with one row per treatment combination, and in each row put in a new Plot number (e.g. 100) and a new Year number (e.g. 2018). Then simulate 1000 times or so, setting re.form to NULL and allow.new.levels to TRUE. Each simulation will produce a new random value for each treatment combination. Those simulated values can then become the inputs to calculate 1000 lambdas per treatment combination. Then those lambdas can be analyzed without regard to the random effects structure (since they are all for unknown Plots and Years). 



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