# this script tests a hierarchical bayes model for Marina's experiment data

library(brms)
library(rjags)

load("../data/m2b.rds")

traits <- read.csv("../data/trait_w.csv")
dem <- read.csv("../data/dem-070518.csv")
head(dem)
levels(dem$Treat.Code)
species <- levels(dem$Species)
dem$species_id <- match(dem$Species, species)
dem$plot_id <- match(dem$Plot, unique(dem$Plot))
#dem$year_id <- match(dem$Year, unique(dem$Year))

jags_data <- list(N = nrow(dem), N_species = nrow(traits), N_plots = length(unique(dem$Plot)), N_years = length(unique(dem$Year)), y = dem$mort, n = dem$viable, grass = as.integer(dem$Subplot=="Grass"), drought = as.integer(dem$Treat.Code=="D"), water = as.integer(dem$Treat.Code=="W"), species = dem$species_id, year = as.integer(dem$year_id==2), plot = dem$plot_id, pc1 = traits$PC1)

# Run jags equivalent of m2b model 

jags_inits_m2b <- list(alpha = 0, a_species = rep(0, 6), a_plot = rep(0, 30), b_year = 0,  b_drought = rep(0, 6), b_water = rep(0, 6), b_grass = rep(0, 6), b_grassxdrought = rep(0, 6), b_grassxwater = rep(0, 6), sd_a_species = 1, sd_a_plot = 1, sd_a_year = 1, sd_b_drought = 1,  sd_b_water = 1, sd_b_grass = 1, sd_b_grassxdrought = 1,  sd_b_grassxwater = 1, mu_b_drought = 0, mu_b_water = 0,  mu_b_grass = 0, mu_b_grassxdrought = 0, mu_b_grassxwater = 0)

jags_model <- "../scripts/m2b_jags.txt"

m2b_fit <- jags.model(file=jags_model, data=jags_data, inits=jags_inits_m2b, n.chains=3, n.adapt=1000)

m2b_samp <- coda.samples(m2b_fit, c("alpha", "a_species","sd_a_species", "a_plot", "sd_a_plot", "mu_b_drought", "mu_b_water", "mu_b_grass", "mu_b_grassxdrought", "mu_b_grassxwater", "b_drought", "b_grass", "b_water", "b_year", "sd_b_drought",  "sd_b_water", "sd_b_grass", "sd_b_grassxdrought", "sd_b_grassxwater"), n.iter=5000)
(m2b_samp)
gelman.diag(m2b_samp)
par(mar=rep(2, 4))
plot(m2b_samp)
summary(m2b_samp)

m2b_dic <- dic.samples(m2b_fit, n.iter=2000)
m2b_dic
# DIC: 
#Mean deviance:  4893 
#penalty 53.67 
#Penalized deviance: 4946 



test <- coda.samples(m2b_fit, "b_year", n.iter=100)
plot(test)

# Run jags equivalent of m2b model plus traits level

jags_inits_m2b_plus_traits <- list(a_species = rep(0, 6), a_plot = rep(0, 30), b_year = 0,  b_drought = rep(0, 6), b_water = rep(0, 6), b_grass = rep(0, 6), b_grassxdrought = rep(0, 6), b_grassxwater = rep(0, 6), sd_a_species = 1, sd_a_plot = 1, sd_a_year = 1, sd_b_drought = 1,  sd_b_water = 1, sd_b_grass = 1, sd_b_grassxdrought = 1,  sd_b_grassxwater = 1, g_drought = rep(0, 2), g_water = rep(0, 2), g_grass = rep(0, 2), g_grassxdrought = rep(0, 2), g_grassxwater = rep(0, 2))

jags_model <- "../scripts/m2b_plus_traits_jags.txt"

m2b_traits_fit <- jags.model(file=jags_model, data=jags_data, inits = jags_inits_m2b_plus_traits, n.chains=3, n.adapt=1000)

m2b_traits_samp <- coda.samples(m2b_traits_fit, c("a_species", "g_drought", "sd_b_drought", "g_water", "sd_b_water",  "g_grass", "sd_b_grass"), n.iter=2000)
gelman.diag(m2b_traits_samp)
summary(m2b_traits_samp)

m2b_traits_dic <- dic.samples(m2b_traits_fit, n.iter=2000)
m2b_traits_dic
#Mean deviance:  4891 
#penalty 55.46 
#Penalized deviance: 4947


# Making model file for the equivalent of Marina's m2b model. 
m2b_jags <-"model{
  for (i in 1:N) {
    y[i] ~ dbinom(p[i], n[i])
    logit(p[i]) <- a_species[species[i]] + b_drought[species[i]]*drought[i] + b_water[species[i]]*water[i] + b_grass[species[i]]*grass[i] + b_grassxdrought[species[i]]*grass[i]*drought[i] + b_grassxwater[species[i]]*grass[i]*water[i]  + a_plot[plot[i]] + b_year * year[i]
  }
  
  # Priors
  alpha ~ dnorm(0, 0.1)
  #a_species[1] <- 0
  for (j in 1:N_species) { a_species[j] ~ dnorm(0, tau_a_species) }
  #a_plot[1] <- 0
  for (j in 1:N_plots) { a_plot[j] ~ dnorm(0, tau_a_plot)  }
  #a_year[1] <- 0
  #a_year[2] ~ dnorm(0, 0.1)
  b_year ~ dnorm(0, 0.1)
  for (j in 1:N_species) { 
    b_drought[j] ~ dnorm(mu_b_drought, tau_b_drought) 
    b_water[j] ~ dnorm(mu_b_water, tau_b_water) 
    b_grass[j] ~ dnorm(mu_b_grass, tau_b_grass) 
    b_grassxdrought[j] ~ dnorm(mu_b_grassxdrought, tau_b_grassxdrought) 
    b_grassxwater[j] ~ dnorm(mu_b_grassxwater, tau_b_grassxwater) 
  }
  mu_b_water ~ dnorm(0, 0.01)
  mu_b_drought ~ dnorm(0, 0.01)
  mu_b_grass ~ dnorm(0, 0.01)
  mu_b_grassxwater ~ dnorm(0, 0.01)
  mu_b_grassxdrought ~ dnorm(0, 0.01)
  tau_a_species <- pow(sd_a_species, -2)
  sd_a_species ~ dnorm(0, 0.1)I(0,)
  tau_a_year <- pow(sd_a_year, -2)
  sd_a_year ~ dnorm(0, 0.1)I(0,)
  tau_a_plot <- pow(sd_a_plot, -2)
  sd_a_plot ~ dnorm(0, 0.1)I(0,)
  tau_b_drought <- pow(sd_b_drought, -2)
  sd_b_drought ~ dnorm(0, 0.1)I(0,)
  tau_b_water <- pow(sd_b_water, -2)
  sd_b_water ~ dnorm(0, 0.1)I(0,)
  tau_b_grass <- pow(sd_b_grass, -2)
  sd_b_grass ~ dnorm(0, 0.1)I(0,)
  tau_b_grassxdrought <- pow(sd_b_grassxdrought, -2)
  sd_b_grassxdrought ~ dnorm(0, 0.1)I(0,)
  tau_b_grassxwater <- pow(sd_b_grassxwater, -2)
  sd_b_grassxwater ~ dnorm(0, 0.1)I(0,)
}"

sink("../scripts/m2b_jags.txt")
cat(m2b_jags)
sink()

# Model file for the equivalent of Marina's m2b model plus trait regression for the treatment effects. 
m2b_plus_traits_jags <-"model{
  for (i in 1:N) {
y[i] ~ dbinom(p[i], n[i])
logit(p[i]) <- a_species[species[i]] + b_drought[species[i]]*drought[i] + b_water[species[i]]*water[i] + b_grass[species[i]]*grass[i] + b_grassxdrought[species[i]]*grass[i]*drought[i] + b_grassxwater[species[i]]*grass[i]*water[i] + a_plot[plot[i]] + b_year * year[i]
}

# Species-level regression of coefficients on traits
for (j in 1:N_species) {
b_drought[j] ~ dnorm(mu_b_drought[j], tau_b_drought)
mu_b_drought[j] <- g_drought[1] + g_drought[2] * pc1[j]
b_water[j] ~ dnorm(mu_b_water[j], tau_b_water)
mu_b_water[j] <- g_water[1] + g_water[2] * pc1[j]
b_grass[j] ~ dnorm(mu_b_grass[j], tau_b_grass)
mu_b_grass[j] <- g_grass[1] + g_grass[2] * pc1[j]
b_grassxdrought[j] ~ dnorm(mu_b_grassxdrought[j], tau_b_grassxdrought)
mu_b_grassxdrought[j] <- g_grassxdrought[1] + g_grassxdrought[2] * pc1[j]
b_grassxwater[j] ~ dnorm(mu_b_grassxwater[j], tau_b_grassxwater)
mu_b_grassxwater[j] <- g_grassxwater[1] + g_grassxwater[2] * pc1[j]
}

# Priors
for (j in 1:N_species) { a_species[j] ~ dnorm(0, tau_a_species) }
for (l in 1:N_plots) { a_plot[l] ~ dnorm(0, tau_a_plot) }
b_year ~ dnorm(0, 0.1)
tau_a_species <- pow(sd_a_species, -2)
sd_a_species ~ dnorm(0, 0.1)I(0,)
#tau_a_year <- pow(sd_a_year, -2)
#sd_a_year ~ dnorm(0, 0.1)I(0,)
tau_a_plot <- pow(sd_a_plot, -2)
sd_a_plot ~ dnorm(0, 0.1)I(0,)
tau_b_drought <- pow(sd_b_drought, -2)
sd_b_drought ~ dnorm(0, 0.1)I(0,)
tau_b_water <- pow(sd_b_water, -2)
sd_b_water ~ dnorm(0, 0.1)I(0,)
tau_b_grass <- pow(sd_b_grass, -2)
sd_b_grass ~ dnorm(0, 0.1)I(0,)
tau_b_grassxdrought <- pow(sd_b_grassxdrought, -2)
sd_b_grassxdrought ~ dnorm(0, 0.1)I(0,)
tau_b_grassxwater <- pow(sd_b_grassxwater, -2)
sd_b_grassxwater ~ dnorm(0, 0.1)I(0,)
for (k in 1:2) {
g_drought[k] ~ dnorm(0, 0.1)
g_water[k] ~ dnorm(0, 0.1)
g_grass[k] ~ dnorm(0, 0.1)
g_grassxwater[k] ~ dnorm(0, 0.1)
g_grassxdrought[k] ~ dnorm(0, 0.1)
}
}"
sink("../scripts/m2b_plus_traits_jags.txt")
cat(m2b_plus_traits_jags)
sink()
