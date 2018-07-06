# this script tests a hierarchical bayes model for Marina's experiment data

library(brms)
library(rjags)

load("../data/m2b.rds")

traits <- read.csv("../data/trait_w.csv")
dem <- read.csv("../data/dem-070518.csv")
head(dem)
levels(dem$Treat.Code)
levels(dem$Species)
dem$species_id <- match(dem$Species, unique(dem$Species))
dem$plot_id <- match(dem$Plot, unique(dem$Plot))
dem$year_id <- match(dem$Year, unique(dem$Year))

jags_data <- list(N = nrow(dem), N_species = nrow(traits), N_plots = length(unique(dem$Plot)), N_years = length(unique(dem$Year)), y = dem$mort, n = dem$viable, grass = as.integer(dem$Subplot=="Grass"), drought = as.integer(dem$Treat.Code=="D"), water = as.integer(dem$Treat.Code=="W"), species = dem$species_id, year = dem$year_id, plot = dem$plot_id, pc1 = traits$PC1)

jags_inits <- list(a_species = rep(0, 6), a_plot = rep(0, 30), a_year = rep(0, 2),  b_drought = rep(0, 6), b_water = rep(0, 6), b_grass = rep(0, 6), b_grassxdrought = rep(0, 6), b_grassxwater = rep(0, 6), sd_a_species = 1, sd_a_plot = 1, sd_a_year = 1, sd_b_drought = 1,  sd_b_water = 1, sd_b_grass = 1, sd_b_grassxdrought = 1,  sd_b_grassxwater = 1, g_drought = rep(0, 2), g_water = rep(0, 2), g_grass = rep(0, 2), g_grassxdrought = rep(0, 2), g_grassxwater = rep(0, 2))

jags_model <- "../scripts/m2b_jags.txt"

m1 <- jags.model(file=jags_model, data=jags_data, inits=jags_inits, n.chains=3, n.adapt=1000)

m1_samp <- jags.samples(m1, c("a_species", "g_drought", "sd_b_drought", "g_water", "sd_b_water",  "g_grass", "sd_b_grass"), n.iter=2000)
(m1_samp)


jags.samples()


# Making model file for the equivalent of Marina's m2b model. 
m2b_jags <-"model{
  for (i in 1:N) {
    y[i] ~ dbinom(p[i], n[i])
    logit(p[i]) <- a_species[species[i]] + b_drought[species[i]]*drought[i] + b_water[species[i]]*water[i] + b_grass[species[i]]*grass[i] + b_grassxdrought[species[i]]*grass[i]*drought[i] + b_grassxwater[species[i]]*grass[i]*water[i] + a_plot[plot[i]] + a_year[year[i]]
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
  for (k in 1:N_years) { a_year[k] ~ dnorm(0, tau_a_year) }
  for (l in 1:N_plots) { a_plot[l] ~ dnorm(0, tau_a_plot) }
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
  for (k in 1:2) {
    g_drought[k] ~ dnorm(0, 0.1)
    g_water[k] ~ dnorm(0, 0.1)
    g_grass[k] ~ dnorm(0, 0.1)
    g_grassxwater[k] ~ dnorm(0, 0.1)
    g_grassxdrought[k] ~ dnorm(0, 0.1)
  }
}"

sink("../scripts/m2b_jags.txt")
cat(m2b_jags)
sink()

# Model file for the equivalent of Marina's m2b model plus trait regression for the treatment effects. 
m2b_plus_traits_jags <-"model{
for (i in 1:N) {
y[i] ~ dbinom(p[i], n[i])
logit(p[i]) <- alpha[species_id[i]] + beta_drought[species_id[i]]*drought[i]
}

for (j in 1:N_species) {
alpha[j] ~ dnorm(0, tau_alpha)
beta_drought[j] ~ dnorm(mu_beta[j], tau_beta)
mu_alpha[j] <- delta[1] + delta[2] * pc1[j]
}

# Priors
tau_alpha <- pow(sd_alpha, -2)
sd_alpha ~ dnorm(0, 0.1)I(0,)
tau_beta <- pow(sd_beta, -2)
sd_alpha ~ dnorm(0, 0.1)I(0,)
for (k in 1:2) {
delta[k] ~ dnorm(0, 0.1)
}
}"

sink("./scripts/jags_m1.txt")
cat("jags_m1")
sink()
