# this script tests a hierarchical bayes model for Marina's experiment data

library(brms)
library(rjags)

load("./data/m2b.rds")

traits <- read.csv("./data/trait_w.csv")
dem <- read.csv("./data/dem-070518.csv")
head(dem)
levels(dem$Treat.Code)
levels(dem$Species)
dem$species_id <- match(dem$Species, unique(dem$Species))


jags_data <- list(N = nrow(dem), N_species = nrow(traits), y = dem$mort, n = dem$viable, grass = as.integer(dem$Subplot=="Grass"), drought = as.integer(dem$Treat.Code=="D"), water = as.integer(dem$Treat.Code=="W"), species = dem$species_id, pc1 = traits$PC1)

jags_inits <- list(alpha = rep(0, 6), beta_drought = rep(0, 6), sd_alpha = 1, sd_beta = 1, gamma = rep(0, 2), delta = rep(0, 2))

jags_model <- "./scripts/jags_m1.txt"

m1 <- jags.model(file=jags_model, data=jags_data, inits=jags_inits, n.chains=3, n.adapt=1000)

m1_samp <- jags.samples(m1, c("alpha", "beta_drought", "delta", "sd_beta"), n.iter=2000)
(m1_samp)

jags.samples()


# Making model file 
jags_m1 <-"model{
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
