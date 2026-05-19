targets.utils::tl()
library(greta)

# Modelling problem
# Have two snapshots rasters of proportion of improved housing ranging from 0—1
# Want to understand rate of change in order to predict over time.
# Assumptions
# - only increase in proportion of improved housing
# - assume logistic growth pattern
#
# Know the difference in time between two points but not where these points
# fall on curve

# simulation study

# Can I recapture the steepness from values where I don't know


logistic_function <- function(x, x_0 = 0, k = 1, L = 1){
  1 / (1 + exp(-k*(x - x_0)))
}


# example illustration
x_eg <- seq(
  from = -10,
  to = 10,
  by = 0.1
)

epsilon_eg <- rnorm(
  n = length(x_eg),
  mean = 0,
  sd = 1
)

k_eg <- 0.5

y_eg1 <- logistic_function(
  x = x_eg,
  k = k_eg
)


y_eg2 <- logistic_function(
  x = x_eg + epsilon_eg,
  k = k_eg
)


plot(
  x = x_eg,
  y = y_eg1
)
points(
  x = x_eg,
  y = y_eg2,
  col = "blue",
  add = TRUE
)


# simulated data

# number of sims
nsim <- 100

# simulate x1 data, i.e. first time point of observation - this will be unknown
# x1_sim <- runif(
#   n = nsim,
#   min = -10,
#   max = 5
# )
x1_sim <- rnorm(
  n = nsim,
  mean = 0,
  sd = 5
)


# simulate difference between successive x points - this will be known
# x_delta_sim <- runif(
#   n = nsim,
#   min = 1,
#   max = 10
# )

x_delta_sim <- 5

# simulate value of x2; second time point - this will be unknown
x2_sim <- x1_sim + x_delta_sim

plot(x1_sim, x2_sim)


# simulate steepness of curve - this will be unknown
k_sim <- 0.5


# simulation with noise on x

# # simulate noise
# epsilon_sd <- 1
#
# epsilon1_sim <- rnorm(
#   n = nsim,
#   mean = 0,
#   sd = epsilon_sd
# )
#
# epsilon2_sim <- rnorm(
#   n = nsim,
#   mean = 0,
#   sd = epsilon_sd
# )
#
#
# # simulate observed data
# # NB adding noise to x rather than y is necessary as it keeps it within range 0-1
# y1_sim <- logistic_function(
#   x = x1_sim + epsilon1_sim,
#   k = k_sim
# )
#
# y2_sim_raw <- logistic_function(
#   x = x2_sim + epsilon2_sim,
#   k = k_sim
# )
#
# # ensure y2 > y1 as noise might have introduced error
# y2_sim <- ifelse(
#   test = y2_sim_raw < y1_sim,
#   yes = y1_sim,
#   no = y2_sim_raw
# )

# # check out relationship, also highlight how many points shifted due to noise
# plot(y1_sim, y2_sim_raw, col = "red")
# points(y1_sim, y2_sim, col = "black")

# simulation with noise on y - better
y1_sim_pred <- logistic_function(
  x = x1_sim,
  k = k_sim
)

y2_sim_pred <- logistic_function(
  x = x2_sim,
  k = k_sim
)

# noise params
# mean of beta dist is alpha/(alpha + beta)
# this formulation means that the mean of y_sim = y_sim_pred
# variance is maximised in the middle and approaches 0 at the tails
# matching the variation we expect at the tails of the logistic function

phi_sim <- 30

alpha1_sim <- y1_sim_pred * phi_sim
alpha2_sim <- y2_sim_pred * phi_sim

beta1_sim <- (1 - y1_sim_pred) * phi_sim
beta2_sim <- (1 - y2_sim_pred) * phi_sim

y1_sim <- rbeta(
  n = nsim,
  shape1 = alpha1_sim,
  shape2 = beta1_sim
)

y2_sim <- rbeta(
  n = nsim,
  shape1 = alpha2_sim,
  shape2 = beta2_sim
)

plot(x = y1_sim_pred, y = y2_sim_pred)
points(x = y1_sim, y = y2_sim, col = "red")

########## model relationship

nobs <- nsim

x_delta_observed <- as_data(x_delta_sim)
y1_observed  <- as_data(y1_sim)
y2_observed  <- as_data(y2_sim)

k_latent <- exponential(1)

x1_latent <- normal(0, 5, dim = nobs)

x2_latent <- x1_latent + x_delta_observed

y1_pred <- 1/(1 + exp(-k_latent * (x1_latent)))
y2_pred <- 1/(1 + exp(-k_latent * (x2_latent)))

phi <- exponential(1)

alpha_1 <-y1_pred * phi
alpha_2 <- y2_pred * phi
beta_1  <- (1 - y1_pred) * phi
beta_2  <- (1 - y2_pred) * phi

distribution(y1_observed) <- beta(alpha_1, beta_1)
distribution(y2_observed) <- beta(alpha_2, beta_2)


m <- model(k_latent, phi)



draws <- mcmc(
  m,
  warmup = 1000,
  n_samples = 1000,
  chains = 4
)

summary(draws)

library(bayesplot)

mcmc_trace(draws)

coda::gelman.diag(draws, autoburnin = FALSE)

# ppc

# prior
prior <- calculate(
  y1_observed,
  nsim = 100
)

y1_prior <- prior$y1_observed[,,1] |>
  as.matrix()

ppc_dens_overlay(
  y = y1_sim,
  yrep = y1_prior
)

plot(
  x1_sim,
  x1_prior |>
    apply(
      X = _,
      MARGIN = 2,
      FUN = mean
    )
)

# posterior

posterior <- calculate(
  x1_latent,
  y1_observed,
  values = draws,
  nsim = 100
)

x1_posterior <- posterior$x1_latent[,,1] |>
  as.matrix()

ppc_dens_overlay(
  y = x1_sim,
  yrep = x1_posterior
)

plot(
  x1_sim,
  x1_posterior |>
    apply(
      X = _,
      MARGIN = 2,
      FUN = mean
    )
)

y1_posterior <- posterior$y1_observed[,,1] |>
  as.matrix()

ppc_dens_overlay(
  y = y1_sim,
  yrep = y1_posterior
)

plot(
  y1_sim,
  y1_posterior |>
    apply(
      X = _,
      MARGIN = 2,
      FUN = mean
    )
)

library(DHARMa)

y1_resids <- DHARMa::createDHARMa(
  simulatedResponse = t(y1_posterior),
  observedResponse = y1_sim
)

plot(y1_resids)

rdat <- tibble(
  y1 = y1_sim,
  residual = y1_resids$scaledResiduals
)

rdat  |>
  ggplot(
    aes(
      x = residual
    )
  ) +
  geom_histogram(
    bins = 10
  )


################################################################################
################################################################################


# h2000 <- values(housing_2000)
# h2015 <- values(housing_2015)



h2000 <- values(
  housing_2000 #|>
    #terra::aggregate(fact = 5, na.rm = TRUE, fun = "max")
)



h2015 <- values(
  housing_2015 #|>
    #terra::aggregate(fact = 5, na.rm = TRUE, fun = "max")
)

notnaidx <- which(!is.na(h2000))

hdat_all <- tibble(
  h0 = h2000[notnaidx],
  h1 = h2015[notnaidx]
) |>
  arrange(h0)

total_non_na <- length(notnaidx)

sampleidx <- seq(
  from = 1,
  to = total_non_na,
  length.out = 2000
) |>
  round()

hdat <- hdat_all[sampleidx,] |>
  mutate(
    t = 15,
    hdelta = h1 - h0
  )

#
# h0 <- as_data(hdat$h0)
# h1 <- as_data(hdat$h1)
#
#
# hdelta <- h1 - h0
#
#
# h1 <- logistic()
#
# t <- ilogit(h0)
#
# t + 15 <- ilogit(h1)


###
nsim <- 1000
tmid_sim <- 50
k_sim <- 0.1

t0_sim <- runif(nsim, 0, 100)

h0_sim <- 1/(1 + exp(-k_sim*(t0_sim - tmid_sim)))

plot(t0_sim, h0_sim)

t_delta_sim <- runif(nsim, 5, 20)
#t_delta_sim <- rep(15, nsim)

h1_sim <- 1 /(1 + exp(-k_sim*(t0_sim + t_delta_sim - tmid_sim)))

par(mfrow = c(1,2))
plot(h0_sim, h1_sim)
plot(hdat$h0, hdat$h1)
par(mfrow = c(1,1))

# ###
#
# h0 <- as_data(hdat$h0)
# h1 <- as_data(hdat$h1)
#
# t_delta <- 15
#
# k <- normal(0, 1)
# #t_mid <- normal(50, 5, truncation(0, Inf))
# #h0_mean <- 1/(1 + exp(-k*(t0 - tmid_sim)))
# #h1_mean <- 1/(1 + exp(-k*(t0 + t_delta - tmid_sim)))
#
#
# h1_mean <- h0 * exp(k * t_delta) / (1 - h0 + h0 * exp(k * t_delta))
#
# sigma <- exponential(1)
#
# distribution(h1) <- normal(h1_mean, sigma)
#
# m <- model(k, sigma)
#
# h1_prior <- calculate(h1)
#
# plot(h1, h1_prior$h1)
#
# draws <- mcmc(
#   m,
#   n_samples = 1000,
#   warmup = 500,
#   chains = 4
# )
#
# summary(draws)
#


#######



h_delta <- as_data(hdat$hdelta)

t_delta <- 15

k <- normal(0, 1)
t_mid <- 0
t0 <- normal(0, 5, dim = dim(h_delta))
# difference between h1 and h0 depends on t0
h0_mean <- 1/(1 + exp(-k*(t0 - t_mid)))
h1_mean <- 1/(1 + exp(-k*(t0 + t_delta - t_mid)))

h_delta_mean <-  h1_mean - h0_mean

sd <- lognormal(0, 3)

distribution(h_delta) <- normal(h_delta_mean, sd)

m <- model(k)


draws <- mcmc(
  m,
  warmup = 500,
  n_samples = 1000,
  chains = 4
)

summary(draws)

h0_pred <- calculate(
  h0_mean,
  values = draws,
  nsim = 100
)






