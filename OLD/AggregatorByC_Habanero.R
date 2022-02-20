# Intentions de votes pour la présidentielle française 2022
# Bayesian Model - Habenero Cluster
# Flavien Ganter

# Created on October 27, 2021
# Modified on February 20, 2022



#### PRELIMINARIES ####

# Clear working space
rm(list = ls())

# Set working directory
.libPaths("/rigel/sscc/users/fg2465/rpackages")

# Packages
library(rstan)
library(splines)
library(tidyverse)

# Commands
"%nin%" <- Negate("%in%")
options(mc.cores = parallel::detectCores())

# Array argument
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
args <- as.numeric(slurm_arrayid)

# Data
load("PollsData.RData")
data <- data %>% 
  filter(id_candidate == args[1]) %>% 
  group_by(id_poll) %>% mutate(id_poll = cur_group_id()) %>% ungroup() %>% 
  group_by(id_house) %>% mutate(id_house = cur_group_id()) %>% ungroup()




#### MODEL ####

# Define splines 
num_knots     <- 10
spline_degree <- 3
num_basis     <- num_knots + spline_degree - 1
B             <- t(bs(1:max(data$id_date), df = num_basis, degree = spline_degree, intercept = TRUE))

# Gather data to feed the model
data_spline_model <- list(N             = nrow(data),
                          tot_eff       = data$tot_eff,
                          vote_eff      = data$vote_eff,
                          id_cand       = args[1],
                          id_date       = round(data$id_date),
                          id_month      = data$id_month,
                          M             = length(unique(data$id_month)),
                          id_poll       = data$id_poll,
                          P             = length(unique(data$id_poll)),
                          id_house       = data$id_house,
                          F             = length(unique(data$id_house)),
                          X             = data[, c("unsure_1", "unsure_2", "rolling_yes")],
                          isn_z         = data$isnot_zemmour,
                          isn_t         = data$isnot_taubira,
                          num_knots     = num_knots,
                          knots         = unname(quantile(1:max(data$id_date), probs = seq(from = 0, to = 1, length.out = num_knots))),
                          spline_degree = spline_degree,
                          num_basis     = num_basis,
                          D             = ncol(B),
                          B             = as.matrix(B))

# Run model
aggregator_model <- stan(file = "ModelSplinesByC.stan",
                         data = data_spline_model,
                         iter = 4000, chains = 4,
                         seed = 94838,
                         include = TRUE,
                         pars = c("prob"),
                         save_warmup = FALSE,
                         control = list(adapt_delta = .99, max_treedepth = 13))

save(aggregator_model, file = paste0("/rigel/sscc/users/fg2465/model_aggregator_", args[1], ".RData"))
