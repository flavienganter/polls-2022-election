# Intentions de votes pour la présidentielle française 2022
# Flavien Ganter

# Created on October 4, 2021
# Last modified on October 12, 2021




#### PRELIMINARIES ####


## Set environment

# Clear working space
rm(list = ls())

# Set working directory
setwd("~/Dropbox/PollsFrance2022/")
options(mc.cores = parallel::detectCores())

# Packages
library(tidyverse)
library(tidylog)
library(cmdstanr)
library(readxl)
library(HDInterval)
library(splines)

# Commands & Functions
"%nin%" <- Negate("%in%")
logit <- function(x) log(x/(1-x))


## Data

# Scaling functions
scale_cont <- function(variable) (variable - mean(variable)) / (2 * weighted.sd(variable))
scale_factor <- function(variable) {
  scale_dummy <- function(dummy) ifelse(dummy == 1, 1 - mean(dummy == 1), - mean(dummy == 1))
  L <- length(unique(variable))

  if (L == 2) {
    out_variable <- scale_dummy(variable)
  } else {
    out_variable <- vector(mode = "list", length = L-1)
    for (l in 1:(L-1)) {
      dummy <- ifelse(variable == unique(variable)[l], 1, 0)
      out_variable[[l]] <- scale_dummy(dummy)
    }
  }
  return(out_variable)
}

# Prepare data
data <- read_excel("PollsData.xlsx") %>% 
  mutate(c_repub = ifelse(!is.na(c_bertrand), c_bertrand, ifelse(!is.na(c_pecresse), c_pecresse, c_barnier))) %>% 
  filter(!is.na(c_zemmour) | id %in% c(10, 5, 13)) %>%
  mutate(isnot_zemmour = ifelse(is.na(c_zemmour), 1, 0)) %>% 
  gather(candidate, share, c_poutou:c_repub) %>% 
  filter(candidate %nin% c("c_bertrand", "c_pecresse", "c_barnier") &
           !is.na(share)) %>% 
  mutate(share = case_when(share == "T_0.5" ~ 0.0025,
                           share == "T_1.5" ~ 0.0075,
                           share == "T_1" ~ 0.005,
                           TRUE ~ as.numeric(share) / 100)) %>% 
  group_by(id, candidate) %>% 
  mutate(share = mean(share, na.rm = TRUE),
         n = mean(n, na.rm = TRUE),
         id_temp = 1:n()) %>%
  filter(id_temp == 1) %>%
  group_by(id) %>% 
  mutate(share = share / sum(share)) %>%
  filter(candidate %nin% c("c_asselineau", "c_lagarde", "c_lassalle", "c_poisson", "c_philippot")) %>% 
  ungroup() %>% 
  mutate(eff = round(n * share),
         log_n = log(n) - mean(log(n)),
         n = round(n),
         unsure_1 = scale_factor(variable = unsure)[[1]],
         unsure_2 = scale_factor(variable = unsure)[[2]]) %>% 
  group_by(candidate) %>%
  mutate(id_candidate = cur_group_id()) %>% 
  group_by(firm) %>% 
  mutate(id_firm = cur_group_id(),
         id_date = as.numeric(as.Date(paste(year, month, day, sep = "-"))) - 18870)

#### CORRECT FOR ZEMMOUR OMMISSION ####
data_zem <- read_excel("PollsData.xlsx") %>% 
  mutate(c_repub = ifelse(!is.na(c_bertrand), c_bertrand, ifelse(!is.na(c_pecresse), c_pecresse, c_barnier))) %>% 
  filter(id %in% c(1, 2, 5, 6, 9, 13, 14)) %>%
  group_by(id) %>% 
  mutate(id_new = cur_group_id()) %>% 
  ungroup() %>% 
  mutate(isn_z = ifelse(is.na(c_zemmour), 1, 0)) %>% 
  gather(candidate, share, c_poutou:c_repub) %>% 
  filter(candidate %nin% c("c_bertrand", "c_pecresse", "c_barnier") &
           !is.na(share)) %>% 
  mutate(share = case_when(share == "T_0.5" ~ 0.0025,
                           share == "T_1.5" ~ 0.0075,
                           share == "T_1" ~ 0.005,
                           TRUE ~ as.numeric(share) / 100)) %>% 
  group_by(id, candidate, isn_z) %>% 
  mutate(share = mean(share, na.rm = TRUE),
         n = mean(n, na.rm = TRUE),
         id_temp = 1:n()) %>%
  filter(id_temp == 1) %>%
  group_by(id, isn_z) %>% 
  mutate(share = share / sum(share)) %>%
  filter(candidate %nin% c("c_asselineau", "c_lagarde", "c_lassalle", "c_poisson", "c_philippot", "c_zemmour")) %>% 
  ungroup() %>% 
  mutate(eff = round(n * share),
         log_n = log(n) - mean(log(n)),
         n = round(n)) %>% 
  group_by(candidate) %>%
  mutate(id_candidate = cur_group_id())
model_zemmour_code <- cmdstan_model("ModelZemmour.stan")
data_zemmour_model <- list(N = nrow(data_zem),
                           tot_eff = data_zem$n,
                           vote_eff = data_zem$eff,
                           id_cand = data_zem$id_candidate,
                           C = length(unique(data_zem$id_candidate)),
                           id_poll = data_zem$id_new,
                           P = length(unique(data_zem$id_new)),
                           isn_z = data_zem$isn_z,
                           prior_mu = logit(c(.01, .04, .09, .11, .2175, .24, .0875, .05, .01, .145, .02)))
estimated_zemmour_model <- model_zemmour_code$sample(data = data_zemmour_model,
                                                     seed = 94836,
                                                     chains = 4,
                                                     parallel_chains = 4,
                                                     iter_warmup = 2000,
                                                     iter_sampling = 2000,
                                                     refresh = 1000,
                                                     save_warmup = FALSE)
zemmour_draws <- estimated_zemmour_model$draws(variables = "kappa", format = "draws_df")
zemmour_correct <- apply(zemmour_draws[,1:11], 2, function(x) c(mean(x), sd(x))) %>%
  t() %>% as.data.frame() %>% 
  rename(zemcorr_mean = V1,
         zemcorr_sd = V2) %>% 
  mutate(id_c = 1:11)
  
  


#### MODEL ####
model_code <- cmdstan_model("ModelSplines.stan")
num_knots <- 3
spline_degree <- 3
num_basis <- num_knots + spline_degree - 1
B <- t(bs(1:max(data$id_date), df = num_basis, degree = spline_degree, intercept = TRUE))
data_spline_model <- list(N = nrow(data),
                          tot_eff = data$n,
                          vote_eff = data$eff,
                          id_cand = data$id_candidate,
                          C = length(unique(data$id_candidate)),
                          id_date = data$id_date,
                          id_firm = data$id_firm,
                          F = length(unique(data$id_firm)),
                          X = data[, c("log_n", "unsure_1", "unsure_2")],
                          isn_z = data$isnot_zemmour,
                          zemcorr_mean = zemmour_correct$zemcorr_mean,
                          zemcorr_sd = zemmour_correct$zemcorr_sd,
                          prior_mu = logit(c(.01, .04, .09, .11, .2175, .24, .0875, .05, .01, .145, .02, .07)),
                          num_knots = num_knots,
                          knots = unname(quantile(1:max(data$id_date), probs = seq(from = 0, to = 1, length.out = num_knots))),
                          spline_degree = spline_degree,
                          num_basis = num_basis,
                          D = ncol(B),
                          B = as.matrix(B))
estimated_spline_model <- model_code$sample(data = data_spline_model,
                                            seed = 94836,
                                            chains = 4,
                                            parallel_chains = 4,
                                            iter_warmup = 2000,
                                            iter_sampling = 2000,
                                            refresh = 1000,
                                            save_warmup = FALSE)
spline_draws <- estimated_spline_model$draws(variables = "prob", format = "draws_df")



#### PLOT ####

## Prepare table
plot_spline_estimates <- apply(spline_draws[,1:432], 2, function(x) c(hdi(x), hdi(x, .5), median(x))) %>%
  t() %>% as.data.frame()
names(plot_spline_estimates) <- c("lower95", "upper95", "lower50", "upper50", "median")
plot_spline_estimates$coef <- row.names(plot_spline_estimates)
plot_spline_estimates <- plot_spline_estimates %>%
  mutate(date = substr(coef, 6, 8),
         date = case_when(substr(date, 3, 3) == "," ~ substr(date, 1, 2),
                          substr(date, 2, 2) == "," ~ substr(date, 1, 1),
                          TRUE ~ substr(date, 1, 3)),
         date = as.Date(as.numeric(date), origin = as.Date("2021-08-31")),
         candidate = substr(coef, 7, 11),
         candidate = case_when(substr(candidate, 1, 1) == "," ~ substr(candidate, 2, 3),
                               substr(candidate, 2, 2) == "," ~ substr(candidate, 3, 4),
                               TRUE ~ substr(candidate, 4, 5)),
         candidate = ifelse(substr(candidate, 2, 2) == "]", substr(candidate, 1, 1), candidate),
         candidate = as.factor(case_when(candidate == 1 ~ "Nathalie Arthaud",
                                         candidate == 2 ~ "Nicolas Dupont-Aignan",
                                         candidate == 3 ~ "Anne Hidalgo",
                                         candidate == 4 ~ "Yannick Jadot",
                                         candidate == 5 ~ "Marine Le Pen",
                                         candidate == 6 ~ "Emmanuel Macron",
                                         candidate == 7 ~ "Jean-Luc Mélenchon",
                                         candidate == 8 ~ "Arnaud Montebourg",
                                         candidate == 9 ~ "Philippe Poutou",
                                         candidate == 10 ~ "Candidat Les Républicain",
                                         candidate == 11 ~ "Fabien Roussel",
                                         candidate == 12 ~ "Eric Zemmour")))

## Plot
base_breaks_y <- function(x, add = 0){
  b <- pretty(x)
  d <- data.frame(x = as.Date(-Inf), xend = as.Date(-Inf), y = min(b), yend = max(b) + add)
  list(geom_segment(data = d, aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE, size = .7),
       scale_y_continuous(breaks = b))
}

candidate_colors <- c("#f7b4b4", "#af8080", "#0070c0", "#ff6600", "black", "#ff1300",
                      "#b30d00", "#002060", "#c80589", "#7030a0", "#8fa02a", "#00b050")
polls <- plot_spline_estimates %>% 
  ggplot(aes(x = date, group = candidate, color = candidate)) +
  geom_line(aes(y = median * 100)) +
  geom_ribbon(aes(ymin = lower50 * 100, ymax = upper50 * 100, fill = candidate), alpha = .1, size = 0) +
  geom_ribbon(aes(ymin = lower95 * 100, ymax = upper95 * 100, fill = candidate), alpha = .1, size = 0) +
  geom_vline(xintercept = as.Date("2022-04-10"), color = "gray", size = 2) +
  annotate(geom = "text", x = as.Date("2022-03-01"), y = 25, family = "HelveticaNeueCond",
           label = "10 avril 2022 – Premier tour de l'élection présidentielle") +
  annotate("segment", x = as.Date("2022-03-29"), y = 24.9, xend = as.Date("2022-04-09"), yend = 24,
           size = .4, arrow = arrow(angle = 30, length = unit(2.5, "mm"))) +
  labs(x = "", y = "Intentions de votes (% votes exprimés)",
       title = "Estimation des intentions de vote au 1er tour de l'élection présidentielle",
       subtitle = paste0("Dernière mise à jour: ", Sys.Date()),
       caption = "Estimations obtenues à partir des enquêtes d'opinion réalisées par IPSOS, IFOP, Harris Interactive, Elabe, Odoxa et OpinionWay depuis septembre 2021 (sur la base des rapports d'enquête publics), et agrégées à l'aide d'une régression locale bayésienne tenant compte des principales caractéristiques des \nenquêtes. Le modèle prend en compte le fait que la candidature de Zemmour n'a pas été testée dans toutes les enquêtes début septembre. Les intentions de vote en faveur du candidat des Républicains aggrège les intentions de vote en faveur de Xavier Bertrand, Valérie Pécresse et Michel Barnier, \nen donnant un poids identique aux trois candidats. Les lignes relient les médianes des distributions a posteriori, et les zones colorées représentent l'étendue des 95% les plus dense de les distributions a posteriori (50%, pour la partie la plus sombre). D'après le modèle estimé, à un moment donné, \nles intentions de votes ont donc 95% de chances de se trouver dans l'intervalle le plus clair, et 50% de chances se trouver dans le plus sombre.") +
  theme_bw() +
  geom_hline(yintercept = 0) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(.3, "cm"),
        axis.text.x = element_text(vjust = 15),
        axis.text = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        text = element_text(family = "HelveticaNeueCond"),
        plot.title = element_text(size = 20),
        legend.box = "vertical",
        legend.title = element_blank(),
        legend.position = c(.79,1)) +
  base_breaks_y(c(0, 30), add = .01) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
  scale_colour_manual(values = candidate_colors) +
  scale_fill_manual(values = candidate_colors) +
  scale_x_date(expand = c(.01,0), date_breaks = "1 month", date_labels = "%b",
               limits = c(as.Date("2021-09-01"), as.Date("2022-04-10")))


ggsave(polls, filename = "PollsFrance2022_latest.pdf",
       height = 10, width = 14, device = cairo_pdf)

