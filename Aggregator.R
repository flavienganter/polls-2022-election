# Intentions de votes pour la présidentielle française 2022
# Flavien Ganter

# Created on October 4, 2021
# Last modified on November 10, 2021




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
extrafont::loadfonts()

# Commands & Functions
"%nin%" <- Negate("%in%")
logit <- function(x) log(x/(1-x))
  
  # Function to center and standardize continuous and categorical variables
  # (helpful for Bayesian prior, see Gelman et al. 2008)
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
  

  
#### POLL AGGREGATION MODEL ####


## Prepare data
data <- read_excel("PollsData.xlsx") %>% 
  
  # Create hypothesis ID and identify the Les Républicains candidate in each hypothesis
  mutate(id_hyp = 1:n(),
         c_repub = ifelse(!is.na(c_bertrand), c_bertrand, ifelse(!is.na(c_pecresse), c_pecresse, c_barnier))) %>% 
  
  # Identify scenarios w/o EZ
  mutate(isnot_zemmour = ifelse(is.na(c_zemmour), 1, 0)) %>%
    
  # Remove scenarios w/o EZ after September
  filter((isnot_zemmour == 1 & month == 9) |
           isnot_zemmour == 0) %>% 
  
  # Wide to long
  gather(candidate, share, c_poutou:c_lagarde, c_repub) %>% 
  
  # Remove rows corresponding to untested candidates
  filter(candidate %nin% c("c_bertrand", "c_pecresse", "c_barnier") &
           !is.na(share)) %>% 
  
  # Recode truncated values and switch to share
    ## The uncertainty due to truncated values could be modeled directly,
    ## but I simplify things here, and adopt a deterministic imputation
    ## approach (namely, I assign a value that is equal to half of the
    ## truncation value).
  mutate(share = case_when(share == "T_0.5" ~ 0.0025,
                           share == "T_1.5" ~ 0.0075,
                           share == "T_1" ~ 0.005,
                           TRUE ~ as.numeric(share) / 100)) %>% 
  
  # Standardize distributions to account for the fact that
  # the number of candidates tested varies between and within polls
  group_by(id_hyp) %>% 
  mutate(share = share / sum(share)) %>%
  
  # Remove candidates who are not consistently tested
    ## Estimates for these candidates are very noisy, and not
    ## necessarily relevant. Omitting then does not affect the
    ## estimation for other candidates
  filter(candidate %nin% c("c_asselineau", "c_lagarde", "c_lassalle", "c_poisson", "c_philippot")) %>% 
  ungroup() %>% 
  
  # Switch from share for numbers, create a logged sample size variable,
  # and center and standardize the variable that indicates what respondents
  # are included in the estimates (all of them, only those who are certain to
  # vote, or a mix of both)
  mutate(eff = round(n * share),
         log_n = log(n) - mean(log(n)),
         unsure_1 = scale_factor(variable = unsure)[[1]],
         unsure_2 = scale_factor(variable = unsure)[[2]]) %>% 
  
  # Create a candidate ID
  group_by(candidate) %>%
  mutate(id_candidate = cur_group_id()) %>% 
  
  # Create a firm ID and recode poll dates
    ## Poll dates are the median dates of data collection, or the first day after
    ## the median, if the median is not properly defined.
  group_by(firm) %>% 
  mutate(id_firm = cur_group_id(),
         id_date = as.numeric(as.Date(paste(year, month, day, sep = "-"))) - 18870) %>% 
  group_by(month) %>% 
  mutate(id_month = cur_group_id())
save(data, file = "PollsData.RData")


## Model

if (0) { # Model run in a separate cluster

# Define splines 
num_knots <- 4
spline_degree <- 3
num_basis <- num_knots + spline_degree - 1
B <- t(bs(1:max(data$id_date), df = num_basis, degree = spline_degree, intercept = TRUE))

# Gather data to feed the model
  ## prior_mu is the prior for each candidate score based on the latest poll in August
data_spline_model <- list(N = nrow(data),
                          tot_eff = data$n,
                          vote_eff = data$eff,
                          id_cand = data$id_candidate,
                          C = length(unique(data$id_candidate)),
                          id_date = data$id_date,
                          id_month = data$id_month,
                          M = length(unique(data$id_month)),
                          id_poll = data$id,
                          P = length(unique(data$id)),
                          id_firm = data$id_firm,
                          F = length(unique(data$id_firm)),
                          X = data[, c("log_n", "unsure_1", "unsure_2")],
                          isn_z = data$isnot_zemmour,
                          num_knots = num_knots,
                          knots = unname(quantile(1:max(data$id_date), probs = seq(from = 0, to = 1, length.out = num_knots))),
                          spline_degree = spline_degree,
                          num_basis = num_basis,
                          D = ncol(B),
                          B = as.matrix(B))
  
# Compile model
model_code <- cmdstan_model("ModelSplines.stan")

# Estimate model
estimated_spline_model <- model_code$sample(data = data_spline_model,
                                            seed = 94836,
                                            chains = 4,
                                            parallel_chains = 4,
                                            iter_warmup = 2000,
                                            iter_sampling = 2000,
                                            adapt_delta = .99,
                                            max_treedepth = 15,
                                            refresh = 1000,
                                            save_warmup = FALSE)

# Get posterior draws
spline_draws <- estimated_spline_model$draws(variables = "prob", format = "draws_df")
save(spline_draws, file = "LatestDraws.RData")

}



#### GET ESTIMATES ####

# Import
load("model_aggregator.RData")

# Prepare draws
spline_draws <- data.frame(`prob[1,1]` = rstan::extract(aggregator_model, pars = "prob[1,1]"))
colnames(spline_draws) <- "prob[1,1]"
for (i in 1:77) {
  for (j in 1:12) {
    if (!(i == 1 & j == 1)) {
      spline_draws <- cbind(spline_draws,
                            `prob[i,j]` = rstan::extract(aggregator_model, pars = paste0("prob[", i, ",", j, "]")))
    }
  }
}




#### GENERATE PLOT ####


## Prepare table for plot

# Calculate statistics of interest
plot_spline_estimates <- apply(spline_draws, 2, function(x) c(hdi(x), hdi(x, .5), median(x))) %>%
  t() %>% as.data.frame()

# Rename names and columns
names(plot_spline_estimates) <- c("lower95", "upper95", "lower50", "upper50", "median")
plot_spline_estimates$coef <- row.names(plot_spline_estimates)

# Name estimates: identify the date and candidate associated with each estimate
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
                                         candidate == 10 ~ "Candidat Les Républicains",
                                         candidate == 11 ~ "Fabien Roussel",
                                         candidate == 12 ~ "Eric Zemmour")))


## Create plot

# Define candidate colors
candidate_colors <- c("#f7b4b4", "#af8080", "#0070c0", "#ff6600", "black", "#ff1300",
                      "#b30d00", "#002060", "#c80589", "#7030a0", "#8fa02a", "#00b050")
  
# Generate plot
poll_plot <- plot_spline_estimates %>% 
  mutate(label = if_else(date == max(date), as.character(candidate), NA_character_),
         median_label = case_when(label == "Arnaud Montebourg" ~ median,
                                  label == "Fabien Roussel" ~ median - .002,
                                  label == "Nicolas Dupont-Aignan" ~ median + .003,
                                  label == "Philippe Poutou" ~ median + .001,
                                  label == "Nathalie Arthaud" ~ median - .002,
                                  label == "Marine Le Pen" ~ median,
                                  label == "Eric Zemmour" ~ median,
                                  !is.na(label) ~ median)) %>% 
  ggplot(aes(x = date, group = candidate, color = candidate)) +
  
  # Plot data
  geom_line(aes(y = median * 100)) +
  geom_ribbon(aes(ymin = lower50 * 100, ymax = upper50 * 100, fill = candidate), alpha = .1, size = 0) +
  geom_ribbon(aes(ymin = lower95 * 100, ymax = upper95 * 100, fill = candidate), alpha = .1, size = 0) +
  
  # Candidate labels
  geom_text(aes(x = date + 1, y = median_label * 100, label = label), na.rm = TRUE,
            hjust = 0, vjust = 0, nudge_y = -.1, family = "Open Sans Condensed") +
  
  # Show 1st round
  geom_vline(xintercept = as.Date("2022-04-10"), color = "gray", size = 2) +
  annotate(geom = "text", x = as.Date("2022-03-04"), y = 24, family = "Open Sans Condensed",
           label = "10 avril 2022 – Premier tour de l'élection présidentielle") +
  annotate("segment", x = as.Date("2022-03-29"), y = 23.9, xend = as.Date("2022-04-09"), yend = 23,
           size = .4, arrow = arrow(angle = 30, length = unit(2.5, "mm"))) +
  
  # Define labs
  labs(x = "", y = "Intentions de votes (% votes exprimés)",
       title = "Intentions de vote au 1er tour de l'élection présidentielle de 2022",
       caption = paste0("Estimations obtenues à partir des enquêtes d'opinion réalisées par IPSOS, IFOP, Harris Interactive, Elabe, Odoxa et OpinionWay depuis septembre 2021 (sur la base des rapports d'enquête publics), et agrégées à l'aide d'une régression locale bayésienne tenant compte des principales caractéristiques des enquêtes. Le \nmodèle prend en compte le fait que la candidature de Zemmour n'a pas été testée dans toutes les enquêtes début septembre. Les intentions de vote en faveur du candidat des Républicains agrègent celles en faveur de Xavier Bertrand, Valérie Pécresse et Michel Barnier, en donnant un poids identique aux trois \ncandidats. Les lignes relient les médianes des distributions a posteriori, et les zones colorées représentent l'étendue des 95% les plus dense des distributions a posteriori (50%, pour la partie la plus sombre). D'après le modèle estimé, à un moment donné, les intentions de votes ont donc 95% de chances de se trouver \ndans l'intervalle le plus clair, et 50% de chances se trouver dans le plus sombre. Dernière mise à jour: ", Sys.Date(), ".")) +
  
  # Specify plot theme
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "#2b2b2b", linetype = "dotted", size = 0.15),
        text = element_text(family = "Open Sans Condensed"),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(hjust = 0),
        axis.text.y = element_text(margin = margin(r = -2)),
        axis.title = element_blank(),
        axis.line.x = element_line(color = "#2b2b2b", size = 0.15),
        axis.ticks.x = element_line(color = "#2b2b2b", size = 0.15),
        axis.ticks.length = unit(.2, "cm"),
        plot.title = element_text(size = 25, family = "Open Sans Condensed ExtraBold", face = "bold"),
        plot.title.position = "plot",
        legend.position = "none",
        plot.caption = element_text(color = "gray30", hjust = 0, margin = margin(t = 15)),
        plot.margin = unit(rep(0.5, 4), "cm")) +
  
  # Candidate colors
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
  scale_colour_manual(values = candidate_colors) +
  scale_fill_manual(values = candidate_colors) +
  
  # Date axis
  scale_x_date(expand = c(.005,1), date_breaks = "1 month",
               date_labels = c("Avril", "Septembre", "Octobre", "Novembre", "Décembre", "Janvier", "Février", "Mars"),
               limits = c(as.Date("2021-09-01"), as.Date("2022-04-10"))) +
  
  # Percent axis
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(.02, 0), breaks = seq(0, 30, 5), lim = c(0, 27.5))


## Export plot
ggsave(poll_plot, filename = "PollsFrance2022_latest.pdf",
       height = 10, width = 14, device = cairo_pdf)

