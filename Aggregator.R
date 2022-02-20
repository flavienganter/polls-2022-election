# Intentions de votes pour la présidentielle française 2022
# Flavien Ganter

# Created on October 4, 2021
# Last modified on February 18, 2022




#### PRELIMINARIES ####


## Set environment

# Clear working space
rm(list = ls())

# Set working directory
setwd("~/Dropbox/PollsFrance2022/")
options(mc.cores = parallel::detectCores())
Sys.setlocale("LC_TIME", "fr_FR.UTF-8")

# Packages
library(tidyverse)
library(tidylog)
library(cmdstanr)
library(readxl)
library(HDInterval)
library(splines)
library(zoo)
extrafont::loadfonts()

# Commands & Functions
"%nin%" <- Negate("%in%")
logit <- function(x) log(x/(1-x))
round2 <- function(x) {
  diff <- c(abs(x-trunc(x)), abs(x-trunc(x)-.5), abs(x-trunc(x)-1))
  if (length(which(min(diff) == diff)) == 1) {
    c(trunc(x), trunc(x)+.5, trunc(x)+1)[which(min(diff) == diff)]
  } else {
    c(trunc(x), trunc(x)+.5, trunc(x)+1)[which(min(diff) == diff)[2]]
  }
}
  
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
  
  # Create hypothesis ID
  mutate(id_hyp = 1:n()) %>% 
  
  # Identify scenarios w/o EZ and those w/o CT
  mutate(isnot_zemmour = ifelse(is.na(c_zemmour), 1, 0),
         isnot_taubira = ifelse(is.na(c_taubira), 1, 0)) %>%
    
  # Remove scenarios w/o EZ after September
  filter((isnot_zemmour == 1 & month_end == 9) |
           isnot_zemmour == 0) %>%
    
  # Remove weird scenarios re: Taubira
  filter(n_wgt > 0) %>% 
    
  # Remove scenarios w/o VP
  filter(!is.na(c_pecresse)) %>% 
  
  # Wide to long
  gather(candidate, share, c_poutou:c_thouy) %>% 
  
  # Remove rows corresponding to untested candidates
  filter(!is.na(share)) %>% 
  
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
  ungroup() %>% 
    
  # Renumber polls (to account for polls that have been removed)
  group_by(id) %>% 
  mutate(id_poll = cur_group_id()) %>% 
  ungroup() %>% 
  
  # Switch from share for numbers, create a logged sample size variable,
  # and center and standardize the variable that indicates what respondents
  # are included in the estimates (all of them, only those who are certain to
  # vote, or a mix of both)
  mutate(vote_eff = round(n_t1 * n_wgt * share),
         unsure_1 = scale_factor(variable = unsure)[[1]],
         unsure_2 = scale_factor(variable = unsure)[[2]],
         rolling_yes = scale_factor(variable = poll_type)) %>% 
    
  # Remove candidates who are not consistently tested
    ## Estimates for these candidates are very noisy, and not
    ## necessarily relevant. Omitting then does not affect the
    ## estimation for other candidates
  filter(candidate %nin% c("c_asselineau", "c_lagarde", "c_lassalle", "c_poisson", "c_montebourg",
                           "c_philippot", "c_thouy")) %>% 
  
  # Create a candidate ID
  mutate(id_candidate = case_when(candidate == "c_arthaud" ~ 1,
                                  candidate == "c_daignant" ~ 2,
                                  candidate == "c_hidalgo" ~ 3,
                                  candidate == "c_jadot" ~ 4,
                                  candidate == "c_lepen" ~ 5,
                                  candidate == "c_macron" ~ 6,
                                  candidate == "c_melenchon" ~ 7,
                                  candidate == "c_pecresse" ~ 8,
                                  candidate == "c_poutou" ~ 9,
                                  candidate == "c_roussel" ~ 10,
                                  candidate == "c_zemmour" ~ 11,
                                  candidate == "c_taubira" ~ 12)) %>% 
  
  # Create a house ID
  group_by(house) %>% 
  mutate(id_house = cur_group_id()) %>% 
  ungroup() %>% 
    
  # Create date IDs
  mutate(id_date_start = as.numeric(as.Date(paste(year, month_start, day_start, sep = "-"))) - 18869,
         id_date_end = as.numeric(as.Date(paste(year, month_end, day_end, sep = "-"))) - 18869,
         id_date = id_date_start + (id_date_end - id_date_start) / 2,
         id_month = as.numeric(format(as.Date(as.numeric(round(id_date)), origin = as.Date("2021-08-31")), "%m")),
         id_month = case_when(id_month > 8 ~ id_month - 8,  TRUE ~ id_month + 4)) %>% 
    
  # Correct for rounding differences in Ns
  group_by(id_hyp) %>%
  mutate(tot_eff = sum(vote_eff)) %>% 
  ungroup()
  
  
# Export
save(data, file = "PollsData.RData")


## Model

if (0) { # Model run in a separate cluster

# Define splines 
num_knots     <- 10
spline_degree <- 3
num_basis     <- num_knots + spline_degree - 1
B             <- t(bs(1:max(data$id_date), df = num_basis, degree = spline_degree, intercept = TRUE))

# Gather data to feed the model
data_spline_model <- list(N             = nrow(data),
                          tot_eff       = data$tot_eff,
                          vote_eff      = data$vote_eff,
                          id_cand       = data$id_candidate,
                          C             = length(unique(data$id_candidate)),
                          id_date       = data$id_date,
                          id_month      = data$id_month,
                          M             = length(unique(data$id_month)),
                          id_poll       = data$id,
                          P             = length(unique(data$id)),
                          id_firm       = data$id_firm,
                          F             = length(unique(data$id_firm)),
                          X             = data[, c("unsure_1", "unsure_2", "rolling_yes")],
                          isn_z         = data$isnot_zemmour,
                          isn_t         = data$isnot_taubira,
                          num_knots     = num_knots,
                          knots         = unname(quantile(1:max(data$id_date),
                                                          probs = seq(from = 0, to = 1,
                                                                      length.out = num_knots))),
                          spline_degree = spline_degree,
                          num_basis     = num_basis,
                          D             = ncol(B),
                          B             = as.matrix(B))
  
# Compile model
model_code <- cmdstan_model("ModelSplines.stan")

# Estimate model
estimated_spline_model <- model_code$sample(data            = data_spline_model,
                                            seed            = 94836,
                                            chains          = 4,
                                            parallel_chains = 4,
                                            iter_warmup     = 2000,
                                            iter_sampling   = 2000,
                                            adapt_delta     = .99,
                                            max_treedepth   = 15,
                                            refresh         = 1000,
                                            save_warmup     = FALSE)

# Get posterior draws
spline_draws <- estimated_spline_model$draws(variables = "prob", format = "draws_df")
save(spline_draws, file = "LatestDraws.RData")

}



#### GET ESTIMATES ####

# Import
load("model_aggregator.RData")

# Prepare draws

  # Splines
  spline_draws <- data.frame(`prob[1,1]` = rstan::extract(aggregator_model, pars = "prob[1,1]"))
  colnames(spline_draws) <- "prob[1,1]"
  for (i in 1:aggregator_model@par_dims[["prob"]][1]) {
    for (j in 1:12) {
      if (!(i == 1 & j == 1)) {
        spline_draws <- cbind(spline_draws,
                              `prob[i,j]` = rstan::extract(aggregator_model, pars = paste0("prob[", i, ",", j, "]")))
      }
    }
  }
  
  # Instantaneous estimates
  inst_draws <- data.frame(`prob[1]` = rstan::extract(aggregator_model, pars = paste0("prob[", aggregator_model@par_dims[["prob"]][1], ",1]")))
  colnames(inst_draws) <- "prob[111,1]"
  for (j in 2:12) inst_draws <- cbind(inst_draws, `prob[j]` = rstan::extract(aggregator_model, pars = paste0("prob[", aggregator_model@par_dims[["prob"]][1], ",", j, "]")))
  


#### SPLINE PLOT ####


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
                                         candidate == 8 ~ "Valérie Pécresse",
                                         candidate == 9 ~ "Philippe Poutou",
                                         candidate == 10 ~ "Fabien Roussel",
                                         candidate == 11 ~ "Éric Zemmour",
                                         candidate == 12 ~ "Christiane Taubira"))) %>% 
  filter(candidate != "Christiane Taubira" | (candidate == "Christiane Taubira" & date >= "2021-12-15"))


## Create plot

# Define candidate colors
candidate_colors <- c("#f7b4b4", "#FFCC33", "#ff6600", "black", "#ff1300", "#b30d00",
                      "#002060", "#8fa02a", "#7030a0", "#c80589", "#0070c0", "#00b050")

# Generate plot
poll_plot <- plot_spline_estimates %>% 
  mutate(label = if_else(date == max(date), paste0(as.character(candidate), " (", unlist(lapply(median*100, round2)), "%)"), NA_character_),
         median_label = case_when(candidate == "Anne Hidalgo" ~ median - .00,
                                  candidate == "Fabien Roussel" ~ median + .00,
                                  candidate == "Nicolas Dupont-Aignan" ~ median - .00,
                                  candidate == "Philippe Poutou" ~ median + .0015,
                                  candidate == "Nathalie Arthaud" ~ median - .0015,
                                  candidate == "Jean-Luc Mélenchon" ~ median,
                                  candidate == "Yannick Jadot" ~ median + .00,
                                  candidate == "Christiane Taubira" ~ median - .00,
                                  candidate == "Marine Le Pen" ~ median + .00,
                                  candidate == "Valérie Pécresse" ~ median + .00,
                                  candidate == "Éric Zemmour" ~ median - .00,
                                  !is.na(candidate) ~ median)) %>% 
  group_by(candidate) %>% 
  mutate(lower50_l = zoo::rollmean(lower50, k = 2, align = "left", fill = NA),
         lower50_r = zoo::rollmean(lower50, k = 2, align = "right", fill = NA),
         upper50_l = zoo::rollmean(upper50, k = 2, align = "left", fill = NA),
         upper50_r = zoo::rollmean(upper50, k = 2, align = "right", fill = NA),
         lower95_l = zoo::rollmean(lower95, k = 2, align = "left", fill = NA),
         lower95_r = zoo::rollmean(lower95, k = 2, align = "right", fill = NA),
         upper95_l = zoo::rollmean(upper95, k = 2, align = "left", fill = NA),
         upper95_r = zoo::rollmean(upper95, k = 2, align = "right", fill = NA),
         lower50_s = zoo::rollmean(lower50, k = 3, fill = NA),
         upper50_s = zoo::rollmean(upper50, k = 3, fill = NA),
         lower95_s = zoo::rollmean(lower95, k = 3, fill = NA),
         upper95_s = zoo::rollmean(upper95, k = 3, fill = NA),
         lower50_s = ifelse(is.na(lower50_s), lower50_l, lower50_s),
         lower50_s = ifelse(is.na(lower50_s), lower50_r, lower50_s),
         upper50_s = ifelse(is.na(upper50_s), upper50_l, upper50_s),
         upper50_s = ifelse(is.na(upper50_s), upper50_r, upper50_s),
         lower95_s = ifelse(is.na(lower95_s), lower95_l, lower95_s),
         lower95_s = ifelse(is.na(lower95_s), lower95_r, lower95_s),
         upper95_s = ifelse(is.na(upper95_s), upper95_l, upper95_s),
         upper95_s = ifelse(is.na(upper95_s), upper95_r, upper95_s)) %>% 
  ggplot(aes(x = date, group = candidate, color = candidate)) +
  
  # Plot data
  geom_line(aes(y = median * 100)) +
  geom_ribbon(aes(ymin = lower50_s * 100, ymax = upper50_s * 100, fill = candidate), alpha = .1, size = 0) +
  geom_ribbon(aes(ymin = lower95_s * 100, ymax = upper95_s * 100, fill = candidate), alpha = .1, size = 0) +
  
  # Candidate labels
  geom_text(aes(x = date + 1, y = median_label * 100, label = label), na.rm = TRUE,
            hjust = 0, vjust = 0, nudge_y = -.1, family = "Open Sans Condensed", size = 3) +
  
  # Show latest poll's date
  annotate("segment", x = max(plot_spline_estimates$date), y = 0, xend = max(plot_spline_estimates$date), yend = 28,
           size = .4) +
  annotate(geom = "text", x = max(plot_spline_estimates$date), y = 28.5, family = "Open Sans Condensed",
           label = format(max(plot_spline_estimates$date), "%d %B %Y"), size = 3) +
  
  # Show 1st round
  annotate("segment", x = as.Date("2022-04-10"), y = 0, xend = as.Date("2022-04-10"), yend = 27,
           size = .4) +
  annotate(geom = "text", x = as.Date("2022-04-10"), y = 28, family = "Open Sans Condensed",
           label = "Premier tour \n10 avril 2022", size = 3) +
  
  # Define labs
  labs(x = "", y = "Intentions de votes (% votes exprimés)",
       title = "Intentions de vote au 1er tour de l'élection présidentielle de 2022",
       subtitle = "Depuis septembre 2021",
       caption = paste0("Estimations obtenues à partir des enquêtes d'opinion réalisées par BVA, Cluster17, Elabe, Harris Interactive, IFOP, IPSOS, Odoxa, et OpinionWay depuis septembre 2021 sur la base des rapports d'enquête publiés sur le site de la Commission des sondages, et agrégées à l'aide d'un modèle bayésien tenant compte \ndes principales caractéristiques des enquêtes. Le graphique présente les médianes et intervalles de crédibilité (95% / 50%)Pour chaque candidat, la ligne solide relie les médianes des distributions a posteriori à chaque date, et la zone colorée représente la partie la plus dense de la distribution a posteriori (95% / 50%) des \ndistributions a posteriori. Dernière mise à jour: ", format(Sys.time(), "%d %B %Y"), ".")) +
  
  # Specify plot theme
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "#2b2b2b", linetype = "dotted", size = 0.05),
        text = element_text(family = "Open Sans Condensed", size = 7.5),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(hjust = 0),
        axis.text.y = element_text(margin = margin(r = -2)),
        axis.title = element_blank(),
        axis.line.x = element_line(color = "#2b2b2b", size = 0.15),
        axis.ticks.x = element_line(color = "#2b2b2b", size = 0.15),
        axis.ticks.length = unit(.2, "cm"),
        plot.title = element_text(size = 20, family = "Open Sans Condensed ExtraBold", face = "bold"),
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
               limits = c(as.Date("2021-09-01"), as.Date("2022-04-15"))) +
  
  # Percent axis
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0), breaks = seq(0, 30, 5), lim = c(0, 29))


## Export plot
ggsave(poll_plot, filename = "PollsFrance2022_evolution.pdf",
       height = 6, width = 10, device = cairo_pdf)
ggsave(poll_plot, filename = "PollsFrance2022_evolution.png",
       height = 6, width = 10, device = "png", bg = "white")




#### INSTANTANEOUS PLOT ####


## Prepare table for plot

# Calculate statistics of interest
plot_inst_estimates <- apply(inst_draws, 2, function(x) c(hdi(x), hdi(x, .9), hdi(x, .8), hdi(x, .5), median(x))) %>%
  t() %>% as.data.frame()

# Rename names and columns
names(plot_inst_estimates) <- c("lower95", "upper95", "lower90", "upper90", "lower80", "upper80", "lower50", "upper50", "median")
plot_inst_estimates$coef <- row.names(plot_inst_estimates)

# Name estimates: identify the date and candidate associated with each estimate
plot_inst_estimates <- plot_inst_estimates %>%
  mutate(label = paste0(unlist(lapply(median*100, round2)), "%"),
         label = ifelse(label == "0.5%", "   0.5%", label),
         candidate = substr(coef, 10, 11),
         candidate = case_when(substr(candidate, 2, 2) == "]" ~ substr(candidate, 1, 1),
                               TRUE ~ candidate),
         candidate = as.factor(case_when(candidate == 1 ~ "Nathalie Arthaud",
                                         candidate == 2 ~ "Nicolas Dupont-Aignan",
                                         candidate == 3 ~ "Anne Hidalgo",
                                         candidate == 4 ~ "Yannick Jadot",
                                         candidate == 5 ~ "Marine Le Pen",
                                         candidate == 6 ~ "Emmanuel Macron",
                                         candidate == 7 ~ "Jean-Luc Mélenchon",
                                         candidate == 8 ~ "Valérie Pécresse",
                                         candidate == 9 ~ "Philippe Poutou",
                                         candidate == 10 ~ "Fabien Roussel",
                                         candidate == 11 ~ "Eric Zemmour",
                                         candidate == 12 ~ "Christiane Taubira")))
plot_inst_estimates$candidate <- factor(plot_inst_estimates$candidate,
                                             levels = c("Anne Hidalgo", "Christiane Taubira", "Emmanuel Macron", "Eric Zemmour", 
                                                        "Fabien Roussel", "Jean-Luc Mélenchon", "Marine Le Pen", "Nathalie Arthaud", 
                                                        "Nicolas Dupont-Aignan", "Philippe Poutou", "Valérie Pécresse", 
                                                        "Yannick Jadot")[c(3, 7, 4, 11, 6, 12, 5, 2, 1, 9, 10, 8)])


## Create plot

# Define candidate colors
candidate_colors <- c("#f7b4b4", "#FFCC33", "#ff6600", "black", "#ff1300", "#b30d00",
                      "#002060", "#8fa02a", "#7030a0", "#c80589", "#0070c0", "#00b050")[c(3, 7, 4, 11, 6, 12, 5, 2, 1, 9, 10, 8)]

# Generate plot
inst_plot <- plot_inst_estimates %>% 
  ggplot(aes(x = candidate, color = candidate)) +
  
  # Plot data
  geom_linerange(aes(ymin = lower50 * 100, ymax = upper50 * 100), alpha = .2, size = 10) +
  geom_linerange(aes(ymin = lower80 * 100, ymax = upper80 * 100), alpha = .2, size = 10) +
  geom_linerange(aes(ymin = lower90 * 100, ymax = upper90 * 100), alpha = .2, size = 10) +
  geom_linerange(aes(ymin = lower95 * 100, ymax = upper95 * 100), alpha = .2, size = 10) +
  geom_linerange(aes(ymin = median * 100 - .05, ymax = median * 100 + .05), size = 12, color = "white") +
  geom_linerange(aes(ymin = median * 100 - .02, ymax = median * 100 + .02), size = 11) +
  geom_text(aes(label = label, y = median*100), family = "Open Sans Condensed SemiBold", vjust = -1.65) +
  
  # Define labs
  labs(x = "", y = "Intentions de votes (% votes exprimés)",
       title = "Intentions de vote au 1er tour de l'élection présidentielle de 2022",
       subtitle = paste("Au", format(Sys.time(), "%d %B %Y")),
       caption = paste0("Estimations obtenues à partir des enquêtes d'opinion réalisées par BVA, Cluster17, Elabe, Harris Interactive, IFOP, IPSOS, Odoxa, et OpinionWay depuis septembre 2021 sur la base des rapports d'enquête publiés sur le site de la Commission des sondages, et agrégées à l'aide d'un \nmodèle bayésien tenant compte des principales caractéristiques des enquêtes. Le graphique présente les médianes et intervalles de crédibilité (95% / 90% / 80% / 50%) des distributions a posteriori.")) +
  
  # Specify plot theme
  coord_flip() +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "#2b2b2b", linetype = "dotted", size = 0.05),
        panel.grid.major.y = element_line(color = "gray99", size = 10),
        text = element_text(family = "Open Sans Condensed", size = 7.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12, family = "Open Sans Condensed", face = "bold"),
        axis.title = element_blank(),
        axis.line.x = element_line(color = "#2b2b2b", size = 0.15),
        axis.ticks.x = element_line(color = "#2b2b2b", size = 0.15),
        axis.ticks.length = unit(.2, "cm"),
        plot.title = element_text(size = 20, family = "Open Sans Condensed ExtraBold", face = "bold"),
        plot.subtitle = element_text(size = 15),
        plot.title.position = "plot",
        legend.position = "none",
        plot.caption = element_text(color = "gray30", hjust = 0, margin = margin(t = 15)),
        plot.margin = unit(rep(0.5, 4), "cm")) +
  
  # Candidate colors
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
  scale_colour_manual(values = candidate_colors) +
  scale_fill_manual(values = candidate_colors) +
  
  # Percent axis
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0), breaks = seq(0, 30, 5), lim = c(0, 27))


## Export plot
ggsave(inst_plot, filename = "PollsFrance2022_latest.pdf",
       height = 7.5, width = 10, device = cairo_pdf)
ggsave(inst_plot, filename = "PollsFrance2022_latest.png",
       height = 7.5, width = 10, device = "png", bg = "white")


