#################################################################
# This script is needed to prepare COVID data and              #
# run model v1.1 or above                                       #
# RULES: playing as if it was February 2020                     #
#                                                               # 
#  Author: Leonardo Gada                                        #   
#                                                               # 
################################################################


############################ Load Libraries ####################################

library(readxl)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
# library(shinystan) at present not possible to use with CmdStanR due to output file issue
library(priorsense)
library(cowplot)
library(coda)

# Stan set up
options(mc.cores = parallel::detectCores()) # if running on multicore machine for parallel computing

############################# Load dataset #####################################

confirmed_cases <- read_excel("/Datasets/time_series_covid19_confirmed_global.xlsx")

confirmed_deaths <- read_excel("/Datasets/time_series_covid19_deaths_global.xlsx")

cluster_lookup <- read_excel("/Datasets/CART_Clusters_1.xlsx", sheet = "Sheet2")

pop_data_2019 <- read_excel("/Datasets/All Countries Population 2019.xlsx", sheet = "Data")

pop_data_2019 <- pop_data_2019 %>% select(`Country Name`, `2019 [YR2019]`) # keep only what we need
pop_data_2019$`2019 [YR2019]` <- as.integer(pop_data_2019$`2019 [YR2019]`) # make into integer or it gives issues later

# We want to take only up to end of February 2020, before measures in place


########################### Clean dataset ######################################

# The data is cumulative, so we take the column at end of February (29/02/2020)

EndOfFeb_Cases <- confirmed_cases %>% 
  select(`Country/Region`,`2/29/20`)%>%
  group_by(`Country/Region`) %>%
  summarise(TotCases = sum(`2/29/20`))

EndOfFeb_Deaths <- confirmed_deaths %>% 
  select(`Country/Region`,`2/29/20`)%>%
  group_by(`Country/Region`) %>%
  summarise(TotDeaths = sum(`2/29/20`))

# Removing non-country entities 

EndOfFeb_Cases <- EndOfFeb_Cases %>%
  filter(!(`Country/Region` %in% c("Diamond Princess", "MS Zaandam", 
                                   "Winter Olympics 2022", "Summer Olympics 2020",
                                   "Holy See", "Antarctica")))

EndOfFeb_Deaths <- EndOfFeb_Deaths %>%
  filter(!(`Country/Region` %in% c("Diamond Princess", "MS Zaandam", 
                                   "Winter Olympics 2022", "Summer Olympics 2020",
                                   "Holy See", "Antarctica")))

# Join

EoFeb_Cases_Deaths <- left_join(EndOfFeb_Cases, EndOfFeb_Deaths, by = "Country/Region")


############################ Cluster Grouping ##################################

# Sort out naming inconsistencies

EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Bahamas"] <- "Bahamas, The"
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Brunei"] <- "Brunei Darussalam"
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Burma"] <- "Myanmar"
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Congo (Brazzaville)"] <- "Congo, Rep."
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Congo (Kinshasa)"] <- "Congo, Dem. Rep."
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Egypt"] <- "Egypt, Arab Rep."
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Gambia"] <- "Gambia, The"
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Iran"] <- "Iran, Islamic Rep."
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Korea, North"] <- "Korea, Dem. People's Rep."
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Korea, South"] <- "Korea, Rep."
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Kyrgyzstan"] <- "Kyrgyz Republic"
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Laos"] <- "Lao PDR"
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Micronesia"] <- "Micronesia, Fed. Sts."
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Russia"] <- "Russian Federation"
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Saint Kitts and Nevis"] <- "St. Kitts and Nevis"
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Saint Lucia"] <- "St. Lucia"
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Saint Vincent and the Grenadines"] <- "St. Vincent and the Grenadines"
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Slovakia"] <- "Slovak Republic"
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Syria"] <- "Syrian Arab Republic"
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Taiwan*"] <- "China"
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Turkey"] <- "Turkiye"
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "US"] <- "United States"
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Venezuela"] <- "Venezuela, RB"
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Vietnam"] <- "Viet Nam"
EoFeb_Cases_Deaths$`Country/Region`[EoFeb_Cases_Deaths$`Country/Region` == "Yemen"] <- "Yemen, Rep."

# Join to add population data from 2020

EoFeb_Cases_Deaths <- left_join(EoFeb_Cases_Deaths, pop_data_2019, by = c("Country/Region" = "Country Name"))

# Join to add cluster memberships

EoFeb_Cases_Deaths_Clust <- left_join(EoFeb_Cases_Deaths, cluster_lookup, by = c("Country/Region" = "Country"))
EoFeb_Cases_Deaths_Clust <- EoFeb_Cases_Deaths_Clust %>% filter(TotCases > 0) # if we decide to use the tot pop as all countries you can remove this line

# Group by cluster 

covid_stan_data <- EoFeb_Cases_Deaths_Clust %>%
  filter(is.na(Cluster) == FALSE) %>% 
  select(Cluster, TotCases, TotDeaths,`2019 [YR2019]`) %>%
  group_by(Cluster) %>%
  summarise(TotClustCases = sum(TotCases), TotClustDeaths = sum(TotDeaths), TotClustPop = sum(`2019 [YR2019]`))

# Add fake year variable to make it fit the stan model structure

covid_stan_data$year <- 1

############################ Prepare Stan Input ################################

covid_list <- list(
  N = nrow(covid_stan_data),                                # Number of observations
  C = length(unique(covid_stan_data$Cluster)),              # Number of countries
  #  Ye = length(unique(covid_stan_data$year)),                # Number of years
  cases = covid_stan_data$TotClustCases,                    # Observed cases by country by year
  deaths = covid_stan_data$TotClustDeaths,                  # Observed deaths by country by year
  country = as.integer(as.factor(covid_stan_data$Cluster)), #  grouping variable, encoding which country the observation comes from
  population_size = covid_stan_data$TotClustPop#,            # Number of people by cluster (repeated avg per year as well, can change to be year specific pop)
  # year = covid_stan_data$year                               # grouping variable, encoding which year the observation comes from
)


############################### Stan Model #####################################

covid_mod_file <- file.path("/ModelsFolder","covid_model_noyear_v1.3.stan") # find model file
covid_mod <- cmdstan_model(covid_mod_file) # compile model file

n_warmup = 5000
n_iter = 10000
n_chains = 4

# Fit with CmdStanR "Sample" method - runs NUTS MCMC algorithm

covid_model_fit <- covid_mod$sample(data = covid_list,
                                    seed = 281295,
                                    init = NULL,
                                    chains = n_chains,
                                    parallel_chains = n_chains,
                                    iter_warmup = n_warmup,
                                    iter_sampling = n_iter,
                                    max_treedepth = NULL,
                                    adapt_delta = 0.99,
                                    step_size = 0.01,
                                    diagnostics = c("divergences", "treedepth", "ebfmi"))

covid_model_summary <- covid_model_fit$summary() # Use 'posterior' package for summaries

covid_model_draws <- covid_model_fit$draws()        # Get posterior draws
covid_model_draws <- as_draws_df(covid_model_draws) # Convert to data frame

############################## Routine Checks ##################################

covid_model_fit$diagnostic_summary() # check sampling diagnostics

rhats <- covid_model_fit$summary(NULL, "rhat") # check to see if Rhat is pathological

names(covid_model_draws) # check column indexes of parameters

traces <- list() # set empty list to save trace plots
acf <- list()    # set empty list to save autocorrelation plots

color_scheme_set("mix-brightblue-gray")
traces[[1]] <- mcmc_trace(covid_model_draws, pars = names(covid_model_draws)[c(2:5,24:27)]) + 
  xlab("Post-warmup iteration") # check convergence diagnostics (Infections and InfProb)
traces[[2]] <- mcmc_trace(covid_model_draws, pars = names(covid_model_draws)[6:9]) + 
  xlab("Post-warmup iteration") # check convergence diagnostics (IFR)
traces[[3]] <- mcmc_trace(covid_model_draws, pars = names(covid_model_draws)[c(10:11,14:23,48:55)]) +
  xlab("Post-warmup iteration") # check convergence diagnostics (alphas and betas)
traces[[4]] <- mcmc_trace(covid_model_draws, pars = names(covid_model_draws)[c(12,13,28:35)]) + 
  xlab("Post-warmup iteration") # check convergence diagnostics (MUs and PHIs)
traces[[5]] <- mcmc_trace(covid_model_draws, pars = names(covid_model_draws)[36:47]) + 
  xlab("Post-warmup iteration") # check convergence diagnostics (ReportingProbs and true deaths)

pdf("traces_covid.pdf") # save traces to pdf
invisible(lapply(traces,print))
dev.off()

acf[[1]] <- mcmc_acf(covid_model_draws, pars = names(covid_model_draws)[c(2:5,24:27)], lags = 10)
acf[[2]] <- mcmc_acf(covid_model_draws, pars = names(covid_model_draws)[c(6:9)], lags = 10)
acf[[3]] <- mcmc_acf(covid_model_draws, pars = names(covid_model_draws)[c(10:11,48:55)], lags = 10)
acf[[4]] <- mcmc_acf(covid_model_draws, pars = names(covid_model_draws)[c(14:23)], lags = 10)
acf[[4]] <- mcmc_acf(covid_model_draws, pars = names(covid_model_draws)[c(12,13,28:35)], lags = 10)
acf[[5]] <- mcmc_acf(covid_model_draws, pars = names(covid_model_draws)[c(36:47)], lags = 10)

pdf("autocorr_covid.pdf") # save traces to pdf
invisible(lapply(acf,print))
dev.off()

color_scheme_set("mix-brightblue-gray")
mcmc_nuts_energy(nuts_params(covid_model_fit)) # check energy

# https://mc-stan.org/bayesplot/reference/MCMC-parcoord.html 
color_scheme_set("darkgray")
mcmc_parcoord(covid_model_fit$draws(), pars = names(covid_model_draws)[6:9])

# Geweke Plots
covid_iter_mcmc <- as_mcmc.list(covid_model_fit) # convert CmdStanMCMC to MCMC list obj
geweke.plot(covid_iter_mcmc, frac1 = 0.50, frac2 = 0.50, ask = TRUE)   # mostly behave well

########################### Display Key Results ################################

color_scheme_set(scheme = "blue")

# Use 'bayesplot' package for histograms
mcmc_hist(covid_model_fit$draws("true_ifr"), bins = 50)

# Use 'bayesplot' package for densities
plot_title <- ggtitle("True IFR Posterior distributions",
                      "medians and 95% intervals")
mcmc_areas(covid_model_fit$draws(format = "array"),
           pars = c("true_ifr[1]","true_ifr[2]","true_ifr[3]","true_ifr[4]"),
           prob = 0.95) + plot_title

mcmc_intervals(covid_model_fit$draws(format = "array"), pars = c("true_ifr[1]", "true_ifr[2]", "true_ifr[3]", "true_ifr[4]"))

# Overall - ALL CLUSTER YEAR COMBOS

covid_overall_ifr_wide <- covid_model_draws[c(6:9,71)] # all ifrs, all chains
covid_overall_ifr_long <- reshape2::melt(covid_overall_ifr_wide,
                                        id.vars = ".draw", 
                                        variable.name = "ClusterYearID",
                                        value.name = "True IFR")
# key overall IFR stats
quantile(covid_overall_ifr_long$`True IFR`, probs = c(0.025,0.50,0.975), digits = 2)

# plot it 
ggplot(covid_overall_ifr_long) +
  annotate('rect', 
           xmin=quantile(covid_overall_ifr_long$`True IFR`, probs = 0.025), 
           xmax=quantile(covid_overall_ifr_long$`True IFR`, probs = 0.975), 
           ymin=0, ymax=1.16e4, alpha=.6, fill='darkgray') + # shade credible interval
  geom_histogram(aes(x=`True IFR`), color ="#00AB8E", fill = "#00AB8E", alpha=0.9, bins = 50) +
  scale_x_continuous(breaks = seq(0,0.05,0.0025), labels = scales::percent_format()) +
  geom_vline(xintercept = median(covid_overall_ifr_long$`True IFR`), linetype = "longdash") +  # median line
  geom_vline(xintercept = quantile(covid_overall_ifr_long$`True IFR`, probs = 0.025), linetype = "dotted") + # LB line CrI
  geom_vline(xintercept = quantile(covid_overall_ifr_long$`True IFR`, probs = 0.975), linetype = "dotted") + # UB line CrI
  ylab("COVID-19") +
  xlab("True IFR Estimate") +
  theme_minimal() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# density plot for all IFRs overlapped 
# Use the white space to create a view of the overall IFR  posterior distribution

p <- covid_overall_ifr_long %>%
  ggplot(aes(x =`True IFR`, color=ClusterYearID, fill=ClusterYearID)) +
  geom_density(alpha=0.00) +
  viridis::scale_fill_viridis(discrete=TRUE) +
  viridis::scale_color_viridis(discrete=TRUE) +
  scale_x_continuous(breaks = seq(0,0.04,0.004), labels = scales::percent_format()) +
  xlab("True IFR Estimate") +
  ylab("Cluster Posterior Density") +
  theme_minimal() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none") 

inset <- ggplot(covid_overall_ifr_long) +
  annotate('rect', 
           xmin=quantile(covid_overall_ifr_long$`True IFR`, probs = 0.025), 
           xmax=quantile(covid_overall_ifr_long$`True IFR`, probs = 0.975), 
           ymin=0, ymax=250, alpha=.5, fill='darkgray') + # shade credible interval
  scale_x_continuous(breaks = seq(0,0.04,0.004), labels = scales::percent_format()) +
  geom_density(aes(x=`True IFR`), colour = "#00AB8E", fill = "#00AB8E", alpha=0.8) +
  geom_vline(xintercept = median(covid_overall_ifr_long$`True IFR`), linetype = "longdash") +  # median line
  geom_vline(xintercept = quantile(covid_overall_ifr_long$`True IFR`, probs = 0.025), linetype = "dotted") + # LB line CrI
  geom_vline(xintercept = quantile(covid_overall_ifr_long$`True IFR`, probs = 0.975), linetype = "dotted") + # UB line CrI
  ylab("Overall Posterior Density") +
  xlab("True IFR Estimate") +
  theme_minimal() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = 0.4))

ggdraw(p) +
  draw_plot(inset, x = .60, y = .60, width = .38, height = .38) +
  draw_plot_label(
    c("a", "b"),
    c(0, 0.55),
    c(1, 0.98),
    size = 12
  )


# Customise mcmc intervals for publication to match the spaghetti

d <- mcmc_intervals_data(covid_model_fit$draws(format = "array"), 
                         pars = c("true_ifr[1]","true_ifr[2]","true_ifr[3]","true_ifr[4]"),
                         prob_outer = 0.95)

# Generate viridis colors for each parameter
n_pars <- length(c("true_ifr[1]","true_ifr[2]","true_ifr[3]","true_ifr[4]"))
param_colors <- viridis::viridis(n_pars)
names(param_colors) <- c("true_ifr[1]","true_ifr[2]","true_ifr[3]","true_ifr[4]")

ggplot(d, aes(x = m, xmin = ll, xmax = hh, y = parameter, color = parameter)) +
  geom_point() +
  geom_errorbar(lineend = "butt", width = 0.5, orientation = "y", ) +
  scale_color_manual(values = param_colors) +
  scale_x_continuous(breaks = seq(0,0.04,0.004), limits = c(0.018,0.04), labels = scales::percent_format()) +
  xlab("True IFR Estimate") +
  ylab("") +
  theme_minimal() +
  theme(legend.position = "none")

########################### Posterior Predictive ###############################

cases_rep <- covid_model_fit$draws("cases_rep", format = "matrix") # workaround CmdStan output Env
dim(cases_rep)
# note that each row is a sample (predicted set of obs for each line in the dataset)

rep_select <- sample(1:40000, 10, replace = FALSE)  # sample IDs of random iterations to plot
ppc_hist(covid_stan_data$TotClustCases, cases_rep[rep_select, ], binwidth = 10)

graph_cases_all <- ppc_dens_overlay(y = covid_stan_data$TotClustCases, yrep = cases_rep) + xlim(0,2e3) # limited x axis loads to see real obs
ggsave("PostPredCases_Density_All_COVID.png", plot = graph_cases_all, bg = "white",
       path = "/FinalGraphs")

deaths_rep <- covid_model_fit$draws("deaths_rep", format = "matrix") # workaround CmdStan output Env
dim(deaths_rep)

graph_deaths_all <- ppc_dens_overlay(y = covid_stan_data$TotClustDeaths, yrep = deaths_rep) + xlim(c(0,100))
ggsave("PostPredDeaths_Density_All_COVID.png", plot = graph_deaths_all, bg = "white",
       path = "/FinalGraphs")

rep_select <- sample(1:40000, 10, replace = FALSE)  # sample IDs of random iterations to plot
ppc_hist(covid_stan_data$TotClustDeaths, deaths_rep[200:208, ], binwidth = 10)

# We want to plot the histograms of Cases in transparency

cases_rep1 <- covid_model_fit$draws("cases_rep", format = "df") %>%
  rename("Cluster 0" = `cases_rep[1]`,
         "Cluster 1" = `cases_rep[2]`,
         "Cluster 2" = `cases_rep[3]`,
         "Cluster 3" = `cases_rep[4]`)

ggplot() +
  geom_histogram(aes(x = cases_rep1$`Cluster 0`, fill = "Cluster 0"), bins = 100, alpha = 0.2) +
  geom_histogram(aes(x = cases_rep1$`Cluster 1`, fill = "Cluster 1"), bins = 100, alpha = 0.2) +
  geom_histogram(aes(x = cases_rep1$`Cluster 2`, fill = "Cluster 2"), bins = 100, alpha = 0.2) +
  geom_histogram(aes(x = cases_rep1$`Cluster 3`, fill = "Cluster 3"), bins = 100, alpha = 0.2) +
  scale_fill_manual(values = c("Cluster 0" = "#440154", 
                               "Cluster 1" = "#31688e", 
                               "Cluster 2" = "#35b779",
                               "Cluster 3" = "#fde725")) +
  theme_minimal() +
  ggtitle("Cases Replicates from Model",
          subtitle = "by cluster") +
  xlab("Cases") +
  ylab("Frequency")

# Removing high highs

ggplot() +
  geom_histogram(aes(x = cases_rep1$`Cluster 0`/1000, fill = "Cluster 0"), bins = 100, alpha = 0.3) +
  geom_histogram(aes(x = cases_rep1$`Cluster 1`/1000, fill = "Cluster 1"), bins = 100, alpha = 0.3) +
  geom_histogram(aes(x = cases_rep1$`Cluster 2`/1000, fill = "Cluster 2"), bins = 100, alpha = 0.3) +
  geom_histogram(aes(x = cases_rep1$`Cluster 3`/1000, fill = "Cluster 3"), bins = 100, alpha = 0.3) +
  scale_fill_manual(values = c("Cluster 0" = "#440154", 
                               "Cluster 1" = "#31688e", 
                               "Cluster 2" = "#35b779",
                               "Cluster 3" = "#fde725")) +
  theme_minimal() +
  ylim(c(0,1.5e4)) + # removing very high highs, only makes sense paired with xlim 0,5e5
  xlim(c(0,4e2)) + # removing very high highs, I think it was just a couple
  ggtitle("Cases Replicates from Model",
          subtitle = "by cluster") +
  xlab("Cases (thousands)") +  
  ylab("Frequency")

ggsave("Cases Replicates from Model by Cluster.png",
       dpi = 300, bg = "white",
       path = "/FinalGraphs")

# We want to plot the histograms of Deaths in transparency

deaths_rep1 <- covid_model_fit$draws("deaths_rep", format = "df") %>%
  rename("Cluster 0" = `deaths_rep[1]`,
         "Cluster 1" = `deaths_rep[2]`,
         "Cluster 2" = `deaths_rep[3]`,
         "Cluster 3" = `deaths_rep[4]`)

ggplot() +
  geom_histogram(aes(x = deaths_rep1$`Cluster 0`, fill = "Cluster 0"), bins = 1000, alpha = 0.3) +
  geom_histogram(aes(x = deaths_rep1$`Cluster 1`, fill = "Cluster 1"), bins = 1000, alpha = 0.3) +
  geom_histogram(aes(x = deaths_rep1$`Cluster 2`, fill = "Cluster 2"), bins = 1000, alpha = 0.3) +
  geom_histogram(aes(x = deaths_rep1$`Cluster 3`, fill = "Cluster 3"), bins = 1000, alpha = 0.3) +
  scale_fill_manual(values = c("Cluster 0" = "#440154", 
                               "Cluster 1" = "#31688e", 
                               "Cluster 2" = "#35b779",
                               "Cluster 3" = "#fde725")) +
  theme_minimal() +
 #xlim(c(0,25000)) +
  #ylim(c(0,7500)) +
  ggtitle("Deaths Replicates from Model",
          subtitle = "by cluster") +
  xlab("Deaths") +  
  ylab("Frequency")

ggsave("Deaths Replicates from Model by Cluster.png",
       dpi = 300, bg = "white",
       path = "/FinalGraphs")


########################## Sensitivity to Priors ###############################

# Check documentation of function
# https://n-kall.github.io/priorsense/reference/powerscale-sensitivity.html 
covid_sens <-powerscale_sensitivity(covid_model_fit,
                          div_measure = "cjs_dist",
                          component = c("prior", "likelihood"),
                          sensitivity_threshold = 0.05,
                          moment_match = FALSE,
                          k_threshold = 0.5,
                          transform = NULL,
                          prediction = NULL)

View(covid_sens)

pws <- list()

pws[[1]] <- powerscale_plot_dens(covid_model_fit, variable = "true_ifr")
pws[[2]] <- powerscale_plot_dens(covid_model_fit, variable = "infection_prob")
pws[[3]] <- powerscale_plot_dens(covid_model_fit, variable = "infections")
pws[[4]] <- powerscale_plot_dens(covid_model_fit, variable = "mu_cases")
pws[[5]] <- powerscale_plot_dens(covid_model_fit, variable = "mu_deaths")
pws[[6]] <- powerscale_plot_dens(covid_model_fit, variable = "true_deaths")

# pdf("powerscale_covid.pdf") # save traces to pdf - might not be necessary if no prior data conflict
# invisible(lapply(pws,print))
# dev.off()


###################### Run double iteration for local conv. ####################

n_warmup = 5000 * 2
n_iter = 10000 * 2
n_chains = 4

# Re run the code above
# Fit with CmdStanR "Sample" method - runs NUTS MCMC algorithm



