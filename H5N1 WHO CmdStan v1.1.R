#################################################################
#  This script is needed to run the model H5N1 V1.1             #
#  with clusters instead of countries                           #
#                                                               #
#                                                               # 
#  Author: Leonardo Gada                                        #   
#                                                               # 
#################################################################


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
  
  ####################### Read Data ##############################################
  
  # Read in data prepared by Joe Gibson by taking info from the WHO H5N1 cases and deaths tracker (2003-2024)
  
  who_0324 <- read.csv("/Datasets/WHO 2003-2024 Clean Dataset.csv")
  
  unique(who_0324$country) # 24 different countries all over the world
  
  
  # Read in World Bank Population Data for the countries in the WHO Tracker
  
  pop_data_0324 <- read_excel("/Datasets/TotPop 2003 2024 WHO Tracker Countries.xlsx", sheet = "Data")
  pop_data_0324$MeanPop <- as.numeric(round(pop_data_0324$MeanPop, digits = 0))   # Round to integer
  
  pop_data_0324$`Country Name`[pop_data_0324$`Country Name` == "Egypt, Arab Rep."] <- "Egypt"                     # consistency changes
  pop_data_0324$`Country Name`[pop_data_0324$`Country Name` == "Turkiye"] <- "Turkey"                             # consistency changes
  pop_data_0324$`Country Name`[pop_data_0324$`Country Name` == "United States"] <- "United States of America"     # consistency changes
  pop_data_0324$`Country Name`[pop_data_0324$`Country Name` == "Lao PDR"] <- "Lao People's Democratic Republic"   # consistency changes
  
  sum(unique(pop_data_0324$`Country Name`) == unique(who_0324$country))  # Check all countries match now (they do)
  
  pop_data_0324 <- pop_data_0324 %>% 
    select(- c(`Series Name`,`Series Code`, `Country Code`)) 
  
  for (i in 2:(length(names(pop_data_0324))-1)){ # we don't want to rename the mean column so -1
    if (i < 9){
      names(pop_data_0324)[i] <- paste0("200",i+1) # rename
    } else {
      names(pop_data_0324)[i] <- paste0("20",i+1) # rename
    }
  }
  # Change format of pop size dataset to match format of cases deaths data
  pop_data_0324_long <- reshape2::melt(pop_data_0324, id.vars = "Country Name", variable.name = "Year", value.name = "PopSize")
  
  # Join pop size info to the cases deaths tracker data
  
  who_0324$year <- as.factor(who_0324$year) # consistency change to make the join work, note no cases in a few years
  who_popsize <- left_join(who_0324, pop_data_0324_long, by = c("country" = "Country Name",
                                                                "year" = "Year"))
  
  ############################ Create Clusters  ##################################
  
  # Instead of working with individual countries, we work with clusters.
  # Clustering was achieved with hierarchical clustering based on a set of 
  # World Development Indicators relevant to characterisation for H5N1 response
  
  # CLUSTER SUBDIVISION v1 (Manhattan + CART + dynamic tree cutting)
  # As per clustering paper, outlined in the excel file: CART_Clusters_1 from Mwandida Afuleni
  
  clustered_data <- who_popsize %>% 
    mutate(
      cluster = case_when(country %in% c("Australia", "Canada", "Spain", "United Kingdom", "United States of America") ~ "Cluster0",
                          country %in% c("Bangladesh", "Cambodia", "China", "Djibouti", 
                                         "India", "Indonesia", "Lao People's Democratic Republic", "Myanmar", "Nepal", 
                                         "Nigeria", "Pakistan", "Thailand", "Viet Nam") ~ "Cluster1",
                          country %in% c("Azerbaijan", "Ecuador", "Egypt", "Iraq", "Turkey") ~ "Cluster2",
                          country %in% c("Chile") ~ "Cluster3"))
  
  # Let's save the mapping country to cluster
  clustered_data_PreAggregation <- clustered_data
  
  # Let's sort the population size, cases and deaths by cluster rather than by country
  
  clustered_data <- clustered_data %>% 
    group_by(cluster, year) %>% 
    summarise(ClusterPop = sum(PopSize), ClusterCases = sum(cases), ClusterDeaths = sum(deaths))
  
  
  ############################## Prepare Stan List ###############################
  
  h5n1_data <- list(N = nrow(clustered_data),                                # Number of observations in the dataset
                    C = length(unique(clustered_data$cluster)),              # Number of distinct countries in the dataset
                    Ye = length(unique(clustered_data$year)),                # Number of distinct years in the dataset
                    cases = clustered_data$ClusterCases,                     # Recorded cases observations by country by year
                    deaths = clustered_data$ClusterDeaths,                   # Recorded deaths observations by country by year
                    country = as.integer(as.factor(clustered_data$cluster)), # Vector of matching cluster per case/death pair row
                    population_size = clustered_data$ClusterPop,             # Vector of pop size  by country by year
                    year =as.integer(as.factor(clustered_data$year)))        # Vector of matching year per case/death pair row
  
  
  ############################### Stan Model #####################################
  
  h5n1_mod_file <- file.path("/ModelsFolder","H5N1_model_v1.1.stan") # find model file
  h5n1_mod <- cmdstan_model(h5n1_mod_file) # compile model file
  
  n_warmup = 5000
  n_iter = 10000
  n_chains = 4
  
  # Fit with CmdStanR "Sample" method - runs NUTS MCMC algorithm
  
  h5n1_model_fit <- h5n1_mod$sample(data = h5n1_data,
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

h5n1_mod_summary <- h5n1_model_fit$summary() # Use 'posterior' package for summaries

h5n1_model_draws <- h5n1_model_fit$draws()        # Get posterior draws
h5n1_model_draws <- as_draws_df(h5n1_model_draws) # Convert to data frame

############################## Routine Checks ##################################

h5n1_model_fit$diagnostic_summary() # check sampling diagnostics

rhats <- h5n1_model_fit$summary(NULL, "rhat") # check to see if Rhat is pathological

names(h5n1_model_draws) # check column indexes of parameters

traces_h5n1 <- list()   # empty list to save trace plots
autocorr_h5n1 <- list() # empty list to save autocorrelation plots

color_scheme_set("mix-brightblue-gray")
traces_h5n1[[1]] <- mcmc_trace(h5n1_model_draws, pars = names(h5n1_model_draws)[41:79]) + 
  xlab("Post-warmup iteration") # check convergence diagnostics (IFR)
traces_h5n1[[2]] <- mcmc_trace(h5n1_model_draws, pars = names(h5n1_model_draws)[2:40]) + 
  xlab("Post-warmup iteration") # check convergence diagnostics (InfectionProb)
traces_h5n1[[3]] <- mcmc_trace(h5n1_model_draws, pars = names(h5n1_model_draws)[c(80,81,84:93)]) + 
  xlab("Post-warmup iteration") # check convergence diagnostics (alphas and betas std)
traces_h5n1[[4]] <- mcmc_trace(h5n1_model_draws, pars = names(h5n1_model_draws)[c(94:115)]) + 
  xlab("Post-warmup iteration") # check convergence diagnostics (gammas 1)
traces_h5n1[[5]] <- mcmc_trace(h5n1_model_draws, pars = names(h5n1_model_draws)[c(116:137)]) + 
  xlab("Post-warmup iteration") # check convergence diagnostics (gammas 2)
traces_h5n1[[6]] <- mcmc_trace(h5n1_model_draws, pars = names(h5n1_model_draws)[c(138:176)]) + 
  xlab("Post-warmup iteration") # check convergence diagnostics (infections)
traces_h5n1[[7]] <- mcmc_trace(h5n1_model_draws, pars = names(h5n1_model_draws)[c(177:215,82)]) + 
  xlab("Post-warmup iteration") # check convergence diagnostics (mu cases and phi cases)
traces_h5n1[[8]] <- mcmc_trace(h5n1_model_draws, pars = names(h5n1_model_draws)[c(216:254,83)]) + 
  xlab("Post-warmup iteration") # check convergence diagnostics (mu deaths and phi deaths)
traces_h5n1[[9]] <- mcmc_trace(h5n1_model_draws, pars = names(h5n1_model_draws)[c(255:293)]) + 
  xlab("Post-warmup iteration") # check convergence diagnostics (reporting prob cases)
traces_h5n1[[10]] <- mcmc_trace(h5n1_model_draws, pars = names(h5n1_model_draws)[c(294:332)]) + 
  xlab("Post-warmup iteration") # check convergence diagnostics (reporting prob deaths)
traces_h5n1[[11]] <- mcmc_trace(h5n1_model_draws, pars = names(h5n1_model_draws)[c(333:371)]) + 
  xlab("Post-warmup iteration") # check convergence diagnostics (true deaths)
traces_h5n1[[12]] <- mcmc_trace(h5n1_model_draws, pars = names(h5n1_model_draws)[c(372:421)]) + 
  xlab("Post-warmup iteration") # check convergence diagnostics (beta and gammas cases and deaths)

# Loop through the plots and save them all as PNG

for (i in seq_along(traces_h5n1)) {
  filename <- paste0("traces_h5n1_", i, ".png")
  
  # Save using ggsave
  ggsave(filename, plot = traces_h5n1[[i]], bg = "white", width = 20, height = 15, units = "in", dpi = 300)
}

autocorr_h5n1[[1]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[41:51], lags = 10) # autocorrelation true ifr
autocorr_h5n1[[2]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[52:65], lags = 10) # autocorrelation true ifr
autocorr_h5n1[[3]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[66:79], lags = 10) # autocorrelation true ifr
autocorr_h5n1[[4]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[2:14], lags = 10)  # autocorrelation infection prob
autocorr_h5n1[[5]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[15:28], lags = 10) # autocorrelation infection prob
autocorr_h5n1[[6]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[29:40], lags = 10) # autocorrelation infection prob
autocorr_h5n1[[7]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[c(80,81,84:93)], lags = 10)  # autocorrelation alphas and betas std
autocorr_h5n1[[8]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[94:104], lags = 10) # autocorrelation gammas
autocorr_h5n1[[9]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[105:115], lags = 10) # autocorrelation gammas
autocorr_h5n1[[10]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[116:126], lags = 10) # autocorrelation gammas
autocorr_h5n1[[11]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[127:137], lags = 10) # autocorrelation gammas
autocorr_h5n1[[12]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[138:149], lags = 10) # autocorrelation infections
autocorr_h5n1[[13]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[150:162], lags = 10) # autocorrelation infections
autocorr_h5n1[[14]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[163:176], lags = 10) # autocorrelation infections
autocorr_h5n1[[15]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[c(82,177:189)], lags = 10) # autocorrelation mu and phi cases
autocorr_h5n1[[16]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[190:202], lags = 10) # autocorrelation mu and phi cases
autocorr_h5n1[[17]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[203:215], lags = 10) # autocorrelation mu and phi cases
autocorr_h5n1[[18]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[c(83,216:226)], lags = 10) # autocorrelation mu deaths and phi deaths
autocorr_h5n1[[19]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[227:240], lags = 10) # autocorrelation mu deaths and phi deaths
autocorr_h5n1[[20]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[241:254], lags = 10) # autocorrelation mu deaths and phi deaths
autocorr_h5n1[[21]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[255:268], lags = 10) # autocorrelation reporting prob cases
autocorr_h5n1[[22]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[269:282], lags = 10) # autocorrelation reporting prob cases
autocorr_h5n1[[23]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[283:293], lags = 10) # autocorrelation reporting prob cases
autocorr_h5n1[[24]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[294:307], lags = 10) # autocorrelation reporting prob deaths
autocorr_h5n1[[25]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[308:320], lags = 10) # autocorrelation reporting prob deaths
autocorr_h5n1[[26]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[321:332], lags = 10) # autocorrelation reporting prob deaths
autocorr_h5n1[[27]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[333:346], lags = 10) # autocorrelation true deaths
autocorr_h5n1[[28]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[347:359], lags = 10) # autocorrelation true deaths
autocorr_h5n1[[29]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[360:371], lags = 10) # autocorrelation true deaths
autocorr_h5n1[[30]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[372:385], lags = 10) # autocorrelation beta and gammas cases and deaths
autocorr_h5n1[[31]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[386:400], lags = 10) # autocorrelation beta and gammas cases and deaths
autocorr_h5n1[[32]] <- mcmc_acf(h5n1_model_draws, pars = names(h5n1_model_draws)[401:421], lags = 10) # autocorrelation beta and gammas cases and deaths


color_scheme_set("blue")
mcmc_nuts_energy(nuts_params(h5n1_model_fit), bins = 50) # check energy, https://arxiv.org/pdf/1701.02434 

color_scheme_set("darkgray")
mcmc_parcoord(h5n1_model_fit$draws(), pars = names(h5n1_model_draws)[41:51]) # divergences, parallel coordinates
mcmc_parcoord(h5n1_model_fit$draws(), pars = names(h5n1_model_draws)[52:65]) # divergences, parallel coordinates
mcmc_parcoord(h5n1_model_fit$draws(), pars = names(h5n1_model_draws)[66:79]) # divergences, parallel coordinates

# Geweke Plots
h5n1_iter_mcmc <- as_mcmc.list(h5n1_model_fit) # convert CmdStanMCMC to MCMC list obj
geweke.plot(h5n1_iter_mcmc, frac1 = 0.50, frac2 = 0.50, ask = TRUE)   # mostly behave well

########################### Display Key Results ################################

param <- names(h5n1_model_draws) # check column indexes of parameters

color_scheme_set(scheme = "blue")

# Use 'bayesplot' package for histograms
mcmc_hist(h5n1_model_fit$draws("true_ifr"), bins = 50)

# Use 'bayesplot' package for densities - not best plot here bcause too many
plot_title <- ggtitle("True IFR Posterior distributions",
                      "medians and 95% intervals")
mcmc_areas(h5n1_model_fit$draws(format = "array"),
           pars = param[41:50],
           prob = 0.95) + plot_title

# Use 'bayesplot' package for CrIs - better than above

mcmc_intervals(h5n1_model_fit$draws(format = "array"), 
               pars = param[41:79],
               prob = 0.50,
               prob_outer = 0.95)
  # Should also be customisable further with ggplot

# Overall - ALL CLUSTER YEAR COMBOS

h5n1_overall_ifr_wide <- h5n1_model_draws[c(41:79,542)] # all ifrs, all chains
h5n1_overall_ifr_long <- reshape2::melt(h5n1_overall_ifr_wide,
                                        id.vars = ".draw", 
                                        variable.name = "ClusterYearID",
                                        value.name = "True IFR")
# key overall IFR stats
quantile(h5n1_overall_ifr_long$`True IFR`, probs = c(0.025,0.50,0.975), digits = 2)

# plot it 
ggplot(h5n1_overall_ifr_long) +
  annotate('rect', 
           xmin=quantile(h5n1_overall_ifr_long$`True IFR`, probs = 0.025), 
           xmax=quantile(h5n1_overall_ifr_long$`True IFR`, probs = 0.975), 
           ymin=0, ymax=1.1e5, alpha=.6, fill='darkgray') + # shade credible interval
  geom_histogram(aes(x=`True IFR`), color ="#1D57A5", fill = "#1D57A5", alpha=0.9, bins = 50) +
  scale_x_continuous(breaks = seq(0,1,0.1), labels = scales::percent_format()) +
  geom_vline(xintercept = median(h5n1_overall_ifr_long$`True IFR`), linetype = "longdash") +  # median line
  geom_vline(xintercept = quantile(h5n1_overall_ifr_long$`True IFR`, probs = 0.025), linetype = "dotted") + # LB line CrI
  geom_vline(xintercept = quantile(h5n1_overall_ifr_long$`True IFR`, probs = 0.975), linetype = "dotted") + # UB line CrI
  ylab("H5N1") +
  xlab("True IFR Estimate") +
  theme_minimal() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# density plot for all IFRs overlapped 
# Use the white space to create a view of the overall IFR  posterior distribution

p <- h5n1_overall_ifr_long %>%
  ggplot(aes(x =`True IFR`, color=ClusterYearID, fill=ClusterYearID)) +
  geom_density(alpha=0.00) +
  viridis::scale_fill_viridis(discrete=TRUE) +
  viridis::scale_color_viridis(discrete=TRUE) +
  scale_x_continuous(breaks = seq(0,1,0.1), labels = scales::percent_format()) +
  xlab("True IFR Estimate") +
  ylab("Cluster-Year Posterior Density") +
  theme_minimal() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none") 

inset <- ggplot(h5n1_overall_ifr_long) +
  annotate('rect', 
           xmin=quantile(h5n1_overall_ifr_long$`True IFR`, probs = 0.025), 
           xmax=quantile(h5n1_overall_ifr_long$`True IFR`, probs = 0.975), 
           ymin=0, ymax=3.8, alpha=.5, fill='darkgray') + # shade credible interval
  scale_x_continuous(breaks = seq(0,1,0.1), labels = scales::percent_format()) +
  geom_density(aes(x=`True IFR`), colour = "#1D57A5", fill = "#1D57A5", alpha=0.8) +
  geom_vline(xintercept = median(h5n1_overall_ifr_long$`True IFR`), linetype = "longdash") +  # median line
  geom_vline(xintercept = quantile(h5n1_overall_ifr_long$`True IFR`, probs = 0.025), linetype = "dotted") + # LB line CrI
  geom_vline(xintercept = quantile(h5n1_overall_ifr_long$`True IFR`, probs = 0.975), linetype = "dotted") + # UB line CrI
  ylab("Overall Posterior Density") +
  xlab("True IFR Estimate") +
  theme_minimal() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = 0.4))

ggdraw(p) +
  draw_plot(inset, x = .40, y = .45, width = .55, height = .55) +
  draw_plot_label(
    c("a", "b"),
    c(0, 0.40),
    c(1, 0.98),
    size = 12
  )

# Customise mcmc intervals for publication to match the spaghetti

d <- mcmc_intervals_data(h5n1_model_fit$draws(format = "array"), 
                         pars = param[41:79],
                         prob_outer = 0.95)

# Generate viridis colors for each parameter
n_pars <- length(param[41:79])
param_colors <- viridis::viridis(n_pars)
names(param_colors) <- param[41:79]

ggplot(d, aes(x = m, xmin = ll, xmax = hh, y = parameter, color = parameter)) +
  geom_point() +
  geom_errorbar(lineend = "butt", width = 0.5, orientation = "y", ) +
  scale_color_manual(values = param_colors) +
  scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0,1), labels = scales::percent_format()) +
  xlab("True IFR Estimate") +
  ylab("") +
  theme_minimal() +
  theme(legend.position = "none")

########################### Posterior Predictive ###############################

cases_rep <- h5n1_model_fit$draws("cases_rep", format = "matrix") # workaround CmdStan output Env
dim(cases_rep)
# note that each row is a sample (predicted set of obs for each line in the dataset)

graph_cases_all <- ppc_dens_overlay(y =clustered_data$ClusterCases, yrep = cases_rep) + xlim(c(0,1.5e3))
ggsave("PostPredCases_Density_All_H5N1.png", plot = graph_cases_all,
       dpi = 300, bg = "white",
       path = "/FinalGraphs")

rep_select <- sample(1:40000, 10, replace = FALSE) # sample IDs of random iterations to plot
ppc_hist(clustered_data$ClusterCases, cases_rep[rep_select, ], binwidth = 10)

deaths_rep <- h5n1_model_fit$draws("deaths_rep", format = "matrix") # workaround CmdStan output Env
dim(deaths_rep)

graph_deaths_all <- ppc_dens_overlay(y = clustered_data$ClusterDeaths, yrep = deaths_rep) + xlim(c(0,100))
ggsave("PostPredDeaths_Density_All_H5N1.png", plot = graph_deaths_all,
       dpi = 300, bg = "white",
       path = "/FinalGraphs")

rep_select <- sample(1:40000, 10, replace = FALSE) # sample IDs of random iterations to plot
ppc_hist(clustered_data$ClusterDeaths, deaths_rep[rep_select,], binwidth = 10)

########################## Sensitivity to Priors ###############################

# Check documentation of function
# https://n-kall.github.io/priorsense/reference/powerscale-sensitivity.html 

h5n1_sens <-powerscale_sensitivity(h5n1_model_fit,
                                   div_measure = "cjs_dist",
                                   component = c("prior", "likelihood"),
                                   sensitivity_threshold = 0.05,
                                   moment_match = FALSE,
                                   k_threshold = 0.5,
                                   transform = NULL,
                                   prediction = NULL)

View(h5n1_sens)

h5n1_sens <- h5n1_sens %>% filter(diagnosis != "-")

pws <- list()

# Go through and save one by one + the other ones that have potential prior-data conflict

# pws[[1]] <- powerscale_plot_dens(h5n1_model_fit, variable = c("true_ifr[1]", "true_ifr[2]","true_ifr[3]","true_ifr[4]","true_ifr[5]"))
# pws[[2]] <- powerscale_plot_dens(h5n1_model_fit, variable = c("true_ifr[6]", "true_ifr[7]","true_ifr[9]","true_ifr[10]"))
# pws[[3]] <- powerscale_plot_dens(h5n1_model_fit, variable = c("true_ifr[11]", "true_ifr[12]","true_ifr[13]","true_ifr[14]","true_ifr[15]"))
# pws[[4]] <- powerscale_plot_dens(h5n1_model_fit, variable = c("true_ifr[17]","true_ifr[18]","true_ifr[19]","true_ifr[20]"))
# pws[[5]] <- powerscale_plot_dens(h5n1_model_fit, variable = c("true_ifr[21]", "true_ifr[22]","true_ifr[23]","true_ifr[24]"))
# pws[[6]] <- powerscale_plot_dens(h5n1_model_fit, variable = c("true_ifr[28]","true_ifr[32]","true_ifr[33]"))
# pws[[7]] <- powerscale_plot_dens(h5n1_model_fit, variable = c("infection_prob[1]","true_ifr[32]","true_ifr[33]"))

pws[[1]] <- powerscale_plot_dens(h5n1_model_fit, variable = h5n1_sens$variable[1:6])   # infection prob
pws[[2]] <- powerscale_plot_dens(h5n1_model_fit, variable = h5n1_sens$variable[7:15])  # true ifr
pws[[3]] <- powerscale_plot_dens(h5n1_model_fit, variable = h5n1_sens$variable[16:23]) # true ifr
pws[[4]] <- powerscale_plot_dens(h5n1_model_fit, variable = h5n1_sens$variable[24:30]) # true ifr
pws[[5]] <- powerscale_plot_dens(h5n1_model_fit, variable = h5n1_sens$variable[31:35]) # true ifr
pws[[6]] <- powerscale_plot_dens(h5n1_model_fit, variable = h5n1_sens$variable[36:41]) # infections
pws[[7]] <- powerscale_plot_dens(h5n1_model_fit, variable = h5n1_sens$variable[42:44]) # mu deaths
pws[[8]] <- powerscale_plot_dens(h5n1_model_fit, variable = h5n1_sens$variable[45:47]) # true deaths

###################### Run double iteration for local conv. ####################

n_warmup = 5000 * 2
n_iter = 10000 * 2
n_chains = 4

# Re run the code above
# Fit with CmdStanR "Sample" method - runs NUTS MCMC algorithm

# https://stackoverflow.com/questions/59608999/cmdstanr-extracting-draws-from-stan-model-fit 




