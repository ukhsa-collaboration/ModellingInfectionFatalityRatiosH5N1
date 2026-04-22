data {
  int<lower=0> N;                // Number of observations
  int<lower=0> C;                // Number of countries
  int<lower=0> Ye;               // Number of years
  array[N] int<lower=0> cases;   // Observed cases by country by year
  array[N] int<lower=0> deaths;  // Observed deaths by country by year
  array[N] int<lower=1,upper=C> country;    // grouping variable, encoding which country the observation comes from
  array[N] int<lower=0> population_size;    // Number of people by country (repeated avg per year as well, can change to be year specific pop)
  array[N] int<lower=1,upper=Ye> year;       // grouping variable, encoding which year the observation comes from
}

parameters {
  
  vector<lower=0,upper=1>[N] infection_prob;  // Probability of infections per country per year 
  vector<lower=0,upper=1>[N] true_ifr;    // Real infection fatality ratio 
  
  real<lower=0,upper=1> alpha_cases;   // Global intercept, baseline reporting rate for cases (from seasonal flu?)
  real<lower=0,upper=1> alpha_deaths;  // Global intercept, baseline reporting rate for deaths (from seasonal flu?)
    
  real<lower=0> phi_cases;             // Dispersion parameter for cases
  real<lower=0> phi_deaths;            // Dispersion parameter for deaths
  
  vector[C] beta_cases_std;       // Define (non-centered param) beta coefficients  (country random effect for cases) 
  real<lower=0> b_sigma_cases;    // Define standard deviation of country coefficient beta cases
  vector[C] beta_deaths_std;      // Define (non-centered param) beta coefficients  (country random effect for deaths)
  real<lower=0> b_sigma_deaths;   // Define standard deviation of country coefficient beta deaths
  
  vector[Ye] gamma_cases_std;      // Define (non-centered param) gamma coefficients (year random effect for cases)
  real<lower=0> g_tau_cases;       // Define standard deviation of year coefficient gamma cases
  vector[Ye] gamma_deaths_std;     // Define (non-centered param) gamma coefficients (year random effect for deaths)
  real<lower=0> g_tau_deaths;      // Define standard deviation of year coefficient gamma deaths
}

transformed parameters{
  
  vector<lower=0>[N] infections;      // Total number of infections per country per year (hidden layer expressed by cases)
  vector<lower=0>[N] mu_cases;        // Negative Binomial parameter for cases
  vector<lower=0>[N] mu_deaths;       // Negative Binomial parameter for deaths
  vector<lower=0,upper=1>[N] p;       // Reporting probability for cases (hidden layer parameter) 
  vector<lower=0,upper=1>[N] q;       // Reporting probability for deaths (hidden layer parameter)
  vector<lower=0>[N]  true_deaths;    // Total number of deaths per country per year (hidden layer expressed by deaths and derived by ifr and infections)
  vector[C] beta_cases;               // Param to host non standardised beta cases (non centred param)
  vector[C] beta_deaths;              // Param to host non standardised beta deaths (non centred param)
  vector[Ye] gamma_cases;             // Param to host non standardised gamma cases (non centred param)
  vector[Ye] gamma_deaths;            // Param to host non standardised gamma deaths (non centred param)
  
  

  
  for (c in 1:C){
    beta_cases[c] = beta_cases_std[c] * b_sigma_cases; // non-centered param for cases by country
    beta_deaths[c] = beta_deaths_std[c] * b_sigma_deaths; // non-centered param for deaths by country
  }
  for (t in 1:Ye){
    gamma_cases[t] = gamma_cases_std[t] * g_tau_cases; // non-centred param for cases by year
    gamma_deaths[t] = gamma_deaths_std[t] * g_tau_deaths; // non-centred param for deaths by year
  }
  
  for (n in 1:N){ 
    
  infections[n] = infection_prob[n] * population_size[n]; // true infections taking into account country size   
    
  p[n] = inv_logit(alpha_cases + beta_cases[country[n]] + gamma_cases[year[n]]);  // p = logit^-1(alpha + beta + gamma)
  mu_cases[n] = infections[n]*p[n]; 
  
  true_deaths[n] = infections[n] * true_ifr[n];    // relationship between infections and true deaths 
  
  q[n] = inv_logit(alpha_deaths + beta_deaths[country[n]] + gamma_deaths[year[n]]); // q = logit^-1(alpha + beta + gamma)
  mu_deaths[n] = true_deaths[n]*q[n]; 
  }
}

model {
  
  // Prior for infections probability and true deaths
  infection_prob ~ beta(0.5,100); // uninformative prior (seasonal flu infection probability)
  true_ifr ~ beta(1.25,2);        // uninformative prior centred on estimate from Ward et al. 2024
  infections ~ gamma(1.85,0.01);  // uninformative prior with mean = max of cases and stdev = 100
  true_deaths ~ gamma(1,0.06);    // uninformative prior with mean = mean deaths + 5 = 17
  
  p ~ beta(1.25,3); // Beta prior compatible with Epi Review of Surveillance during COVID-19 first wave
  q ~ beta(3,1.50); // Beta prior compatible with assumption that deaths are more evident than infections
  
  // Priors for case and death parameters
  alpha_cases ~ normal(0,100);     // uninformative prior
  alpha_deaths ~ normal(0,100);    // uninformative prior
  
  // Priors for 
  mu_cases ~ gamma(0.54,0.02);  // solving for alpha*theta = mean and alpha*(theta)^2 = var = 1000
  mu_deaths ~ gamma(0.47,0.04); // solving for alpha*theta = mean and alpha*(theta)^2 = var = 300
  
  // Priors for dispersion parameters
  phi_cases ~ gamma(1,1);        // Gamma prior for case dispersion 
  phi_deaths ~ gamma(1,1);       // Gamma prior for death dispersion 

  // Model the parameters for cases and deaths (betas, gammas)
  
  b_sigma_cases ~ cauchy(0,2.5);         // Half Cauchy due to truncation at 0 in parameter definition    
  b_sigma_deaths ~ cauchy(0,2.5);        // Half Cauchy due to truncation at 0 in parameter definition
  
  for (c in 1:C){
  beta_cases_std[c] ~ std_normal();      // Standard normal prior for each country cases parameter
  beta_deaths_std[c] ~ std_normal();     // Standard normal prior for each country deaths parameter
  }
  
  g_tau_cases ~ cauchy(0,2.5);           // Half Cauchy due to truncation at 0 in parameter definition
  g_tau_deaths ~ cauchy(0,2.5);          // Half Cauchy due to truncation at 0 in parameter definition
  
  for (t in 1:Ye){
  gamma_cases_std[t] ~ std_normal();     // Standard normal prior for each year cases parameter
  gamma_deaths_std[t] ~ std_normal();    // Standard normal prior for each year deaths parameter
  }

  // Likelihood for cases and deaths using Negative Binomial distribution
  for (n in 1:N) {
    cases[n] ~ neg_binomial_2(mu_cases[n], phi_cases);  // Negative Binomial for cases
    deaths[n] ~ neg_binomial_2(mu_deaths[n], phi_deaths); // Negative Binomial for deaths
  }
}

generated quantities{
  array[N] int<lower=0> cases_rep;
  array[N] int<lower=0> deaths_rep;
  vector[N] log_lik;
  real lprior;
  
  for (n in 1:N){
   cases_rep[n] = neg_binomial_2_rng(mu_cases[n], phi_cases); // Posterior Predictive
   deaths_rep[n] = neg_binomial_2_rng(mu_deaths[n], phi_deaths); // Posterior Predictive
   log_lik[n] = beta_lpdf(true_ifr[n] | 1.25, 2); // prior sensitivity check
  }
  
  lprior =  beta_lpdf(true_ifr | 1.25, 2); // prior sensitivity check 
}



