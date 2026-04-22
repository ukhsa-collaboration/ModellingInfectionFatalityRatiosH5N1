data {
  int<lower=0> N;                // Number of observations
  int<lower=0> C;                // Number of countries
  array[N] int<lower=0> cases;   // Observed cases by country by year
  array[N] int<lower=0> deaths;  // Observed deaths by country by year
  array[N] int<lower=1,upper=C> country;    // grouping variable, encoding which country the observation comes from
  array[N] real<lower=0> population_size;    // Number of people by country (repeated avg per year as well, can change to be year specific pop)
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
  
  for (c in 1:C){
    beta_cases[c] = beta_cases_std[c] * b_sigma_cases; // non-centered param for cases by country
    beta_deaths[c] = beta_deaths_std[c] * b_sigma_deaths; // non-centered param for deaths by country
  }
  
  for (n in 1:N){ 
    
  infections[n] = infection_prob[n] * population_size[n]; // true infections taking into account country size   
    
  p[n] = inv_logit(alpha_cases + beta_cases[country[n]]);  // p = logit^-1(alpha + beta)
  mu_cases[n] = infections[n]*p[n]; 
  
  true_deaths[n] = infections[n] * true_ifr[n];    // relationship between infections and true deaths 
  
  q[n] = inv_logit(alpha_deaths + beta_deaths[country[n]]); // q = logit^-1(alpha + beta)
  mu_deaths[n] = true_deaths[n]*q[n]; 
  }
}

model {
  
  // Prior for infections probability and true deaths
  infection_prob ~ beta(0.5,100); // uninformative prior (seasonal flu infection probability)
  infections ~ gamma(324,0.003);    // from SitRep Modelling 02/03/2020 Word doc - John Hopkins Dashboard Y:\Sitrep Modelling\2020\Sit Rep Word
  true_deaths ~ gamma(3600,1.20); // from SitRep Modelling 02/03/2020 Word doc - John Hopkins Dashboard Y:\Sitrep Modelling\2020\Sit Rep Word
  true_ifr ~  beta(0.96,47); // prior as if we were back at Feb 2020
                             // https://assets.publishing.service.gov.uk/media/5eccfb2fe90e0754d96bb6fd/18-spi-m-o-consensus-statement-17022020.pdf
                             // https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30244-9/fulltext
                            // https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30245-0/fulltext
                            // solved as https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
  
  p ~ beta(2,2); // weakly informative Beta prior
  q ~ beta(2,2); // weakly informative Beta prior 
  
  // Priors for case and death parameters
  alpha_cases ~ normal(0,100);     // uninformative prior
  alpha_deaths ~ normal(0,100);    // uninformative prior
  
  // Priors for 
  mu_cases ~ gamma(0.91,4.27e-5);  // solving for alpha*theta = mean and alpha*(theta)^2 = var = 1000^2
  mu_deaths ~ gamma(1.10,1.50e-3); // solving for alpha*theta = mean and alpha*(theta)^2 = var = 700^2
  
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

  // Likelihood for cases and deaths using Negative Binomial distribution
  for (n in 1:N) {
    cases[n] ~ neg_binomial_2(mu_cases[n], phi_cases);  // Negative Binomial for cases
    deaths[n] ~ neg_binomial_2(mu_deaths[n], phi_deaths); // Negative Binomial for deaths
  }
}

generated quantities{
 vector[N] log_lik;
 real lprior;
 array[N] int cases_rep;
 array[N] int deaths_rep;
 
 // Posterior predictive check
 // Reference: https://avehtari.github.io/BDA_R_demos/demos_rstan/ppc/poisson-ppc.html
 
 for (n in 1:N){
   cases_rep[n] = neg_binomial_2_rng(mu_cases[n], phi_cases);
   deaths_rep[n] = neg_binomial_2_rng(mu_deaths[n], phi_deaths);
   log_lik[n] = beta_lpdf(true_ifr[n] | 0.96, 47);
 }
  
  lprior =  beta_lpdf(true_ifr | 0.96, 47);
 
 // Priorsense sensitivity checks 
 //for (n in 1:N) log_lik[n] = beta_lpdf(true_ifr[n] | 0.96, 47);
 //lprior =  beta_lpdf(true_ifr | 0.96, 47);
}


