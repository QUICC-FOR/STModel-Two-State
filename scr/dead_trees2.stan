data
{
	int<lower=1> num_data_points;
	int<lower=1> num_species;
	int<lower=1> num_plots;
	int<lower=0> dead[num_data_points];
	int<lower=0> alive[num_data_points];
	int<lower=1,upper=num_species> species[num_data_points];
	int<lower=1,upper=num_plots> plot[num_data_points];
	int<lower=0,upper=1> type[num_data_points];
}
parameters
{
	// hyperparameters
	real alpha_mu;
	real beta_mu;
	real<lower=0> alpha_sigma;
	real<lower=0> beta_sigma;
	
	// individual-level pars
	real alpha [num_species, num_plots];
	vector[num_species] beta;
		
	// plot-level
	real alphaPl_mu[num_plots];
	real<lower=0> alphaPl_sigma[num_plots];
	
}
model
{
	alpha_mu ~ normal(0,1000);
	alpha_sigma ~ gamma(2,0.1);
	beta_mu ~ normal(0,1000);
	beta_sigma ~ gamma(2,0.1);

	for(pl in 1:num_plots)
	{
		alphaPl_mu[pl] ~ normal(alpha_mu, alpha_sigma);
		alphaPl_sigma[pl] ~ gamma(2,0.1);
	}

	for(sp in 1:num_species)
	{
		// slopes are ONLY nested within species
		beta[sp] ~ normal(beta_mu, beta_sigma);

		// intercept is nested within species within  plot
		for(pl in 1:num_plots)
		{
			alpha[sp,pl] ~ normal(alphaPl_mu[pl], alphaPl_sigma[pl]);
		}
	}
	
	for(i in 1:num_data_points)
	{
		dead[i] ~ binomial_logit(dead[i]+alive[i], alpha[species[i], plot[i]] + beta[species[i]] * type[i]);
	}
}