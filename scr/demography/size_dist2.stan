data
{
	int<lower=1> num_data_points;
	int<lower=1> num_species;
	int<lower=1> num_years;
	int<lower=1> num_plots;
	vector<lower=0>[num_data_points] dbh;
	int<lower=1,upper=num_species> species[num_data_points];
	int<lower=1,upper=num_years> year[num_data_points];
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
	real alpha [num_species, num_years, num_plots];
	vector[num_species] beta;
	real<lower=0> sigma [num_species, num_years, num_plots];
	
	// year-level pars
	real alphaYr_mu [num_years, num_plots];
	real<lower=0> alphaYr_sigma[num_years, num_plots];
	
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

	for(yr in 1:num_years)
	{
		for(pl in 1:num_plots)
		{
			alphaYr_mu[yr,pl] ~ normal(alphaPl_mu[pl], alphaPl_sigma[pl]);
			alphaYr_sigma[yr,pl] ~ gamma(2,0.1);
		}
	}

	for(sp in 1:num_species)
	{
		// slopes are ONLY nested within species
		beta[sp] ~ normal(beta_mu, beta_sigma);

		// intercept is nested within species within year within plot
		// variance is unpooled across species, year, and plot
		for(yr in 1:num_years)
		{
			for(pl in 1:num_plots)
			{
				alpha[sp,yr,pl] ~ normal(alphaYr_mu[yr,pl], alphaYr_sigma[yr,pl]);
				sigma[sp,yr,pl] ~ gamma(2,0.1);
			}
		}
	}
	
	for(i in 1:num_data_points)
	{
		dbh[i] ~ normal(alpha[species[i], year[i], plot[i]] + beta[species[i]] * type[i], sigma[species[i], year[i], plot[i]]);
	}
}