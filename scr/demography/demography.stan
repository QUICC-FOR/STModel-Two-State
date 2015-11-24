data
{
	int<lower=1> num_data_points;
	int<lower=1> num_species;
	int<lower=0> recruit[num_data_points];
	int<lower=0> N[num_data_points];
	int<lower=1,upper=num_species> species[num_data_points];
	int<lower=0,upper=1> type[num_data_points];
	int<lower=1,upper=15> interval[num_data_points];
}
parameters
{
	// hyperparameters
	real alpha_mu;
	real beta_mu;
	real<lower=0> alpha_sigma;
	real<lower=0> beta_sigma;
	
	// individual-level pars
	real alpha [num_species];
	vector[num_species] beta;
}
model
{
	real phi [num_data_points];
	real phi_annual[num_data_points];

	alpha_mu ~ normal(0,1000);
	alpha_sigma ~ gamma(2,0.1);
	beta_mu ~ normal(0,1000);
	beta_sigma ~ gamma(2,0.1);

	for(sp in 1:num_species)
	{
		beta[sp] ~ normal(beta_mu, beta_sigma);
		alpha[sp] ~ normal(alpha_mu, alpha_sigma);
	}
	
	for(i in 1:num_data_points)
	{
		phi[i] <- inv_logit(alpha[species[i]] + beta[species[i]] * type[i]);
		phi_annual[i] <- 1-((1 - phi[i])^interval[i]);
		recruit[i] ~ binomial(N[i], phi_annual[i]);
	}
}