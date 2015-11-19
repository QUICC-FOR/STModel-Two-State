data
{
	int<lower=1> num_data_points;
	int<lower=1> num_species;
	int<lower=0> dead[num_data_points];
	int<lower=0> alive[num_data_points];
	int<lower=1,upper=num_species> species[num_data_points];
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
	real alpha [num_species];
	vector[num_species] beta;	
}
model
{
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
		dead[i] ~ binomial_logit(dead[i]+alive[i], alpha[species[i]] + beta[species[i]] * type[i]);
	}
}