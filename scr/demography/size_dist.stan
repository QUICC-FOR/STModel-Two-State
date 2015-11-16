data
{
	int<lower=1> num_data_points;
	int<lower=1> num_species;
	vector<lower=0>[num_data_points] dbh;
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
	
	// group pars
	vector[num_species] alpha;
	vector[num_species] beta;
	vector<lower=0>[num_species] sigma;
	
}
model
{
	alpha_mu ~ normal(0,1000);
	alpha_sigma ~ gamma(2,0.1);
	beta_mu ~ normal(0,1000);
	beta_sigma ~ gamma(2,0.1);

// 	for(sp in 1:num_species)
// 	{
// 		alpha[sp] ~ normal(alpha_mu, alpha_sigma);
// 		beta[sp] ~ normal(beta_mu, beta_sigma);
// 		sigma[sp] ~ gamma(2,0.1);
// 	}
	alpha ~ normal(alpha_mu, alpha_sigma);
	beta ~ normal(beta_mu, beta_sigma);
	sigma ~ gamma(2,0.1);
	
	for(i in 1:num_data_points)
	{
		dbh[i] ~ normal(alpha[species[i]] + beta[species[i]] * type[i], sigma[species[i]]);
	}
}