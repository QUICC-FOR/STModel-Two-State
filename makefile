# set the species to be processed manually here

# MCMC to-do

# MCMC done
#SPECIES=19462-FAG-GRA
#SPECIES=183295-PIC-GLA
#SPECIES=195773-POP-TRE

# posterior processing done
#SPECIES=19049-ULM-AME
#SPECIES=19290-QUE-ALB
#SPECIES=183319-PIN-BAN
#SPECIES=183302-PIC-MAR
#SPECIES=19481-BET-ALL
#SPECIES=19489-BET-PAP
#SPECIES=28731-ACE-SAC
#SPECIES=28731-ACE-SAC2
#SPECIES=18032-ABI-BAL
#SPECIES=28728-ACE-RUB

# figures done
#SPECIES=28731-ACE-SAC3 -- running
#SPECIES=32931-FRA-AME -- running


#all: results/$(SPECIES)/posterior.csv

# step 1: data preparation
s1: dat/$(SPECIES)/$(SPECIES)_processed.rdata
dat/$(SPECIES)/$(SPECIES)_processed.rdata: scr/1_prep_data.r \
			dat/$(SPECIES)/transition_twostate_$(SPECIES).rdata dat/SDMClimate_grid.csv
	./scr/1_prep_data.r -s $(SPECIES)

dat/$(SPECIES)/$(SPECIES)_climGrid_scaled.rds: dat/$(SPECIES)/$(SPECIES)_processed.rdata



# step 2: sdm
s2: results/$(SPECIES)/$(SPECIES)_sdm_models.rdata
results/$(SPECIES)/$(SPECIES)_sdm_models.rdata: scr/2_fit_sdm.r \
			dat/$(SPECIES)/$(SPECIES)_processed.rdata
	./scr/2_fit_sdm.r -s $(SPECIES) -p 3


# step 3 - project the SDMs
sdm: dat/$(SPECIES)/$(SPECIES)_transitions_projected.rds
s3: dat/$(SPECIES)/$(SPECIES)_transitions_projected.rds
dat/$(SPECIES)/$(SPECIES)_transitions_projected.rds: scr/3_project_sdm.r \
			dat/$(SPECIES)/$(SPECIES)_processed.rdata \
			dat/$(SPECIES)/$(SPECIES)_climGrid_scaled.rds \
			results/$(SPECIES)/$(SPECIES)_sdm_models.rdata
	./scr/3_project_sdm.r -s $(SPECIES)

dat/$(SPECIES)/$(SPECIES)_climGrid_projected.rds: dat/$(SPECIES)/$(SPECIES)_transitions_projected.rds
img/$(SPECIES)/$(SPECIES)_sdm_maps.pdf: dat/$(SPECIES)/$(SPECIES)_transitions_projected.rds
img/$(SPECIES)/$(SPECIES)_response_curves.pdf: dat/$(SPECIES)/$(SPECIES)_transitions_projected.rds


## step 4 - estimate some params with simulated annealing
anneal: results/$(SPECIES)/$(SPECIES)_anneal_parameters.rds
results/$(SPECIES)/$(SPECIES)_anneal_parameters.rds: scr/4_fit_stm.r \
			dat/$(SPECIES)/$(SPECIES)_transitions_projected.rds
	# for glm
#	./scr/4_fit_stm.r -s $(SPECIES)
	# for gam
#	./scr/4_fit_stm.r -s $(SPECIES) -m
	# for random forest
	./scr/4_fit_stm.r -s $(SPECIES) -r


## step 5 - make MCMC input files
mcmc: run2/$(SPECIES)/trans.txt run2/$(SPECIES)/inits.txt
run2/$(SPECIES)/trans.txt: scr/5_prep_mcmc.r \
			dat/$(SPECIES)/$(SPECIES)_transitions_projected.rds \
			results/$(SPECIES)/$(SPECIES)_anneal_parameters.rds
	# for glm
#	./scr/5_prep_mcmc.r -s $(SPECIES)
	# for gam
#	./scr/5_prep_mcmc.r -s $(SPECIES) -m
	# for random forest
	./scr/5_prep_mcmc.r -s $(SPECIES) -r

run2/$(SPECIES)/inits.txt: run2/$(SPECIES)/trans.txt


## run the mcmc
## we don't actually do it using make, since we run it by hand using qsub scripts
## this is just here to give an idea how it should be lauched
## the program logs extensively to std err, so it might be a good idea to redirect using 
## 2>log.txt at the end of the command
##
## ./stm2_mcmc -v 2 -n 25 -i 10000 -b 5000 -c 20 -p run/$(SPECIES)/inits.txt -t \
## 				run/$(SPECIES)/trans.txt -o ./ -l 5
##
## for resuming an interrupted session - note that after -i the number gives the number 
## of ADDITIONAL samples to take
## ./stm2_mcmc -r run/$(SPECIES)/resumeData.txt -t -t run/$(SPECIES)/trans.txt -i 10000

## step 6 - process the posterior distribution
s6: results/$(SPECIES)/$(SPECIES)_posterior.rds
results/$(SPECIES)/$(SPECIES)_rangeRaster.rds: results/$(SPECIES)/$(SPECIES)_posterior.rds
results/$(SPECIES)/$(SPECIES)_posteriorSummary.rds: results/$(SPECIES)/$(SPECIES)_posterior.rds
results/$(SPECIES)/$(SPECIES)_posterior.rds: scr/6_process_posterior.r \
			results/$(SPECIES)/posterior.csv \
			dat/$(SPECIES)/$(SPECIES)_climGrid_projected.rds
	./scr/6_process_posterior.r -s $(SPECIES) -r 4
	
	
## step 7 - make figures
figures: img/$(SPECIES)/$(SPECIES)_posterior_maps.pdf
img/$(SPECIES)/$(SPECIES)_posterior_maps.pdf: scr/7_make_figs.r \
			results/$(SPECIES)/$(SPECIES)_rangeRaster.rds \
			results/$(SPECIES)/$(SPECIES)_posteriorSummary.rds
	./scr/7_make_figs.r -s $(SPECIES)

