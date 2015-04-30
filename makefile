# set the species to be processed manually here
#SPECIES=28731-ACE-SAC
#SPECIES=18032-ABI-BAL
SPECIES=28728-ACE-RUB

all: results/$(SPECIES)_mcmc_out.txt

## step 1
s1: dat/$(SPECIES)_processed.rdata
dat/$(SPECIES)_processed.rdata: scr/1_prep_data.r dat/transition_twostate_$(SPECIES).rdata
	./scr/1_prep_data.r -s "$(SPECIES)"
	
dat/$(SPECIES)_climGrid_scaled.rds: dat/$(SPECIES)_processed.rdata

## step 2 - this must be done MANUALLY for each species in an interactive session
#dat/$(SPECIES)_sdmData.rdata: scr/2_interactive_select_sdm_vars.r dat/$(SPECIES)_processed.rdata
#	./scr/2_interactive_select_sdm_vars.r


## step 3
s3: results/$(SPECIES)_sdm_models.rdata
results/$(SPECIES)_sdm_models.rdata: scr/3_fit_sdm.r dat/$(SPECIES)_sdmData.rdata
	./scr/3_fit_sdm.r -s "$(SPECIES)" -p 2
	
## step 4
s4: dat/$(SPECIES)_climGrid_projected.rds
dat/$(SPECIES)_climGrid_projected.rds: scr/4_project_sdms.r \
dat/$(SPECIES)_climGrid_scaled.rds results/$(SPECIES)_sdm_models.rdata \
dat/$(SPECIES)_processed.rdata results/$(SPECIES)_sdm_models.rdata
	./scr/4_project_sdms.r -s "$(SPECIES)"

dat/$(SPECIES)_transitions_projected.rds: dat/$(SPECIES)_climGrid_projected.rds

## step 5
s5: img/$(SPECIES)_sdm_maps.pdf
img/$(SPECIES)_sdm_maps.pdf: scr/5_plot_sdm.r results/$(SPECIES)_sdm_models.rdata \
dat/$(SPECIES)_sdmData.rdata dat/$(SPECIES)_climGrid_projected.rds
	./scr/5_plot_sdm.r -s "$(SPECIES)"

img/$(SPECIES)_response_curves.pdf: img/$(SPECIES)_sdm_maps.pdf

# step 6
anneal: results/$(SPECIES)_anneal_parameters_0.1.rds
results/$(SPECIES)_anneal_parameters_0.1.rds: scr/6_fit_stm.r \
dat/$(SPECIES)_transitions_projected.rds
	./scr/6_fit_stm.r -s "$(SPECIES)" -g -f 0.1


# step 7 - prep mcmc
s7: dat/mcmc_trans_$(SPECIES).txt
dat/mcmc_trans_$(SPECIES).txt: scr/7_prep_mcmc_dat.r \
dat/$(SPECIES)_transitions_projected.rds results/$(SPECIES)_anneal_parameters_0.1.rds
	./scr/7_prep_mcmc_dat.r -s "$(SPECIES)" -f 0.25

dat/mcmc_pars_$(SPECIES).txt: dat/mcmc_trans_$(SPECIES).txt


#./scr/7_prep_mcmc_dat.r -s "18032-ABI-BAL" -f 0.25
#./scr/7_prep_mcmc_dat.r -s "18032-ABI-BAL" -f 0.001

# step 8 - run mcmc
mcmc: results/$(SPECIES)_mcmc_out.txt
results/$(SPECIES)_mcmc_out.txt: stm2_mcmc dat/mcmc_trans_$(SPECIES).txt \
dat/mcmc_pars_$(SPECIES).txt
	./stm2_mcmc -v 2 -n 25 -i 10000 -b 5000 -c 10 -p dat/mcmc_pars_$(SPECIES).txt \
	-t dat/mcmc_trans_$(SPECIES).txt -o results/$(SPECIES)/

#./stm2_mcmc -v 2 -n 25 -i 10000 -b 5000 -c 10 -p dat/mcmc_pars_18032-ABI-BAL.txt -t dat/mcmc_trans_18032-ABI-BAL.txt -o results/18032-ABI-BAL/
#./stm2_mcmc -v 2 -i 500 -c 10 -p dat/mcmc_pars_18032-ABI-BAL.txt -t dat/mcmc_trans_18032-ABI-BAL.txt -o results/18032-ABI-BAL