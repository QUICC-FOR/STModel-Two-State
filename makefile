## step 1
dat/28731-ACE-SAC_processed.rdata: scr/1_prep_data.r dat/transition_twostate_28731-ACE-SAC.rdata
	./scr/1_prep_data.r -s "28731-ACE-SAC"
	
dat/18032-ABI-BAL_processed.rdata: scr/1_prep_data.r dat/transition_twostate_18032-ABI-BAL.rdata
	./scr/1_prep_data.r -s "18032-ABI-BAL"

dat/28731-ACE-SAC_climGrid_scaled.rds: dat/28731-ACE-SAC_processed.rdata

dat/18032-ABI-BAL_climGrid_scaled.rds: dat/18032-ABI-BAL_processed.rdata

## step 2
#dat/28731-ACE-SAC_sdmData.rdata: scr/2_interactive_select_sdm_vars.r dat/28731-ACE-SAC_processed.rdata
#	./scr/2_interactive_select_sdm_vars.r

#dat/18032-ABI-BAL_sdmData.rdata: scr/2_interactive_select_sdm_vars.r dat/18032-ABI-BAL_processed.rdata
# 	./scr/2_interactive_select_sdm_vars.r


## step 3
results/28731-ACE-SAC_sdm_models.rdata: scr/3_fit_sdm.r dat/28731-ACE-SAC_sdmData.rdata
	./scr/3_fit_sdm.r -s "28731-ACE-SAC" -p 2
	
results/18032-ABI-BAL_sdm_models.rdata: scr/3_fit_sdm.r dat/18032-ABI-BAL_sdmData.rdata
	./scr/3_fit_sdm.r -s "18032-ABI-BAL" -p 2

## step 4
dat/28731-ACE-SAC_climGrid_projected.rds: scr/4_project_sdms.r \
dat/28731-ACE-SAC_climGrid_scaled.rds results/28731-ACE-SAC_sdm_models.rdata \
dat/28731-ACE-SAC_processed.rdata results/28731-ACE-SAC_sdm_models.rdata
	./scr/4_project_sdms.r -s "28731-ACE-SAC"

dat/28731-ACE-SAC_transitions_projected.rds: dat/28731-ACE-SAC_climGrid_projected.rds



#./scr/4_project_sdms.r -s "18032-ABI-BAL"


## step 5
img/28731-ACE-SAC_sdm_maps.pdf: scr/5_plot_sdm.r results/28731-ACE-SAC_sdm_models.rdata \
dat/28731-ACE-SAC_sdmData.rdata dat/28731-ACE-SAC_climGrid_projected.rds
	./scr/5_plot_sdm.r -s "28731-ACE-SAC"

img/28731-ACE-SAC_response_curves.pdf: img/28731-ACE-SAC_sdm_maps.pdf

# step 6
results/28731-ACE-SAC_anneal_parameters_0.1.rds: scr/6_fit_stm.r \
dat/28731-ACE-SAC_transitions_projected.rds
	./scr/6_fit_stm.r -s "28731-ACE-SAC" -g -f 0.1



# step 7
## requires
annealResults = readRDS(paste("results/", argList$species, "_anneal_parameters_", argList$fraction, ".rds", sep=""))
load(paste("dat/", argList$species, "_sdmData.rdata", sep=""))
load(paste("results/", argList$species, "_sdm_models.rdata", sep=""))
	climGrid = readRDS(paste("dat/", argList$species, "_climGrid_projected.rds", sep=''))


# produces
		pdf(paste("img/", argList$species, "_stm_response_curves.pdf", sep=""), w=width, h=height)
