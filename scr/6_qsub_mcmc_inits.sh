#!/bin/sh
# should be run from the species/ directory
for SPECIES in *; do
	echo "#!/bin/sh" >../scr/mcmc_scripts/$SPECIES-0.sh
	echo "#PBS -q default" >>../scr/mcmc_scripts/$SPECIES-0.sh
	echo "#PBS -l walltime=24:00:00" >>../scr/mcmc_scripts/$SPECIES-0.sh
	echo "#PBS -l nodes=1:ppn=40" >>../scr/mcmc_scripts/$SPECIES-0.sh
	echo "#PBS -r n" >>../scr/mcmc_scripts/$SPECIES-0.sh
	echo "#PBS -N ${SPECIES}-0" >>../scr/mcmc_scripts/$SPECIES-0.sh
	echo "SPECIES=${SPECIES}" >>../scr/mcmc_scripts/$SPECIES-0.sh
	cat ../scr/mcmc_init_template.sh >>../scr/mcmc_scripts/$SPECIES-0.sh
	qsub ../scr/mcmc_scripts/$SPECIES-0.sh
done



