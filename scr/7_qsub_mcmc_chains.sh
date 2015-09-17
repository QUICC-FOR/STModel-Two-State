#!/bin/sh
# should be run from the species/ directory
nums=(1 2 3)
flags=("" "-a" "-g")

for FL in "${flags[@]}"; do
	for SPECIES in *; do
		for N in "${nums[@]}"; do
			echo "#!/bin/sh" >../scr/mcmc_scripts/"$SPECIES"_"$N"_"$FL".sh
			echo "#PBS -q default" >>../scr/mcmc_scripts/"$SPECIES"_"$N"_"$FL".sh
			echo "#PBS -l walltime=72:00:00" >>../scr/mcmc_scripts/"$SPECIES"_"$N"_"$FL".sh
			echo "#PBS -l nodes=1:ppn=40" >>../scr/mcmc_scripts/"$SPECIES"_"$N"_"$FL".sh
			echo "#PBS -r n" >>../scr/mcmc_scripts/"$SPECIES"_"$N"_"$FL".sh
			echo "#PBS -N ${SPECIES}_${N}_${FL}" >>../scr/mcmc_scripts/"$SPECIES"_"$N"_"$FL".sh
			echo "SPECIES=${SPECIES}" >>../scr/mcmc_scripts/"$SPECIES"_"$N"_"$FL".sh
			echo "N=${N}" >>../scr/mcmc_scripts/"$SPECIES"_"$N"_"$FL".sh
			echo 'FLAG=''"'$FL'"' >>../scr/mcmc_scripts/"$SPECIES"_"$N"_"$FL".sh
			cat ../scr/mcmc_chain_template.sh >>../scr/mcmc_scripts/"$SPECIES"_"$N"_"$FL".sh
			qsub ../scr/mcmc_scripts/"$SPECIES"_"$N"_"$FL".sh
		done
	done
done



