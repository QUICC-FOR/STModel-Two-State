if [ -z "$SPECIES" ]; then
    echo "Error: species not set for model selection"
    exit 1
fi

DIR=~/STModel-Two-State

declare -a TempVars=(annual_mean_temp mean_diurnal_range mean_temp_wettest_quarter)
declare -a PrecipVars=(tot_annual_pp pp_seasonality pp_warmest_quarter)
declare -a AllVars=( ${TempVars[@]} ${PrecipVars[@]} )
declare -a DesignStrings=(111111 111110 111100 110111 110110 110100 100111 100110)
declare -a SingleVar=(110000 111000)

RUNNING=0
DONE=0
MAXCORES=30

# start with an intercept-only model
TEMP=NA
PRECIP=NA
COL=000000
cd $DIR; Rscript ./scr/3b_fit_stm.r $SPECIES $TEMP $PRECIP $COL $COL &
RUNNING=$((RUNNING + 1)) 

# loop through all single variable models
for TEMP in ${AllVars[@]}
do
    for COL in ${SingleVar[@]}
    do
        cd $DIR; Rscript ./scr/3b_fit_stm.r $SPECIES $TEMP $PRECIP $COL $COL &
        RUNNING=$((RUNNING + 1))
    done
done

for TEMP in ${TempVars[@]}
do
    for PRECIP in ${PrecipVars[@]}
    do
        for COL in ${DesignStrings[@]}
        do
            cd $DIR; Rscript ./scr/3b_fit_stm.r $SPECIES $TEMP $PRECIP $COL $COL &
            RUNNING=$((RUNNING + 1))
            if [ "$RUNNING" -ge "$MAXCORES" ] ; then
                wait
                DONE=$((DONE + RUNNING))
                RUNNING=0
            fi
        done
    done
done

wait