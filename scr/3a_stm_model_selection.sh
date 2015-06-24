if [ -z "$SPECIES" ]; then
  echo "Error: species not set for model selection"
  exit 1
fi

DIR=~/STModel-Two-State

declare -a TempVars=(annual_mean_temp mean_diurnal_range mean_temp_wettest_quarter)
declare -a PrecipVars=(tot_annual_pp pp_seasonality pp_warmest_quarter)
declare -a AllVars=( ${TempVars[@]} ${PrecipVars[@]} )
declare -a DesignStrings=(110110 110100 100110)
SingleVar=(110000 100000)

RUNNING=0
DONE=0
MAXCORES=40

# start with an intercept-only model - 1 model
TEMP=NA
PRECIP=NA
COL=000000
cd $DIR; Rscript ./scr/3b_fit_stm.r $SPECIES $TEMP $PRECIP $COL $COL &
RUNNING=$((RUNNING + 1)) 

# loop through all single variable models
# 6 variables * 2 col * 2 ext models each = 24 models
for TEMP in ${AllVars[@]}
do
  for COL in ${SingleVar[@]}
  do
    for COL in ${SingleVar[@]}
    do
      cd $DIR; Rscript ./scr/3b_fit_stm.r $SPECIES $TEMP $PRECIP $COL $EXT &
      RUNNING=$((RUNNING + 1))
  done
done

# 3 temp * 3 precip * 3 col models * 3 ext models = 81 models
# with previous, we have 106 total models
for TEMP in ${TempVars[@]}
do
  for PRECIP in ${PrecipVars[@]}
  do
    for COL in ${DesignStrings[@]}
    do
      for EXT in ${DesignStrings[@]}
      do
        cd $DIR; Rscript ./scr/3b_fit_stm.r $SPECIES $TEMP $PRECIP $COL $EXT &
        RUNNING=$((RUNNING + 1))
        if [ "$RUNNING" -ge "$MAXCORES" ] ; then
          wait
          DONE=$((DONE + RUNNING))
          RUNNING=0
        fi
      done
    done
  done
done

wait