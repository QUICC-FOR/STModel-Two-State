if [ -z "$SPECIES" ]; then
  echo "Error: species not set for model selection"
  exit 1
fi

DIR=~/STModel-Two-State

declare -a TempVars=(annual_mean_temp mean_diurnal_range mean_temp_wettest_quarter)
declare -a PrecipVars=(tot_annual_pp pp_seasonality pp_warmest_quarter)
declare -a AllVars=( ${TempVars[@]} ${PrecipVars[@]} )

# design strings are six characters, each a 1 or 0
# 1 indicates the parameter will be fit
# they alternate env1 then env2, in order of increasing power
# so 110000 fits a first order eqn for env1 and env2
# 111010 first third-order for env1 and first order for env2
# the intercept is always included and added automatically
declare -a DesignStrings=(110000 111000 111100 110100)
SingleVar=(101000 100000)

RUNNING=0
DONE=0
MAXCORES=40

# start with an intercept-only model - 1 model
TEMP=NA
PRECIP=NA
COL=000000
cd $DIR; Rscript ./scr/fit_stm_4.r $SPECIES $TEMP $PRECIP $COL $COL &
RUNNING=$((RUNNING + 1)) 

# loop through all single variable models
# 6 variables * 2 col * 2 ext models each = 24 models
for TEMP in ${AllVars[@]}
do
  for COL in ${SingleVar[@]}
  do
    for EXT in ${SingleVar[@]}
    do
      cd $DIR; Rscript ./scr/fit_stm_4.r $SPECIES $TEMP $PRECIP $COL $EXT &
      RUNNING=$((RUNNING + 1))
    done
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
        cd $DIR; Rscript ./scr/fit_stm_4.r $SPECIES $TEMP $PRECIP $COL $EXT &
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