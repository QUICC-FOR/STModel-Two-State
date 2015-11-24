#!/usr/bin/env Rscript

# depends:
#    res/eval/SPNAME_MOD__rdeAreas.rds

# produces:
#    res/si_tab_eval.tex
library(xtable)
library(foreach)

spList = readRDS('dat/speciesList.rds')
spInfo = read.csv('dat/speciesInfo.csv')
get.info = function(sp)  spInfo[spInfo$spName == sp,]

# get the datasets
evalDat = foreach(spName = spList, .final=function(x) {names(x) = spList; x}) %do% {
	readRDS(paste0('res/eval/', spName, '_stm_eval_posterior.rds'))
}


spTab = data.frame(
	"Scientific name"=foreach(sp = spList, .combine=c) %do% paste(get.info(sp)$genus, get.info(sp)$species),
	"ROC" = foreach(sp = spList, .combine=c) %do% median(evalDat[[sp]][,'roc']),
	"ROC Lower" = foreach(sp = spList, .combine=c) %do% quantile(evalDat[[sp]][,'roc'], 0.05),
	"ROC Upper" = foreach(sp = spList, .combine=c) %do% quantile(evalDat[[sp]][,'roc'], 0.95),
	check.names=F
)

spXTab = xtable(spTab, label="tab:eval", align='llccc',
	caption="Evaluation statistics for predictions of species ranges. Values given are the median area under the receiver operating curve (ROC) with 90\\% credible intervals")
print(spXTab, file="res/table/si_tab_s2_eval.tex", caption.placement="top", booktabs=TRUE, include.rownames=FALSE)
