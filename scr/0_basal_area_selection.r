trees = read.csv("/Volumes/Data/matt/Dropbox/work/projects_git/STModel-Data/out_files/treeData.csv")
plots = read.csv("/Volumes/Data/matt/Dropbox/work/projects_git/STModel-Data/out_files/plotInfoData.csv")
species = read.csv("/Volumes/Data/matt/Dropbox/work/projects_git/STModel-Data/out_files/speciesCode.csv")
species = species[,c('id_spe', 'genus', 'species', 'en_common_name')]

speciesList = readRDS("dat/speciesList.rds")


# roughly the transition zone
lat.lims = c(42,50)

trees.ll = merge(trees, plots)

trees.transition = trees.ll[which(trees.ll$lat <= max(lat.lims) & trees.ll$lat >= min(lat.lims)),]
## trees.ba.plot = aggregate(basal_area ~ plot_id + id_spe, trees.transition, FUN=sum)
trees.ba = aggregate(basal_area ~ id_spe, trees.transition, FUN=sum)
trees.ba$ba_prop = trees.ba$basal_area / sum(trees.ba$basal_area)
trees.ba = merge(trees.ba, species)
trees.ba = trees.ba[rev(order(trees.ba$ba_prop)),]
trees.ba$ba_cum = cumsum(trees.ba$ba_prop)
trees.ba$rank = 1:nrow(trees.ba)

trees.ba.selected = trees.ba[trees.ba$id_spe %in% speciesList,]
#the basal area of selected species in the transition zone
sum(trees.ba.selected$ba_prop)

# used to get the plots used for the stm; we just select one particular species
calibDat = readRDS(paste0("dat/stm_calib/", speciesList[1], "_stm_calib.rds"))
presDat = readRDS(paste0("dat/presence/", speciesList[1], "_presence.rds"))
trees.ll.stmPlots = trees.ll[trees.ll$plot_id %in% calibDat$plot,]
trees.ll.presPlots = trees.ll[trees.ll$plot_id %in% presDat$plot_id,]

trees.ba.stmPlots = aggregate(basal_area ~ id_spe, trees.ll.stmPlots, FUN=sum)
trees.ba.stmPlots$ba_prop = trees.ba.stmPlots$basal_area / sum(trees.ba.stmPlots$basal_area)
trees.ba.stmPlots = merge(trees.ba.stmPlots, species)
trees.ba.stmPlots = trees.ba.stmPlots[rev(order(trees.ba.stmPlots$ba_prop)),]
trees.ba.stmPlots$ba_cum = cumsum(trees.ba.stmPlots$ba_prop)

#the basal area of selected species in all stm plots
trees.ba.stmPlots.selected = trees.ba.stmPlots[trees.ba.stmPlots$id_spe %in% speciesList,]
sum(trees.ba.stmPlots.selected$ba_prop)




par(cex=0.8)
plot(trees.ba$rank, trees.ba$ba_cum, type='l', xlab="Species Rank Abundance", ylab="Cumulative basal area", lwd=1.5)
polygon(c(-10, 21, 21, -10), c(0.8907451, 0.8907451, -.1, -0.1), lwd=0.7)
text(0, 0.8907451+.005, "89% of species", pos=4, cex=0.7)

polygon(c(-10, 40, 40, -10), c(0.9860828, 0.9860828, -.1, -0.1), lwd=0.7)
text(0, 0.9860828+.005, "98.5% of species", pos=4, cex=0.7)


trees.ba[trees.ba$ba_cum < 0.9, c('genus', 'species')]