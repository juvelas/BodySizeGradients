
library(mvMORPH)
library(RPANDA)
library(OUwie)
library(geiger)
library(phytools)
library(plyr)
library(ggplot2)
library(emdbook)
library(phangorn)


# clibrated tree
tree<-read.nexus("mcc_thinned_allruns.trees")
# # # fixing tree
# tree$edge.length <- tree$edge.length * 600
# tree<-nnls.tree(cophenetic(tree),tree,rooted=TRUE)
# is.ultrametric(tree)


# dataset
datos<-read.csv("Datos_SVL_SSD_377spp_full.csv")
regions <- read.csv("regiones_spp.csv")
datos.regions <- merge(datos, regions, "species")
datos <- datos.regions[,c(1,8,15)]
rownames(datos) <- datos$species


# match tree with data
TreeOnly <- setdiff(tree$tip.label,rownames(datos))
TreeOnly # Enter the name of the object we just created to see what's in it.
DataOnly <- setdiff(rownames(datos), tree$tip.label)
DataOnly # Enter to see what species are in the data set but not the tree

# tree pruned
pruned <- drop.tip(tree, tip=TreeOnly)

# eliminate species without phylo and data info
datos <- datos[-match(DataOnly,rownames(datos)),]

# select areas
reg <- datos$code
names(reg)<-row.names(datos)

# fit data to other models
svl<-datos$mean_SVL
names(svl) <- datos$species
sd<-sd(svl)

# fit data to OUWie
datos.ouwie <- as.data.frame(datos)
rownames(datos.ouwie) <- c()
datos.ouwie <- datos.ouwie[,c(1,3,2)]
datos.ouwie["mserr"]<-sd

# SIMMAP
simmap.reg <- make.simmap(pruned, reg, model="ER", nsim=100)
summary.simmap <-summary(simmap.reg)


# Model fitting
# Geographical using OUwie

# ou1.svl <- OUwie(simmap.reg[[1]], datos.ouwie, model="OU1", simmap.tree=TRUE, scaleHeight = 1)
# oum.svl <- OUwie(simmap.reg[[1]], datos.ouwie, model="OUM", simmap.tree=TRUE, scaleHeight = 1)
# oumv.svl <- OUwie(simmap.reg[[1]], datos.ouwie, model="OUMV", simmap.tree=TRUE, scaleHeight = 1)
# ouma.svl <- OUwie(simmap.reg[[1]], datos.ouwie, model="OUMA", simmap.tree=TRUE, scaleHeight = 1)
# oumva.svl <- OUwie(simmap.reg[[1]], datos.ouwie, model="OUMVA", simmap.tree=TRUE, scaleHeight = 1)

ou.ouwie <- OUwie(simmap.reg[[1]], datos.ouwie, model="OU1", simmap.tree=TRUE, mserr = "known")
oum.ouwie <- OUwie(simmap.reg[[1]], datos.ouwie, model="OUM", simmap.tree=TRUE, mserr = "known")
oumv.ouwie <- OUwie(simmap.reg[[1]], datos.ouwie, model="OUMV", simmap.tree=TRUE, mserr = "known")
ouma.ouwie <- OUwie(simmap.reg[[1]], datos.ouwie, model="OUMA", simmap.tree=TRUE, mserr = "known")
oumva.ouwie <- OUwie(simmap.reg[[1]], datos.ouwie, model="OUMVA", simmap.tree=TRUE, mserr = "known")

# Geographical using mvMORPH
oum.mvmorph <- mvOU(simmap.reg[[1]], svl, model="OUM", error = sd)
ou.mvmorph <- mvOU(simmap.reg[[1]], svl, model="OU1", error = sd)

# Environmental
temp.data<-read.csv("temperature_age_Morlon.csv", row.names=1)
smooth.spline(temp.data[,1], temp.data[,2])$df
env.exp.1<-fit_t_env(simmap.reg[[1]], svl, env_data=temp.data, model="EnvExp", scale=TRUE)
env.exp.2<-fit_t_env(simmap.reg[[1]], svl, env_data=temp.data, model="EnvExp", df=10, scale=TRUE)
env.exp.3<-fit_t_env(simmap.reg[[1]], svl, env_data=temp.data, model="EnvExp", df=50, scale=TRUE)
env.lin.4<-fit_t_env(simmap.reg[[1]], svl, env_data=temp.data, model="EnvLin", scale=TRUE)
env.lin.5<-fit_t_env(simmap.reg[[1]], svl, env_data=temp.data, model="EnvLin", df=10, scale=TRUE)
env.lin.6<-fit_t_env(simmap.reg[[1]], svl, env_data=temp.data, model="EnvLin", df=50, scale=TRUE)


# Paleo-inspired and other models
etm<-52 # Eocene Thermal Maximum ~52 ma
eogm<-33 # Eocene-Oligocene Glacial Maximum ~ 33 ma
mmco<-14 # Middle Miocene Climatic Optimum ~ 14 ma

source("fitContinuous_paleoModels.R") # Graham Slater (MEE 2012 paper)

bm.fit <- fitContinuous_paleo(pruned, svl, model="BM", meserr = sd)
white.fit <- fitContinuous_paleo(pruned, svl, model="white", meserr = sd)
trend.fit <- fitContinuous_paleo(pruned, svl, model="trend", meserr = sd)
kappa.fit <- fitContinuous_paleo(pruned, svl, model="kappa", meserr = sd)
delta.fit <- fitContinuous_paleo(pruned, svl, model="delta", meserr = sd)
eb.fit <- fitContinuous_paleo(pruned, svl, model="EB", meserr = sd)

shift.fit.ETM <- fitContinuous_paleo(pruned, svl, model="timeshift", shift.time=etm, meserr=sd)
shift.fit.EOGM <- fitContinuous_paleo(pruned, svl, model="timeshift", shift.time=eogm, meserr=sd)
shift.fit.MMCO <- fitContinuous_paleo(pruned, svl, model="timeshift", shift.time=mmco, meserr=sd)
release.fit.ETM <- fitContinuous_paleo(pruned, svl, model="release",bounds=list(alpha=c(0,10)),shift.time=etm, meserr=sd)
release.fit.EOGM <- fitContinuous_paleo(pruned, svl, model="release",bounds=list(alpha=c(0,10)),shift.time=eogm, meserr=sd)
release.fit.MMCO <- fitContinuous_paleo(pruned, svl, model="release",bounds=list(alpha=c(0,10)),shift.time=mmco, meserr=sd)
releaseradiate.fit.ETM  <- fitContinuous_paleo(pruned, svl, model="releaseradiate", bounds=list(alpha=c(0,10)),shift.time=etm, meserr=sd)
releaseradiate.fit.EOGM <- fitContinuous_paleo(pruned, svl, model="releaseradiate", bounds=list(alpha=c(0,10)),shift.time=eogm, meserr=sd)
releaseradiate.fit.MMCO <- fitContinuous_paleo(pruned, svl, model="releaseradiate", bounds=list(alpha=c(0,10)),shift.time=mmco, meserr=sd)

save.image("FittingModels_CalibratedTree_Ultrametric.RData")

# Model comparisons
aicc.all <- c(ou.ouwie$AICc, oum.ouwie$AICc, oumv.ouwie$AICc, ouma.ouwie$AICc, oumva.ouwie$AICc,
                 oum.mvmorph$AICc, ou.mvmorph$AICc,
                 env.exp.1$aicc, env.exp.2$aicc, env.exp.3$aicc, env.lin.4$aicc, env.lin.5$aicc, env.lin.6$aicc,
                 bm.fit$Trait1$aicc, white.fit$Trait1$aicc, trend.fit$Trait1$aicc, kappa.fit$Trait1$aicc, delta.fit$Trait1$aicc, eb.fit$Trait1$aicc,
                 shift.fit.ETM$Trait1$aicc, shift.fit.EOGM$Trait1$aicc, shift.fit.MMCO$Trait1$aicc,
                 release.fit.ETM$Trait1$aicc, release.fit.EOGM$Trait1$aicc, release.fit.MMCO$Trait1$aicc,
                 releaseradiate.fit.ETM$Trait1$aicc, releaseradiate.fit.EOGM$Trait1$aicc, releaseradiate.fit.MMCO$Trait1$aicc)

lik.all <- c(ou.ouwie$loglik, oum.ouwie$loglik, oumv.ouwie$loglik, ouma.ouwie$loglik, oumva.ouwie$loglik,
                 oum.mvmorph$LogLik, ou.mvmorph$LogLik,
                 env.exp.1$LH, env.exp.2$LH, env.exp.3$LH, env.lin.4$LH, env.lin.5$LH, env.lin.6$LH,
                 bm.fit$Trait1$lnl, white.fit$Trait1$lnl, trend.fit$Trait1$lnl, kappa.fit$Trait1$lnl, delta.fit$Trait1$lnl, eb.fit$Trait1$lnl,
                 shift.fit.ETM$Trait1$lnl, shift.fit.EOGM$Trait1$lnl, shift.fit.MMCO$Trait1$lnl,
                 release.fit.ETM$Trait1$lnl, release.fit.EOGM$Trait1$lnl, release.fit.MMCO$Trait1$lnl,
                 releaseradiate.fit.ETM$Trait1$lnl, releaseradiate.fit.EOGM$Trait1$lnl, releaseradiate.fit.MMCO$Trait1$lnl)

k.all <- c(ou.ouwie$param.count, oum.ouwie$param.count, oumv.ouwie$param.count, ouma.ouwie$param.count, oumva.ouwie$param.count,
           oum.mvmorph$param$nparam, ou.mvmorph$param$nparam,
           env.exp.1$free.parameters, env.exp.2$free.parameters, env.exp.3$free.parameters, env.lin.4$free.parameters,
           env.lin.5$free.parameters, env.lin.6$free.parameters,
           bm.fit$Trait1$k, white.fit$Trait1$k, trend.fit$Trait1$k, kappa.fit$Trait1$k, delta.fit$Trait1$k, eb.fit$Trait1$k,
           shift.fit.ETM$Trait1$k, shift.fit.EOGM$Trait1$k, shift.fit.MMCO$Trait1$k,
           release.fit.ETM$Trait1$k, release.fit.EOGM$Trait1$k, release.fit.MMCO$Trait1$k,
           releaseradiate.fit.ETM$Trait1$k, releaseradiate.fit.EOGM$Trait1$k, releaseradiate.fit.MMCO$Trait1$k)

models <-cbind(lik.all, k.all, aicc.all)
rownames(models)<-c("OU", "OUM", "OUMV", "OUMA", "OUMVA", "OUM_mvMorph", "OU_mvMorph",
                     "Env-Depend_Exp1", "Env-Depend_Exp2", "Env-Depende_Exp3", "Env-Depend_Lin1", "Env-Depend_Lin2", "Env-Depend_Lin3",
                      "BM single-rate", "White-noise", "Trend", "Kappa", "Delta", "Early-burst", "Shift ETM", "Shift EOGM", "Shift MMCO",
                      "Release ETM", "Release EOGM", "Release MMCO", "Release & Radiate ETM", "Release & Radiate EOGM",
                      "Release & Radiate MMCO")

colnames(models)<-c("Lik", "parameters", "AICc")
models <- as.data.frame(models)
models["AIC_New"] <- -2*(models$Lik) + 2*(models$parameters)
weights <- aicw(models$AICc)
models.weights <- cbind(models, weights)
write.csv(models.weights, file="AkaikeWeights_AllModels_CalibratedTree.csv")



###########################################################
###########################################################


# full tree (a subset of species)

# clibrated tree
tree<-read.nexus("AnolisFullyResolvedCloseToConsensus.tre")
# # # fixing tree
tree$edge.length <- tree$edge.length * 600
tree<-nnls.tree(cophenetic(tree),tree,rooted=TRUE)
is.ultrametric(tree)


# dataset
datos<-read.csv("Datos_SVL_SSD_377spp_full.csv")
regions <- read.csv("regiones_spp.csv")
datos.regions <- merge(datos, regions, "species")
datos <- datos.regions[,c(1,8,15)]
rownames(datos) <- datos$species


# match tree with data
TreeOnly <- setdiff(tree$tip.label,rownames(datos))
TreeOnly # Enter the name of the object we just created to see what's in it.
DataOnly <- setdiff(rownames(datos), tree$tip.label)
DataOnly # Enter to see what species are in the data set but not the tree

# tree pruned
pruned <- drop.tip(tree, tip=TreeOnly)

# select areas
reg <- datos$code
names(reg)<-row.names(datos)

# fit data to other models
svl<-datos$mean_SVL
names(svl) <- datos$species
sd<-sd(svl)

# fit data to OUWie
datos.ouwie <- as.data.frame(datos)
rownames(datos.ouwie) <- c()
datos.ouwie <- datos.ouwie[,c(1,3,2)]
datos.ouwie["mserr"]<-sd

# SIMMAP
simmap.reg <- make.simmap(pruned, reg, model="ER", nsim=100)
summary.simmap <-summary(simmap.reg)

# Model fitting
# Geographical using OUwie

# ou1.svl <- OUwie(simmap.reg[[1]], datos.ouwie, model="OU1", simmap.tree=TRUE, scaleHeight = 1)
# oum.svl <- OUwie(simmap.reg[[1]], datos.ouwie, model="OUM", simmap.tree=TRUE, scaleHeight = 1)
# oumv.svl <- OUwie(simmap.reg[[1]], datos.ouwie, model="OUMV", simmap.tree=TRUE, scaleHeight = 1)
# ouma.svl <- OUwie(simmap.reg[[1]], datos.ouwie, model="OUMA", simmap.tree=TRUE, scaleHeight = 1)
# oumva.svl <- OUwie(simmap.reg[[1]], datos.ouwie, model="OUMVA", simmap.tree=TRUE, scaleHeight = 1)

ou.ouwie <- OUwie(simmap.reg[[1]], datos.ouwie, model="OU1", simmap.tree=TRUE, mserr = "known")
oum.ouwie <- OUwie(simmap.reg[[1]], datos.ouwie, model="OUM", simmap.tree=TRUE, mserr = "known")
oumv.ouwie <- OUwie(simmap.reg[[1]], datos.ouwie, model="OUMV", simmap.tree=TRUE, mserr = "known")
ouma.ouwie <- OUwie(simmap.reg[[1]], datos.ouwie, model="OUMA", simmap.tree=TRUE, mserr = "known")
oumva.ouwie <- OUwie(simmap.reg[[1]], datos.ouwie, model="OUMVA", simmap.tree=TRUE, mserr = "known")

# Geographical using mvMORPH
oum.mvmorph <- mvOU(simmap.reg[[1]], svl, model="OUM", error = sd)
ou.mvmorph <- mvOU(simmap.reg[[1]], svl, model="OU1", error = sd)

# Environmental
temp.data<-read.csv("temperature_age_Morlon.csv", row.names=1)
smooth.spline(temp.data[,1], temp.data[,2])$df
env.exp.1<-fit_t_env(simmap.reg[[1]], svl, env_data=temp.data, model="EnvExp", scale=TRUE)
env.exp.2<-fit_t_env(simmap.reg[[1]], svl, env_data=temp.data, model="EnvExp", df=10, scale=TRUE)
env.exp.3<-fit_t_env(simmap.reg[[1]], svl, env_data=temp.data, model="EnvExp", df=50, scale=TRUE)
env.lin.4<-fit_t_env(simmap.reg[[1]], svl, env_data=temp.data, model="EnvLin", scale=TRUE)
env.lin.5<-fit_t_env(simmap.reg[[1]], svl, env_data=temp.data, model="EnvLin", df=10, scale=TRUE)
env.lin.6<-fit_t_env(simmap.reg[[1]], svl, env_data=temp.data, model="EnvLin", df=50, scale=TRUE)


# Paleo-inspired and other models
etm<-52 # Eocene Thermal Maximum ~52 ma
eogm<-33 # Eocene-Oligocene Glacial Maximum ~ 33 ma
mmco<-14 # Middle Miocene Climatic Optimum ~ 14 ma

source("fitContinuous_paleoModels.R")

bm.fit <- fitContinuous_paleo(pruned, svl, model="BM", meserr = sd)
white.fit <- fitContinuous_paleo(pruned, svl, model="white", meserr = sd)
trend.fit <- fitContinuous_paleo(pruned, svl, model="trend", meserr = sd)
kappa.fit <- fitContinuous_paleo(pruned, svl, model="kappa", meserr = sd)
delta.fit <- fitContinuous_paleo(pruned, svl, model="delta", meserr = sd)
eb.fit <- fitContinuous_paleo(pruned, svl, model="EB", meserr = sd)

shift.fit.ETM <- fitContinuous_paleo(pruned, svl, model="timeshift", shift.time=etm, meserr=sd)
shift.fit.EOGM <- fitContinuous_paleo(pruned, svl, model="timeshift", shift.time=eogm, meserr=sd)
shift.fit.MMCO <- fitContinuous_paleo(pruned, svl, model="timeshift", shift.time=mmco, meserr=sd)
release.fit.ETM <- fitContinuous_paleo(pruned, svl, model="release",bounds=list(alpha=c(0,10)),shift.time=etm, meserr=sd)
release.fit.EOGM <- fitContinuous_paleo(pruned, svl, model="release",bounds=list(alpha=c(0,10)),shift.time=eogm, meserr=sd)
release.fit.MMCO <- fitContinuous_paleo(pruned, svl, model="release",bounds=list(alpha=c(0,10)),shift.time=mmco, meserr=sd)
releaseradiate.fit.ETM  <- fitContinuous_paleo(pruned, svl, model="releaseradiate", bounds=list(alpha=c(0,10)),shift.time=etm, meserr=sd)
releaseradiate.fit.EOGM <- fitContinuous_paleo(pruned, svl, model="releaseradiate", bounds=list(alpha=c(0,10)),shift.time=eogm, meserr=sd)
releaseradiate.fit.MMCO <- fitContinuous_paleo(pruned, svl, model="releaseradiate", bounds=list(alpha=c(0,10)),shift.time=mmco, meserr=sd)

save.image("FittingModels_FullTree_Ultrametric.RData")

# Model comparisons
aicc.all <- c(ou.ouwie$AICc, oum.ouwie$AICc, oumv.ouwie$AICc, ouma.ouwie$AICc, oumva.ouwie$AICc,
              oum.mvmorph$AICc, ou.mvmorph$AICc,
              env.exp.1$aicc, env.exp.2$aicc, env.exp.3$aicc, env.lin.4$aicc, env.lin.5$aicc, env.lin.6$aicc,
              bm.fit$Trait1$aicc, white.fit$Trait1$aicc, trend.fit$Trait1$aicc, kappa.fit$Trait1$aicc, delta.fit$Trait1$aicc, eb.fit$Trait1$aicc,
              shift.fit.ETM$Trait1$aicc, shift.fit.EOGM$Trait1$aicc, shift.fit.MMCO$Trait1$aicc,
              release.fit.ETM$Trait1$aicc, release.fit.EOGM$Trait1$aicc, release.fit.MMCO$Trait1$aicc,
              releaseradiate.fit.ETM$Trait1$aicc, releaseradiate.fit.EOGM$Trait1$aicc, releaseradiate.fit.MMCO$Trait1$aicc)

lik.all <- c(ou.ouwie$loglik, oum.ouwie$loglik, oumv.ouwie$loglik, ouma.ouwie$loglik, oumva.ouwie$loglik,
             oum.mvmorph$LogLik, ou.mvmorph$LogLik,
             env.exp.1$LH, env.exp.2$LH, env.exp.3$LH, env.lin.4$LH, env.lin.5$LH, env.lin.6$LH,
             bm.fit$Trait1$lnl, white.fit$Trait1$lnl, trend.fit$Trait1$lnl, kappa.fit$Trait1$lnl, delta.fit$Trait1$lnl, eb.fit$Trait1$lnl,
             shift.fit.ETM$Trait1$lnl, shift.fit.EOGM$Trait1$lnl, shift.fit.MMCO$Trait1$lnl,
             release.fit.ETM$Trait1$lnl, release.fit.EOGM$Trait1$lnl, release.fit.MMCO$Trait1$lnl,
             releaseradiate.fit.ETM$Trait1$lnl, releaseradiate.fit.EOGM$Trait1$lnl, releaseradiate.fit.MMCO$Trait1$lnl)

k.all <- c(ou.ouwie$param.count, oum.ouwie$param.count, oumv.ouwie$param.count, ouma.ouwie$param.count, oumva.ouwie$param.count,
           oum.mvmorph$param$nparam, ou.mvmorph$param$nparam,
           env.exp.1$free.parameters, env.exp.2$free.parameters, env.exp.3$free.parameters, env.lin.4$free.parameters,
           env.lin.5$free.parameters, env.lin.6$free.parameters,
           bm.fit$Trait1$k, white.fit$Trait1$k, trend.fit$Trait1$k, kappa.fit$Trait1$k, delta.fit$Trait1$k, eb.fit$Trait1$k,
           shift.fit.ETM$Trait1$k, shift.fit.EOGM$Trait1$k, shift.fit.MMCO$Trait1$k,
           release.fit.ETM$Trait1$k, release.fit.EOGM$Trait1$k, release.fit.MMCO$Trait1$k,
           releaseradiate.fit.ETM$Trait1$k, releaseradiate.fit.EOGM$Trait1$k, releaseradiate.fit.MMCO$Trait1$k)

models <-cbind(lik.all, k.all, aicc.all)
rownames(models)<-c("OU", "OUM", "OUMV", "OUMA", "OUMVA", "OUM_mvMorph", "OU_mvMorph",
                    "Env-Depend_Exp1", "Env-Depend_Exp2", "Env-Depende_Exp3", "Env-Depend_Lin1", "Env-Depend_Lin2", "Env-Depend_Lin3",
                    "BM single-rate", "White-noise", "Trend", "Kappa", "Delta", "Early-burst", "Shift ETM", "Shift EOGM", "Shift MMCO",
                    "Release ETM", "Release EOGM", "Release MMCO", "Release & Radiate ETM", "Release & Radiate EOGM",
                    "Release & Radiate MMCO")

colnames(models)<-c("Lik", "parameters", "AICc")
models <- as.data.frame(models)
models["AIC_New"] <- -2*(models$Lik) + 2*(models$parameters)
weights <- aicw(models$AICc)
models.weights <- cbind(models, weights)
write.csv(models.weights, file="AkaikeWeights_AllModels_FullTree.csv")
