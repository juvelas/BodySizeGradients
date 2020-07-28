library(raster)
library(sp)
library(maptools)
library(rgdal)
library(dismo)
library(XML)
library(maps)
library(letsR)
library(plyr)
library(geiger)
library(phytools)
library(vegan)
library(PVR)
library(picante)
library(rgeos)
library(spgwr)
library(classInt)
library(RColorBrewer)
library(rcompanion)
library(lmtest)
library(spdep)
library(spatialreg)


# Velasco etal's dataset points
pts <- as.data.frame(read.csv("ALL_RECORDS_COMBINED_31052017.csv"))
recs <- count(pts, "species")
pts <- merge(pts, recs, by="species")

# Points for PAM
xy<-pts[,c(4:5)]
species<-pts$species
pts.anole<-pts[,c(-2,-3)]

# Velasco etal's dataset convex hulls
ranges <-readOGR("./ranges/ANOLE_RANGE_MAPS_POE_etal_TAXONOMY_27032017.shp")
proj4string(ranges) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
projection(ranges)<-crs.geo

ranges <- ranges[,c(1,3)]
names(ranges)
ranges_spp <-rename(ranges, replace=c("Species"="sciname"))
names(ranges_spp)

# SDM's
models <-readOGR("./ranges/merge_of_merges.shp")
proj4string(models) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
projection(models)<-crs.geo


PAM_convex<-lets.presab(ranges, xmn=-115, xmx=-36, ymn=-23, ymx=36, resol=1,
                    remove.cells=TRUE, remove.sp=TRUE, show.matrix=FALSE,
                    crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"), cover=0, presence=NULL,
                    origin=NULL, seasonal=NULL, count=FALSE)

PAM_points<-lets.presab.points(xy, species, xmn=-115, xmx=-36, ymn=-23, ymx=36, resol=1,
                               remove.cells=TRUE, remove.sp=TRUE, show.matrix=F,
                               crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"), count=FALSE)

PAM_models<-lets.presab(models, xmn=-115, xmx=-36, ymn=-23, ymx=36, resol=1,
                        remove.cells=TRUE, remove.sp=TRUE, show.matrix=FALSE,
                        crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"), cover=0, presence=NULL,
                        origin=NULL, seasonal=NULL, count=FALSE)



plot(PAM_convex, main="Species richness")
plot(PAM_points, main="Species richness")

#adicionar los datos de tama?o corporal
body_size<-read.csv("./SVL_Datos/svl_anoles_SPoe_matched_maps.csv", row.names = 1)
# extraer mean Body size
max.svl<-body_size[,8]
names(max.svl)<-row.names(body_size)

svl_pam <-lets.maplizer(PAM_models, max.svl, PAM_models$Species_name, func=median, ras=T)
svl_pam.obs <- lets.addvar(x=PAM_models, y=svl_pam$Raster, fun=median)
svl_pam.obs <- as.data.frame(svl_pam.obs[,c(1,2,388)])
svl_pam.obs <- rename(svl_pam.obs, replace=c("layer_median"="median_log_svl"))


# randomize svl column 
rand.svl <- lapply(1:100, function(x) sample(body_size$Log10_meanSVL))
svl_pam_rand <-lapply(1:100, function(x) lets.maplizer(PAM_models, rand.svl[[x]], PAM_models$Species_name, func=median, ras=T))
svl_pam_rand2 <- lapply(1:100, function(x) lets.addvar(x=PAM_models, y=svl_pam_rand[[x]]$Raster, fun=median))
svl.y <- lapply(1:100, function(x) as.data.frame(svl_pam_rand2[[x]][,c(1,2,388)]))
svl.y <- lapply(1:100, function(x) plyr::rename(svl.y[[x]], replace=c("layer_median"="median_log_svl")))



# variables
# Heat balance
varclim <- stack(list.files(path="./Current_5km/", pattern='asc',full.names=T))
tpa_h1<-varclim[[1]] # variable for heat balance hypothesis (H1)
projection(tpa_h1)<-crs.geo

# Primary productivity 
aet<-raster("./Envir_vars/aet_yr1.bil", pattern='bil', full.names=T)
projection(aet)<-crs.geo

# Seasonality
seas.temp<-varclim[[4]]
seas.prec<-varclim[[3]]
seas.r<-stack(seas.temp, seas.prec)
projection(seas.r)<-crs.geo

# PAM quarter
#interpolar PAM (presence-absence matrix) con rasters
heatbal1 <-lets.addvar(x=PAM_models, y=tpa_h1, fun=median)
primpro <-lets.addvar(x=PAM_models, y=aet, fun=median)
seas <-lets.addvar(x=PAM_models, y=seas.r, fun=median)

# dataframes
heatbal.t <-as.data.frame(heatbal1[,c(1,2,388)])
heatbal.t<-rename(heatbal.t, replace=c("bio1_median"="anntemp"))
prim_prod<-as.data.frame(primpro[,c(1,2,388)])
prim_prod<-rename(prim_prod, replace=c("aet_yr1_median"="aet"))
seas.tp<-as.data.frame(seas[,c(1,2,388,389)])
seas.tp<-rename(seas.tp, replace=c("bio4_median"="season_temp", "bio15_median"="season_prec"))

# juntar todo
pam_vars<-cbind(heatbal.t, prim_prod, seas.tp)
pam_vars<-pam_vars[,c(3,6,9,10)]
pam_vars_rand <- lapply(1:100, function(x) cbind(svl.y[[x]], pam_vars))
pam_vars_rand[[1]]

pam_vars_svl <- cbind(svl_pam.obs, pam_vars)

# OLS randomizing svl value

heatbal.obs <- lm(median_log_svl  ~ anntemp, data=pam_vars_svl)
primprod.obs <- lm(median_log_svl  ~ aet, data=pam_vars_svl)
seas.obs <- lm(median_log_svl  ~ season_temp + season_prec, data=pam_vars_svl)

heatbal.rand <- lapply(1:100, function(x) lm(median_log_svl  ~ anntemp, data=pam_vars_rand[[x]]))
summary(heatbal.rand[[1]])
primprod.rand <- lapply(1:100, function(x) lm(median_log_svl  ~ aet, data=pam_vars_rand[[x]]))
summary(primprod.rand[[1]])
seas.rand <- lapply(1:100, function(x) lm(median_log_svl  ~ season_temp + season_prec, data=pam_vars_rand[[x]]))
summary(seas.rand[[1]])


heatbal.r2.rand <- lapply(1:100, function(x) summary(heatbal.rand[[x]])$r.squared)
heatbal.r2.rand <- unlist(heatbal.r2.rand)
primprod.r2.rand <- lapply(1:100, function(x) summary(primprod.rand[[x]])$r.squared)
primprod.r2.rand <- unlist(primprod.r2.rand)
seas.r2.rand <- lapply(1:100, function(x) summary(seas.rand[[x]])$r.squared)
seas.r2.rand <- unlist(seas.r2.rand)


# histograms
pdf("Histograms_NullModel_RandomizationSVL_SDMs.pdf", height = 3, width = 7)
par(mfrow=c(1,3))
hist(heatbal.r2.rand, xlab="random r2 values", main="heat balance")
abline(v=summary(heatbal.obs)$r.squared, col="red")
hist(primprod.r2.rand, xlab="random r2 values", main="primary productivity")
abline(v=summary(primprod.obs)$r.squared, col="red")
hist(seas.r2.rand, xlab="random r2 values", main="seasonality")
abline(v=summary(seas.obs)$r.squared, col="red")
dev.off()


save.image("Assemblage_Approach_NullModel_SDMs.RData")




# Spatial Autoregressive models
library(spdep)
library(ncf)

#------------------   Define coordinates, neighbourhoods, and spatial weights ----------------------------#
#Make a matrix of coordinates (X and Y coordinates)
coords<-cbind(svl_pam.obs$`Longitude(x)`,svl_pam.obs$`Latitude(y)`)
coords<-as.matrix(coords)

#Define neighbourhood (here distance in kms, min=0, max=500)
nb1.5<-dnearneigh(coords,0,1000,longlat = TRUE)

#Spatial weights, illustrated with coding style "W" (row standardized)
nb1.5.w<-nb2listw(nb1.5, glist=NULL, style="W", zero.policy=F)


#----------------------------------------   SARerr model ---------------------------------------------------------#

#SARerr model with neighbourhood distance 1.5 and coding style "W"

pam.heatbal <- pam_vars_svl[,c(1,2,3,4)]
pam.heatbal <- na.omit(pam.heatbal)
pred.ols.heatbal <- predict.lm(heatbal.obs)

# SARerr model observed
sar.heatbal <-errorsarlm(heatbal.obs, listw=nb1.5.w, data=pam_vars_svl, na.action = "na.omit")
summary(sar.heatbal)
res.sar.heatbal <- residuals(sar.heatbal)
sem.heatbal <- predict.sarlm(sar.heatbal, newdata=NULL,listw=nb1.5.w, pred.type="TS", zero.policy=TRUE)
heatbal.sum <- spatialreg::summary.sarlm(sar.heatbal, correlation=TRUE, Nagelkerke = TRUE)

pred.heatbal <- as.data.frame(sem.heatbal)
plot(pam.heatbal$median_log_svl, pred.heatbal$fit)
plot(pam.heatbal$median_log_svl, pred.ols.heatbal)

sar.primprod <-errorsarlm(primprod.obs, listw=nb1.5.w, data=pam_vars_svl, na.action = "na.omit")
summary(sar.primprod)
res.sar.primprod <- residuals(sar.primprod)
sem.primprod <- predict.sarlm(sar.primprod, newdata=NULL,listw=nb1.5.w, pred.type="TS", zero.policy=TRUE)
primprod.sum <- summary.sarlm(sar.primprod, correlation=TRUE, Nagelkerke = TRUE)

sar.seas <-errorsarlm(seas.obs, listw=nb1.5.w, data=pam_vars_svl, na.action = "na.omit")
summary(sar.seas)
res.sar.seas <- residuals(sar.seas)
sem.seas <- predict.sarlm(sar.seas, newdata=NULL,listw=nb1.5.w, pred.type="TS", zero.policy=TRUE)
seas.sum <- summary.sarlm(sar.seas, correlation=TRUE, Nagelkerke = TRUE)


# null model
sar.heatbal.rand <- lapply(1:100, function(x) errorsarlm(heatbal.obs, listw=nb1.5.w, data=pam_vars_rand[[x]], na.action = "na.omit"))
heatbal.sum.rand <- lapply(1:100, function(x) summary.sarlm(sar.heatbal.rand[[x]], correlation=TRUE, Nagelkerke = TRUE))
nak.heatbal.rand <- lapply(1:100, function(x) heatbal.sum.rand[[x]]$NK)
nak.heatbal.rand <- unlist(nak.heatbal.rand)

sar.heatbal.rand <- lapply(1:100, function(x) errorsarlm(heatbal.obs, listw=nb1.5.w, data=pam_vars_rand[[x]], na.action = "na.omit"))
heatbal.sum.rand <- lapply(1:100, function(x) summary.sarlm(sar.heatbal.rand[[x]], correlation=TRUE, Nagelkerke = TRUE))
nak.primprod.rand <- lapply(1:100, function(x) heatbal.sum.rand[[x]]$NK)
nak.primprod.rand <- unlist(nak.primprod.rand)

sar.seas.rand <- lapply(1:100, function(x) errorsarlm(seas.obs, listw=nb1.5.w, data=pam_vars_rand[[x]], na.action = "na.omit"))
seas.sum.rand <- lapply(1:100, function(x) summary.sarlm(sar.seas.rand[[x]], correlation=TRUE, Nagelkerke = TRUE))
nak.seas.rand <- lapply(1:100, function(x) seas.sum.rand[[x]]$NK)
nak.seas.rand <- unlist(nak.seas.rand)


library(DescTools)
library(fmsb)

# histograms
pdf("Histograms_NullModel_RandomizationSVL_SARerr_SDMs.pdf", height = 6, width = 7)
par(mfrow=c(2,3))
hist(heatbal.r2.rand, xlab="random r2 values", main="heat balance")
abline(v=summary(heatbal.obs)$r.squared, col="red")
hist(primprod.r2.rand, xlab="random r2 values", main="primary productivity")
abline(v=summary(primprod.obs)$r.squared, col="red")
hist(seas.r2.rand, xlab="random r2 values", main="seasonality")
abline(v=summary(seas.obs)$r.squared, col="red")

hist(nak.seas.rand, xlab="random pseudo-r2 Nagelkerke values", main="heat balance")
abline(v=heatbal.sum$NK, col="red")
hist(nak.primprod.rand, xlab="random pseudo-r2 Nagelkerke values", main="primary productivity")
abline(v=primprod.sum$NK, col="red")
hist(nak.seas.rand, xlab="random pseudo-r2 Nagelkerke values", main="seasonality")
abline(v=seas.sum$NK, col="red")
dev.off()


save.image("Assemblage_Approach_NullModel_SDMs.RData")


# individual hypothesis

pdf("univariate_regressions_SMDs.pdf", height = 6, width = 5)
par(mfrow=c(2,2))
plot(pam_vars_svl$anntemp, pam_vars_svl$median_log_svl, xlab="Annual mean temperature", ylab="median log10 SVL", pch=1, col="blue")
abline(sar.heatbal$coefficients[c(1,2)], col="red")
plot(pam_vars_svl$aet, pam_vars_svl$median_log_svl, xlab="Actual evapotranspiration", ylab="median log10 SVL", pch=1, col="blue")
abline(sar.primprod$coefficients[c(1,2)], col="red")
plot(pam_vars_svl$season_temp/100, pam_vars_svl$median_log_svl, xlab="Temperature seasonality", ylab="median log10 SVL", pch=1, col="blue")
abline(sar.seas$coefficients[1:2], col="red")
plot(pam_vars_svl$season_prec, pam_vars_svl$median_log_svl, xlab="Precipitation seasonality", ylab="median log10 SVL", pch=1, col="blue")
abline(sar.seas$coefficients[c(1,3)], col="red")
dev.off()

###############################


# FIGURAS NULL MODEL ALL
# body size gradients maps

svl_pam_convex <-lets.maplizer(PAM_convex, max.svl, PAM_convex$Species_name, func=median, ras=T)
svl_pam_roll <-lets.maplizer(PAM_roll, max.svl, PAM_roll$Species_name, func=median, ras=T)
svl_pam_pts <-lets.maplizer(PAM_points, max.svl, PAM_points$Species_name, func=median, ras=T)
svl_pam_models <-lets.maplizer(PAM_models, max.svl, PAM_points$Species_name, func=median, ras=T)

cols<-brewer.pal(n=9, name="RdYlBu")
data(wrld_simpl)

pdf("BodySizeGradientes_DifferentDatasets.pdf", width = 5, height = 5)
par(mfrow=c(2,2))
plot(svl_pam_convex$Raster, main="-Convex hulls-", col=cols)
plot(wrld_simpl, add = TRUE,  xlim=c(-115,-36), ylim=c(-20,36))

plot(svl_pam_roll$Raster, main="Roll et al's ranges", col=cols)
plot(wrld_simpl, add = TRUE,  xlim=c(-115,-36), ylim=c(-20,36))

plot(svl_pam_models$Raster, main="SDM's", col=cols)
plot(wrld_simpl, add = TRUE,  xlim=c(-115,-36), ylim=c(-20,36))

plot(svl_pam_pts$Raster, main="Occurrence records", col=cols)
plot(wrld_simpl, add = TRUE,  xlim=c(-115,-36), ylim=c(-20,36))
dev.off()



# 

heatbal.r2.rand.models <- heatbal.r2.rand
primprod.r2.rand.models <- primprod.r2.rand
seas.r2.rand.models <- seas.r2.rand
nak.heatbal.rand.models <- nak.heatbal.rand 
nak.primprod.rand.models <- nak.primprod.rand
nak.seas.rand.models <- nak.seas.rand

heatbal.r2.rand.points <- heatbal.r2.rand
primprod.r2.rand.points <- primprod.r2.rand
seas.r2.rand.points <- seas.r2.rand
nak.heatbal.rand.points <- nak.heatbal.rand 
nak.primprod.rand.points <- nak.primprod.rand
nak.seas.rand.points <- nak.seas.rand

heatbal.r2.rand.roll <- heatbal.r2.rand
primprod.r2.rand.roll <- primprod.r2.rand
seas.r2.rand.roll <- seas.r2.rand
nak.heatbal.rand.roll <- nak.heatbal.rand 
nak.primprod.rand.roll <- nak.primprod.rand
nak.seas.rand.roll <- nak.seas.rand

heatbal.r2.rand.convex <- heatbal.r2.rand
primprod.r2.rand.convex <- primprod.r2.rand
seas.r2.rand.convex <- seas.r2.rand
nak.heatbal.rand.convex <- nak.heatbal.rand 
nak.primprod.rand.convex <- nak.primprod.rand
nak.seas.rand.convex <- nak.seas.rand
save.image("NULL_MODELS_OLS_SAR.RData")




load("NULL_MODELS_OLS_SAR.RData")

# histograms
pdf("Histograms_NullModels_SARs.pdf", height = 10, width = 8)
par(mfrow=c(4,3))
hist(nak.seas.rand.convex, xlab="Nagelkerke r2 values", main="HeatBal -Convex hulls-")
abline(v=0.51, col="red")
hist(nak.primprod.rand.convex, xlab="Nagelkerke r2 values", main="PrimPro -Convex hulls-")
abline(v=0.51, col="red")
hist(nak.seas.rand.convex, xlab="Nagelkerke r2 values", main="SEAS -Convex hulls-")
abline(v=0.50, col="red")

hist(nak.seas.rand.roll, xlab="Nagelkerke r2 values", main="HeatBal -Roll etal's-")
abline(v=0.29, col="red")
hist(nak.primprod.rand.roll, xlab="Nagelkerke r2 values", main="PrimPro -Roll etal's-")
abline(v=0.37, col="red")
hist(nak.seas.rand.roll, xlab="Nagelkerke r2 values", main="SEAS -Roll etal's-")
abline(v=0.34, col="red")

hist(nak.seas.rand.models, xlab="Nagelkerke r2 values", main="HeatBal -SDM's-")
abline(v=0.64, col="red")
hist(nak.primprod.rand.models, xlab="Nagelkerke r2 values", main="PrimPro -SDM's-")
abline(v=0.64, col="red")
hist(nak.seas.rand.models, xlab="Nagelkerke r2 values", main="SEAS -SDM's-")
abline(v=0.68, col="red")

hist(nak.seas.rand.points, xlab="Nagelkerke r2 values", main="HeatBal -points-")
abline(v=0.11, col="red")
hist(nak.primprod.rand.points, xlab="Nagelkerke r2 values", main="PrimPro -points-")
abline(v=0.12, col="red")
hist(nak.seas.rand.points, xlab="Nagelkerke r2 values", main="SEAS -points-")
abline(v=0.12, col="red")
dev.off()















