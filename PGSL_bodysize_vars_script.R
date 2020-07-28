

library(phytools)
library(caper)
library(geiger)
library(ape)
library(picante)
library(nlme)
library(rcompanion)

datos <- read.csv("Cross_Species_BodySize_Vars.csv")
tree<-read.nexus("AnolisFullyResolvedCloseToConsensus.tre")

# match
match.all<-match.phylo.data(tree, datos)
pruned<-match.all$phy
datos.pruned<-match.all$data

#pruned.scaled <- pruned # this is to fit the MrBayes tree to accomodate branch lenghts in a proper way
#pruned.scaled$edge.length <- pruned.scaled$edge.length*1000



# PGLS

comp.data<-comparative.data(pruned, datos, names.col="species",vcv.dim=2,warn.dropped=TRUE)

model.full.v1 <- pgls(Log10_meanSVL ~ Anntemp + AET + GVI + NDVI + season_temp + season_prec, data=comp.data, lambda="ML")
summary(model.full.v1)

model.full.v2 <- pgls(Log10_meanSVL ~ Anntemp + AET + GVI + NDVI + season_temp + season_prec, data=comp.data)
summary(model.full.v2)

par(mfrow=c(2,2))
plot(model.full.v1)

par(mfrow=c(2,2))
plot(model.full.v2)

# PGLS models with caper
model.heatbal.v1 <- pgls(Log10_meanSVL ~ Anntemp, data=comp.data, lambda="ML")
summary(model.heatbal.v1)
AIC(model.heatbal.v1)

model.primpro.v1 <- pgls(Log10_meanSVL ~ AET, data=comp.data, lambda="ML")
summary(model.primpro.v1)
AIC(model.primpro.v1)

model.seas.v1 <- pgls(Log10_meanSVL ~ season_temp + season_prec, data=comp.data, lambda="ML")
summary(model.seas.v1)
AIC(model.seas.v1)


pdf("pgls_regressions.pdf", height = 6, width = 5)
par(mfrow=c(2,2))
plot(datos$Anntemp, datos$Log10_meanSVL, xlab="Annual mean temperature", ylab="log10 SVL", col="blue")
abline(model.heatbal.v1, col="red")
plot(datos$AET, datos$Log10_meanSVL, xlab="Actual evapotranspiration", ylab="log10 SVL", col="blue")
abline(model.primpro.v1, col="red")
plot(datos$season_temp/100, datos$Log10_meanSVL, xlab="Temperature seasonality", ylab="log10 SVL", col="blue")
abline(model.seas.v1$model$coef[c(1,2)], col="red")
plot(datos$season_prec, datos$Log10_meanSVL, xlab="Precipitation seasonality", ylab="log10 SVL", col="blue")
abline(model.seas.v1$model$coef[c(1,3)], col="red")
dev.off()


save.image("PGLS_BodySizeClim.RData")


pdf("MacroecologicalModels.pdf", height = 6, width = 10)
par(mfrow=c(2,4))
plot(pam_vars_svl$anntemp, pam_vars_svl$median_log_svl, xlab="AMT", ylab="median log10 SVL", pch=1, col="blue", main="Assemblage-based", cex.main=0.8)
abline(sar.heatbal, col="red")
plot(pam_vars_svl$aet, pam_vars_svl$median_log_svl, xlab="AET", ylab="median log10 SVL", pch=1, col="blue", main="Assemblage-based", cex.main=0.8)
abline(sar.primprod, col="red")
plot(pam_vars_svl$season_temp/100, pam_vars_svl$median_log_svl, xlab="TS", ylab="median log10 SVL", pch=1, col="blue", main="Assemblage-based", cex.main=0.8)
abline(sar.seas$coefficients[1:2], col="red")
plot(pam_vars_svl$season_prec, pam_vars_svl$median_log_svl, xlab="PS", ylab="median log10 SVL", pch=1, col="blue", main="Assemblage-based", cex.main=0.8)
abline(sar.seas$coefficients[c(1,3)], col="red")

plot(datos$Anntemp, datos$Log10_meanSVL, xlab="AMT", ylab="log10 SVL", col="blue", main="Cross-Species", cex.main=0.8)
abline(model.heatbal.v1, col="red")
plot(datos$AET, datos$Log10_meanSVL, xlab="AET", ylab="log10 SVL", col="blue", main="Cross-Species", cex.main=0.8)
abline(model.primpro.v1, col="red")
plot(datos$season_temp/100, datos$Log10_meanSVL, xlab="TS", ylab="log10 SVL", col="blue", main="Cross-Species", cex.main=0.8)
abline(model.seas.v1$model$coef[c(1,2)], col="red")
plot(datos$season_prec, datos$Log10_meanSVL, xlab="PS", ylab="log10 SVL", col="blue", main="Cross-Species", cex.main=0.8)
abline(model.seas.v1$model$coef[c(1,3)], col="red")
dev.off()














# PGLS BM model
bm<-corBrownian(1,pruned)
bm

model1 <- gls(Log10_meanSVL ~ Anntemp + AET + GVI + NDVI + season_temp + season_prec, data=datos, correlation=bm)
summary(model1)
phylosig(pruned,residuals(model1))

model2<-gls(Log10_meanSVL ~ Anntemp + AET + GVI + NDVI + season_temp + season_prec, data=datos, correlation=corPagel(1,pruned))
summary(model2)
phylosig(pruned,residuals(model2))

model2.heatbal <- gls(Log10_meanSVL ~ Anntemp, data=datos, correlation=corPagel(1,pruned))
model2.primpro <- gls(Log10_meanSVL ~ AET + GVI + NDVI, data=datos, correlation=corPagel(1,pruned))
model2.seas <- gls(Log10_meanSVL ~ season_temp + season_prec, data=datos, correlation=corPagel(1,pruned))

model2.heatbal <- gls(Log10_meanSVL ~ Anntemp, data=datos, correlation=bm)
model2.primpro <- gls(Log10_meanSVL ~ AET + GVI + NDVI, data=datos, correlation=bm)
model2.seas <- gls(Log10_meanSVL ~ season_temp + season_prec, data=datos, correlation=bm)


summary(model2.heatbal)
AIC(model2.heatbal)

summary(model2.primpro)
AIC(model2.primpro)

summary(model2.seas)
AIC(model2.seas)

nagelkerke(model2.heatbal)
nagelkerke(model2.primpro)
nagelkerke(model2.seas)

# fabro
source("C:/Users/juvel_000/Dropbox/PosDoc_DGAPA/p_s_components_fabro.R")

#correr la funci?n
anole.pscomps <- p.s.comps(phy=pruned, sppnames=pruned$tip.label, data = datos.pruned, colnum=4)
pscomps.fabro <- as.data.frame(anole.pscomps[[2]])
rn<-row.names(pscomps.fabro)
pscomps.fabro <-pscomps.fabro[order(as.numeric(rn)), ]


# combine
datos.pruned2 <- cbind(datos.pruned, pscomps.fabro)
plot(datos.pruned2$P.comp, datos.pruned2$Phy_component)
plot(datos.pruned2$S.comp, datos.pruned2$Spec_component)

# linear regresions

heatbal.p <- lm(datos.pruned2$P.comp~datos.pruned2$Anntemp + datos.pruned2$Midpoint_Lat)
heatbal.s <- lm(datos.pruned2$S.comp~datos.pruned2$Anntemp + datos.pruned2$Midpoint_Lat)

primpro.p <- lm(datos.pruned2$P.comp~datos.pruned2$AET+datos.pruned2$GVI+datos.pruned2$NDVI + datos.pruned2$Midpoint_Lat)
primpro.s <- lm(datos.pruned2$S.comp~datos.pruned2$AET+datos.pruned2$GVI+datos.pruned2$NDVI + datos.pruned2$Midpoint_Lat)

seas.p <- lm(datos.pruned2$P.comp~datos.pruned2$season_temp+datos.pruned2$season_prec + datos.pruned2$Midpoint_Lat)
seas.s <- lm(datos.pruned2$S.comp~datos.pruned2$season_temp+datos.pruned2$season_prec + datos.pruned2$Midpoint_Lat)

library(spdep)
library(ncf)

# pseudo R2 for OLS
heatbal.pcomp.nagelkerke <- nagelkerke(heatbal.p)
heatbal.scomp.nagelkerke <- nagelkerke(heatbal.s)

primpro.pcomp.nagelkerke <- nagelkerke(primpro.p)
primpro.scomp.nagelkerke <- nagelkerke(primpro.s)

seas.nagelkerke <- nagelkerke(seas.p)
seas.scomp.nagelkerke <- nagelkerke(seas.s)


AIC(heatbal.p)
AIC(heatbal.s)
AIC(primpro.p)
AIC(primpro.s)
AIC(seas.p)
AIC(seas.s)


# partial regression
library(plsdepot)
pls1 = plsreg1(datos.pruned2[,7:12], datos.pruned2[,5], comps = 6)

pls1$R2
pls1$resid
plot(pls1)

#plot each observation predicted versus actual
plot(datos.pruned2$P.comp, pls1$y.pred, type ="n", xlab="Original", ylab = "Predicted")
title("Comparison of responses", cex.main = 0.9) 
abline(a = 0, b = 1, col = "gray85", lwd = 2) 
text(datos.pruned2$P.comp, pls1$y.pred, col = "#5592e3") 
       
       


# scaling
ou.pruned<-corMartins(1,phy=pruned)

ou.gls1<-gls(Log10_meanSVL ~ Anntemp + Midpoint_Lat,correlation=ou.pruned,data=datos.pruned)
summary(ou.gls1)
AIC(ou.gls1)

ou.gls2<-gls(Log10_meanSVL ~ AET + GVI + NDVI + Midpoint_Lat,correlation=ou.pruned,data=datos.pruned)
summary(ou.gls2)
AIC(ou.gls2)

ou.gls3<-gls(Log10_meanSVL ~ season_temp + season_prec + Midpoint_Lat,correlation=ou.pruned,data=datos.pruned)
summary(ou.gls3)
AIC(ou.gls3)




bm<-corBrownian(1,pruned)
bm


pgls.model1 <- gls(Log10_meanSVL ~ Anntemp + Midpoint_Lat, correlation=corBrownian(phy=pruned), data=datos.pruned, method = "ML")
summary(pgls.model1)
AIC(pgls.model1)

pgls.model2 <- gls(Log10_meanSVL ~ AET + GVI + NDVI + Midpoint_Lat, correlation=corBrownian(phy=pruned), data=datos.pruned, method = "ML")
summary(pgls.model2)
AIC(pgls.model2)

pgls.model3 <- gls(Log10_meanSVL ~ season_temp + season_prec + Midpoint_Lat, correlation=corBrownian(phy=pruned), data=datos.pruned, method = "ML")
summary(pgls.model3)
AIC(pgls.model3)





# effect size
library(yhat)
effect.size(model1)

phylosig(pruned,residuals(pgls.model1))


datos <- read.csv("C:/Users/juvel_000/Dropbox/PosDoc_DGAPA/Cross_Species_BodySize_Vars.csv")
comp.data<-comparative.data(pruned, datos, names.col="species",vcv.dim=2,warn.dropped=TRUE)


model.full.v1 <- pgls(Log10_meanSVL ~ Anntemp + AET + GVI + NDVI + season_temp + season_prec, data=comp.data, lambda="ML")
summary(model.full.v1)


model1<-pgls(Log10_meanSVL ~ Anntemp + Midpoint_Lat, data=comp.data,lambda="ML")
summary(model1)
nagelkerke(model1)

model2<-pgls(Log10_meanSVL ~ AET + GVI + NDVI + Midpoint_Lat, data=comp.data,lambda="ML")
summary(model2)


model3<-pgls(Log10_meanSVL ~ season_temp + season_prec + Midpoint_Lat, data=comp.data,lambda="ML")
summary(model3)

# Cohen's f2 = r2 / (1  - r2)
# and is interpreted as  0.02=small, 0.15=medium, 0.35=large





# without geography
pgls.model4 <- gls(Log10_meanSVL ~ Anntemp, correlation=corBrownian(phy=pruned), data=datos.pruned, method = "ML")
summary(pgls.model4)

pgls.model5 <- gls(Log10_meanSVL ~ AET + GVI + NDVI, correlation=corBrownian(phy=pruned), data=datos.pruned, method = "ML")
summary(pgls.model5)

pgls.model6 <- gls(Log10_meanSVL ~ season_temp + season_prec, correlation=corBrownian(phy=pruned), data=datos.pruned, method = "ML")
summary(pgls.model6)





# geography
pgls.model7 <- gls(P.comp ~ Midpoint_Lat, correlation=corBrownian(phy=pruned), data=datos.pruned, method = "ML")
summary(pgls.model7)

pgls.model8 <- gls(S.comp ~ Midpoint_Lat, correlation=corBrownian(phy=pruned), data=datos.pruned, method = "ML")
summary(pgls.model8)



# Brownian full
pgls.model1.BM <- gls(P.comp ~ Anntemp + AET + GVI + NDVI + season_temp + season_prec + Midpoint_Lat, 
                       correlation=corBrownian(phy=pruned), data=datos.pruned, method = "ML")
summary(pgls.model1.BM)

pgls.model2.BM <- gls(S.comp ~ Anntemp + AET + GVI + NDVI + season_temp + season_prec + Midpoint_Lat, 
                      correlation=corBrownian(phy=pruned), data=datos.pruned, method = "ML")
summary(pgls.model2.BM)


# OU full
pgls.model1.OU <- gls(P.comp ~ Anntemp + AET + GVI + NDVI + season_temp + season_prec + Midpoint_Lat, 
                          correlation=corPagel(1, phy=pruned, fixed=FALSE), data=datos.pruned, method = "ML")
summary(pgls.model1.OU)

pgls.model2.OU <- gls(S.comp ~ Anntemp + Midpoint_Lat, 
                      correlation=corMartins(1, phy=pruned, fixed=FALSE), data=datos.pruned, method = "ML")
summary(pgls.model2.OU)


windows()
par(mfrow=c(1,3))
plot(P.comp ~ season_temp, data=datos.pruned)
abline(a=coef(pgls.model5)[1], b=coef(pgls.model5)[2])

plot(P.comp ~ season_prec, data=datos.pruned)
abline(a=coef(pgls.model5)[1], b=coef(pgls.model5)[3])

plot(S.comp ~ Anntemp, data=datos.pruned)
abline(a=coef(pgls.model2)[1], b=coef(pgls.model2)[2])


