##### AREAL DATA

rm(list=ls())
packages <- c("maps", "maptools", "mapview", "osmar", "osmdata", 
              "raster", "rdgal", "sf", "spdep", "gstat","dplyr", 
              "splancs","tidyr", "ggplot", "reshape2", "sp", "readxl",
              "elasticnet","pgirmess", "RColorBrewer", "classInt", "spgwr")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())),
                   dependencies = TRUE)  
}

library(maps)         ## Projections
library(maptools)     ## Data management
library(sp)           ## Data management
library(spdep)        ## Spatial autocorrelation
library(gstat)        ## Geostatistics
library(splancs)      ## Kernel Density
library(spatstat)     ## Geostatistics
library(pgirmess)     ## Spatial autocorrelation
library(RColorBrewer) ## Visualization
library(rgdal)
library(classInt)     ## Class intervals
library(spgwr)        ## GWR
library(raster)


guerry <- readOGR("guerry/Guerry.shp")
summary(guerry)
names(guerry)

guerry$Suicids <- as.integer(guerry$Suicids)
guerry$Wealth <- as.integer(guerry$Wealth)
guerry$Clergy <- as.integer(guerry$Clergy)

p1 <- spplot(guerry, "Wealth", main="Wealth")
p2 <- spplot(guerry, "Clergy", main="Clergy")
p3 <- spplot(guerry, "Suicids", main="Suicides")

require(gridExtra)
grid.arrange(p1, p2, p3, ncol = 3)

# regression Y
## Linear Model
reg <- lm(Suicids ~ Wealth + Clergy, data=guerry)
summary(reg)

## Plot residuals
## fix color palette once for all plots of this type, for better comparison
res.max = 98661.66  # maximum residual of (lm, car, sar, sdm)
res.palette <- colorRampPalette(c("red","orange","white", "lightgreen","green"), space = "rgb")
pal <- res.palette(5)
m <- res.max
breaks = round(c(-m, -0.4*m, -0.1*m, 0.1*m, 0.4*m, m), digits=0)

res <- reg$residuals
classes_fx <- classIntervals(res, n=5, style="fixed", fixedBreaks=breaks, rtimes = 1)
cols <- findColours(classes_fx,pal)
cols <- findColours(classes_fx,pal)
plot(guerry,col=cols, border="grey",pretty=T)
legend(x="bottom",cex=1,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),title="Residuals from SDM Model",ncol=5)

## Neighbors
xy <- coordinates(guerry)
W_cont_el <- poly2nb(guerry, queen=FALSE)
W_cont_el_mat <- nb2listw(W_cont_el, style="B", zero.policy=TRUE)
plot(W_cont_el_mat, coords=xy, cex=0.1, col="gray")

moran.test(res, listw=W_cont_el_mat, zero.policy=T)      # TODO: whats the difference? to moran below?
moran.test(guerry$Suicids, listw=W_cont_el_mat, zero.policy=T)  # TODO: Interpretation
geary.test(guerry$Suicids, listw=W_cont_el_mat, zero.policy=T)


# CAR - SAR - SDM

## CAR
car.out <- spdep::spautolm(Suicids ~ Wealth + Clergy, data=guerry,
                           listw=W_cont_el_mat, family="CAR")
mod.car <- fitted(car.out)
print(car.out)
summary(car.out)

res <- car.out$fit$residuals
classes_fx <- classIntervals(res, n=5, style="fixed", fixedBreaks=breaks, rtimes = 1)
cols <- findColours(classes_fx,pal)
plot(guerry,col=cols, border="grey",pretty=T)
legend(x="bottom",cex=1,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),
       title="Residuals from CAR Model",ncol=5)


## SAR
mod.sar <- lagsarlm(Suicids ~ Wealth + Clergy, data=guerry,
                    listw=W_cont_el_mat, zero.policy=T, tol.solve=1e-12)
summary(mod.sar)

res <- mod.sar$residuals
classes_fx <- classIntervals(res, n=5, style="fixed", fixedBreaks=breaks, rtimes = 1)
cols <- findColours(classes_fx,pal)
plot(guerry,col=cols, border="grey",pretty=T)
legend(x="bottom",cex=1,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),
       title="Residuals from SAR Model",ncol=5)


## SDM
mod.sdm <- lagsarlm(Suicids ~ Wealth + Clergy, data=guerry, listw=W_cont_el_mat, 
                    zero.policy=T, type="mixed", tol.solve=1e-12)
summary(mod.sdm)

res <- mod.sdm$residuals
classes_fx <- classIntervals(res, n=5, style="fixed", fixedBreaks=breaks, rtimes = 1)
cols <- findColours(classes_fx,pal)
plot(guerry,col=cols, border="grey")
legend(x="bottom",cex=1,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),title="Residuals from SDM Model",ncol=5)

## Plot values
Original.plot.values <- spplot(guerry, "Suicids", main="Guerry Data")
guerry$CAR <- car.out$fit$fitted.values
CAR.plot.values <- spplot(guerry, "CAR", main="CAR")
guerry$SAR <- mod.sar$fitted.values
SAR.plot.values <- spplot(guerry, "SAR", main="SAR")
guerry$SDM <- mod.sdm$fitted.values
SDM.plot.values <- spplot(guerry, "SDM", main="SDM")

grid.arrange(Original.plot.values, 
             CAR.plot.values, 
             SAR.plot.values, 
             SDM.plot.values, ncol = 2)



res <- res.sdm
print(max(-min(res), max(res)))
  
min(res)

res.max = 98661.66  # maximum residual in (lm, car, sar, sdm)

## Versuch die residuals hübscher zu plotten und auch in einem Grid dar zu stellen. 
## p1 <- plot(...residuals...) gibt leider ein NULL-element zurück, somit kann 
## p1 nicht in grid.arrange() gegeben werden.
## also versuch mit spplot zu plotten, aber dann weiß ich nicht, wie man die colormap 
## an die selbst gewählten FixedBreaks anpasst. 


library(RColorBrewer)
display.brewer.all()
my.palette <- brewer.pal(n = 7, name = "RdBu")

res <- mod.sdm$residuals
m = max(-min(res), max(res))
breaks = round(c(-m, -0.4*m, -0.1*m, 0.1*m, 0.4*m, m), digits=2)
classes_fx <- classIntervals(res, n=5, style="fixed", fixedBreaks=breaks, rtimes = 1)
cols <- findColours(classes_fx,pal)
SDM.residuals <- plot(guerry,col=cols, border="grey",pretty=T)
legend(x="bottom",cex=1,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),title="Residuals from SDM Model",ncol=5)

guerry$SDM.res <- mod.sdm$residuals
spplot(guerry, "SDM.res", col.regions=my.palette, zlim=c(-2,-1), cuts=6, main="SDM.res")

spplot(guerry, "SDM.res", colorkey = colorkey, col.regions=my.palette, cuts=6, zlim=c(0,5))

list(at = breaks)

colorkey = list(height = 1, labels = list(at = breaks), labels = at)
at <- breaks
seq(0.5, length(at) -0.5)

?seq
