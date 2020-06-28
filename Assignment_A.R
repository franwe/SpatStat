#####################################################################################
#                                                                                   #
#     Assignemnt Block A : Bianca Neubert, Franziska Wehrmann                       #
#                                                                                   #
#####################################################################################

################### THE JURA DATA ###################################################

setwd("~/Uni/Semester_10/Spatial_Statistics/Assignment_A")
source("build.convex.grid.R")

library(gstat)
library(sp)
library(lattice)

# Import the Jura file, inspect data and set up an object cj 
jura <- read.table("jura.txt", header = TRUE)
cj <- read.table("jura.txt", header = TRUE)
coordinates(cj) = ~x+y
class(cj)
summary(cj)
plot(cj)

# Generate a convex grid
jg <- data.frame(build.convex.grid(jura[,1], jura[,2],10000))
names(jg) <- c("x","y")
jg <- SpatialPoints(jg)

# Check jg to be spatial  
class(jg)
plot(jg)


#### For this task we only consider concentration of Ni
# Compute a map of Ni
# Distance matrix could be chosen to be Euclidean, which is the default of dist
d <- dist(jura[1:2])

spplot(cj, zcol = "Ni")
bubble(cj, "Ni", col="blue", main = "Ni concentrations")

xyplot(log(Ni) ~ sqrt(d), as.data.frame(jura))


# Compute (semi-)variograms from the data using the variogram function
### Constant trend for variable log(Ni)

v1  = variogram(log(Ni)~1, cj)
plot(v1)
v = variogram(log(Ni)~Landuse, cj)
plot(v)
# not much change

# alternative variograms
# eg change cutoff so that decreasing part is not considered
va1  = variogram(log(Ni)~1, cj, cutoff = 1.5)
plot(va1)
# change width
va2 = variogram(log(Ni)~1, cj,cutoff = 1.5, width =  0.05)
plot(va2)

# look at semivariogram cloud
va2 = variogram(log(Ni)~1, cj,cutoff = 1.5, width =  0.05, cloud = T)
plot(va2)

# fitting variogram models
v11.fit <- fit.variogram(v1, model = vgm(1, "Sph", 10, 1))
v12.fit <- fit.variogram(v1, model = vgm(1, "Exp", 10, 1))
v13.fit <- fit.variogram(v1, model = vgm(1, "Gau", 10, 1))
v14.fit <- fit.variogram(v1, model = vgm(1, "Mat", 10, 1))
v15.fit <- fit.variogram(v1, model = vgm("Nug"))
v16.fit <- fit.variogram(v1, model = vgm("Lin"))

# Produce plots
plot(v1, v1.fit)
plot(v2, v2.fit)
plot(v1, v11.fit)
plot(v1, v12.fit)
plot(v1, v13.fit)
plot(v1, v14.fit)
plot(v1, v15.fit, cutoff = 3)
plot(v1, v16.fit)


# Determine how to implement Cressie's version and compute a 
# robustified version
v3 = variogram(log(Ni)~1, cressie  = TRUE, cj)
v3.fit <- fit.variogram(v3, model = vgm(1, "Sph", 10, 1 ))
plot(v3, v3.fit)

# Compute OK and UK from specified model
# Ordinary kriging (with arbitrary variogram model and fitted variogram)
x <- krige(log(Ni)~1, cj, jg, model = vgm(1, "Sph", 5,1))
x2 <- krige(log(Ni)~1, cj, jg, model = v3.fit)

spplot(x["var1.pred"], main = "ordinary kriging predictions")
spplot(x["var1.var"],  main = "ordinary kriging variance")
spplot(x2["var1.pred"], main = "ordinary kriging predictions")
spplot(x2["var1.var"],  main = "ordinary kriging variance")
# fitted variogram yields better results which was expected

# Universal Kriging
# with cressie vario
y2 <- krige(log(Ni)~ x+y, cj, jg, model = v3.fit)
spplot(y2["var1.pred"], main = "universal kriging predictions")
spplot(y2["var1.var"],  main = "universal kriging variance")

# not much difference to ordinary

###### Now: geoR
library(geoR)

# Transform data into geodata object
# Inspect the geodata object, and compute and plot a variogram
# from the jura data only for Ni using
# - variofit  (ordinary and Cressie's specification)
ju.geo <- as.geodata(jura, data.col = 9)
vario <- variog(ju.geo)

# Compute Empirical Variogram, once classical once Cressie
vario1 <- variog(ju.geo)
vario2 <- variog(ju.geo, estimator.type = "modulus")
plot(vario1)
plot(vario2)

f1 <- variofit(vario1)
f2 <- variofit(vario2, cov.model = "matern", weights = "cressie")
plot(vario1)
lines(f1)
plot(vario2)
lines(f2)
lines(f3)
# Warum so komisch?

# - ML, REML (likfit)
ml <- likfit(ju.geo, ini.cov.pars = c(0.5, 0.5))
reml <- likfit(ju.geo, ini.cov.pars = c(0.5,0.5), lik.method = "REML")
plot(vario1)
lines(ml)
lines(reml, col = "red")

# - Profile LL (proflik)
# Computes profile likelihoods for model parameters 
pl  <- proflik(ml, ju.geo,ill.values=seq(0.5, 1.5, l=4),
               range.val=seq(0.1, .5, l=4))

plot(pl, nlevels = 16)

# - REML (gstat::fit.variogram.reml)
reml.fit <- fit.variogram.reml(log(Ni)~1, cj, jg, model = vgm(1, "Sph", range=5))
plot(reml.fit, cutoff=8)
va3 <- variogram(log(Ni)~1, cj)
plot(va3)
# Why are do they have totally different domain and values?

# Perform a visual modelling using the eyefit funciton
eyefit(vario1, silent = F)
# Compare all results
#### ??

################### THE GAMBIA MALARIA DATA #########################################

# Usa data on malaria prevalence in children obtained at 65 villages in Gambia
# Inspect the data and check wheter this data meets assumptions of geostatistical data
# If not: chose a reasonable modification
# Recall that we do not allow for multiple conicident locations
library(geoR)
library(leaflet)
library(viridis)
library(dplyr)
library(sp)
library(rgdal)
library(raster)

data(gambia)
str(gambia)
dim(unique(gambia[,c("x","y")]))

# Compute the prevalence per lovation defined as the number of positively tested
# children divided by the total number of tested children. To this end, define a
# dataframe d. Finally, add the longitude and latitude variables to d
gambia$count <- rep(1, 2035)
d <- gambia %>% 
      group_by(x,y) %>% 
      summarize(positive = sum(pos), total = sum(count)) %>%
      mutate(prev = positive / total) %>% ungroup

# Determine the CRS and specify the projection to CRS("+proj=longlat +datum=WGS84")
# Data is in UTM format (Easting/Northing)
spd <- SpatialPoints(d[, c("x", "y")],
                     proj4string = CRS("+proj=utm +zone=28 +datum=WGS84")
)
spdt <- spTransform(spd, CRS("+proj=longlat +datum=WGS84"))

d[,c("long", "lat")] <- coordinates(spdt)
head(d)

# Construct a map with the locations of the villages and the malaria prevalence
pal <- colorBin("viridis", bins = c(0, 0.25, 0.5, 0.75, 1))
leaflet(d) %>% 
  addProviderTiles(providers$CartoDB.Positron) %>% 
  addCircles(lng = ~long, lat = ~lat, color = ~pal(prev)) %>% 
  addLegend("bottomright", pal = pal, values = ~prev, title = "Prevalence") %>%
  addScaleBar(position = c("bottomleft"))

# To join covariate information on the evelation in Gambia to our model, we will use the
# getData() function from the raster library

# Define an R object called r through the getData() function where name is set to alt
# country is set to GMB and mask is set to TRUE
r <- getData(name = "alt", country = "GMB", mask = T)

# Compute a map of the relevant raster using the capacities of the leaflet through
# addRasterImage(r, \textit{colour .... }) where colour are defined by the code below
pal <- colorNumeric("viridis", values(r), na.color = "transparent")
leaflet(d) %>% 
  addProviderTiles(providers$CartoDB.Positron) %>% 
  addRasterImage(r, color = pal) %>%
  addLegend("bottomright",
            pal = pal, values = values(r),
            title = "Altitude") %>%
  addScaleBar(position = c("bottomleft"))


# Extract the altitude values at the villages locations using the extract function of
# the raster package
d$alt <- raster::extract(r, d[, c("long", "lat")])
head(d)

# Now, we turn to fitting the prevalence model using the INLA and SPDE approach where
# we consider the following specifications of the prevalences
# Y_i | P(s_i) ~ Bin(N_i, P(s_i))
# where P(s_i) is the true prevalence at location s_i, i = 1,...,n, and Y_i is the
# number of positive results out of N_i people sampled at s_i such that
# logit(P(s_i)) = beta0 + beta1 * altitude + f(s_i)
# where f(s_i) is a spatial random effect following a zero-mean Gaussian process
# with Mat?rn covariance function
library(INLA)

# Define a mesh setting max.edge to c(0.1, 5) and cutoff to 0.01 using the INLA 
# library and plot the mesh
coo <- cbind(d$long, d$lat)
mesh <- inla.mesh.2d(
  loc = coo, max.edge = c(0.1, 5),
  cutoff = 0.01
)
mesh$n
plot(mesh)

# Look at different specifications of max.edge and cutoff and discuss the effect
mesh1 <- inla.mesh.2d(
  loc = coo, max.edge = c(0.1, 10),
  cutoff = 0.01
)
plot(mesh1)

spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)

# Generate an index set and projection matrix
A <- inla.spde.make.A(mesh = mesh, loc = coo)
indexs <- inla.spde.make.index("s", spde$n.spde)

# Compute the predicting location through rasterToPoints(r) and inspect the new object
# eg dim(s0)
s0 <- rasterToPoints(r)
dim(s0)

# Compute a lower number of predicting locations by aggregation from r specifying
# fact = 5 in the aggregate function and extract the coordinates from this reduced
# set of prediction points and compute a prediction matrix Ap
ag <- aggregate(r, fact = 5)
c <- rasterToPoints(ag)
dim(c)
coop <- c[, c("x", "y")]

Ap <- inla.spde.make.A(mesh, loc = coop)

# Set up the stack data for the estimation and prediction
stk.e <- inla.stack(tag = "est",
                    data = list(y = d$positive, numtrials = d$total),
                    A = list(1, A),
                    effects = list(data.frame(b0 = 1, cov = d$alt),  s = indexs)
                    )

stk.p <- inla.stack(tag = "pred",
                    data = list(y = NA, numtrials = NA),
                    A = list(1, Ap),
                    effects = list(data.frame(b0 = 1, cov = c[, 3]), s = indexs)
                    )

stk.full <- inla.stack(stk.e, stk.p)

# For the model, we specify the following formula
formula <- y ~ 0 + b0 + cov + f(s, model = spde)

# What is the effect of 0 in the above formula?
# 0 removes intercept and adding a covariate term b0 so that all coovariate
# terms can be captured in the projection matrix

# Finally, we can call inla
res <- inla(formula, family = "binomial", Ntrials = numtrials,
            control.family = list(link = "logit"),
            data = inla.stack.data(stk.full),
            control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.full)))

# Inspect the results

# Create a map of the posteriori means from the results using an index obtaind from
index <- inla.stack.index(stack = stk.full, tag = "pred")$data

prev_mean <- res$summary.fitted.values[index, "mean"]

pal <- colorNumeric("viridis", c(0, 1), na.color = "transparent")

leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircles(
    lng = coop[, 1], lat = coop[, 2],
    color = pal(prev_mean)
  ) %>%
  addLegend("bottomright",
            pal = pal, values = prev_mean,
            title = "Prev."
  ) %>%
  addScaleBar(position = c("bottomleft"))

r_prev_mean <- rasterize(
  x = coop, y = ag, field = prev_mean,
  fun = mean
)

pal <- colorNumeric("viridis", c(0, 1), na.color = "transparent")

leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(r_prev_mean, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright",
            pal = pal,
            values = values(r_prev_mean), title = "Prev."
  ) %>%
  addScaleBar(position = c("bottomleft"))
# Hint: the predicted values correspond to a set of points. We can create a raster with 
# the predicted values using rasterize()
 



