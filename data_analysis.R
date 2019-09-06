setwd("C:/Users/gjbowen/Dropbox/Hypomirror/Chile/chileLakeWater")
setwd("C:/Users/u0133977/Dropbox/Hypomirror/Chile/chileLakeWater")

#####
#Model selection for Chilean lake data
#####

#Read data
l = read.csv("lakes.csv")
#Give lat and lon correct sign
l$Latitude = -l$Latitude
l$Longitude = -l$Longitude
#Clean up vars
l = l[,!names(l) %in% c("Salp.ID", "Lake.name")]
l = l[!is.na(l$Latitude),]
#plot Dex
plot(l$d2H - 8 * l$d18O, col = "white")
text(row.names(l), l$d2H - 8 * l$d18O, row.names(l))
#A couple sites are evaporated
l = rbind(l[1:20,], l[23:32,])
#Split the data
l.o = l[,!names(l) == "d2H"]
l.h = l[,!names(l) == "d18O"]

#Fit base model using only non-seasonal terms (seasonal P not significant)
mo.all = lm(l$d18O ~ ., data = l.o[,1:5])
summary(mo.all)

#Find and add interactions
add1(mo.all, .~. + .^2, test = "F")
#Add strongest term
add1(update(mo.all, .~. + Latitude:Modeled.MAP), .~. + .^2, test = "F")
#No more significant interactions
mo.1 = update(mo.all, .~. + Latitude:Modeled.MAP)
summary(mo.1)

#Start dropping variables
drop1(mo.1, test = "F")
#Remove weakest term
drop1(update(mo.1, .~. - Elevation), test = "F")
#Next
drop1(update(mo.1, .~. - Elevation - Modeled.MAT), test = "F")
#That's it
mo.2 = update(mo.1, .~. - Elevation - Modeled.MAT)
summary(mo.2)

#Base model
mh.all = lm(l$d2H ~ ., data = l.h[,1:5])
summary(mh.all)

#Find and add interactions
add1(mh.all, .~. + .^2, test = "F")
#Add strongest term
add1(update(mh.all, .~. + Latitude:Elevation), .~. + .^2, test = "F")
#No more significant interactions
mh.1 = update(mh.all, .~. + Latitude:Elevation)
summary(mh.1)

#Start dropping variables
drop1(mh.1, test = "F")
#Remove weakest term
drop1(update(mh.1, .~. - Modeled.MAT), test = "F")
#That's all
mh.2 = update(mh.1, .~. - Modeled.MAT)
summary(mh.2)

#Lets look at using the H form for O model
mo.h = lm(d18O ~ Latitude + Longitude + Elevation + Modeled.MAP + Latitude:Elevation, data = l)
summary(mo.h)

#Multivariate model
m.all = lm(cbind(l$d18O, l$d2H) ~ ., data = l[,1:5])
summary(m.all)
library(car)
Anova(m.all)

m.all.1 = update(m.all, .~. - Elevation - Modeled.MAT)
summary(m.all.1)
Anova(m.all.1)
anova(m.all, m.all.1)
linearHypothesis(m.all, hypothesis.matrix = c("Elevation = 0", "Modeled.MAT = 0"))

#Check for spatial autocorrelation
library(ape)
r = data.frame(l, m.all.1$residuals)
ro = r$X1
rh = r$X2
#weights
w = matrix(double(), nrow = nrow(r), ncol = nrow(r))
for(i in 1:nrow(r)){
  for(j in 1:nrow(r)){
    w[i, j] = 1 / sqrt((r$Latitude[i] - r$Latitude[j])^2 + (r$Longitude[i] - r$Longitude[j])^2)
  }  
}
diag(w) = 0

#none
Moran.I(ro, w)
Moran.I(rh, w)

#####
#Make isoscape
#####
library(raster)

#Import DEM
dem = raster("../chile.tif")
str(dem)
plot(dem)
l.spdf = SpatialPointsDataFrame(data.frame(l$Longitude, l$Latitude), data = l, proj4string = CRS(proj4string(dem)))
points(l.spdf)
#grid resolution
res(dem)

#Set extent for analysis
b = bbox(l.spdf)
bb = b + matrix(c(-2, 2, -2, 2), nrow = 2, byrow = TRUE) 

#Import monthly P
library(ncdf4)
p = brick("../CR2MET_p.nc")

#Procees to MAP
map = sum(p) / nlayers(p) * 12
plot(map)
str(map)
#grid resolution
res(map)

#Rescale and rectify grids
dem2 = resample(dem, map)
bbox(dem2) - bbox(map)
res(dem2) - res(map)
dem3 = crop(dem2, bb)
map2 = crop(map, bb)

#Prepare dataframe of predictors
idv = data.frame("Modeled.MAP" = values(map2), "Elevation" = values(dem3))

#Add lat and long values
xmn = bb[1,1]
xmx = bb[1,2]
ymn = bb[2,1]
ymx = bb[2,2]
xlen = ncol(map2)
ylen = nrow(map2)
xres = res(map2)[1]
yres = res(map2)[2]
lons = rep(seq(xmn + xres/2, xmx - xres/2, by = xres), ylen)
lat.bits = seq(ymx - yres/2, ymn + yres/2, by = -yres)
lats = double()
for(i in 1:ylen){
  lats = c(lats, rep(lat.bits[i], xlen))
}
idv$Latitude = lats
idv$Longitude = lons

#Check it!
tst = map2
tst = setValues(tst, idv$Elevation)
plot(tst)
points(l.spdf)

#Apply H and O model equations
d2H = predict(mh.2, idv)
d18O = predict(mo.2, idv)
iso = predict(m.all.1, idv)

#Package as rasters
d2H.r = map2
d18O.r = map2
d2H.r = setValues(d2H.r, iso[,2])
d18O.r = setValues(d18O.r, iso[,1])
plot(d2H.r)
plot(d18O.r)
plot(d2H.r - 8*d18O.r)

#Save
writeRaster(d2H.r, "d2H.tif")
writeRaster(d18O.r, "d18O.tif")

#Load precip data
library(RODBC)
ch = odbcConnect("WIDB")
p = sqlQuery(ch, "SELECT Sites.Site_ID, Sites.Latitude, Sites.Longitude, Sites.Elevation_mabsl, Samples.Sample_ID,
              Samples.Collection_Date, AVG(Water_Isotope_Data.d2H), AVG(Water_Isotope_Data.d18O) FROM 
              Water_Isotope_Data INNER JOIN Samples ON Water_Isotope_Data.Sample_ID = Samples.Sample_ID
             INNER JOIN Sites ON Samples.Site_ID = Sites.Site_ID WHERE Samples.Type = 'Precipitation' AND 
             Sites.Latitude < -35.9 AND Sites.Latitude > -46.01 AND Sites.Longitude > -75.9 AND Sites.Longitude < -69
             GROUP BY Samples.Sample_ID")
View(p)
pamt = sqlQuery(ch, "SELECT Samples.Sample_ID, Climate_Data.Precipitation_mm FROM 
              Climate_Data INNER JOIN Samples ON Climate_Data.Sample_ID = Samples.Sample_ID
             INNER JOIN Sites ON Samples.Site_ID = Sites.Site_ID WHERE Samples.Type = 'Precipitation' AND 
             Sites.Latitude < -35.9 AND Sites.Latitude > -46.01 AND Sites.Longitude > -75.9 AND Sites.Longitude < -69
             GROUP BY Samples.Sample_ID")
pp = merge(p, pamt)
names(pp)[7:8] = c("d2H", "d18O")
View(pp)
po = data.frame(Site_ID = pp$Site_ID,  pp$d18O * pp$Precipitation_mm, pp$Precipitation_mm)
ph = data.frame(Site_ID = pp$Site_ID, pp$d2H * pp$Precipitation_mm, pp$Precipitation_mm)
po = po[!is.na(pp$d18O),]
ph = ph[!is.na(pp$d2H),]

sites = unique(p[,1:3])
sites$Count = integer(nrow(sites))
sites$d2H = sites$d18O = double(nrow(sites))
for(i in 1:nrow(sites)){
  sites$Count[i] = length(p$Sample_ID[p$Site_ID == sites$Site_ID[i]])
  sites$d18O[i] = sum(po[po$Site_ID == sites$Site_ID[i], 2]) / sum(po[po$Site_ID == sites$Site_ID[i], 3])
  sites$d2H[i] = sum(ph[ph$Site_ID == sites$Site_ID[i], 2]) / sum(ph[ph$Site_ID == sites$Site_ID[i], 3])
}
View(sites)
sites = sites[sites$Site_ID != "isomap_8776501",]
sites.spdf = SpatialPointsDataFrame(data.frame(sites$Longitude, sites$Latitude), data = sites, 
                                    proj4string = CRS(proj4string(dem)))
write.csv(pp, "all_precip.csv")

#Plots
library(rgdal)
bnd = readOGR("South_America.shp")
bnd = bnd[bnd$COUNTRY == "Chile",]

png("../map.png", res=600, units = "in", width = 8, height = 6)
layout(matrix(c(1,2), nrow = 1))

par(mai = c(0.6,0.6,0.2,1.25))

plot(d2H.r, xlim = c(-76, -70),
     legend.args = list(text = expression("Lake isoscape "*delta^{2}*"H (\u2030)"),
                               side = 4, line = 2.8))
lines(bnd, col="dark grey")
points(l.spdf, col = "blue")

plot(d18O.r, xlim = c(-76, -70), 
     legend.args = list(text = expression("Lake isoscape "*delta^{18}*"O (\u2030)"),
                              side = 4, line = 2.8))
lines(bnd, col="dark grey")
points(sites.spdf, col = "blue")

dev.off()

sites.spdf$mapH = extract(d2H.r, sites.spdf)
sites.spdf$mapO = extract(d18O.r, sites.spdf)

hmod = lm(sites.spdf$d2H ~ sites.spdf$mapH)
omod = lm(sites.spdf$d18O ~ sites.spdf$mapO)

equation <- function(mod) {
  lm_coef <- list(a = as.numeric(round(coef(mod)[1], digits = 2)),
                  b = as.numeric(round(coef(mod)[2], digits = 2)), r2 = round(summary(mod)$r.squared,
                                                                              digits = 2))
  lm_eq <- substitute(italic(y) == a + b %.% italic(x) *
                        "," ~ ~italic(R)^2 ~ "=" ~ r2, lm_coef)
  as.expression(lm_eq)
}

png("../pred.png", res=600, units = "in", width = 10, height = 5)
layout(matrix(c(1,2), nrow = 1))

par(mai = c(1, 1, 0.2, 0.2))

plot(sites.spdf$mapH, sites.spdf$d2H, pch = 21, bg = "grey", 
     xlab = expression("Lake isoscape "*delta^{2}*"H (\u2030)"),
     ylab = expression("Mean annual precipitation  "*delta^{2}*"H (\u2030)"))
abline(hmod)
points(sites.spdf$mapH, sites.spdf$d2H, pch = 21, bg = "grey")
xl = max(sites.spdf$mapH)
yl = min(sites.spdf$d2H) + 0.05 * diff(range(sites.spdf$d2H))
text(xl, yl, equation(hmod), pos = 2)

plot(sites.spdf$mapO, sites.spdf$d18O, pch = 21, bg = "grey", 
     xlab = expression("Lake isoscape "*delta^{18}*"O (\u2030)"),
     ylab = expression("Mean annual precipitation  "*delta^{18}*"O (\u2030)"))
abline(omod)
points(sites.spdf$mapO, sites.spdf$d18O, pch = 21, bg = "grey")
xl = max(sites.spdf$mapO)
yl = min(sites.spdf$d18O) + 0.05 * diff(range(sites.spdf$d18O))
text(xl, yl, equation(omod), pos = 2)

dev.off()