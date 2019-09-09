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

png("../map.png", res=600, units = "in", width = 9, height = 4)
layout(matrix(c(1,2,3), nrow = 1))

par(mai = c(0.5,0.5,0.1,0.9), cex = 0.8)

plot(d2H.r, xlim = c(-76, -70), legend = FALSE)
addRasterLegend(d2H.r, side = 4, location = c(-69.6, -69.3, -43.5, -38.5), cex = 1)
axis(4, at = -41, labels = expression("Lake isoscape "*delta^{2}*"H (\u2030)"), tick = FALSE, line=3.8, cex.axis = 1.2)
lines(bnd, col="dark grey")
points(l.spdf, col = "blue")
text(-75.2, -36.2, "a", pos = 2)
box()

plot(d18O.r, xlim = c(-76, -70), legend = FALSE)
addRasterLegend(d2H.r, side = 4, location = c(-69.6, -69.3, -43.5, -38.5), cex = 1)
axis(4, at = -41, labels = expression("Lake isoscape "*delta^{18}*"O (\u2030)"), tick = FALSE, line=3.8, cex.axis = 1.2)
lines(bnd, col="dark grey")
points(sites.spdf, col = "blue")
text(-75.2, -36.2, "b", pos = 2)
box()

plot(D.r, xlim = c(-76, -70), legend = FALSE)
addRasterLegend(d2H.r, side = 4, location = c(-69.6, -69.3, -43.5, -38.5), cex = 1)
axis(4, at = -41, labels = "Lake isoscape D-excess (\u2030)", tick = FALSE, line=3.8, cex.axis = 1.2)
lines(bnd, col="dark grey")
text(-75.2, -36.2, "c", pos = 2)
box()

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


addRasterLegend <- function(r, direction, side, location = 'right', nTicks = 2, shortFrac = 0.02, longFrac = 0.3, axisOffset = 0, border = TRUE, ramp = "terrain", isInteger = 'auto', ncolors = 64, breaks = NULL, minmax = NULL, locs = NULL, cex.axis = 0.8, labelDist = 0.7, digits = 2, ...) {
  
  if (class(r) != 'RasterLayer') {
    stop("r must be a RasterLayer.");
  }
  
  if(!hasArg('direction')) {
    direction <- 'auto'
  }
  
  if (!direction %in% c('auto', 'vertical', 'horizontal')) {
    stop("direction must be auto, vertical or horizontal.");
  }
  
  if (is.character(location)) {
    if (!location %in% c('bottomleft','bottomright','topleft','topright','bottom','top','left','right')) {
      stop('location is not recognized.');
    }
  }
  
  if (!isInteger %in% c('auto', TRUE, FALSE)) {
    stop('isInteger must be "auto", TRUE or FALSE.')
  }
  
  if (length(ramp) == 1) {
    if (ramp == 'terrain') {
      pal <- rev(terrain.colors(ncolors));
    }
  } else if (length(ramp) > 1) {
    pal <- colorRampPalette(ramp)(ncolors);
  }
  
  #if minmax provided, use to generate linear color breaks
  if (is.null(minmax)) {
    colorbreaks <- seq(minValue(r), maxValue(r), length.out = (ncolors+1));
  } else {
    colorbreaks <- seq(minmax[1], minmax[2], length.out = (ncolors+1));
  }
  
  #if supplied, use custom set of breaks
  if (!is.null(breaks)) {
    colorbreaks <- breaks;
    ncolors <- length(breaks) + 1;
  }
  
  n <- length(colorbreaks);
  
  #return plot region extremes and define outer coordinates
  minX <- grconvertX(par('fig')[1], from = 'ndc', to = 'user') 
  maxX <- grconvertX(par('fig')[2], from = 'ndc', to = 'user')
  minY <- grconvertY(par('fig')[3], from = 'ndc', to = 'user')
  maxY <- grconvertY(par('fig')[4], from = 'ndc', to = 'user')
  
  xrange <- maxX - minX
  yrange <- maxY - minY
  minX <- minX + xrange * 0.05
  maxX <- maxX - xrange * 0.05
  minY <- minY + yrange * 0.05
  maxY <- maxY - yrange * 0.05
  
  if (is.character(location)) {
    
    if (location == 'topleft' & direction %in% c('auto', 'vertical')) {
      location <- vector('numeric', length = 4);
      location[1] <- minX
      location[2] <- minX + (maxX - minX) * shortFrac
      location[3] <- maxY - (maxY - minY) * longFrac
      location[4] <- maxY
    } else
      
      if (location == 'topleft' & direction == 'horizontal') {
        location <- vector('numeric', length = 4);
        location[1] <- minX
        location[2] <- minX + (maxX - minX) * longFrac
        location[3] <- maxY - (maxY - minY) * shortFrac
        location[4] <- maxY
      } else
        
        if (location == 'topright' & direction %in% c('auto', 'vertical')) {
          location <- vector('numeric', length = 4);
          location[1] <- maxX - (maxX - minX) * shortFrac
          location[2] <- maxX
          location[3] <- maxY - (maxY - minY) * longFrac
          location[4] <- maxY
        } else
          
          if (location == 'topright' & direction == 'horizontal') {
            location <- vector('numeric', length = 4);
            location[1] <- maxX - (maxX - minX) * longFrac
            location[2] <- maxX
            location[3] <- maxY - (maxY - minY) * shortFrac
            location[4] <- maxY
          } else
            
            if (location == 'bottomleft' & direction %in% c('auto', 'vertical')) {
              location <- vector('numeric', length = 4);
              location[1] <- minX
              location[2] <- minX + (maxX - minX) * shortFrac
              location[3] <- minY
              location[4] <- minY + (maxY - minY) * longFrac
            } else
              
              if (location == 'bottomleft' & direction == 'horizontal') {
                location <- vector('numeric', length = 4);
                location[1] <- minX
                location[2] <- minX + (maxX - minX) * longFrac
                location[3] <- minY
                location[4] <- minY + (maxY - minY) * shortFrac
              } else
                
                if (location == 'bottomright' & direction %in% c('auto', 'vertical')) {
                  location <- vector('numeric', length = 4);
                  location[1] <- maxX - (maxX - minX) * shortFrac
                  location[2] <- maxX
                  location[3] <- minY
                  location[4] <- minY + (maxY - minY) * longFrac
                } else
                  
                  if (location == 'bottomright' & direction == 'horizontal') {
                    location <- vector('numeric', length = 4);
                    location[1] <- maxX - (maxX - minX) * longFrac
                    location[2] <- maxX
                    location[3] <- minY
                    location[4] <- minY + (maxY - minY) * shortFrac 
                  } else
                    
                    if (location == 'left') {
                      location <- vector('numeric', length = 4);
                      location[1] <- minX
                      location[2] <- minX + (maxX - minX) * shortFrac
                      location[3] <- mean(par('usr')[3:4]) - ((maxY - minY) * longFrac)/2
                      location[4] <- mean(par('usr')[3:4]) + ((maxY - minY) * longFrac)/2
                      direction <- 'vertical'
                    } else
                      
                      if (location == 'right') {
                        location <- vector('numeric', length = 4);
                        location[1] <- maxX - (maxX - minX) * shortFrac
                        location[2] <- maxX
                        location[3] <- mean(par('usr')[3:4]) - ((maxY - minY) * longFrac)/2
                        location[4] <- mean(par('usr')[3:4]) + ((maxY - minY) * longFrac)/2
                        direction <- 'vertical'
                      } else
                        
                        if (location == 'top') {
                          location <- vector('numeric', length = 4);
                          location[1] <- mean(par('usr')[1:2]) - ((maxX - minX) * longFrac)/2
                          location[2] <- mean(par('usr')[1:2]) + ((maxX - minX) * longFrac)/2
                          location[3] <- maxY - (maxY - minY) * shortFrac
                          location[4] <- maxY
                          direction <- 'horizontal'
                        } else
                          
                          if (location == 'bottom') {
                            location <- vector('numeric', length = 4);
                            location[1] <- mean(par('usr')[1:2]) - ((maxX - minX) * longFrac)/2
                            location[2] <- mean(par('usr')[1:2]) + ((maxX - minX) * longFrac)/2
                            location[3] <- minY
                            location[4] <- minY + (maxY - minY) * shortFrac
                            direction <- 'horizontal'
                          }
  }
  
  # infer direction based on dimensions of legend box
  if (direction == 'auto') {
    if (((location[2] - location[1]) / (par('usr')[2] - par('usr')[1])) >= ((location[4] - location[3]) / (par('usr')[4] - par('usr')[3]))) {
      direction <- 'horizontal';
    } else {
      direction <- 'vertical';
    }
  }
  
  if (direction == 'horizontal') {
    axisOffset <- axisOffset * (par('usr')[4] - par('usr')[3]);
  } else if (direction == 'vertical') {
    axisOffset <- axisOffset * (par('usr')[2] - par('usr')[1]);
  }
  
  #determine side for labels based on location in plot and direction
  if (!hasArg('side')) {
    if (direction == 'vertical') { #side = 1 or 4
      if (mean(location[1:2]) <= mean(par('usr')[1:2])) {
        side <- 4;
      } else {
        side <- 2;
      }
    }
    if (direction == 'horizontal') { #side = 2 or 3
      if (mean(location[3:4]) > mean(par('usr')[3:4])) {
        side <- 1;
      } else {
        side <- 3;
      }
    }
  }
  
  if (direction == 'horizontal') {
    x <- seq(from = location[1], to = location[2], length.out = n);
    width <- location[3:4];
  } else {
    x <- seq(from = location[3], to = location[4], length.out = n);
    width <- location[1:2];
  }
  
  #get bin coordinates
  x <- rep(x,each = 2);
  x <- x[-c(1,length(x))];
  x <- matrix(x, ncol = 2, byrow = TRUE);
  
  #find tick locations
  #get equivalent color bins
  z <- rep(colorbreaks,each = 2);
  z <- z[-c(1,length(z))];
  z <- matrix(z, ncol = 2, byrow = TRUE);
  
  #if tick locations are supplied, use them, otherwise generate regularly spaced tick locations
  if (!is.null(locs)) {
    # if max value is included, add in manually
    tol <- 1e-10
    maxValueIncluded <- FALSE
    if (any(abs(maxValue(r) - locs) < tol)) {
      locs <- locs[which(abs(maxValue(r) - locs) >= tol)]
      maxValueIncluded <- TRUE
    }
    tickLocs <- x[sapply(locs, function(x) which((x >= z[,1] & x < z[,2]) == TRUE)),1];
    if (maxValueIncluded) {
      tickLocs <- c(tickLocs, max(x[,2]));
      locs <- c(locs, max(colorbreaks))
    }
    tx <- locs;
  } else {
    tx <- trunc(seq(from = 1, to = nrow(x), length.out = nTicks + 2));
    tickLocs <- x[tx,1];
    tx <- z[tx,1];
    tickLocs[length(tickLocs)] <- max(x[,2]);
    tx[length(tx)] <- max(z[,2]);
  }
  
  #if raster values are integer, then make legend have integer values
  if (isInteger == 'auto') {
    randomSample <- sample(values(r)[which(!is.na(values(r)))], size=1000, replace = TRUE)
    if (identical(randomSample, trunc(randomSample))) {
      tx <- round(tx, 0)
    }
  } else if (isInteger) {
    tx <- round(tx, 0)
  }
  
  #plot bar
  if (direction == 'horizontal') {
    rect(xleft = x[,1], ybottom = width[1], xright = x[,2], ytop = width[2], border = pal, col = pal, xpd = NA);
  } else {
    rect(xleft = width[1], ybottom = x[,1], xright = width[2], ytop = x[,2], border = pal, col = pal, xpd = NA);
  }
  
  if (border) {
    rect(location[1], location[3], location[2], location[4], border='black', xpd = NA);
  }
  
  #add tickmarks
  if (side == 1) { #bottom
    axis(side, at = tickLocs, pos = location[3] - axisOffset, labels = round(tx, digits), xpd = NA, las = 1, cex.axis = cex.axis, mgp = c(3, labelDist, 0), ...);
  } 
  if (side == 3) { #top
    axis(side, at = tickLocs, pos = location[4] + axisOffset, labels = round(tx, digits), xpd = NA, las = 1, cex.axis = cex.axis, mgp = c(3, labelDist, 0), ...);
  }
  if (side == 2) { #left
    axis(side, at = tickLocs, pos = location[1] - axisOffset, labels = round(tx, digits), xpd = NA, las = 1, cex.axis = cex.axis, mgp = c(3, labelDist, 0), ...);
  }
  if (side == 4) { #right
    axis(side, at = tickLocs, pos = location[2] + axisOffset, labels = round(tx, digits), xpd = NA, las = 1, cex.axis = cex.axis, mgp = c(3, labelDist, 0), ...);
  }
}

