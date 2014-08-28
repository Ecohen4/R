setwd("/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013")

library(fields)
library(MASS)

# Define custom functions...
#### Cross-Validated and Fitted Estimates ####
crossval.Krig<-function(X, Y, rho, theta, sigma){
  
  # fitted values... 
  X<-xs # assign subset of predictor variables (lon, lat) to X
  # Kriging is an "exact esimator" so this returns the observed values...
  zz=Krig(X,Y,rho=rho,theta=theta,sigma2=sigma,m=1)
  Fitted = predict.Krig(zz,x=X,drop.Z=TRUE) # predict at the same lat-long coords  
  par(mfrow=c(1,1))
  plot(Y, Fitted, pch=20, col="black", xlab="Observed Precip (mm)",ylab="Modeled Precip (mm)", main="Observed vs Modeled Precipitation")
  abline(a=0, b=1)
  
  # Now check the cross validated estimates (drop observations one at a time, refit the model to the remaining data (N-1) and predict at the dropped point; repeat for all observations)
  par(mfrow=c(1,1))
  n=length(Y) # number of observations
  yest=1:n # number of estimates to be made
  nvar=dim(X)[2] # number of regressor variables (lon, lat)
  nobs=length(Y) # Number of observations
  nobs1=nobs-1 # number of observations less the one being cross-validated
  yhatcross=1:nobs # define structure for assigning yhat (yest)
  ksecross=1:nobs # define structure for assigning kse
  yest=1:nobs
  se=1:nobs
  
  index=1:n
  for(i in 1:n){
    index1=index[index != i] #drop one observation at a time
    xval=X[index1,] #X data less the dropped observation
    yval=Y[index1]  #Y data less the dropped observation
    xp<-X[i,] # the dropped point becomes the predictant
    
    # fit the model without the dropped observation...
    
    # Doug's method....
    #K1<- Matern( rdist( xval, xval)/.3, smoothness=1.5) # 490x490
    #Kstar1<- Matern(rdist( xp, xval)/.3, smoothness=1.5) # 1x490
    
    # Balaji's method...
    K <- sigma^2 + rho*(1 - exp(-1*rdist(xval,xval)/theta)) #dim(K)=490 x 490
    Kstar <- sigma^2 + rho*(1 - exp(-1*rdist(as.matrix(xp),xval)/theta)) # 1 x 490
    ghat<- Kstar %*% solve( K + diag( sigma^2, nobs1)) %*% yval # 1x1 --> 1 number
    kse <- sqrt(var(yval) - diag((Kstar) %*% solve( K + diag( sigma^2, nobs1)) %*% t(Kstar))) # original from Balaji's code... # 1 number
    yhatcross[i] <- ghat
    ksecross[i] <- kse
    
    # Compare with the Krig package and predict function....
    zz=Krig(xval,yval,rho=rho,theta=theta,sigma2=sigma,m=1) #fit the model to all observations less the dropped point.
    yest[i] = predict.Krig(zz,x=as.matrix(xp),drop.Z=TRUE) #predict at the dropped point....
    ## Compute the standard error
    # se[i] = predict.se(zz,x=as.matrix(xp))
    
  }  ### close for-loop
  
  print(cbind(Y,yhatcross, yest)) # compare Y observed and y cross-validated from two different methods.
  # cbind(ksecross, se) # compare the standard erros from the two methods.
  
  #now surface the x-validated estimates..
  points(Y, yest, col="green", pch=20) #green filled-in bullets
  points(Y, yhatcross, col="blue", pch=0) #blue squares
  legend("topleft", legend=c("Krig Package", "First Principles"), col=c("green","blue"), pch=c(20,0), bty="n")
  
}  ############# Close function

### Doug's Lecture Notes 1 through 4
### Chapter 8 of Springer online book
#### Applied Spatial Data Analysis with R
#Authors:
#    Roger S. Bivand,
#    Edzer Pebesma,
#    Virgilio GÃ³mez-Rubio 

### Three Models 
#### Fit a Kriging (i.e. spatial) Model to the precipitation data
####  Predict on the DEM grid

### Model 2  ### Spatial Additive Model - Linear term for the mean function
### plus the spatial model on the residual

#### Model 3 - Fit a GLM to lat, lon and elevation and fit a spatial model to
#### the residuals

### Model 2 estimates the mean component and residual spatial model parameters
### simultaneously, while Model 3 does it in two steps - i.e. Hierarchy!.

#### Spatial trend, mean function all refer to the same..

#### Kriging is an exact estimator - i.e., at the observation locations Yhat will
### be exactly equal to Y without nugget effect, but with one it will be close.

### Therefore you have to do cross-validation to evaluate the performance.

## You can write a loop to perform cross validation for all these models

### Problem 6 - Model 1

### Problem 7 - Models 2 and 3

##### Model 1 - Lectures 2B and 3 
## Kriging on the precipitation..
#Import data: Mean annual precipitation at 491 locations in Colorado based on data for the period 1980-2002
data<-read.table(file="/Users/elliotcohen/Dropbox/Data/Climate/Rainfall/Colorado_Annual_Precip_Gridded.csv", sep=",", header=TRUE)
data<-as.data.frame(data)

# data$Lat<-scale(data$Lat)
# data$Long<-scale(data$Long)
# data$Elev<-scale(data$Elev)
# data$Precip<-scale(data$Precip)

X<-data[,-4]  # df of Independent variable
Y<-data$Precip #vector of dependent variable (annual rainfall at lat-long position)
n<-length(Y)
nvar<-dim(X)[2]

# colopcp=read.table("http://cires.colorado.edu/~aslater/CVEN_6833/colo_precip.dat")
# lat=colopcp$V1
# lon=colopcp$V2
# elev=colopcp$V3
# precip=colopcp$V4
# vars = as.data.frame(cbind(lon,lat,elev)) #df of independent variables
# varsdata<-vars
# nobs<-length(precip)

fitall<-glm(Y~., data=X)
xs<-X[,1:2] # lat, long
xse<-X[,] # lat, long and elev

### Also If Kriging with elevation then replace xps with xpse below

## Compute empirical variogram and plot ...
par(mfrow=c(1,1))
par(mar=c(4,4,1,1))
look<- vgram(xs, Y, N=15, lon.lat=FALSE)
bplot.xy( look$d, look$vgram, breaks= look$breaks,
          outline=FALSE,
          xlab="distance (degrees)", ylab="variogram")
points( look$centers, look$stats[2,], col="blue")

# fit of exponential by nonlinear least squares
xd<- look$centers
ym<- look$stats[2,]
sigma<- 1.0 # 
nls.fit<- nls( ym ~ sigma^2 + rho*( 1- exp(-xd/theta)),
               start= list(rho=max(look$stats[2,]), theta=quantile(xd,0.1)), control=list( maxiter=200) )
pars<- coefficients(nls.fit)
rho<-pars[1]
theta<-pars[2]
xr = round(max(xd)+sd(xd))
dgrid<- seq(0,xr,,400)
lines(dgrid, sigma^2 + rho*(1 - exp(-1*dgrid/theta)), col="blue", lwd=3)

### Predict at observation points.. is sigma = 0 then it will be exact.
K = sigma^2 + rho*(1 - exp(-1*rdist(xs,xs)/theta))
#Kstar<- Matern(rdist( xgrid, x)/.3, smoothness=1.5)
Kstar = K
nobs=length(Y)
ghat<- Kstar %*% solve( K + diag( sigma^2, nobs)) %*% Y

## Check with Krig and predict.krig
zz=Krig(xs,Y,rho=rho,theta=theta,m=1,sigma2=sigma)
yp = predict.Krig(zz)

cbind(ghat,yp) #visually inspect the two approaches above --> highly similar

## Predict on the grid.. and standard error..
#predpts=read.table("http://cires.colorado.edu/~aslater/CVEN_6833/colo_dem.dat")
predpts<-read.table(file="/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013/colo_dem.txt")
names(predpts)<-c("Lat","Long","Elev") #gridded lat, long and elev
xps<-predpts[,1:2] #lat, long
xpse<-predpts[,] #lat, long, elev

Kstar = sigma^2 + rho*(1 - exp(-1*rdist(xps,xs)/theta)) #xps and xs can have different dimensions (differing nrow)
ghat<- Kstar %*% solve( K + diag( sigma^2, nobs)) %*% Y
kse = sqrt(var(Y) - diag((Kstar) %*% solve( K + diag( sigma^2, nobs)) %*% t(Kstar)))

## Check with Krig and predict.Krig
zz=Krig(xs,Y,rho=rho,theta=theta,sigma2=sigma,m=1)
yp = predict.Krig(zz,x=xps,drop.Z=TRUE)

### spatial map of estimates and errors
par(mar=c(4,3,4,2)+0.1)

xlon = sort(unique(predpts$Long))
nrs = length(xlon)
ylat = sort(unique(predpts$Lat))
ncs = length(ylat)
zmat1a = matrix(ghat, nrow=ncs, ncol= nrs) #zmat=predicted precip from computation
zmat1 = matrix(yp, nrow=ncs, ncol=nrs) #zmat=predicted precip from kriging package

image.plot(xlon,ylat,t(zmat1),zlim=range(Y), col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Precipitation (mm) predicted on DEM Grid\nModel 1 w. Kriging package")
world(add=TRUE, lwd=4)
contour(xlon,ylat,t(zmat1),add=T) #contour plot of model 1 precip
dev.copy(png,"Precipitation (mm) predicted on DEM Grid - Model 1 w. Kriging package.png")
dev.off()

image.plot(xlon,ylat,t(zmat1a),zlim=range(Y), col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Precipitation (mm) predicted on DEM\nModel 1 w. manual computation of ghat")
world(add=TRUE, lwd=4)
contour(xlon,ylat,t(zmat1a),add=T) #contour plot of model 1 precip
dev.copy(png,"Precipitation (mm) predicted on DEM - Model 1 w. manual computation.png")
dev.off()

zmat1b = matrix(kse, nrow=ncs, ncol= nrs)  # zmat = standard errors on the DEM
image.plot(xlon,ylat,t(zmat1b), col=topo.colors(n=20, alpha=0.5), ylab="latitude", xlab="longitude", main="Model 1 errors")
contour(xlon,ylat,t(zmat1b),add=T)
dev.copy(png,"Model 1 errors.png")
dev.off()


## Cross-validated Y estimates...
crossval.Krig(X=xs, Y=Y, rho=rho, theta=theta, sigma=sigma)
dev.copy(png,paste(zz$call[1],"Cross-Validation.png", sep="_"))
dev.off()

#############################################################
### Model 2  - Lecture 3
#### Spatial trend with a linear model + Kriging residuals #####
### Estimated together..

#xx = cbind(rep(1,nobs),colopcp$V2,colopcp$V1)
#xx = cbind(rep(1,nobs),colopcp$V2,colopcp$V1,colopcp$V3)
#xe = cbind(colopcp$V2,colopcp$V1,colopcp$V3)
# Why is a column of 1's added to some example df's, but not others????

xx<-X[,1:2] # lat, long
xe<-X[,] # lat, long and elev

## select the desired xx
xs<-xx

nobs = length(Y)
sigma=1

par(mfrow=c(1,1))
look<- vgram(xs, Y, N=15, lon.lat=FALSE)
bplot.xy(look$d, look$vgram, breaks= look$breaks,
         outline=FALSE,
         xlab="distance (degrees)", ylab="variogram")
points( look$centers, look$stats[2,], col="blue")

# fit of exponential by nonlinear least squares
xd<- look$centers
ym<- look$stats[2,]
sigma<- 1.0 # 
nls.fit<- nls( ym ~ sigma^2 + rho*( 1- exp(-xd/theta)),
               start= list(rho=max(look$stats[2,]), theta=quantile(xd,0.1)), control=list( maxiter=200) )
pars<- coefficients(nls.fit)
rho<-pars[1]
theta<-pars[2]
xr=round(max(xd)+sd(xd))
dgrid<- seq(0,xr,,400)
lines(dgrid, sigma^2 + rho*(1 - exp(-1*dgrid/theta)), col="blue", lwd=3)


K = sigma^2 + rho*(1 - exp(-1*rdist(xs,xs)/theta))
MM = K + diag(sigma^2, nobs)
xx<-as.matrix(xx)
dhat = solve(t(xx) %*% solve(MM) %*% xx) %*% t(xx) %*% solve(MM) %*% Y

Kstar = K
nobs=length(Y)
ghat<- xx%*%dhat + Kstar %*% solve( K + diag( sigma^2, nobs)) %*% (Y-xx%*%dhat)
# Note: ghat is computed differently for Models 1 and 2.  First term of ghat in model 2 (xx%*%dhat) is *not* present in model 1, and last term (Y-xx%*%dhta) is simply (Y).


## Predict on the grid.. and standard error..
npred<-dim(predpts)[1] # number of prediction nodes
#xps = cbind(lon,lat)
#xp=cbind(rep(1,npred),lon,lat,elev)
xp<-cbind(rep(1,npred), predpts) # augmented X matrix
xp<-as.matrix(xp)
Kstar = sigma^2 + rho*(1 - exp(-1*rdist(xps,xs)/theta)) # 

#ghat<- xp%*%dhat + Kstar %*% solve( K + diag( sigma^2, nobs)) %*% (Y-xx%*%dhat)
# problem computing ghat... xp %*% dhat : non-conformable arguments
# try xps instead of xp...
xps<-as.matrix(xps)
ghat<- xps%*%dhat + Kstar %*% solve( K + diag( sigma^2, nobs)) %*% (Y-xx%*%dhat)

kse = sqrt(var(Y) - diag((Kstar) %*% solve( K + diag( sigma^2, nobs)) %*% t(Kstar)))


## Check with Krig and predict.Krig
zz=Krig(xe[,1:2],Y,Z=xe[,3],rho=rho,theta=theta,sigma2=sigma,m=2)
yp = predict.Krig(zz,x=predpts[,1:2],Z=predpts[,3],drop.Z=FALSE) #similar to ghat..

### spatial map of estimates and errors - modify to add beautify..
xlon = sort(unique(predpts$Long))
nrs = length(xlon)
ylat = sort(unique(predpts$Lat))
ncs = length(ylat)

zmat2a = matrix(ghat, nrow=ncs, ncol= nrs) # matrix of estimated precip on a rectangular grid of width longitude, and height latitude   
zmat2<-matrix(yp, nrow=ncs, ncol=nrs) #using output from krig package

image.plot(xlon,ylat,t(zmat2a),zlim=range(Y), col=topo.colors(20,0.5), ylab="latitude", xlab="longitude", main="Precipitation (mm) predicted on DEM grid\nModel 2: Kriging + Linear Model")
contour(xlon,ylat,t(zmat2a),add=T)
dev.copy(png,paste(zz$call[1],"Model_2_DEM_Predictions.png", sep="_"))
dev.off()

zmat2b = matrix(kse, nrow=ncs, ncol= nrs) 
image.plot(xlon,ylat,t(zmat2b), col=topo.colors(20,0.5), ylab="latitude", xlab="longitude", main="Model 2 standard errors")
contour(xlon,ylat,t(zmat2b),add=T)
dev.copy(png,paste(zz$call[1],"Model_2_Standard_Errors.png", sep="_"))
dev.off()

# Cross-validated Y estimates....
crossval.Krig(X=xs, Y=Y, rho=rho, theta=theta, sigma=sigma)
dev.copy(png,paste(zz$call[1],"Cross-Validation_Model_2.png", sep="_"))
dev.off()

############## Model 3
############## GLM to get the mean function + Kriging the residuals
#######################

### fit GLM and get residuals..
zglm = glm(Y ~ ., data=X)
Yresid = residuals(zglm) # let Y be the residuals of the glm model....

#xs = cbind(colopcp$V2,colopcp$V1)
#Y = yres # set Y to the residuals of the glm model....
nobs = length(Yresid)

sigma=1
look<- vgram(xs, Yresid, N=15, lon.lat=FALSE)
bplot.xy(look$d, look$vgram, breaks= look$breaks,
         outline=FALSE,
         xlab="distance (degrees)", ylab="variogram")
points( look$centers, look$stats[2,], col="blue")

# fit of exponential by nonlinear least squares
xd<- look$centers
ym<- look$stats[2,]
sigma<- 1.0 # 
nls.fit<- nls( ym ~ sigma^2 + rho*( 1- exp(-xd/theta)),
               start= list(rho=max(look$stats[2,]), theta=quantile(xd,0.1)), control=list( maxiter=200) )
pars<- coefficients(nls.fit)
rho<-pars[1]
theta<-pars[2]
xr=round(max(xd)+sd(xd))
dgrid<- seq(0,xr,,400)
lines(dgrid, sigma^2 + rho*(1 - exp(-1*dgrid/theta)), col="blue", lwd=3)


### Predict at observation points.. is sigma = 0 then it will be exact.
K = sigma^2 + rho*(1 - exp(-1*rdist(xs,xs)/theta))
#Kstar<- Matern(rdist( xgrid, x)/.3, smoothness=1.5)
Kstar = K
nobs=length(Yresid)
ghat<- Kstar %*% solve( K + diag( sigma^2, nobs)) %*% Yresid
ghat = predict(zglm) + ghat

## Check with Krig and predict.krig
zz=Krig(xs,Yresid,rho=rho,theta=theta,m=1,sigma2=sigma)
yp = predict(zglm) + predict.Krig(zz)

## Predict on the grid.. and standard error..
# xe = cbind(lon,lat,elev)
# xps=cbind(lon,lat) 

Kstar = sigma^2 + rho*(1 - exp(-1*rdist(xps,xs)/theta))
ghat<- Kstar %*% solve( K + diag( sigma^2, nobs)) %*% Yresid
kse = sqrt(var(Y) - diag((Kstar) %*% solve( K + diag( sigma^2, nobs)) %*% t(Kstar)))
ghat = ghat + predict(zglm,newdata=predpts)

## Check with Krig and predict.Krig
zz=Krig(xs,Yresid,rho=rho,theta=theta,sigma2=sigma,m=1)
yp = predict(zglm,newdata=predpts)+ predict.Krig(zz,x=xps,drop.Z=TRUE)

## Now check the cross-validated estimates of Y...
crossval.Krig(X=xs, Y=Y, rho=rho, theta=theta, sigma=sigma)
dev.copy(png,paste(zz$call[1],"Cross-Validation_Model_3.png", sep="_"))
dev.off()

### spatial map of estimates and errors
# xlon = sort(unique(data$Lon))
# nrs = length(xlon)
# ylat = sort(unique(data$Lat))
# ncs = length(ylat)
zmat3a = matrix(ghat, nrow=ncs, ncol= nrs) 
image.plot(xlon,ylat,t(zmat3a),zlim=range(Y),col=topo.colors(20,0.5), ylab="latitude", xlab="longitude", main="Precipitation (mm) predicted on DEM grid\nModel 3: Kriging + GLM")
contour(xlon,ylat,t(zmat3a),add=T)
dev.copy(png,paste(zz$call[1],"Model_3_DEM_Predictions.png", sep="_"))
dev.off()

zmat3b = matrix(kse, nrow=ncs, ncol= nrs) 
image.plot(xlon,ylat,t(zmat3b), col=topo.colors(20,0.5), ylab="latitude", xlab="longitude", main="Model 3 standard errors")
contour(xlon,ylat,t(zmat3b),add=T)
dev.copy(png,paste(zz$call[1],"Model_3_Standard_Errors.png", sep="_"))
dev.off()


###############################################################
##############################################################
#visually compare all the Kriging precip models....
#v. Make a 3-D plot (or spatial colored/contour map) of model estimates (i.e. latitude, longitude and the model estimates) and model error. 
# (v) Spatial map of actual precipitation, model estimates, and model error
plot.new()
par(mfrow=c(2,2))
par(mar=c(4, 3, 3, 2) + 0.1) # numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot
library(akima)
library(fields)

cex.main=0.9
zz0<-interp(x=data$Long, y=data$Lat, z=data$Precip, duplicate="mean")
image.plot(zz0, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Observed Precipitation (mm)")
contour(zz0, add=T)
world(add=TRUE, lwd=4)

image.plot(xlon,ylat,t(zmat1a), col=topo.colors(20,0.5), ylab="latitude", xlab="longitude", main="Model 1 precip predictions (mm)\nKriging")
contour(xlon,ylat,t(zmat1a),add=T)

image.plot(xlon,ylat,t(zmat2a), col=topo.colors(20,0.5), ylab="latitude", xlab="longitude", main="Model 2 precip predictions (mm)\nKriging + Mean Trend")
contour(xlon,ylat,t(zmat1a),add=T)

image.plot(xlon,ylat,t(zmat3a), col=topo.colors(20,0.5), ylab="latitude", xlab="longitude", main="Model 3 precip predictions (mm)\nKriging + GLM")
contour(xlon,ylat,t(zmat1a),add=T)

dev.copy(png,"Observed vs. Modeled Precip - 3 models.png")
dev.off()
######
###### Conclusions from Q6 and Q7...
# As expected, Kriging is effective at caputuring localized effects, and provides  much more fidelty to actual observed precipitation patterns than do linear models or even local polynomial methods (from HW questions 1, 2, 4 and 5.)
# By visual inspeciton alone, it is difficult to discern between the three spatial models... 
# Refer to the cross-validated estimates (above) for a quantitative assessment

