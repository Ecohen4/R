################# Start Q3

#Repeat 1. by fitting a GLM with appropriate link function.

#Fit the GLM - using Gamma PDF and appropriate canonical link function.
zz=glm(Y ~ ., data=data,family=Gamma(link="inverse"))
summary(zz)  ## you can issue this command to see the model

## Now we wish to estimate/predict at the points on the DEM
predpoints = read.table("http://cires.colorado.edu/~aslater/CVEN_6833/colo_dem.dat",sep=",")

Xpred = cbind(predpoints$V1,predpoints$V2,predpoints$V3,predpoints$V4)

## Now Xpred should be renamed as X (the same name of the data frame that went into fitting the GLm model)

X=as.data.frame(Xpred)

# now predict
ypred = predict.glm(zz,newdata=X, type="response",se.fit=TRUE)

#We now wish to do a spatial plot - Latitude, Longitude and the predicted value..

# since the predpoints are on a unifrom spatial grid but not on 
# rectangular grid we can do the following..

library(akima)  	# we need this library
zz = interp(X[,2],X[,1],ypred$fit)

library(fields)
surface(zz,xlab="Longitude",ylab ="Latitude")

# you can do the same for the ypred$se.fit  - the error of the predictions.

## you can also use the image command and add contours on top using the contour command

# Fields package does spline smoothing as well - another nonparametric method. You can give it a tri
## at the end just for comparison.
X = cbind(test$V3,test$V2,test$V4,test$V18)	# put longitude first and then latitude
zspline = Tps(X,Y)		# fit the spline on the original data..
surface(zspline,xlab="Longitude",ylab ="Latitude")

# note that the surface command in fields plots a 4-d model - i.e., 4 pred variables.
# It takes the first two to do the plotting.

# you can also use the predict command, se.fit etc. etc. - please see help.

#################################### subset selection ##############
library(leaps)		# to provide combinations
library(MPV)		# to help estimate PRESS and consequently, GCV

X = cbind(test$V3,test$V2,test$V4,test$V18)		# all the predictor set.
#X = as.data.frame(X)

N = length(Y)

combs = leaps(X,Y, nbest=25)     #  GEt upto 25 combinations for each
# number of predictors
combos = combs$which

ncombos = length(combos[,1])

xpress=1:ncombos
xmse = 1:ncombos

for(i in 1:ncombos){
  xx = X[,combos[i,]]
  xx=as.data.frame(xx)
  
  zz=glm(Y ~ ., data=xx,family=Gamma(link="inverse"))
  
  xpress[i]=PRESS(zz)
  
  xmse[i] = sum((zz$res)^2) / (N - length(zz$coef))
  
}

## you can compute GCV using the PRESS values

# best model based on MSE
zc=order(xmse)[1]   #best model for PRESS
bestmodelm = lsfit(X[,combos[zc,]], Y)

#best model based on PRESS
zc=order(xpress)[1]   #best model for PRESS
bestmodelp = lsfit(X[,combos[zc,]], Y)

## You can also compute GCV

# using AIC
X1 = as.data.frame(X)
## Fit the model with all the predictors
zz=glm(Y ~ ., data=X1,family=Gamma(link="inverse"))
zms=stepAIC(zz, scope=list(upper = ~., lower = ~1), trace=5)

#arguement trace =5 (or a high number spits out the details on the screen at each step)
