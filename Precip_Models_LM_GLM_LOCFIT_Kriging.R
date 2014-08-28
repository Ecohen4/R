# September 17, 2013
# CVEN6833 - Advance Data Analysis
################# 
################# 
################# Start Q1

setwd("/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013")
#Import data
#Mean annual precipitation at 491 locations in Colorado based on data for the period 1980-2002
data<-read.table(file="/Users/elliotcohen/Dropbox/Data/Climate/Rainfall/Colorado_Annual_Precip_Gridded.csv", sep=",", header=TRUE)

## create a data frame to be used in GLM
data<-as.data.frame(data)
X<-data[,-4]  # df of Independent variable
Y<-data$Precip #vector of dependent variable (annual rainfall at lat-long position)

# i. Fit a ‘best’ linear regression model (use one of the objective functions - GCV, AIC, SBC or PRESS; you can also try a couple of them to see any differences). This entails fitting the model with all possible combinations of covariates and selecting the model with the minimum objective function.
lm<-lm(Precip ~ Lat + Long + Elev, data=data)
summary(lm)

# Load the MASS package to use the stepAIC function to choose the "best model"
require(MASS)
bestlm<-stepAIC(lm)

#ii. Perform ANOVA (i.e. model significance) and model diagnostics (i.e., check the assumptions of the residuals – Normality, independence,homoskedasticity). 
summary(bestlm) #"bestlm" has longitutde as only predictor variable
#Note: other models have similar AIC values, indicating that "best" model may not be that much better than the next best model.

# Model diagnostics
par(mfrow=c(2,2))
plot(bestlm)

# Results: visual inspection of plot 1 finds that residuals increase for larger values of the response variable yhat, indicating that the model does *not* capture certains attributes of the data (perhaps non-linearities); Plot 2 (Q-Q plot) shows that residuals are *not* distributed normally at the upper tail. In summary, at least 2 of the model diagnostics suggest that the linear model is inadequate for this data set.

# can also use custom-built  model diagnostics...
##### Manual Model Diagnostics (6 visual checks) #####
par(mfrow=c(2,3)) 

#But here is how to do it manually...
# Get the residuals of the model..
modresid=residuals(bestlm)

Y=data$Precip
nX<-length(bestlm$coefficients)
X=subset(data, select=names(bestlm$coefficients)[2:nX])
Yhat=Y-modresid  #residuals = Y - Yestimate ==> Yestimate = Y - residuals
k=dim(X)[2]  	  #number of regressor variables
p=k+1              #number of model parameters 
n=length(Y)

# Compute ANOVA quantities for use down below
SST = sum((Y - mean(Y))^2)   #SST = Total corrected sum of squares, n-1 dof
#Yhat = Y - bestGlm$res 
#(Y - Yhat = residuals), Yhat is the modeled response of Y
SSR = sum((Yhat - mean(Y))^2)  #SSR = Regression Sum of Squares = sum[(yhati-ybar)^2], dof = k (number of predictor  variables)
SSE = sum((modresid)^2)   #SSE = Error Sum of Squares = sum[(yi-yhati)^2], dof = n-p 
MSR = SSR / ncol(X)            
MSE = SSE/(n - length(bestlm$coef))         

# Now start computing diagnostics and plotting them...
par(mfrow=c(2,3))

# (1) Check if residuals fit a normal distribution
qqnorm(modresid)
qqline(modresid)	

# (2-3) Plot the residuals vs X.  Check to make sure there is *no* apparent pattern.  Distribution of residuals should be random.
for(i in 1:k){
  plot(X[,i],modresid,xlab="X",ylab="residuals",main="Residuals vs. X[i]") 
}

# (4) Plot the residuals vs the model estimates of Y. 
#Check to make sure there is *no* apparent pattern or structure.  In other words, the distribution of the residuals should look random.
plot(Yhat,modresid,xlab="estiimate of Y", ylab="residuals",main="Residuals vs model estimates of Y")

# (5) Plot the autocorrelation function - to make sure the residuals are *not* related to each other.  
z1=acf(modresid,main="autocorrelation of residuals")

# (6) Cooks Distance - to make sure outliers do not exert undue influence on the regression model.

# Compute the Hat matrix
hatm= hatvalues(bestlm)

#studentized residuals - ri  - equation 12-42
ri = modresid/sqrt((1 - diag(hatm)) * MSE)

#Compute Cook's Distance Equation 12-44
Di = ri*ri * diag(hatm) / ((1-diag(hatm)) * length(bestlm$coef))
plot(Y, diag(Di), main="Cook's Distance")
#If Dis are greater than 1 then there are points that have undue influence on the fit.

# Comments on Model Diagnostics: 
# The Q-Q plot reveals non-normality at the upper tail
# Plot of residuals vs model estimates of Y reveals heteroscdasticity. 
#### END MODEL DIAGNOSTICS ####
#### END PART (ii) ####

# iii. Compute cross-validated and fitted estimates at each observation points and plot them against the observed values. This is to visually see how the model performs in a fitting and cross-validated mode.

## Observed vs modeled estimates of response variable Y
par(mfrow=c(1,1))
#preds<-data.frame(c(paste("data$",names(bestlm$coefficients)[2:nX],sep="")), quote=FALSE)  #trying to automate regressor variables to include in preds
preds<-data.frame(data$Long)
names(preds)<-c("Long")
Fitted<-predict(bestlm, newdata=preds, type="response") #predict precip using best linear model with Long as only predictor variable.

plot(data$Precip,Fitted,pch=20,col="black",xlab="Observed Precip (mm)",ylab="Modeled Precip (mm)", main="Observed vs Modeled Precipitation")
abline(a=0, b=1)

# (iii) Now compute cross validated estimates
par(mfrow=c(1,1))
n=length(Y)
yest=1:n
nvar=dim(X)[2]

index=1:n
for(i in 1:n){
  index1=index[index != i] #drop one observation at a time
  Xval=X[index1,] #X data less the dropped observation
  Yval=Y[index1]  #Y data less the dropped observation
  
  zz=lsfit(Xval,Yval) #re-fit the model without the dropped observation
  
  xpred=c(1,X[i,1:nvar]) #now estimate at the point that was dropped
  xpred<-as.numeric(xpred)
  yest[i]=sum(zz$coef * xpred)
}
#now surface the x-validated estimates..
points(Y, yest, col="blue", pch=20)

# Comments:  Cross-validated model estimates and estimates modeled on the full data set appear internally consistent to each other, but niether are accurate to observed precipitation values for larger values of precip (e.g. above 800 mm/yr).

# iv. Drop 10% of observations, fit the model (i.e., the ‘best’ model from i. above) to the rest of the data and predict the dropped points. Compute RMSE and R2 and show them as boxplots. 
# Drop some % of points, fit the model and predict the dropped points..
library(arules)
source("/Users/elliotcohen/Dropbox/Advance Data Analysis/R source files/myboxplot-stats.r")
source("/Users/elliotcohen/Dropbox/Advance Data Analysis/R source files/myboxplot.r")
nsim = 500
rmseskill=1:nsim
corskill=1:nsim  
N = length(Y)
N10 = round(0.10*N)  	#choose % of points to drop (e.g. 10%)
index=1:N

for(i in 1:nsim){
  drop=sample(c(1:N),N10)  #sample 10% of the row indices from 1:N at random
  keep=setdiff(index,drop)  #discard values at the intersection of index and drop (e.g. drop 10% of the data, keep the rest)
  x<-as.data.frame(X[keep,])
  y<-Y[keep]
  xx<-as.data.frame(x)
  zz<-lm(y~xx$Long,data=xx) #fit model to remaining data
  xpred<-X[drop,]
  xx<-as.data.frame(xpred) #re-assign xx to the dropped data
  yhat<-predict.lm(zz,newdata=xx) #predict at dropped points using model fit w.out those points
  rmseskill[i]<-sqrt(mean(Y[drop]-yhat)^2)
  corskill[i]<-cor(Y[drop],yhat)
}

par(mfrow=c(1,2))
zz=myboxplot(rmseskill, main="Simulated RMSE skill",plot=FALSE)
zz$names=rep("",length(zz$names))
z1=bxp(zz,xlab="",ylab="RMSE",cex=1.25)
rmse=mean(((modresid)/sd(Y))^2)
points(z1,rmse,col="red",cex=2,pch=19)  #add a point showing the true RMSE of the data.
title(main="RMSE skill")

zz=myboxplot(corskill, main="Simulated Correlation skill",plot=FALSE)
zz$names=rep("",length(zz$names))
z1=bxp(zz,xlab="",ylab="Cor",cex=1.25)
title(main="Cor skill")

#v. Make a 3-D plot (or spatial colored/contour map) of model estimates (i.e. latitude, longitude and the model estimates) and model error. 
# (v) Spatial map of actual precipitation, model estimates, and model error
plot.new()
par(mfrow=c(2,2))
library(akima)
library(fields)
zz0<-interp(x=data$Long, y=data$Lat, z=data$Precip, duplicate="mean")
image.plot(zz0, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Observed Precipitation (mm)")
contour(zz0, add=T)
world(add=TRUE, lwd=4)

zz1<-interp(data$Long,data$Lat,bestlm$fitted.values, duplicate="mean")
image.plot(zz1, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="LM Modeled precipitation (mm)")
contour(zz1, add=T)
world(add=TRUE, lwd=4)

zz2<-interp(data$Long,data$Lat,residuals(bestlm, type="response"), duplicat="mean")
image.plot(zz2, col=topo.colors(20,0.5), ylab="latitude", xlab="longitude", main="LM residuals (mm/yr)")
contour(zz2, add=T)
world(add=TRUE, lwd=4)

lmPrecip<-zz1 #save the linear model contour plot
lmResid<-zz2 #save the linear model residuals

#excellent!

# vi. Estimate the precipitation and the standard errors on the DEM grid and show them as 3-D plots as in (v) above.
# get DEM grid Lat-Long-Elev
#predpoints = read.table("http://cires.colorado.edu/~aslater/CVEN_6833/colo_dem.dat",sep=",")
predpoints<-read.table(file="/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013/colo_dem.txt")
names(predpoints)<-c("Lat","Long","Elev")

Xpred<-predpoints #gridded lat, long and elev

# (vi) Estimate the precipitation and the SEs on the DEM grid and show as 3-D plots as in (v) above.
# Predict response variable (Precip) on the grid of predictor variables (long, Lat, Elev)
ypred = predict.lm(bestlm,newdata=predpoints, type="response",se.fit=TRUE)

# Create spatial plot - Latitude, Longitude and the predicted value..
# since the predpoints are on a unifrom spatial grid but not on 
# rectangular grid we can do the following..

library(akima)
zz = interp(predpoints$Long,predpoints$Lat,ypred$fit)

library(fields)
surface(zz,xlab="Longitude",ylab ="Latitude", main="DEM Grid Predictions") 
# how to adjust color (heat) scale to match other plots???

# Discussion
# LM does not capture non-linearities and complexity inherent in the data
# LM attempts to "spread" precipitation evenly across the state, decreasing from West to East.  As can be seen in the LM residuals plot, the errors explain most of the variance in the data, not the model, indicating that the model is inadequate.  

# See plots
################# End Q1
################# 
################# Start Q2
## repeat Q1 with GLM....
setwd("/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013")
#Import data
#Mean annual precipitation at 491 locations in Colorado based on data for the period 1980-2002
data<-read.table(file="/Users/elliotcohen/Dropbox/Data/Climate/Rainfall/Colorado_Annual_Precip_Gridded.csv", sep=",", header=TRUE)

## create a data frame to be used in GLM
data<-as.data.frame(data)
X<-data[,-4]  # df of Independent variable
Y<-data$Precip #vector of dependent variable (annual rainfall at lat-long position)

# i. Fit a ‘best’ generalized linear model (use one of the objective functions - GCV, AIC, SBC or PRESS; you can also try a couple of them to see any differences). This entails fitting the model with all possible combinations of covariates and selecting the model with the minimum objective function.
#Fit the GLM - using Gamma PDF and appropriate canonical link function.
glm<-glm(Precip ~ Lat + Long + Elev, data=data, family=Gamma(link="inverse"))
summary(glm)

# Load the MASS package to use the stepAIC function to choose the "best model"
require(MASS)
bestglm<-stepAIC(glm)

#ii. Perform ANOVA (i.e. model significance) and model diagnostics (i.e., check the assumptions of the residuals – Normality, independence,homoskedasticity). 
summary(bestglm) #"bestglm" has longitutde as only predictor variable
#Note: other models have similar AIC values, indicating that "best" model may not be that much better than the next best model.

# Model diagnostics
par(mfrow=c(2,2))
plot(bestglm)

# Results: Model diagnostics check if key assumptions are met, namely:
# (1) Residuals are "identically and independently distributed" (iid), e.g. their distribution does not change substantially for different values of x.
# To check iid, we look at Plot 1 (upper left) and plot 3 (lower left).
# --> Plot 1 shows a trend of decreasing variance in the residuals for larger values of the response variable, thus violating IID.
# (2) Residuals are normally distributed, which is what Q-Q plot helps us checks (upper right).
# --> QQ plots looks good
# Plot 4 (lower right) helps identify outliers (their row numbers are shown such that we can go back and find them in the dataframe).  Outliers are not an explicit part of the assumptions required for regression, but are important to be aware of.

# can also use custom-built  model diagnostics...
##### Manual Model Diagnostics (6 visual checks) #####
par(mar=c(4,4,4,2) + 0.1)
par(mfrow=c(2,3)) 

#But here is how to do it manually...
# Get the residuals of the model..
modresid=residuals(bestglm)

Y=data$Precip
nX<-length(bestglm$coefficients)
X=subset(data, select=names(bestglm$coefficients)[2:nX])
Yhat=Y-modresid  #residuals = Y - Yestimate ==> Yestimate = Y - residuals
k=dim(X)[2]      #number of regressor variables
p=k+1              #number of model parameters 
n=length(Y)

# Compute ANOVA quantities for use down below
SST = sum((Y - mean(Y))^2)   #SST = Total corrected sum of squares, n-1 dof
#Yhat = Y - bestGlm$res 
#(Y - Yhat = residuals), Yhat is the modeled response of Y
SSR = sum((Yhat - mean(Y))^2)  #SSR = Regression Sum of Squares = sum[(yhati-ybar)^2], dof = k (number of predictor  variables)
SSE = sum((modresid)^2)   #SSE = Error Sum of Squares = sum[(yi-yhati)^2], dof = n-p 
MSR = SSR / ncol(X)            
MSE = SSE/(n - length(bestglm$coef))         

# Now start computing diagnostics and plotting them...
par(mfrow=c(2,3))

# (1) Check if residuals fit a normal distribution
qqnorm(modresid)
qqline(modresid)	

# (2-3) Plot the residuals vs X.  Check to make sure there is *no* apparent pattern.  Distribution of residuals should be random.
for(i in 1:k){
  plot(X[,i],modresid,xlab="X",ylab="residuals",main="Residuals vs. X[i]") 
}

# (4) Plot the residuals vs the model estimates of Y. 
#Check to make sure there is *no* apparent pattern or structure.  In other words, the distribution of the residuals should look random.
plot(Yhat,modresid,xlab="estiimate of Y", ylab="residuals",main="Residuals vs model estimates of Y")

# (5) Plot the autocorrelation function - to make sure the residuals are *not* related to each other.  
z1=acf(modresid,main="autocorrelation of residuals")

# (6) Cooks Distance - to make sure outliers do not exert undue influence on the regression model.

# Compute the Hat matrix
hatm= hatvalues(bestglm)

#studentized residuals - ri  - equation 12-42
ri = modresid/sqrt((1 - diag(hatm)) * MSE)

#Compute Cook's Distance Equation 12-44
Di = ri*ri * diag(hatm) / ((1-diag(hatm)) * length(bestglm$coef))
plot(Y, diag(Di), main="Cook's Distance")
#If Dis are greater than 1 then there are points that have undue influence on the fit.

# Comments on Model Diagnostics: 
# The Q-Q plot looks good
# Plot of residuals vs model estimates of Y reveals EXTREME heteroscdasticity.
# Autocorrelation appears to be present up to lag 3.
# --> glm with gamma distribution is *not* adequate for precip data.
#### END MODEL DIAGNOSTICS ####

# iii. Compute cross-validated and fitted estimates at each observation points and plot them against the observed values. This is to visually see how the model performs in a fitting and cross-validated mode.

## Observed vs modeled estimates of response variable Y
par(mfrow=c(1,1))
#preds<-data.frame(c(paste("data$",names(bestglm$coefficients)[2:nX],sep="")), quote=FALSE)  #trying to automate regressor variables to include in preds
preds<-data.frame(data$Long)
names(preds)<-c("Long")
Fitted<-predict(bestglm, newdata=preds, type="response") #predict precip using best linear model with Long as only predictor variable.

plot(data$Precip,Fitted,pch=20,col="black",xlab="Observed Precip (mm)",ylab="Modeled Precip (mm)", main="Observed vs Modeled Precipitation")
abline(a=0, b=1)

# (iii) Now compute cross validated estimates
par(mfrow=c(1,1))
n=length(Y)
yest=1:n
nvar=dim(X)[2]

index=1:n
for(i in 1:n){
  index1=index[index != i] #drop one observation at a time
  Xval=X[index1,] #X data less the dropped observation
  Yval=Y[index1]  #Y data less the dropped observation
  
  zz=lsfit(Xval,Yval) #re-fit the model without the dropped observation
  
  xpred=c(1,X[i,1:nvar]) #now estimate at the point that was dropped
  xpred<-as.numeric(xpred)
  yest[i]=sum(zz$coef * xpred)
}
#now surface the x-validated estimates..
points(Y, yest, col="blue", pch=20)

# Comments:  Cross-validated model estimates and estimates modeled on the full data set appear internally consistent to each other, but niether are accurate to observed precipitation values for larger values of precip (e.g. above 800 mm/yr).

# iv. Drop 10% of observations, fit the model (i.e., the ‘best’ model from i. above) to the rest of the data and predict the dropped points. Compute RMSE and R2 and show them as boxplots. 
# Drop some % of points, fit the model and predict the dropped points..
library(arules)
source("/Users/elliotcohen/Dropbox/Advance Data Analysis/R source files/myboxplot-stats.r")
source("/Users/elliotcohen/Dropbox/Advance Data Analysis/R source files/myboxplot.r")
nsim = 500
rmseskill=1:nsim
corskill=1:nsim  
N = length(Y)
N10 = round(0.10*N)  	#choose % of points to drop (e.g. 10%)
index=1:N

#get original observed data
X<-data[,-4]  # df of Independent variable
Y<-data$Precip #vector of dependent variable (annual rainfall)

for(i in 1:nsim){
  drop=sample(c(1:N),N10)  #sample 10% of the row indices from 1:N at random
  keep=setdiff(index,drop)  #discard values at the intersection of index and drop (e.g. drop 10% of the data, keep the rest)
  x<-as.data.frame(X[keep,])
  y<-Y[keep]
  xx<-as.data.frame(x)
  zz<-glm(y~xx$Long,data=xx) #fit model to remaining data
  xpred<-X[drop,]
  xx<-as.data.frame(xpred) #re-assign xx to the dropped data
  yhat<-predict.glm(zz,newdata=xx) #predict at dropped points using model fit w.out those points
  rmseskill[i]<-sqrt(mean(Y[drop]-yhat)^2)
  corskill[i]<-cor(Y[drop],yhat)
}

par(mfrow=c(1,2))
zz=myboxplot(rmseskill, main="Simulated RMSE skill",plot=FALSE)
zz$names=rep("",length(zz$names))
z1=bxp(zz,xlab="",ylab="RMSE",cex=1.25)
rmse=mean(((modresid)/sd(Y))^2)
points(z1,rmse,col="red",cex=2,pch=19)  #add a point showing the true RMSE of the data.
title(main="RMSE skill")

zz=myboxplot(corskill, main="Simulated Correlation skill",plot=FALSE)
zz$names=rep("",length(zz$names))
z1=bxp(zz,xlab="",ylab="Cor",cex=1.25)
title(main="Cor skill")

#v. Make a 3-D plot (or spatial colored/contour map) of model estimates (i.e. latitude, longitude and the model estimates) and model error. 
# (v) Spatial map of actual precipitation, model estimates, and model error
plot.new()
par(mfrow=c(2,2))
library(akima)
library(fields)
zz0<-interp(x=data$Long, y=data$Lat, z=data$Precip, duplicate="mean")
image.plot(zz0, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Observed Precipitation (mm)")
contour(zz0, add=T)
world(add=TRUE, lwd=4)

zz1<-interp(data$Long,data$Lat,bestglm$fitted.values, duplicate="mean")
image.plot(zz1, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="GLM Modeled precipitation (mm)")
contour(zz1, add=T)
world(add=TRUE, lwd=4)
glmPrecip<-zz1

zz2<-interp(data$Long,data$Lat,residuals(bestglm, type="response"), duplicat="mean")
image.plot(zz2, col=topo.colors(20,0.5), ylab="latitude", xlab="longitude", main="GLM residuals (mm/yr)")
contour(zz2, add=T)
world(add=TRUE, lwd=4)
glmResid<-zz2

#excellent! As in the R-commands work...
# GLM is clearly *NOT* adequate for modeling precipitation.

# vi. Estimate the precipitation and the standard errors on the DEM grid and show them as 3-D plots as in (v) above.
# get DEM grid Lat-Long-Elev
#predpoints = read.table("http://cires.colorado.edu/~aslater/CVEN_6833/colo_dem.dat",sep=",")
predpoints<-read.table(file="/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013/colo_dem.txt")
names(predpoints)<-c("Lat","Long","Elev")

Xpred<-predpoints #gridded lat, long and elev

# par(mfrow=c(1,2))
# actual<-Tps(X,Y)
# model.est<-Tps(X,bestglm$fitted.values)
# surface(model.est,zlim=range(Y,yhat),xlab='Longitude',ylab='Latitude',main='Predicted Precipitation (mm)')
# surface(actual,zlim=range(Y,yhat),xlab='Longitude',ylab='Latitude',main='Actual Precipitation (mm)')

# (vi) Estimate the precipitation and the SEs on the DEM grid and show as 3-D plots as in (v) above.
# Predict response variable (Precip) on the grid of predictor variables (long, Lat, Elev)
ypred = predict.glm(bestglm,newdata=predpoints, type="response",se.fit=TRUE)

# Create spatial plot - Latitude, Longitude and the predicted value..
# since the predpoints are on a unifrom spatial grid but not on 
# rectangular grid we can do the following..

library(akima)
zz = interp(predpoints$Long,predpoints$Lat,ypred$fit)

library(fields)
surface(zz,xlab="Longitude",ylab ="Latitude", main="DEM Grid Predictions") 
# how to adjust color (heat) scale to match other plots???

# Discussion
# GLM does not capture non-linearities and complexity inherent in the data
# GLM attempts to "spread" precipitation evenly across the state, decreasing from West to East.  As can be seen in the LM residuals plot, the errors explain most of the variance in the data, not the model, indicating that the model is inadequate.

################# End Q2
################# 
################# Start Q3

#Develop an R code to perform local polynomial method. Test this by selecting an alpha (say alpha = 0.5) and order of polynomial (say p = 1) on the following synthetic data
# Evaluate:
# (i) L matrix; 
# (ii) GCV and 
# (iii) Estimate the value of Y at each X and also the 95% confidence interval (plot the fit and the confidence interval). Check your results with the output from LOCFIT package – they both should match.

myLocfit<- function(y,x, alpha, porder){
  
  #Define model parameters
  x<-as.matrix(x) # more general structure allows for multiple regressor variables 
  N<-dim(x)[1] # number of observations
  np<-dim(x)[2] # number of predictor variables
  nparam<-porder*np+1 #number of parameters to be estimated (e.g. beta coefficients)
  K=round(alpha*N)
  
  #Define structures to be populated in the loop... 
  evect<-matrix(0,ncol=N, nrow=1) # create an Nx1 vector
  L<-matrix(0,ncol=N, nrow=N) # create an NxN matrix
  yhat<-matrix(0,ncol=1, nrow=N) # create an Nx1 vector
  beta<-matrix(0, ncol=nparam, nrow=N) # Nx2 matrix
  
  for (i in 1:N){ 
    xp=x[i,]  #xp is the point at which you want the estimation.
    
    #calculate the distance between xp and all other points
    #if xp is not part of the data then simply append it to the x matrix:
    
    #   N1<-N+1  
    #   xnew=rbind(xp,x)
    #   xx=as.matrix(dist(xnew)) #matrix of distances between each point (N1xN1)
    #   xdist=xx[1,2:N1] #grab the 1st row, columns 2:N1; distance from xp to other pts
    #   length(xdist)
    
    # OR....
    #if xp is already in the data, then xnew = x
    xnew=x
    xx=as.matrix(dist(xnew)) #distances between each point (NxN)
    xdist<-xx[i,] 
    length(xdist)
    distorder<-order(xdist) #ordered index of pts in x from closest to farthest.
    
    #add the first column to be 1 to get the augmented matrix..
    ax = t(t(x)-xp)
    xk=cbind(rep(1,N),ax) # deg = 1
    #xk=cbind(rep(1,N),axm ax^2)  # deg = 2
    yk=y
    
    #weight the distances using a bisquare kernel
    distk=xdist[distorder[K]]   #distance to Kth nearest neighbor
    xdist[xdist >= distk]<-distk   #replace all the distances beyond the Kth neighbor to be the same as the kth so the weights become 0.
    weights=15*((1-((xdist/distk)*(xdist/distk)))^2)/16
    weights=weights/sum(weights)  #normalize the weights
    
    #create a diagonal matrix with these weights..
    weights=diag(weights)
    
    ## solve for beta
    ## betai = (xk^T W xk)^-1  (xk^T W yk)
    betai<-solve(t(xk) %*% weights %*% xk) %*% (t(xk) %*% weights %*% yk) 
    # compare betas to coefficients from lsfit...
    zz = lsfit(xk[,2],yk,wt=diag(weights))
    zz$coef # the coefficients in betai should match. --> TRUE!
    
    #get the vector of the L matrix for the ith observation
    evect<-matrix(0,ncol=N, nrow=1) # create an Nx1 vector
    evect[i]<-1 
    
    # Lx = e^T xk (xk^T W xk)^-1  (xk^T W)
    # Lx is a row vector of length N
    Lx<-evect %*% xk %*% solve(t(xk) %*% weights %*% xk) %*% (t(xk) %*% weights)
    L[i,]<-Lx # assign to L matrix
    
    # There are a couple of ways to estimate Yhat at xp
    #yhat[i]<-sum((c(1,xp)) *betai1)   # for a linear polynomial
    #yhat[i]<-sum((c(1,xp,xp^2)) *betai1) # for second order poly
    #yhat[i]=sum(Lx*yk)
    yhat[i]<-Lx %*% yk #compute yhat value at point xp and assign to yhat matrix
  }
  
  ## Now compare with locfit package...
  zz<-locfit(y~x,alpha=alpha,deg=porder, kern="bisq")
  results<-as.data.frame(cbind(y, yhat, fitted(zz)))
  names(results)<-c("Yobs","my_Yest","locfit_Yest")
  print(results)
  
  # and compare the L matrices from my fn and locfit package....
  L1<-matrix(0,nrow=N, ncol=N) # create an NxN array
  for(i in 1:N){L1[i,]<-locfit(y~x,alpha=alpha,deg=porder,ev=x[i], kern="bisq", geth=1)}
  print(L)
  print(L1)
  identical(round(L, digits=3), round(L1, digits=3))
}
##############################################################
# Now impliment my local polynomial function...
x<-1:10
y<-c(-1.4,2.45,1.45,5.38,5.6,5.99,8.1,7.54,8.24,10.8)
myLocfit(x,y, alpha=0.5, porder=1)

###############################################################
# can also check the GCV if selecting the bestalpha and bestdeg...
# In the example above, we assumed alpha=0.5 and deg=1
# use trace(L), trace(L^T L) to get GCV compare with gcvplot
# compute the GCV for this alpha..
gcvalpha=(N*sum((y-yhat)^2)) / ((N-sum(diag(L)))^2)

# compute gcv from the gcvplot command
zz=gcvplot(y ~ x, alpha= alpha, deg=porder, kern="bisq", ev=dat(), scale=TRUE)
list(alpha=alpha, gcvmanual=gcvalpha, gcvplot=zz$values)


# (iv)
# If we set alpha=1 (fit using all the data) and weights=1 (linear global fit) then we should get the same results as linear regression...
#define model parameters alpha and polynomial order....
alpha=1
porder=1

### Modify myLocfit function to produce a least squares fit by setting alpha=1, porder=1 and weights=1.

myLsfit<- function(y,x, alpha, porder){
  
  #number of parameters to be estimated - i.e. number of beta coefficients
  x<-as.matrix(x) #more general in case x has multiple columns...
  N<-dim(x)[1]
  np<-dim(x)[2] # number of predictor variables
  nparam<-porder*np+1 # number of model parameters
  K=round(alpha*N)
  
  #Define structures to be populated in the loop... 
  evect<-matrix(0,ncol=N, nrow=1) # create an Nx1 vector
  L<-matrix(0,ncol=N, nrow=N) # create an NxN matrix
  yhat<-matrix(0,ncol=1, nrow=N) # create an Nx1 vector
  beta<-matrix(0, ncol=nparam, nrow=N) # Nx2 matrix
  
  
  for (i in 1:N){ 
    xp=x[i,]  #xp is the point at which you want the estimation.
    #xp=x[i]  #if x is a single vector as apposed to a matrix
    
    #calculate the distance between xp and all other points
    #if xp is a point that is not part of the data then you will simply append this to the x matrix as such:
    
    #   N1<-N+1  
    #   xnew=rbind(xp,x)
    #   xx=as.matrix(dist(xnew)) #matrix of distances between each point (N1xN1)
    #   xdist=xx[1,2:N1] #grab the 1st row, columns 2:N1; distance from xp to other pts
    #   length(xdist)
    
    # OR....
    #if xp is already in the data, then xnew = x
    xnew=x
    xx=as.matrix(dist(xnew)) #distances between each point (NxN)
    xdist<-xx[i,] 
    length(xdist)
    distorder<-order(xdist) #ordered index of pts in x from closest to farthest.
    
    #add the first column to be 1 to get the augmented matrix..
    ax = t(t(x)-xp)
    xk=cbind(rep(1,N),ax) # deg = 1
    #xk=cbind(rep(1,N),axm ax^2)  # deg = 2
    yk=y
    
    #weight the distances using a bisquare kernel
    distk=xdist[distorder[K]]   #distance to Kth nearest neighbor
    xdist[xdist >= distk]<-distk   #replace all the distances beyond the Kth neighbor to be the same as the kth so the weights go to 0.
    #weights=15*((1-((xdist/distk)*(xdist/distk)))^2)/16
    weights=rep(1,N)  #coerce weights to 1 instead of using bisquare kernel...
    
    weights=weights/sum(weights)  #normalize the weights
    weights=diag(weights) #create a diagonal matrix with these weights..
    
    ## solve for beta
    ## betai = (xk^T W xk)^-1  (xk^T W yk)
    betai<-solve(t(xk) %*% weights %*% xk) %*% (t(xk) %*% weights %*% yk) 
    # 2x2 %*% 2x1 = 2x1
    # We don't actually need beta to compute the L matrix....
    zz = lsfit(xk[,2],yk,wt=diag(weights))
    zz$coef # the coefficients in betai should match. --> TRUE!
    
    #get the vector of the L matrix for the ith observation
    evect<-matrix(0,ncol=N, nrow=1) # create an Nx1 vector
    evect[i]<-1
    
    # Lx = e^T xk (xk^T W xk)^-1  (xk^T W)
    # Lx is a row vector of length N
    Lx<-evect %*% xk %*% solve(t(xk) %*% weights %*% xk) %*% (t(xk) %*% weights)
    L[i,]<-Lx # assign to L matrix
    
    # There are a couple of ways to estimate Yhat at xp
    #yhat[i]<-sum((c(1,xp)) *betai1)   # for a linear polynomial
    #yhat[i]<-sum((c(1,xp,xp^2)) *betai1) # for second order poly
    #yhat[i] = sum(Lx*yk)
    yhat[i]<-Lx %*% yk #compute yhat value at point xp and assign to yhat matrix
  }
  
  ## Now compare with lsfit package...
  zz<-lm(y~x)
  zz$fitted.values # yhat
  zz$residuals # residuals
  identical(y-zz$residuals, zz$fitted.values) # y[i] - yhat[i] = residual[i]
  #list(yobs=y, my_yhat=as.numeric(yhat), lsfit_yhat=as.numeric(zz$fitted.values))
  results<-as.data.frame(cbind(y, yhat, zz$fitted.values))
  names(results)<-c("Yobs","my_Yhat","lsfit_Yhat")
  print(identical(round(as.numeric(yhat), digits=3), round(as.numeric(zz$fitted.values), digits=3))) 
  print(results)
  
  # and compare the L matrices from my fn and locfit package....
  H1<-matrix(0,nrow=N, ncol=N) # create an Nx1 matrix
  #H1<-hatvalues(zz)
  X<-cbind(rep(1,N),x)
  H1<-X %*% solve(t(X) %*% X) %*% t(X)
  print(L)
  print(H1)
  identical(round(L, digits=3), round(H1, digits=3))
}

myLsfit(x,y,alpha=1,porder=1)
# Yhat estimated from myLsfit (modified myLocfit with alpha=1 and weights=diag(1)) is identical to Yhat estimated from linear regression. 
# L matrix estimated from myLsfit (modified myLocfit) is identical to Yhat estimated from linear regression.

###############

# (v) Replace the last value of y with an outlier.  Set alpha=0.9 and recompute the L matrix from local polynomial regression and the Hat matrix from linear regression.  Compare the weights and comment
# Now impliment my local polynomial function...
x<-1:10
ymod<-c(-1.4,2.45,1.45,5.38,5.6,5.99,8.1,7.54,8.24,10.8*10)
myLocfit(y=ymod, x=x, alpha=0.9, porder=1) # Prints L matrix from myLocfit and locfit package as well as yhat estimates from myLocfit and locfit package...

# Now compare with hat matrix from linear regression...
X<-cbind(rep(1,N),x) # the augmented X matrix
H1<-X %*% solve(t(X) %*% X) %*% t(X) # the Hat matrix
H1


# Comments
# The L matrix from myLocfit gives weight to only 8 of 10 observations (alpha=0.9) compared to all the observations in linear regression. # This alone helps reduce the influence of outliers
# In addition, the non-zero weights farthest from the fitted point xp are smaller in the Lmatrix than in the Hat matrix.
# Finally, if we loook at the diagonal of the L matrix, we see that y(xp) (e.g. L[i,i]) is given the most weight in estimating the function at xp, which is what you want from a local polynomial fit.  In contrast, the value of y(xp) is often not given the most weight in estimating the function at xp in linear regression, thus making the fit more susecptible to outliers.
################# End Q3
################# 
################# Start Q4
# 4. Repeat 1 with Local polynomial method.
setwd("/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013")

#Import data: Mean annual precipitation at 491 locations in Colorado for the period 1980-2002
data<-read.table(file="/Users/elliotcohen/Dropbox/Data/Climate/Rainfall/Colorado_Annual_Precip_Gridded.csv", sep=",", header=TRUE)


# i. Fit a ‘best’ local regression model
################ subset selection ##############
library(leaps)  	# to provide combinations
library(MPV)		# to help estimate PRESS and consequently, GCV
library(locfit)
#Select independent and dependent variables
X<-data[,1:3] # all the predictor set.
Y<-data[,4] # response variable Y
N = length(Y)

combs = leaps(X,Y, nbest=25)  #Get upto 25 combinations for each
combos = combs$which  # number of predictors
ncombos = length(combos[,1])

for(i in 1:ncombos){
  xx<-as.data.frame(X[,combos[i,]])
  names(xx)<-names(X)[combos[i,]]
  
  ## search for best alpha over a range of alpha values between 0,1
  nvar=dim(xx)[2]
  porder=1
  minalpha=3*(nvar*porder+1)/N  #try 3x instead of 2x....
  alpha1=seq(minalpha,1.0,by=0.05)
  n=length(alpha1)
  
  porder=2
  minalpha=3*(nvar*porder+1)/N 
  alpha2=seq(minalpha,1.0,by=0.05)
  alpha=c(alpha1,alpha2)
  
  #gcvplot accepts matrix/vector input....
  y<-as.matrix(Y)
  xx<-as.matrix(xx)
  #dimnames(xx)<-list(rep("obs",N),names(X)[combos[i,]])
  zz<-gcvplot(y ~ xx, alpha=alpha1, deg=1, kern='bisq', ev=dat(), scale=TRUE)
  z1=gcvplot(y ~ xx, alpha=alpha2, deg=2, kern="bisq", ev=dat(), scale=TRUE)
  
  # pick the best alpha and the degree of the polynomial that gives the least GCV
  z2=order(c(zz$values,z1$values)) #order from lowest GCV to highest
  
  deg1=1
  if(z2[1] > n)deg1=2 #choose the degree (1 or 2) that yields the lowest GCV
  bestalpha=alpha[z2[1]] #choose the alpha that yields the lowest GCV
  bestdeg=deg1
  
  ## now apply the best alpha and bestdeg to fit local polynomial
  zz<-locfit(y~xx, alpha=bestalpha, deg=bestdeg)
  if(i==1) {GCVs <- gcv(zz)[[4]]} else {GCVs <- rbind(GCVs,gcv(zz)[[4]])} #create vector of GCV values
  
}

# select the model with the overall lowest GCV and re-fit the "bestmodel"
besti<-which.min(GCVs)
names<-names(X)[combos[besti,]]
X<-as.data.frame(X[,combos[besti,]]) #capital X = df; lowercase x = matrix
names(X)<-names

## search for best alpha over a range of alpha values between 0,1
nvar=dim(X)[2]
porder=1
minalpha=3*(nvar*porder+1)/N  #try 3x instead of 2x....
alpha1=seq(minalpha,1.0,by=0.05)
n=length(alpha1)

porder=2
minalpha=3*(nvar*porder+1)/N 
alpha2=seq(minalpha,1.0,by=0.05)
alpha=c(alpha1,alpha2)

#gcvplot accepts matrix/vector input....
y<-as.matrix(Y)
x<-as.matrix(X)
#dimnames(xx)<-list(rep("obs",N),names(X)[combos[besti,]])
zz<-gcvplot(y ~ x, alpha=alpha1, deg=1, kern='bisq', ev=dat(), scale=TRUE)
z1=gcvplot(y ~ x, alpha=alpha2, deg=2, kern="bisq", ev=dat(), scale=TRUE)

# pick the best alpha and the degree of the polynomial that gives the least GCV
z2=order(c(zz$values,z1$values)) #order from lowest GCV to highest

deg1=1
if(z2[1] > n)deg1=2 #choose the degree (1 or 2) that yields the lowest GCV
bestalpha=alpha[z2[1]] #choose the alpha that yields the lowest GCV
bestdeg=deg1

## now apply the best alpha and bestdeg to fit local polynomial
#bestmodel<-locfit(y~x, alpha=bestalpha, deg=bestdeg) #matrix input
newdf<-as.data.frame(cbind(Y,X))
bestmodel<-locfit(Y~Long, data=newdf, alpha=bestalpha, deg=bestdeg) #df input
summary(bestmodel)
zz<-bestmodel

##### Manual Model Diagnostics (6 visual checks) #####
par(mfrow=c(2,2)) 

# Get the residuals of the model..
modresid=residuals(bestmodel)

#for lm...
#Y=data$Precip
#X<-subset(data, select=names(bestmodel$coefficients)[2:nX])
#nX<-length(bestmodel$coefficients)
#nX<-dim(xx)[2]

# for locfit...
#Y=data$Precip
#X<-data[,1:3] # all the predictor set.
#X<-subset(data, select=names(X)[combos[besti,]]) #best predictor set
nX<-dim(X)[2]

Yhat=Y-modresid  #residuals = Y - Yestimate ==> Yestimate = Y - residuals
k=dim(X)[2]      #number of regressor variables
p=k+1             #number of model parameters 
n=length(Y)

# Compute ANOVA quantities for use down below
SST = sum((Y - mean(Y))^2)   #SST = Total corrected sum of squares, n-1 dof
#(Y - Yhat = residuals), Yhat is the modeled response of Y
SSR = sum((Yhat - mean(Y))^2)  #SSR = Regression Sum of Squares = sum[(yhati-ybar)^2], dof = k (number of predictor  variables)
SSE = sum((modresid)^2)   #SSE = Error Sum of Squares = sum[(yi-yhati)^2], dof = n-p 
MSR = SSR / ncol(X)            
MSE = SSE/nu11 # MSE = (residual sum of squares)/dof. The denominator is the sample size reduced by the number of model parameters estimated from the same data, (n-p) for p regressors or (n-p-1) if an intercept is used (e.g. degrees of freedom)         

# Now start computing diagnostics and plotting them...
par(mfrow=c(2,2))

# (1) Check if residuals fit a normal distribution
qqnorm(modresid)
qqline(modresid)	

# (2-3) Plot the residuals vs X.  Check to make sure there is *no* apparent pattern. Distribution of residuals should be random.
for(i in 1:k){
  plot(X[,i],modresid,xlab="X",ylab="residuals",main="Residuals vs. X[i]") 
}

# (4) Plot the residuals vs the model estimates of Y. 
#Check to make sure there is *no* apparent pattern or structure.  In other words, the distribution of the residuals should look random.
plot(Yhat,modresid,xlab="estiimate of Y", ylab="residuals",main="Residuals vs model estimates of Y")

# (5) Plot the autocorrelation function - to make sure the residuals are *not* related to each other.  
z1=acf(modresid,main="autocorrelation of residuals")

# Comments on Model Diagnostics: 
# Q-Q plot reveals non-normality in the upper tail --> FAIL diagnostic.
# Plot of residuals vs predictor variable X1 reveals heteroscdasticity --> FAIL diagnostic.
# Plot of residuals vs model estimates of Y reveals heteroscdasticity --> FAIL diagnostic.
# Autocorrelation appears to be present up to lag 3 --> FAIL diagnostic
# --> locfit with gaussian distribution is *not* adequate for modeling precip data.
#### END MODEL DIAGNOSTICS ####

#### Ftest comparing locfit to lm
RSS1 = sum(residuals(zz)^2)
nu1 = sum(fitted(zz,what="infl"))     # trace(L)
nu2 = sum(fitted(zz,what="vari"))     # trace(L^T L)
## Or
#nu1 = zz$dp[6]
#nu2 = zz$dp[7]
nu11 = N-2*nu1 + nu2

## Linear regression..
zz=lm(Y ~ xx)

N = length(Y)
XX = cbind(rep(1,N), xx)
# Compute the Hat matrix
hatm = XX %*% solve(t(XX) %*% XX) %*% t(XX)

II = diag(N)
delta0 = t(II-hatm)%*%(II-hatm)    #Equation 9.2
nu00 = sum(diag(delta0))
RSS0 = sum(residuals(zz)^2)

Fdata = (RSS0 - RSS1)/(nu00 - nu11)
Fdata = (Fdata / (RSS1 / nu11))
Ftheory = qf(0.95,(nu00-nu11), nu11)    #95% confidence level..
Fdata>Ftheory # TRUE
## Fdata > Ftheor   - reject null - i.e., reject that the data fits a linear model

# iii. Compute cross-validated and fitted estimates at each observation points and plot them against the observed values. This is to visually see how the model performs in a fitting and cross-validated mode.

## Observed vs modeled estimates of response variable Y
par(mfrow=c(1,1))
#preds<-data.frame(c(paste("data$",names(bestmodel$coefficients)[2:nX],sep="")), quote=FALSE)  #trying to automate regressor variables to include in preds
preds<-data.frame(data$Long)
names(preds)<-c("Long")
Fitted<-predict(bestmodel, newdata=preds, type="response") #predict precip using best local model with Long as only predictor variable.

plot(data$Precip,Fitted,pch=20,col="black",xlab="Observed Precip (mm)",ylab="Modeled Precip (mm)", main="Observed vs Modeled Precipitation")
abline(a=0, b=1)

# (iii) Now compute cross validated estimates
par(mfrow=c(1,1))
n=length(Y)
yest=1:n
nvar=dim(X)[2]

index=1:n
for(i in 1:n){
  index1=index[index != i] #drop one observation at a time
  Xval=X[index1,] #X data less the dropped observation
  Yval=Y[index1]  #Y data less the dropped observation
  
  zz<-locfit(Yval ~ Xval, alpha=bestalpha, deg=bestdeg) #fit the model less the dropped observation
  
  #xpred=c(1,X[i,1:nvar]) #now estimate at the point that was dropped
  xpred<-X[i,1:nvar] #now estimate at the point that was dropped
  xpred<-as.numeric(xpred)
  yest[i]<-predict(zz, newdata=xpred, type="response")
}
#now surface the x-validated estimates..
points(Y, yest, col="blue", pch=20)

modresid<-Y-yest  #model residuals given by actual minus estimated value of response variable

# Comments:  Cross-validated model estimates and estimates modeled on the full data set appear internally consistent to each other, but niether are accurate to observed precipitation values for larger values of precip (e.g. above 1000 mm/yr).


# iv. Drop 10% of observations, fit the model (i.e., the ‘best’ model from i. above) to the rest of the data and predict the dropped points. Compute RMSE and R2 and show them as boxplots. 
# Drop some % of points, fit the model and predict the dropped points..
library(arules)
source("/Users/elliotcohen/Dropbox/Advance Data Analysis/R source files/myboxplot-stats.r")
source("/Users/elliotcohen/Dropbox/Advance Data Analysis/R source files/myboxplot.r")

nsim = 500
rmseskill=1:nsim
corskill=1:nsim  
N = length(Y)
N10 = round(0.10*N)    #choose % of points to drop (e.g. 10%)
index=1:N

for(i in 1:nsim){
  drop=sample(c(1:N),N10)  #sample 10% of the row indices from 1:N at random
  keep=setdiff(index,drop)  #discard values at the intersection of index and drop (e.g. drop 10% of the data, keep the rest)
  xx<-as.data.frame(X[keep,])
  y<-as.data.frame(Y[keep])
  newdf<-as.data.frame(cbind(y,xx))
  names(newdf)<-c(names(data)[4],names(data)[combos[besti,]])
  # For locfit, *either* specify column names of a data frame and supply data=df *or* supply matrices and omit data=df
  zz<-locfit(Precip~Long, data=newdf, alpha=bestalpha, deg=bestdeg) #fit model to remaining data
  xpred<-X[drop,] #the dropped data
  yhat<-predict(zz,newdata=xpred) #predict at dropped points using model fit w.out those points
  rmseskill[i]<-sqrt(mean(Y[drop]-yhat)^2)
  corskill[i]<-cor(Y[drop],yhat)
}

par(mfrow=c(1,2))
zz=myboxplot(rmseskill, main="Simulated RMSE skill",plot=FALSE)
zz$names=rep("",length(zz$names))
z1=bxp(zz,xlab="",ylab="RMSE",cex=1.25)
rmse=mean(((modresid)/sd(Y))^2) 
points(z1,rmse,col="red",cex=2,pch=19)  #add a point showing the true RMSE of the data.
title(main="RMSE skill")

zz=myboxplot(corskill, main="Simulated Correlation skill",plot=FALSE)
zz$names=rep("",length(zz$names))
z1=bxp(zz,xlab="",ylab="Cor",cex=1.25)
title(main="Cor skill")

# How should I interpret the RMSE and cor boxplots?

#v. Make a 3-D plot (or spatial colored/contour map) of model estimates (i.e. latitude, longitude and the model estimates) and model error. 
# (v) Spatial map of actual precipitation, model estimates, and model error
plot.new()
par(mfrow=c(2,2))
library(akima)
library(fields)
zz0<-interp(x=data$Long, y=data$Lat, z=data$Precip, duplicate="mean")
image.plot(zz0, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Observed Precipitation (mm)")
contour(zz0, add=T)
world(add=TRUE, lwd=4)

zz1<-interp(data$Long,data$Lat,yest, duplicate="mean")
image.plot(zz1, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Local Polynomial (family=Gaussian)\nFitted Precipitation (mm)")
contour(zz1, add=T)
world(add=TRUE, lwd=4)

zz2<-interp(data$Long,data$Lat,modresid, duplicate="mean")
image.plot(zz2, col=topo.colors(20,0.5), ylab="latitude", xlab="longitude", main="Local Polynomial (family=Gaussian)\nResiduals (mm/yr)")
contour(zz2, add=T)
world(add=TRUE, lwd=4)

# vi. Estimate the precipitation and the standard errors on the DEM grid and show them as above.
# get DEM grid Lat-Long-Elev
# predpoints = read.table("http://cires.colorado.edu/~aslater/CVEN_6833/colo_dem.dat",sep=",")
predpoints<-read.table(file="/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013/colo_dem.txt")
names(predpoints)<-c("Lat","Long","Elev")
Xpred<-predpoints$Long #gridded lat, long and elev
ypred = predict(bestmodel,newdata=Xpred, type="response",se.fit=TRUE)

zz = interp(predpoints$Long,predpoints$Lat,ypred$fit)
surface(zz,xlab="Longitude",ylab ="Latitude", main="DEM Grid Predictions") 
# how to adjust color (heat) scale to match other plots???

# Discussion:
# The local polynomial regression is better than the linear model and the generalized linear model with appropriate link function, but still not very good as can be seen in the residuals and from the model diagnostics.  We see that most of the variance is still explained by the model residuals rather than the local polynomial model, indicating that the model is insufficienct.  
# Let's repeat one more time, but now with an appropriate link function, e.g. local glm.
################# End Q4
################# 
################# Start Q5
# 5. Repeat 4 with Local Polynomial method but using the appropriate link function (i.e. ‘Local GLM’). [For the Local Polynomial approach the ‘best model’ involves fitting the best subset of predictors and the smoothing parameter, alpha. You can also compare the GCV from these four different methods.]

setwd("/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013")

#Import data: Mean annual precipitation at 491 locations in Colorado for the period 1980-2002
data<-read.table(file="/Users/elliotcohen/Dropbox/Data/Climate/Rainfall/Colorado_Annual_Precip_Gridded.csv", sep=",", header=TRUE)

#Select independent and dependent variables
X<-data[,1:3] # all the predictor set.
Y<-data[,4] # response variable Y
N = length(Y)

# choose an appropriate link function....
#install.packages("fitdistrplus")
library("fitdistrplus")
plotdist(Y) # interpreting this, use a gamma distribution
descdist(Y) # interpreting this, use a gamma distribution
fig1 <- fitdist(Y, "gamma")
plot(fig1)
summary(fig1)  

# i. Fit a ‘best’ local regression model
################ subset selection ##############
library(leaps)    # to provide combinations
library(MPV)		# to help estimate PRESS and consequently, GCV
library(locfit)
#Select independent and dependent variables
# X<-data[,1:3] # all the predictor set.
# Y<-data[,4] # response variable Y
# N = length(Y)

combs = leaps(X,Y, nbest=25)  #Get upto 25 combinations for each
combos = combs$which  # number of predictors
ncombos = length(combos[,1])

for(i in 1:ncombos){
  xx<-as.data.frame(X[,combos[i,]])
  names(xx)<-names(X)[combos[i,]]
  
  ## search for best alpha over a range of alpha values between 0,1
  nvar=dim(xx)[2]
  porder=1
  minalpha=3*(nvar*porder+1)/N  #try 3x instead of 2x....
  alpha1=seq(minalpha,1.0,by=0.05)
  n=length(alpha1)
  
  porder=2
  minalpha=3*(nvar*porder+1)/N 
  alpha2=seq(minalpha,1.0,by=0.05)
  alpha=c(alpha1,alpha2)
  
  #gcvplot accepts matrix/vector input....
  y<-as.matrix(Y)
  xx<-as.matrix(xx)
  #dimnames(xx)<-list(rep("obs",N),names(X)[combos[i,]])
  zz<-gcvplot(y ~ xx, alpha=alpha1, deg=1, kern='bisq', ev=dat(), scale=TRUE, family="gamma")
  z1=gcvplot(y ~ xx, alpha=alpha2, deg=2, kern="bisq", ev=dat(), scale=TRUE, family="gamma")
  
  # pick the best alpha and the degree of the polynomial that gives the least GCV
  z2=order(c(zz$values,z1$values)) #order from lowest GCV to highest
  
  deg1=1
  if(z2[1] > n)deg1=2 #choose the degree (1 or 2) that yields the lowest GCV
  bestalpha=alpha[z2[1]] #choose the alpha that yields the lowest GCV
  bestdeg=deg1
  
  ## now apply the best alpha and bestdeg to fit local polynomial
  zz<-locfit(y~xx, alpha=bestalpha, deg=bestdeg, family="gamma")
  

  if(i==1) {GCVs <- gcv(zz)[[4]]} else {GCVs <- rbind(GCVs,gcv(zz)[[4]])} #create vector of GCV values
  
}

# select the model with the overall lowest GCV and re-fit the "bestmodel"
besti<-which.min(GCVs)
names<-names(X)[combos[besti,]]
X<-as.data.frame(X[,combos[besti,]]) #capital X = df; lowercase x = matrix
names(X)<-names

## search for best alpha over a range of alpha values between 0,1
nvar=dim(X)[2]
porder=1
minalpha=3*(nvar*porder+1)/N  #try 3x instead of 2x....
alpha1=seq(minalpha,1.0,by=0.05)
n=length(alpha1)

porder=2
minalpha=3*(nvar*porder+1)/N 
alpha2=seq(minalpha,1.0,by=0.05)
alpha=c(alpha1,alpha2)

#gcvplot accepts matrix/vector input....
y<-as.matrix(Y)
x<-as.matrix(X)
#dimnames(xx)<-list(rep("obs",N),names(X)[combos[besti,]])
zz<-gcvplot(y ~ x, alpha=alpha1, deg=1, family="gamma", kern='bisq', ev=dat(), scale=TRUE)
z1=gcvplot(y ~ x, alpha=alpha2, deg=2, family="gamma", kern="bisq", ev=dat(), scale=TRUE)

# pick the best alpha and the degree of the polynomial that gives the least GCV
z2=order(c(zz$values,z1$values)) #order from lowest GCV to highest

deg1=1
if(z2[1] > n)deg1=2 #choose the degree (1 or 2) that yields the lowest GCV
bestalpha=alpha[z2[1]] #choose the alpha that yields the lowest GCV
bestdeg=deg1

## now apply the best alpha and bestdeg to fit local polynomial
#bestmodel<-locfit(y~x, alpha=bestalpha, deg=bestdeg) #matrix input
newdf<-as.data.frame(cbind(Y,X))
bestmodel<-locfit(Y~Long, data=newdf, alpha=bestalpha, deg=bestdeg, family="gamma") #df input
summary(bestmodel)
zz<-bestmodel

##### Manual Model Diagnostics (6 visual checks) #####
par(mfrow=c(2,2)) 

# Get the residuals of the model..
modresid=residuals(bestmodel)

#for lm...
#Y=data$Precip
#X<-subset(data, select=names(bestmodel$coefficients)[2:nX])
#nX<-length(bestmodel$coefficients)
#nX<-dim(xx)[2]

# for locfit...
#Y=data$Precip
#X<-data[,1:3] # all the predictor set.
#X<-subset(data, select=names(X)[combos[besti,]]) #best predictor set
nX<-dim(X)[2]

Yhat=Y-modresid  #residuals = Y - Yestimate ==> Yestimate = Y - residuals
k=dim(X)[2]      #number of regressor variables
p=k+1             #number of model parameters 
n=length(Y)

# Compute ANOVA quantities for use down below
SST = sum((Y - mean(Y))^2)   #SST = Total corrected sum of squares, n-1 dof
#(Y - Yhat = residuals), Yhat is the modeled response of Y
SSR = sum((Yhat - mean(Y))^2)  #SSR = Regression Sum of Squares = sum[(yhati-ybar)^2], dof = k (number of predictor  variables)
SSE = sum((modresid)^2)   #SSE = Error Sum of Squares = sum[(yi-yhati)^2], dof = n-p 
MSR = SSR / ncol(X)            
MSE = SSE/nu11 # MSE = (residual sum of squares)/dof. The denominator is the sample size reduced by the number of model parameters estimated from the same data, (n-p) for p regressors or (n-p-1) if an intercept is used (e.g. degrees of freedom)         

# Now start computing diagnostics and plotting them...
par(mfrow=c(2,2))

# (1) Check if residuals fit a normal distribution
qqnorm(modresid)
qqline(modresid)	

# (2-3) Plot the residuals vs X.  Check to make sure there is *no* apparent pattern. Distribution of residuals should be random.
for(i in 1:k){
  plot(X[,i],modresid,xlab="X",ylab="residuals",main="Residuals vs. X[i]") 
}

# (4) Plot the residuals vs the model estimates of Y. 
#Check to make sure there is *no* apparent pattern or structure.  In other words, the distribution of the residuals should look random.
plot(Yhat,modresid,xlab="estiimate of Y", ylab="residuals",main="Residuals vs model estimates of Y")

# (5) Plot the autocorrelation function - to make sure the residuals are *not* related to each other.  
z1=acf(modresid,main="autocorrelation of residuals")

# Comments on Model Diagnostics: 
# Q-Q plot looks normal --> PASS diagnostic.
# Plot of residuals vs predictor variable X1 reveals heteroscdasticity --> FAIL diagnostic.
# Plot of residuals vs model estimates of Y reveals EXTREME heteroscdasticity --> FAIL diagnostic.
# Autocorrelation does *not* seem to be a problem--> PASS diagnostic
# Overall --> FAIL diagnostics
# --> locfit with gamma distribution is *not* adequate for modeling precip data.

#### END MODEL DIAGNOSTICS ####

#### Ftest comparing locfit to lm
RSS1 = sum(residuals(zz)^2)
nu1 = sum(fitted(zz,what="infl"))     # trace(L)
nu2 = sum(fitted(zz,what="vari"))     # trace(L^T L)
## Or
#nu1 = zz$dp[6]
#nu2 = zz$dp[7]
nu11 = N-2*nu1 + nu2

## Linear regression..
zz=lm(y ~ x)

N = length(y)
XX = cbind(rep(1,N), x)
# Compute the Hat matrix
hatm = XX %*% solve(t(XX) %*% XX) %*% t(XX)

II = diag(N)
delta0 = t(II-hatm)%*%(II-hatm)    #Equation 9.2
nu00 = sum(diag(delta0))
RSS0 = sum(residuals(zz)^2)

Fdata = (RSS0 - RSS1)/(nu00 - nu11)
Fdata = (Fdata / (RSS1 / nu11))
Ftheory = qf(0.95,(nu00-nu11), nu11)    #95% confidence level..
Fdata>Ftheory # TRUE
## Fdata > Ftheor   - reject null - i.e., reject that the data fits a linear model

# iii. Compute cross-validated and fitted estimates at each observation points and plot them against the observed values. This is to visually see how the model performs in a fitting and cross-validated mode.

## Observed vs modeled estimates of response variable Y
par(mfrow=c(1,1))
#preds<-data.frame(c(paste("data$",names(bestmodel$coefficients)[2:nX],sep="")), quote=FALSE)  #trying to automate regressor variables to include in preds
preds<-data.frame(data$Long)
names(preds)<-c("Long")
Fitted<-predict(bestmodel, newdata=preds, type="response") #predict precip using best local model with Long as only predictor variable.

plot(data$Precip,Fitted,pch=20,col="black",xlab="Observed Precip (mm)",ylab="Modeled Precip (mm)", main="Observed vs Modeled Precipitation")
abline(a=0, b=1)

# (iii) Now compute cross validated estimates
par(mfrow=c(1,1))
n=length(Y)
yest=1:n
nvar=dim(X)[2]

index=1:n
for(i in 1:n){
  index1=index[index != i] #drop one observation at a time
  Xval=X[index1,] #X data less the dropped observation
  Yval=Y[index1]  #Y data less the dropped observation
  
  zz<-locfit(Yval ~ Xval, alpha=bestalpha, deg=bestdeg, family="gamma") #fit the model less the dropped observation
  
  #xpred=c(1,X[i,1:nvar]) #now estimate at the point that was dropped
  xpred<-X[i,1:nvar] #now estimate at the point that was dropped
  xpred<-as.numeric(xpred)
  yest[i]<-predict(zz, newdata=xpred, type="response")
}
#now surface the x-validated estimates..
points(Y, yest, col="blue", pch=20)

#
modresid<-Y-Fitted  #model residuals given by actual minus estimated value of response variable

# Comments:  Cross-validated model estimates and estimates modeled on the full data set appear internally consistent to each other, but niether are accurate to observed precipitation values for larger values of precip (e.g. above 1000 mm/yr).


# iv. Drop 10% of observations, fit the model (i.e., the ‘best’ model from i. above) to the rest of the data and predict the dropped points. Compute RMSE and R2 and show them as boxplots. 
# Drop some % of points, fit the model and predict the dropped points..
library(arules)
source("/Users/elliotcohen/Dropbox/Advance Data Analysis/R source files/myboxplot-stats.r")
source("/Users/elliotcohen/Dropbox/Advance Data Analysis/R source files/myboxplot.r")

nsim = 500
rmseskill=1:nsim
corskill=1:nsim  
N = length(Y)
N10 = round(0.10*N)    #choose % of points to drop (e.g. 10%)
index=1:N

for(i in 1:nsim){
  drop=sample(c(1:N),N10)  #sample 10% of the row indices from 1:N at random
  keep=setdiff(index,drop)  #discard values at the intersection of index and drop (e.g. drop 10% of the data, keep the rest)
  xx<-as.data.frame(X[keep,])
  y<-as.data.frame(Y[keep])
  newdf<-as.data.frame(cbind(y,xx))
  names(newdf)<-c(names(data)[4],names(data)[combos[besti,]])
  # For locfit, *either* specify column names of a data frame and supply data=df *or* supply matrices and omit data=df
  zz<-locfit(Precip~Long, data=newdf, alpha=bestalpha, deg=bestdeg) #fit model to remaining data
  xpred<-X[drop,] #the dropped data
  yhat<-predict(zz,newdata=xpred) #predict at dropped points using model fit w.out those points
  rmseskill[i]<-sqrt(mean(Y[drop]-yhat)^2)
  corskill[i]<-cor(Y[drop],yhat)
}

par(mfrow=c(1,2))
zz=myboxplot(rmseskill, main="Simulated RMSE skill",plot=FALSE)
zz$names=rep("",length(zz$names))
z1=bxp(zz,xlab="",ylab="RMSE",cex=1.25)
rmse=mean(((modresid)/sd(Y))^2) 
points(z1,rmse,col="red",cex=2,pch=19)  #add a point showing the true RMSE of the data.
title(main="RMSE skill")

zz=myboxplot(corskill, main="Simulated Correlation skill",plot=FALSE)
zz$names=rep("",length(zz$names))
z1=bxp(zz,xlab="",ylab="Cor",cex=1.25)
title(main="Cor skill")

# How should I interpret the RMSE and cor boxplots?

#v. Make a 3-D plot (or spatial colored/contour map) of model estimates (i.e. latitude, longitude and the model estimates) and model error. 
# (v) Spatial map of actual precipitation, model estimates, and model error
plot.new()
par(mfrow=c(2,2))
library(akima)
library(fields)
zz0<-interp(x=data$Long, y=data$Lat, z=data$Precip, duplicate="mean")
image.plot(zz0, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Observed Precipitation (mm)")
contour(zz0, add=T)
world(add=TRUE, lwd=4)

zz1<-interp(data$Long,data$Lat,Fitted, duplicate="mean")
image.plot(zz1, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Local Polynomial (family=Gamma)\nFitted Precipitation (mm)")
contour(zz1, add=T)
world(add=TRUE, lwd=4)

zz2<-interp(data$Long,data$Lat,modresid, duplicate="mean")
image.plot(zz2, col=topo.colors(20,0.5), ylab="latitude", xlab="longitude", main="Local Polynomial (family=Gamma)\nResiduals (mm/yr)")
contour(zz2, add=T)
world(add=TRUE, lwd=4)

# vi. Estimate the precipitation and the standard errors on the DEM grid and show them as above.
# get DEM grid Lat-Long-Elev
# predpoints = read.table("http://cires.colorado.edu/~aslater/CVEN_6833/colo_dem.dat",sep=",")
predpoints<-read.table(file="/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013/colo_dem.txt")
names(predpoints)<-c("Lat","Long","Elev")
Xpred<-predpoints$Long #gridded lat, long and elev
ypred = predict(bestmodel,newdata=Xpred, type="response",se.fit=TRUE)

zz = interp(predpoints$Long,predpoints$Lat,ypred$fit)
surface(zz,xlab="Longitude",ylab ="Latitude", main="DEM Grid Predictions") 
# how to adjust color (heat) scale to match other plots???
###########

# 6. Estimate the spatial surface using Kriging
# i. Fit a variogram using co-kriging (i.e., using latitude, longitude and elevation)
# ii. Repeat iii – vi of problem 1.

###### up to here on Sunday October 13, 2013  ###############
###### up to here on Sunday October 13, 2013  ###############
###### up to here on Sunday October 13, 2013  ###############
###### up to here on Sunday October 13, 2013  ###############
###### up to here on Sunday October 13, 2013  ###############
###### up to here on Sunday October 13, 2013  ###############
###### up to here on Sunday October 13, 2013  ###############
# 7. Repeat 6. with a Hierarchical Spatial Model


################# 
################# 
################# Start Q8
data<-read.table(file="http://cires.colorado.edu/~aslater/CVEN_6833/colo_pcp_daily_1997_01_11.dat")
# Select Jan 11 1997

names(data)<-c("Lat","Long","Elev","Precip")
drop<-which(data$Precip=="-999.999" ) #remove NA values "-999.999"
data<-data[-drop,]
data$Rain<-999
N<-dim(data)[1]
for(i in 1:N){
  if(data$Precip[i]>0) {data$Rain[i]<-1} else
  {data$Rain[i]<-0}}
n<-dim(data)[1]

model1<-glm(Rain~Lat+Long+Elev, family=binomial(link="logit"), data=data)
model2<-glm(Rain~Lat+Long+Elev, family=binomial(link="probit"), data=data)
model3<-glm(Rain~Lat+Long+Elev, family=binomial(link="cloglog"), data=data)
bestmodel1<-stepAIC(model1, k=2) #k=2 for AIC, k=log(n) for BIC or SBC
bestmodel2<-stepAIC(model2, k=2)
bestmodel3<-stepAIC(model3, k=2)
bestmodel1<-stepAIC(model1, k=log(n)) #k=2 for AIC, k=log(n) for BIC or SBC
bestmodel2<-stepAIC(model2, k=log(n))
bestmodel3<-stepAIC(model3, k=log(n))

#--> bestmodel1 is the overall best model based on AIC
bestmodel<-bestmodel1 
#glm(formula = Rain ~ Lat + Elev, family = binomial(link = "logit"), data = data)


# Model diagnostics
par(mfrow=c(2,2))
plot(bestmodel)

#ii. Estimate the function on the DEM grade and plot the surface.
predpoints<-read.table(file="/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013/colo_dem.txt")
names(predpoints)<-c("Lat","Long","Elev")

Xpred<-predpoints #gridded lat, long and elev
ypred = predict.glm(bestmodel,newdata=predpoints, type="response",se.fit=TRUE)

# Create spatial plot - Latitude, Longitude and the predicted value..
# since the predpoints are on a unifrom spatial grid but not on 
# rectangular grid we can do the following..

library(akima)
zz = interp(predpoints$Long,predpoints$Lat,ypred$fit)

library(fields)
surface(zz,xlab="Longitude",ylab ="Latitude", main="DEM Grid Predictions") 

# Compare them with the surface plot of the elevation and also with the results from Slater and Clark (2006) Figure 4

#v. Make a 3-D plot (or spatial colored/contour map) of model estimates (i.e. latitude, longitude and the model estimates) and model error. 
# (v) Spatial map of actual precipitation, model estimates, and model error
par(mfrow=c(2,2))
par(mar=c(4, 3, 3, 2) + 0.1) 
library(akima)
library(fields)
zz0<-interp(x=data$Long, y=data$Lat, z=data$Precip, duplicate="mean")
image.plot(zz0, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Observed Precipitation (mm)")
contour(zz0, add=T)
world(add=TRUE, lwd=4)

zz1<-interp(data$Long,data$Lat,bestmodel$fitted.values, duplicate="mean")
image.plot(zz1, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Logit modeled precipitation (mm)")
contour(zz1, add=T)
world(add=TRUE, lwd=4)

zz2<-interp(data$Long,data$Lat,residuals(bestmodel, type="response"), duplicate="mean")
image.plot(zz2, col=topo.colors(20,0.5), ylab="latitude", xlab="longitude", main="Logit residuals (mm)")
contour(zz2, add=T)
world(add=TRUE, lwd=4)

surface(zz,xlab="Longitude",ylab ="Latitude", main="DEM Grid Predictions") 

# excellent!  
# The plot of observed precipitation looks nearly identical to that of Clark and Slater (2006).  My plot of modeled probability of precipitation (GLM~binomial) shares many features of the corresponding figure 4 in Clark and Slater (2006), although the color scales are opposite, making it hard to inspect visually!  

################# End Q8
################# 
################# Start Q9
## repeat Q8 for LOCAL glm
library(locfit)
library(arules)

# The commands below 
# 1. Fits a best "local polynomial" model with appropriate link function
# 2. F-Test of model significance:
#    a. Compute the RSS from the local model and from linear regression
#    b. Compute the F statistic from data and from theory
#    c. If Fdata > Ftheory - reject null hypothesis (i.e, reject that the linear
#       model is a best fit to the data
# 3. Estimate the funciton on the DEM grid and plot the surface
# 4. Compare with Results from Slater and Clark (2006)

# Select Jan 11 1997 precipitation data
data<-read.table(file="http://cires.colorado.edu/~aslater/CVEN_6833/colo_pcp_daily_1997_01_11.dat")

names(data)<-c("Lat","Long","Elev","Precip")
drop<-which(data$Precip=="-999.999" ) #remove NA values "-999.999"
data<-data[-drop,]
data$Rain<-999
N<-dim(data)[1]

# create binary variable: rain / no rain
for(i in 1:N){
  if(data$Precip[i]>0) {data$Rain[i]<-1} else
  {data$Rain[i]<-0}
}

n<-dim(data)[1]
N<-n
#Select independent and dependent variables
X<-data[,1:3]
Y<-data[,5]

#search for best alpha over a range of alpha values between 0,1
#bestparam(X, Y, family="binomial")

##### try this...
# i. Fit a ‘best’ local regression model
################ subset selection ##############
library(leaps)    # to provide combinations
library(MPV)    # to help estimate PRESS and consequently, GCV
library(locfit)
# #Select independent and dependent variables
# X<-data[,1:3] # all the predictor set.
# Y<-data[,4] # response variable Y
# N = length(Y)

combs = leaps(X,Y, nbest=25)  #Get upto 25 combinations of predictor and predictants
combos = combs$which  # logical combinations of predictor variables
ncombos = length(combos[,1]) # number of combinations of predictors variables

for(i in 1:ncombos){
  xx<-as.data.frame(X[,combos[i,]])
  names(xx)<-names(X)[combos[i,]]
  bestparam(X=xx, Y=Y, family="binomial") # find bestalpha and bestdeg for given set of predictor variables
  
  # apply bestalpha and bestdeg to fit the local polynomial
  xx<-as.matrix(xx)
  y<-as.matrix(Y)
  zz<-locfit(y~xx, alpha=bestalpha, deg=bestdeg, family="binomial") 
  # create vector of GCV values for each model with its own bestalpha and bestdeg
  if(i==1) {GCVs <- gcv(zz)[[4]]} else {GCVs <- rbind(GCVs,gcv(zz)[[4]])} 
}

# select the model with the overall lowest GCV and re-fit the "bestmodel"
besti<-which.min(GCVs)            #best combination of predictor variables based on GCV
names<-names(X)[combos[besti,]]     #predictors
X<-as.data.frame(X[,combos[besti,]]) #capital X = df; lowercase x = matrix
names(X)<-names
x<-as.matrix(X)
bestparam(X=X, Y=Y, family="binomial") # alpha=0.07241379, deg=2
bestmodel<-locfit(y~x, alpha=bestalpha, deg=bestdeg, family="binomial") #fit the best model
plot(bestmodel)

dev.copy(png,"Precipitation Probability Spatial Map - Locfit with binomial Link function.png")
dev.off()

#########################
## test goodness of fit..
## The two ways to compute nu1 and nu2 return different results!
RSS1 = sum(residuals(bestmodel)^2)
#nu1 = sum(fitted(bestmodel,what="infl"))     # trace(L) = 410.5157
#nu2 = sum(fitted(bestmodel,what="vari"))     # trace(L^T L) = 319.3362
#nu11 = N- (2*nu1) + nu2 # nu1 = -211.6951
## -Or-
nu1 = bestmodel$dp[6] # Dof1 = 45.86
nu2 = bestmodel$dp[7] # Dof2 = 46.30
nu11 = N- (2*nu1) + nu2  #nu11 = 244.57

## Question for Balaji: 
# (1) WHAT'S THE DISTINCTION BETWEEN DOF1 AND DOF2? 
# (2) Why are nu1 and nu2 different from the two methods above?

## Linear regression..
zz=lsfit(X,Y)
N = dim(Y)[1]
N=length(Y)
XX = cbind(rep(1,N), X)
XX<-as.matrix(XX)

# Compute the Hat matrix
hatm = XX %*% solve(t(XX) %*% XX) %*% t(XX)

II = diag(N)
delta0 = t(II-hatm)%*%(II-hatm)    #Equation 9.2
nu00 = sum(diag(delta0))
RSS0 = sum(residuals(zz)^2)

# Ftest
# H0: Local Polynomial = Linear Regression
# H1: Local Polynomial != Linear Regression
# if Fdata > Ftheory, reject the null

Fdata = ((RSS0 - RSS1)/(nu00 - nu11))/(RSS1/nu11)
#Fdata = (RSS0 - RSS1)/(nu00 - nu11)
#Fdata = (Fdata / (RSS1 / nu11))
Ftheory = qf(0.95,(nu00-nu11), nu11)    #95% confidence level..
Fdata
Ftheory
Fdata>Ftheory
## IF Fdata > Ftheory - reject null - i.e., reject that the data fits a linear model
## Results: 
## Fdata>Ftheory # FALSE - fail to reject the null
## Therefore no need to use local polynomial instead of linear regression... (this seems wrong!!!)

#ii. Estimate the function on the DEM grade and plot the surface.
predpoints<-read.table(file="/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013/colo_dem.txt")
names(predpoints)<-c("Lat","Long","Elev")

# TO PREDICT... MUST ASSIGN PREDICTION POINTS TO SAME NAME AND STRUCTURE AS USED TO FIT THE MODEL!  
x<-predpoints #gridded lat, long and elev
x<-as.matrix(x)
ypred = predict(bestmodel, newdata=x, type="response",se.fit=TRUE)
ypred2 = predict(bestmodel, newdata=x, type="response", se.fit=FALSE, family="binomial", link="logit")


# Create spatial plot - Latitude, Longitude and the predicted value..
# since the predpoints are on a unifrom spatial grid but not on 
# rectangular grid we can do the following...
par(mfrow=c(2,2))
par(mar=c(3,3,3,2)+0.1)

zz0<-interp(x=data$Long, y=data$Lat, z=data$Precip, duplicate="mean")
image.plot(zz0, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Observed Precipitation (mm)")
contour(zz0, add=T)
world(add=TRUE, lwd=4)

zz = interp(predpoints$Long,predpoints$Lat,as.matrix(ypred$fit))
image.plot(zz, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Probability of Precipitation\nDEM Grid Predictions")
#contour(zz, add=T)
world(add=TRUE, lwd=4)

#test = interp(predpoints$Long,predpoints$Lat,ypred2)
#surface(test,xlab="Longitude",ylab ="Latitude", main="Probability of Precipitation\nDEM Grid Predictions") 

zz1 = interp(predpoints$Long,predpoints$Lat,as.matrix(ypred$se.fit))
image.plot(zz1, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Standard Errors\n on DEM Predictions")
contour(zz1, add=T)
world(add=TRUE, lwd=4)

# Compare them with the surface plot of the elevation and also with the results from Slater and Clark (2006) Figure 4
zz2<-interp(predpoints$Long, predpoints$Lat, predpoints$Elev )
image.plot(zz2, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Elevation Surface")
contour(zz2, add=T)
world(add=TRUE, lwd=4)

dev.copy(png,"Observed vs. Modeled Precip and Standard Errors\nLocfit with binomial link function.png")
dev.off()

################# End Q9
################# 
################# Start Q10
# hierarchical spatial model with locfit.
# repeat problem 7 but using the local polynomial.
# (basically, model 3 but using your local polynomial model from problem 4)

# select a rho and theta from visual inspection of the variogram boxplots.
# Also don't bother doing a 10% drop cross validation (if you have not already done) - leave one out cross validation is fine.
# If you did the 10% cross validation that is great, include them in your submission.

##################################
#Import data: Mean annual precipitation at 491 locations in Colorado based on data for the period 1980-2002
data<-read.table(file="/Users/elliotcohen/Dropbox/Data/Climate/Rainfall/Colorado_Annual_Precip_Gridded.csv", sep=",", header=TRUE)
data<-as.data.frame(data)
X<-data[,-4]  # df of Independent variable
Y<-data$Precip #vector of dependent variable (annual rainfall at lat-long position)
n<-length(Y)
nvar<-dim(X)[2]

xs<-X[,1:2] # lat, long
xse<-X[,] # lat, long and elev

### In Q7 we apply GLM to get the mean function + Kriging on the residuals
## For Q10 we try locfit + kriging on the residuals

#### Find the "best" local model
#### subset selection 
# combs = leaps(X,Y, nbest=25)  #Get upto 25 combinations of predictor and predictants
# combos = combs$which  # logical combinations of predictor variables
# ncombos = length(combos[,1]) # number of combinations of predictors variables
# 
# for(i in 1:ncombos){
#   xx<-as.data.frame(X[,combos[i,]])
#   names(xx)<-names(X)[combos[i,]]
#   bestparam(X=xx, Y=Y, family="gaussian") # find bestalpha and bestdeg for ith combination of predictor variables
#   
#   # apply bestalpha and bestdeg to fit the local polynomial
#   xx<-as.matrix(xx)
#   y<-as.matrix(Y)
#   zz<-locfit(y~xx, alpha=bestalpha, deg=bestdeg) 
#   # create vector of GCV values for each model with its own bestalpha and bestdeg
#   if(i==1) {GCVs <- gcv(zz)[[4]]} else {GCVs <- rbind(GCVs,gcv(zz)[[4]])} 
# }

# # select the model with the overall lowest GCV and re-fit the "bestmodel"
# besti<-which.min(GCVs) #best combination of predictor variables based on GCV

## lowest GCV is for model with Long only but I want to retain Lat and Long for kriging... so choose the second-lowest GCV, which has both Lat and Long....
besti<-4

names<-names(X)[combos[besti,]] #predictors
X<-as.data.frame(X[,combos[besti,]]) #capital X = df; lowercase x = matrix
names(X)<-names
x<-as.matrix(X)
y<-as.matrix(Y)
bestparam(X=X, Y=Y,family="gaussian") #alpha=0.0683299, deg=1
bestmodel<-locfit(y~x, alpha=bestalpha, deg=bestdeg) #fit the best model

# Compute useful quantities for locfit...
modresid=residuals(bestmodel) # Get the residuals of the model..
nX<-dim(X)[2]
Yhat=Y-modresid  #residuals = Y - Yestimate ==> Yestimate = Y - residuals
k=dim(X)[2]      #number of regressor variables
p=k+1             #number of model parameters 
n=length(Y)       #number of observations

#### Now let Y be the residuals of the best locfit model...
Yresid = residuals(bestmodel) # let Y be the residuals of the glm model....
nobs = length(Yresid)

par(mfrow=c(1,1))
par(mar=c(5,4,4,2)+0.1)
sigma=1
look<- vgram(xs, Yresid, N=15, lon.lat=FALSE)
bplot.xy(look$d, look$vgram, breaks= look$breaks,
         outline=FALSE,
         xlab="distance (degrees)", ylab="variogram")
points( look$centers, look$stats[2,], col="blue")

# fit of exponential by nonlinear least squares
xd<- look$centers
ym<- look$stats[2,]
sigma<- 1.0

# nls.fit<- nls( ym ~ sigma^2 + rho*( 1- exp(-xd/theta)),
#                start= list(rho=max(look$stats[2,]), theta=quantile(xd,0.1)), control=list( maxiter=200) )
# pars<- coefficients(nls.fit)
# rho<-pars[1]
# theta<-pars[2]

# Choose rho and theta from visual inspection of the variogram....
rho<-38000
theta<-0.12
xr=round(max(xd)+sd(xd))
dgrid<- seq(0,xr,,400)
lines(dgrid, sigma^2 + rho*(1 - exp(-1*dgrid/theta)), col="blue", lwd=3)


### Predict at observation points.. if sigma = 0 then it will be exact.
K = sigma^2 + rho*(1 - exp(-1*rdist(xs,xs)/theta))
#Kstar<- Matern(rdist( xgrid, x)/.3, smoothness=1.5)
Kstar = K
nobs=length(Yresid)
ghat.krig<- Kstar %*% solve( K + diag( sigma^2, nobs)) %*% Yresid
ghat.locfit<-predict(bestmodel, newdata=x)
ghat = ghat.krig + ghat.locfit

## Check with Krig and predict.krig
zz<-Krig(xs,Yresid,rho=rho,theta=theta,m=1,sigma2=sigma)
yp.krig<-predict.Krig(zz)
yp.locfit<-predict(bestmodel, newdata=x)
yp <- yp.krig + yp.locfit

compare<-cbind(ghat,yp)
head(compare)

# ghat and yp are identical to the tens place.
identical(round(ghat,-1), round(yp,-1))

## Predict on the grid.. and standard error..
# xs contains 491 lat-long coordinates from the observed data
# xps contains 1705 lat-long coordinates from the DEM grid
x<-xps # re-assign xps to x....
Kstar = sigma^2 + rho*(1 - exp(-1*rdist(xps,xs)/theta))
ghat.krig<- Kstar %*% solve( K + diag( sigma^2, nobs)) %*% Yresid
kse = sqrt(var(Y) - diag((Kstar) %*% solve( K + diag( sigma^2, nobs)) %*% t(Kstar)))
ghat.locfit<- predict(bestmodel, newdata=as.matrix(x))
ghat = ghat.krig + ghat.locfit

## Check with Krig and predict.Krig
zz=Krig(xs,Yresid,rho=rho,theta=theta,sigma2=sigma,m=1)
yp.krig<-predict.Krig(zz, x=x)
yp.locfit<-predict(bestmodel, newdata=as.matrix(x))
yp <- yp.krig + yp.locfit
#yp = predict(bestmodel,newdata=predpts) + predict.Krig(zz,x=xps,drop.Z=TRUE)

compare<-cbind(ghat,yp)
head(compare)

# ghat and yp are identical to the hundreds (*not hundredths*) place.
identical(round(ghat,-2), round(yp,-2))

# ## Now check the cross-validated estimates of Y...
# crossval.Krig(X=xs, Y=Y, rho=rho, theta=theta, sigma=sigma)
# dev.copy(png,paste(zz$call[1],"Cross-Validation_Model_4.png", sep="_"))
# dev.off()

### spatial map of estimates and errors
xps<-as.data.frame(xps)
xlon = sort(unique(xps$Long))
nrs = length(xlon)
ylat = sort(unique(xps$Lat))
ncs = length(ylat)

zmat = matrix(ghat, nrow=ncs, ncol= nrs) 
image.plot(xlon,ylat,t(zmat),zlim=range(Y),col=topo.colors(20,0.5), ylab="latitude", xlab="longitude", main="Precip Predictions on DEM (mm/yr)\nHierarchical Spatial Model: Kriging + locfit")
contour(xlon,ylat,t(zmat),add=T)
dev.copy(png,paste(zz$call[1],"Model_4_DEM_Preidction_Plot.png", sep="_"))
dev.off()

zmat1 = matrix(kse, nrow=ncs, ncol= nrs) 
image.plot(xlon,ylat,t(zmat1), col=topo.colors(20,0.5), ylab="latitude", xlab="longitude", main="Model 4 standard errors")
contour(xlon,ylat,t(zmat1),add=T)
dev.copy(png,paste(zz$call[1],"Model_4_Standard_Errors.png", sep="_"))
dev.off()