# September 17, 2013
# CVEN6833 - Advance Data Analysis
################# 
################# 
################# Start Q1

setwd("/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013")
#Import data: Mean annual precipitation at 491 locations in Colorado based on data for the period 1980-2002
data<-read.table(file="/Users/elliotcohen/Dropbox/Data/Climate/Rainfall/Colorado_Annual_Precip_Gridded.csv", sep=",", header=TRUE)
data<-as.data.frame(data)
X<-data[,-4]  # df of Independent variable
Y<-data$Precip #vector of dependent variable (annual rainfall at lat-long position)

# i. Fit a ‘best’ linear regression model (use one of the objective functions - GCV, AIC, SBC or PRESS; you can also try a couple of them to see any differences). This entails fitting the model with all possible combinations of covariates and selecting the model with the minimum objective function.
lm<-lm(Precip ~ Lat + Long + Elev, data=data)
summary(lm)

# Load the MASS package to use the stepAIC function to choose the "best model"
require(MASS)
bestmodel<-stepAIC(lm)
nX<-length(bestmodel$coefficients)
X<-subset(data, select=names(bestmodel$coefficients)[2:nX]) #subset X based on best model regressor variables..

#ii. Perform ANOVA (i.e. model significance) and model diagnostics (i.e., check the assumptions of the residuals – Normality, independence, homoskedasticity). 
summary(bestmodel) #"bestmodel" has longitutde as only predictor variable
#Note: other models have similar AIC values, indicating that "best" model may not be that much better than the next best model.

# Model diagnostics
par(mfrow=c(2,2))
plot(bestmodel)
dev.copy(png,paste(bestmodel$call[1],"model_diagnostics", sep="_"))
dev.off()

# Results: visual inspection of plot 1 finds that residuals increase for larger values of the response variable yhat, indicating that the model does *not* capture certains attributes of the data (perhaps non-linearities); Plot 2 (Q-Q plot) shows that residuals are *not* distributed normally at the upper tail. In summary, at least 2 of the model diagnostics suggest that the linear model is inadequate for this data set.

# can also use custom-built  model diagnostics...
####### Manual Model Diagnostics (6 visual checks) ############
GLMdiagnostics<- function(bestmodel, X, Y){
  
  # Get the residuals of the model..
  modresid=residuals(bestmodel)
  nX<-length(bestmodel$coefficients)
  Yhat=Y-modresid  #residuals = Y - Yestimate ==> Yestimate = Y - residuals
  k=dim(X)[2]      #number of regressor variables
  p=k+1            #number of model parameters 
  n=length(Y)
  
  # Compute ANOVA quantities for use down below
  #SST = Total corrected sum of squares, n-1 dof
  #SSR = Regression Sum of Squares = sum[(yhati-ybar)^2], dof = k (number of predictor  variables)
  #SSE = Error Sum of Squares = sum[(yi-yhati)^2], dof = n-p 
  #Yhat = Y - bestGlm$res 
  #(Y - Yhat = residuals), Yhat is the modeled response of Y
  SST = sum((Y - mean(Y))^2)   
  SSR = sum((Yhat - mean(Y))^2)  
  SSE = sum((modresid)^2)   
  MSR = SSR / ncol(X)            
  MSE = SSE/(n - length(bestmodel$coef))         
  
  # Now start computing diagnostics and plotting them...
  par(mfrow=c(2,p)) 
  
  # (1) Check if residuals fit a normal distribution
  qqnorm(modresid)
  qqline(modresid)	
  
#   jpeg(filename=paste(bestmodel$call[1],"Q-Q plot", sep="_"))
#   plot(bestmodel)
#   dev.off()
  
  # (2-3) Plot the residuals vs X.  Check to make sure there is *no* apparent pattern.  Distribution of residuals should be random.
  for(i in 1:k){
    plot(X[,i],modresid,xlab="X",ylab="residuals",main="Residuals vs. X[,i]") 
  }
  
  # (4) Plot the residuals vs the model estimates of Y. 
  #Check to make sure there is *no* apparent pattern or structure.  In other words, the distribution of the residuals should look random.
  plot(Yhat,modresid,xlab="estiimate of Y", ylab="residuals",main="Residuals vs Fitted Y")
  
  # (5) Plot the autocorrelation function - to make sure the residuals are *not* related to each other.  
  z1=acf(modresid,main="autocorrelation of residuals")
  
#   # (6) Cooks Distance - to make sure outliers do not exert undue influence on the regression model.
#   
#   # Compute the Hat matrix
#   #hatm= hatvalues(bestmodel)
#   XX<-cbind(rep(1,n),X) # augmented X matrix
#   XX<-as.matrix(XX)
#   hatm<-XX %*% solve(t(XX) %*% XX) %*% t(XX)
#   
#   #studentized residuals - ri  - equation 12-42
#   # ri = modresid/sqrt((1 - diag(hatm)) * MSE) #if using hatvalues(bestmodel)
#   ri = modresid/sqrt((1-hatm) * MSE)
#   #Compute Cook's Distance Equation 12-44
#   Di = ri*ri * diag(hatm) / ((1-diag(hatm)) * length(bestmodel$coef))
#   plot(Y, diag(Di), main="Cook's Distance")
#   #If Dis are greater than 1 then there are points that have undue influence on the fit.

}
###############################
# Now impliment my model diagnostics function
GLMdiagnostics(bestmodel, X, Y)
dev.copy(png,paste(bestmodel$call[1],"model_diagnostics", sep="_"))
dev.off()

# Comments on Model Diagnostics: 
# The Q-Q plot reveals non-normality at the upper tail
# Plot of residuals vs model estimates of Y reveals heteroscdasticity. 
#### END MODEL DIAGNOSTICS ####
#### END PART (ii) ####

# iii. Compute cross-validated and fitted estimates at each observation points and plot them against the observed values. This is to visually see how the model performs in a fitting and cross-validated mode.

## Observed vs modeled estimates of response variable Y
par(mfrow=c(1,1))
#preds<-data.frame(c(paste("data$",names(bestmodel$coefficients)[2:nX],sep="")), quote=FALSE)  #trying to automate regressor variables to include in preds
preds<-data.frame(data$Long)
names(preds)<-c("Long")
Fitted<-predict(bestmodel, newdata=preds, type="response") #predict precip using best linear model with Long as only predictor variable.

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

# save plot
dev.copy(png,paste(bestmodel$call[1],"cross-validated estimates", sep="_"))
dev.off()

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
  #x<-as.data.frame(X[keep,])
  x<-X[keep,]
  xx<-as.data.frame(x)
  y<-Y[keep]
  zz<-lm(y~x, data=xx) #fit model to remaining data
  xpred<-as.data.frame(X[drop,])
  xx<-as.data.frame(xpred) #assign dropped data to xx
  yhat<-predict.lm(zz, newdata=xx) #predict at dropped points using model fit w.out those points
  #Warning message:'newdata' had 49 rows but variable(s) found have 442 rows 
  rmseskill[i]<-sqrt(mean(Y[drop]-yhat)^2)
  corskill[i]<-cor(Y[drop],yhat)
}

###
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

# save plot
#dev.copy(png,paste(bestmodel$call[1],"RMSE and COR skill", sep="_"))
#dev.off()

#v. Make a 3-D plot (or spatial colored/contour map) of model estimates (i.e. latitude, longitude and the model estimates) and model error. 
# (v) Spatial map of actual precipitation, model estimates, and model error
plot.new()
par(mfrow=c(2,2))
par(mar=c(5,4,4,2) + 0.1) #A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. The default is c(5, 4, 4, 2) + 0.1.

library(akima)
library(fields)
zz0<-interp(x=data$Long, y=data$Lat, z=data$Precip, duplicate="mean")
image.plot(zz0, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Observed Precipitation (mm)")
contour(zz0, add=T)
world(add=TRUE, lwd=4)

zz1<-interp(data$Long,data$Lat,bestmodel$fitted.values, duplicate="mean")
image.plot(zz1, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="LM Fitted Precipitation (mm)")
contour(zz1, add=T)
world(add=TRUE, lwd=4)

zz2<-interp(data$Long,data$Lat,residuals(bestmodel, type="response"), duplicat="mean")
image.plot(zz2, col=topo.colors(20,0.5), ylab="latitude", xlab="longitude", main="LM Residuals (mm/yr)")
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
ypred = predict.lm(bestmodel,newdata=predpoints, type="response",se.fit=TRUE)

# Create spatial plot - Latitude, Longitude and the predicted value..
# since the predpoints are on a unifrom spatial grid but not on 
# rectangular grid we can do the following..

library(akima)
zz = interp(predpoints$Long,predpoints$Lat,ypred$fit)

library(fields)
surface(zz,xlab="Longitude",ylab ="Latitude", main="DEM Grid Predictions")

# save plot
dev.copy(png,paste(bestmodel$call[1],"Contour Plots", sep="_"))
dev.off()

# Discussion
# LM does not capture non-linearities and complexity inherent in the data
# LM attempts to "spread" precipitation evenly across the state, decreasing from West to East.  As can be seen in the LM residuals plot, the residuals explain most of the variance in the data, not the model, indicating that the model is inadequate.  

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