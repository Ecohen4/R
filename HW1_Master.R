# October 17, 2013
# CVEN6833 - Advance Data Analysis
# HW 1

setwd("/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013/Final")
source("/Users/elliotcohen/Dropbox/Advance Data Analysis/R source files/myboxplot-stats.r")
source("/Users/elliotcohen/Dropbox/Advance Data Analysis/R source files/myboxplot.r")

require(MASS)
library(arules)
library(locfit)   # local polynomial regression 
library(akima)    # for interp function
library(fields)   # for surface function
library(leaps)    # to provide combinations
library(MPV)       # to help estimate PRESS and consequently, GCV
library(fitdistrplus)    # fitting a distribution to data


# Define custom fucntions.....

#### Model Diagnostics (6 visual checks) ####
##### Begin Function #####
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
  if(p>2)par(mfrow=c(3,3)) else par(mfrow=c(2,2))
  
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
  #   #If Dis are greater than 1 then there are points that have undue influence on the fit 
} #### close function ####



#### Best alpha and best polynomial order for locfit ####
##### Begin Function #####
# Custom function to search for the best alpha (between 0 and 1) and best polynomial order (1 or 2) for the local polynomial fit.

bestparam<-function(X, Y, family){
  
  N = length(Y)
  nvar=dim(X)[2]
  porder=1
  minalpha=3*(nvar*porder+1)/N  #try 3x instead of 2x....
  alpha1=seq(minalpha, 1.0, by=0.05)
  n=length(alpha1)
  
  porder=2
  minalpha=3*(nvar*porder+1)/N 
  alpha2=seq(minalpha, 1.0, by=0.05)
  alpha=c(alpha1,alpha2)
  
  #gcvplot accepts matrix/vector input....
  y<-as.matrix(Y)
  x<-as.matrix(X)
  #dimnames(xx)<-list(rep("obs",N),names(X)[combos[besti,]])
  zz<-gcvplot(y ~ x, alpha=alpha1, deg=1, kern='bisq',family=family, ev=dat(), scale=TRUE)
  z1<-gcvplot(y ~ x, alpha=alpha2, deg=2, kern="bisq", family=family, ev=dat(), scale=TRUE)
  
  # pick the best alpha and the degree of the polynomial that gives the least GCV
  z2=order(c(zz$values,z1$values)) #order from lowest GCV to highest
  
  deg1=1
  if(z2[1] > n)deg1=2 #select degree (1 or 2) that yields lowest GCV
  bestalpha<<-alpha[z2[1]] #select alpha that yields lowest GCV and assign to global environment
  bestdeg<<-deg1 #select polynomial degree that yields lowest GCV and assign to global environment
  print(bestalpha)
  print(bestdeg)
}  #### close function ####



#### Cross-Validated and Fitted Estimates ####
##### Begin Function #####
crossval<-function(X, Y, xpred){
  
  Fitted<-predict(bestmodel, newdata=xpred, type="response")
  
  par(mfrow=c(1,1))
  plot(Y, Fitted, pch=20, col="black", xlab="Observed Precip (mm)",ylab="Modeled Precip (mm)", main="Observed vs Modeled Precipitation")
  abline(a=0, b=1)
  
  # cross validated estimates (drop observations one at a time, refit the model to the remaining data (N-1) and predict at the dropped point; repeat for all observations)
  par(mfrow=c(1,1))
  n=length(Y)
  yest=1:n
  test=1:n
  nvar=dim(X)[2]
  index=1:n
  
  for(i in 1:n){
    index1=index[index != i] #drop one observation at a time
    Xval=X[index1,] #X data less the dropped observation
    Yval=Y[index1]  #Y data less the dropped observation
    newdf<-as.data.frame(cbind(Yval, Xval))
    names(newdf)<-c("Yval",names(X))
    
    
    # fit the model without the dropped observation...
    if(bestmodel$call[1]=="lm()"){
      Xval<-as.matrix(Xval) 
      zz=lm(Yval ~ Xval) #re-fit the model without the dropped observation
      xpred=c(1,X[i,1:nvar]) #now estimate at the point that was dropped
      xpred<-as.numeric(xpred)
      yest[i]=sum(zz$coef * xpred)
      #test[i]<-xpred %*% solve(t(xpred) %*% xpred) %*% t(xpred) %*% Yval[i]
      
    } else if (bestmodel$call[1]=="locfit()"){
      zz<-locfit(Yval ~ Xval, alpha=bestalpha, deg=bestdeg) 
      xpred<-X[i,1:nvar] #now estimate at the point that was dropped
      xpred<-as.numeric(xpred)
      yest[i]<-predict(zz, newdata=xpred, type="response")
      
    } else if (bestmodel$call[1]=="glm()") {
      zz<-glm(Yval~Xval, family=bestmodel$family)
      xpred=X[i,1:nvar]
      xpred<-as.numeric(xpred)
      #yest[i]<-predict(zz, newdata=xpred, type="response")
      # or...
      #xpred=c(1,X[i,1:nvar]) #now estimate at the point that was dropped
      #yest[i]=sum(zz$coef * xpred)
      yest[i]<-xpred %*% solve(t(xpred) %*% xpred) %*% t(xpred) %*% Yval[i]
      
    }
  }  ## close for loop
  
  #now surface the x-validated estimates..
  points(Y, yest, col="blue", pch=20)
  
} #### close function ####



#### Cross-Validated and Fitted Estimates for Krigging #####
##### Begin Function #####
crossval.Krig<-function(X, Y, rho, theta, sigma){
  
  # fitted values... 
  X<-xs # assign subset of predictor variables (lon, lat) to X
  zz=Krig(X,Y,rho=rho,theta=theta,sigma2=sigma,m=1)
  Fitted = predict.Krig(zz,x=X,drop.Z=TRUE) # predict at the same lat-long coords  
  par(mfrow=c(1,1))
  plot(Y, Fitted, pch=20, col="black", xlab="Observed Precip (mm)",ylab="Modeled Precip (mm)", main="Observed vs Modeled Precipitation")
  abline(a=0, b=1)
  # Krigin is an "exact esimator" so returns the observed values...
  
  # cross validated estimates (drop observations one at a time, refit the model to the remaining data (N-1) and predict at the dropped point; repeat for all observations)
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
  legend("bottomright", legend=c("Krig Package", "First Principles"), col=c("green","blue"), pch=c(20,0), bty="n")
  
}  ############# Close function  ########



##### Simulated RMSE and Correlation (a.k.a. "droptest") #####
##### Begin Function #####
## Drop some % of points, fit the model and predict the dropped points..
droptest<-function(X, Y, drop, bestmodel){
  
  source("/Users/elliotcohen/Dropbox/Advance Data Analysis/R source files/myboxplot.r")
  source("/Users/elliotcohen/Dropbox/Advance Data Analysis/R source files/myboxplot-stats.r")
  
  nsim = 500
  rmseskill=1:nsim
  corskill=1:nsim  
  N = length(Y)
  N10 = round(drop*N)    #choose % of points to drop (e.g. 10%)
  index=1:N
  
  for(i in 1:nsim){
    drop=sample(c(1:N),N10)  #sample 10% of the row indices from 1:N at random
    keep=setdiff(index,drop)  #discard values at the intersection of index and drop (e.g. drop 10% of the data, keep the rest)
    xx<-X[keep,]
    yy<-Y[keep]
    ydrop<-Y[drop]
    xpred<-X[drop,]    #the dropped data
    #xpred<-as.data.frame(xpred)
    
    if(bestmodel$call[1]=="lm()"){
      zz=lm(yy ~ xx)       #re-fit the model without the dropped points
      yhat=zz$coef[1] + zz$coef[2]*xpred  #estimate at the dropped points
      
    } else if (bestmodel$call[1]=="locfit()"){
      zz<-locfit(yy ~ xx, alpha=bestalpha, deg=bestdeg) 
      yhat<-predict(zz, newdata=xpred, type="response")
      
    } else if (bestmodel$call[1]=="glm()") {
      zz<-glm(yy ~ xx, family=bestmodel$family)
      #xpred<-cbind(rep(1,N10), X[drop,])
      #xpred<-as.numeric(xpred)
      yhat<-xpred %*% solve(t(xpred) %*% xpred) %*% t(xpred) %*% Y[drop]
      #yhat<-predict(zz, newdata=xpred, type="response")
    }
    
    rmseskill[i]<-sqrt(mean(((Y[drop]-yhat)^2)))
    corskill[i]<-cor(Y[drop],yhat)
  }
  
  par(mfrow=c(1,2))
  boxplot(rmseskill, main="Simulated RMSE") #simple version
  boxplot(corskill, main="Simulated Cor." )  #simple version
  
  zz=myboxplot(rmseskill, main="Simulated RMSE skill",plot=FALSE)
  zz$names=rep("",length(zz$names))
  z1=bxp(zz,xlab="",ylab="RMSE",cex=1.25)
  
  if(bestmodel$call[1]=="lm()"){
    modresid = Y - bestmodel$fitted
    rmse<-sqrt(sum(modresid^2)/N)
  }
  if(bestmodel$call[1]=="glm()"){
    modresid = Y - bestmodel$fitted.values
    rmse<-sqrt(sum(modresid^2)/N)
  }
  if(bestmodel$call[1]=="locfit()"){
    modresid = residuals(bestmodel)
    rmse<-sqrt(sum(modresid^2)/N)
    #rmse<-sqrt(sum(residuals(bestmodel)^2)/N)
  }
  
  points(z1,rmse,col="red",cex=2,pch=19)  #add a point showing the true RMSE of the data.
  title(main="RMSE skill")
  
  zz=myboxplot(corskill, main="Simulated Correlation skill",plot=FALSE)
  zz$names=rep("",length(zz$names))
  z1=bxp(zz,xlab="",ylab="Cor",cex=1.25)
  title(main="Cor skill")
  
} #### close function ####



################# 
################# 
################# Start Q1
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

#### Model diagnostics ####
par(mfrow=c(2,2))
plot(bestmodel)
dev.copy(png,paste(bestmodel$call[1],"model_diagnostics.png", sep="_"))
dev.off()

# Results: visual inspection of plot 1 finds that residuals increase for larger values of the response variable yhat, indicating that the model does *not* capture certains attributes of the data (perhaps non-linearities); Plot 2 (Q-Q plot) shows that residuals are *not* distributed normally at the upper tail. In summary, at least 2 of the model diagnostics suggest that the linear model is inadequate for this data set.

# can also use custom-built model diagnostics...
GLMdiagnostics(bestmodel, X, Y)
dev.copy(png,paste(bestmodel$call[1],"model_diagnostics", sep="_"))
dev.off()

# Comments on Model Diagnostics: 
# The Q-Q plot reveals non-normality at the upper tail
# Plot of residuals vs model estimates of Y reveals heteroscdasticity. 
#### END MODEL DIAGNOSTICS ####

# iii. Compute cross-validated and fitted estimates at each observation points and plot them against the observed values. This is to visually see how the model performs in a fitting and cross-validated mode.
crossval(X=X,Y=Y,xpred=X)

# save plot
dev.copy(png,paste(bestmodel$call[1],"cross-validated estimates", sep="_"))
dev.off()

# Comments:  Cross-validated model estimates and estimates modeled on the full data set appear internally consistent to each other, but niether are accurate to observed precipitation values for larger values of precip (e.g. above 800 mm/yr).

# iv. Drop 10% of observations, fit the model (i.e., the ‘best’ model from i. above) to the rest of the data and predict the dropped points. Compute RMSE and R2 and show them as boxplots. 
droptest(X=X, Y=Y, drop=0.1, bestmodel=bestmodel)

# save plot
dev.copy(png,paste(bestmodel$call[1],"Simulated_RMSE_COR.png", sep="_"))
dev.off()

#v. Make a 3-D plot (or spatial colored/contour map) of model estimates (i.e. latitude, longitude and the model estimates) and model error. 
# (v) Spatial map of actual precipitation, model estimates, and model error
png(filename=paste(bestmodel$call[1],"ContourPlots.png" ,sep="_"), width=800, height=600, units="px", pointsize=18)

par(mfrow=c(2,2))
par(mar=c(4,4,4,6)+0.1) 
# numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. The default is c(5, 4, 4, 2) + 0.1.
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
#dev.copy(png,paste(bestmodel$call[1],"Contour Plots", sep="_"))
#dev.off()
dev.off()

# Discussion
# LM does not capture non-linearities and complexity inherent in the data
# LM attempts to "spread" precipitation evenly across the state, decreasing from West to East.  As can be seen in the LM residuals plot, the errors explain most of the variance in the data, not the model, indicating that the model is inadequate.  

# See plots
################# End Q1
################# 
################# Start Q2
## repeat Q1 with GLM....
setwd("/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013/Final")
#Import data: Mean annual precipitation at 491 locations in Colorado based on data for the period 1980-2002
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
bestmodel<-stepAIC(glm)
np<-length(bestmodel$coefficients) # number of model paramenters
# X<-X[,names(bestmodel$coefficients[2:np])]
# X<-as.data.frame(X)
X<-subset(X, select=c(names(bestmodel$coefficients[2:np])))


#ii. Perform ANOVA (i.e. model significance) 
summary(bestmodel) #"bestmodel" has longitutde as only predictor variable
#Note: other models have similar AIC values, indicating that "best" model may not be that much better than the next best model.
#... and model diagnostics (i.e., check the assumptions of the residuals – Normality, independence, homoskedasticity). 

par(mfrow=c(2,2))
plot(bestmodel)
dev.copy(png,paste(bestmodel$call[1],"built-in_model_diagnostics", sep="_"))
dev.off()
# Check if key assumptions are met, namely:
# (1) Residuals are "identically and independently distributed" (iid), e.g. their distribution does not change substantially for different values of x.
# To check iid, we look at Plot 1 (upper left) and plot 3 (lower left).
# --> Plot 1 shows a trend of decreasing variance in the residuals for larger values of the response variable, thus violating IID.
# (2) Residuals are normally distributed (Q-Q plot, upper right)
# --> QQ plots looks good
# Plot 4 (lower right) helps identify outliers (their row numbers are shown such that we can go back and find them in the dataframe).  Outliers are not an explicit part of the assumptions required for regression, but are important to be aware of.

# Can also try custom model diagnostics...
GLMdiagnostics(bestmodel, X, Y)
dev.copy(png,paste(bestmodel$call[1],"my_model_diagnostics", sep="_"))
dev.off()

# Comments on Model Diagnostics: 
# The Q-Q plot looks good
# Plot of residuals vs model estimates of Y reveals EXTREME heteroscdasticity.
# Autocorrelation appears to be present up to lag 3.
# --> glm with gamma distribution is *not* adequate for precip data.
#### END MODEL DIAGNOSTICS ####

# iii. Compute cross-validated and fitted estimates at each observation points and plot them against the observed values. This is to visually see how the model performs in a fitting and cross-validated mode.
crossval(X=X, Y=Y, xpred=X)
dev.copy(png,paste(bestmodel$call[1],"crossval.png", sep="_"))
dev.off()

#dev.copy(png,paste(bestmodel$call[1],"cross-validated estimates", sep="_"))
#dev.off()

# Comments:  Cross-validated model estimates and estimates modeled on the full data set appear inconsistent to one another and niether are accurate to observed precipitation values for larger values of precip (e.g. above 800 mm/yr).

# iv. Drop 10% of observations, fit the model (i.e., the ‘best’ model from i. above) to the rest of the data and predict the dropped points. Compute RMSE and R2 and show them as boxplots. 
droptest(X=X,Y=Y, drop=0.1, bestmodel)

dev.copy(png,paste(bestmodel$call[1],"RMSE_Skill", sep="_"))
dev.off()

#v. Make a 3-D plot (or spatial colored/contour map) of model estimates (i.e. latitude, longitude and the model estimates) and model error. 
# (v) Spatial map of actual precipitation, model estimates, and model error

png(filename=paste(bestmodel$call[1],"ContourPlots.png" ,sep="_"), width=800, height=600, units="px", pointsize=18)
par(mfrow=c(2,2))
par(mar=c(4,4,4,6)+0.1) 

zz0<-interp(x=data$Long, y=data$Lat, z=data$Precip, duplicate="mean")
image.plot(zz0, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Observed Precipitation (mm)")
contour(zz0, add=T)
world(add=TRUE, lwd=4)

zz1<-interp(data$Long,data$Lat,bestmodel$fitted.values, duplicate="mean")
image.plot(zz1, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="GLM Modeled precipitation (mm)")
contour(zz1, add=T)
world(add=TRUE, lwd=4)
glmPrecip<-zz1

zz2<-interp(data$Long,data$Lat,residuals(bestmodel, type="response"), duplicat="mean")
image.plot(zz2, col=topo.colors(20,0.5), ylab="latitude", xlab="longitude", main="GLM residuals (mm/yr)")
contour(zz2, add=T)
world(add=TRUE, lwd=4)
glmResid<-zz2

# GLM is clearly *NOT* adequate for modeling precipitation.

# vi. Estimate the precipitation and the standard errors on the DEM grid
#predpoints = read.table("http://cires.colorado.edu/~aslater/CVEN_6833/colo_dem.dat",sep=",")
predpoints<-read.table(file="/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013/colo_dem.txt")
names(predpoints)<-c("Lat","Long","Elev")
Xpred<-predpoints #gridded lat, long and elev

# (vi) Estimate the precipitation and the SEs on the DEM grid and show as 3-D plots as in (v) above.
# Predict response variable (Precip) on the grid of predictor variables (long, Lat, Elev)
ypred = predict.glm(bestmodel,newdata=predpoints, type="response",se.fit=TRUE)

# Create spatial plot - Latitude, Longitude and the predicted value..
# since the predpoints are on a unifrom spatial grid but not on 
# rectangular grid we can do the following..

library(akima)
zz = interp(predpoints$Long,predpoints$Lat,ypred$fit)
library(fields)
surface(zz,xlab="Longitude",ylab ="Latitude", main="DEM Grid Predictions") 

#dev.copy(png,paste(bestmodel$call[1],"Contour_Plot.png", sep="_"))
dev.off()
# Discussion
# GLM does not capture non-linearities and complexity inherent in the data
# GLM attempts to "spread" precipitation evenly across the state, decreasing from West to East.  As can be seen in the LM residuals plot, the errors explain most of the variance in the data, not the model, indicating that the model is inadequate.

################# End Q2
################# 
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
  
  ####
  par(mfrow=c(1,1))
  plot(y, yhat, pch=20, col="blue", xlab="Observed Y", ylab="Modeled Y", main="Observed vs modeled behavior\nof a generic response variable Y")
  abline(a=0, b=1)
  #now surface the locfit package estimates..
  points(y, fitted(zz), col="green", pch=0) #Locfit package - green squares
  legend("topleft", legend=c("Locfit Package", "First Principles"), col=c("green","blue"), pch=c(0,20), bty="n")
    
  ####
  # and compare the L matrices from my fn and locfit package....
  L1<-matrix(0,nrow=N, ncol=N) # create an NxN array
  for(i in 1:N){L1[i,]<-locfit(y~x,alpha=alpha,deg=porder,ev=x[i], kern="bisq", geth=1)}
  print(round(L,digits=4))
  print(round(L1, digits=4))
  identical(round(L, digits=3), round(L1, digits=3))
}
##############################################################
# Now impliment my local polynomial function...
x<-1:10
y<-c(-1.4,2.45,1.45,5.38,5.6,5.99,8.1,7.54,8.24,10.8)
myLocfit(x,y, alpha=0.5, porder=1)

dev.copy(png, "Model Estimates from First Principles and Locfit Package.png")
dev.off()

#### Results
#### Comparing Y observed, Yhat from myLocfit() and Yhat from locfit().
#### Yhat estimated from myLocfit is identical to Yhat estimated from locfit package. 
# Yobs   my_Yest locfit_Yest
# 1     1  1.124187    1.124187
# 2     2  2.631286    2.577466
# 3     3  2.325444    2.365975
# 4     4  4.384105    4.411680
# 5     5  4.843559    4.836667
# 6     6  5.616142    5.606959
# 7     7  8.087243    8.077347
# 8     8  7.642312    7.646826
# 9     9  8.191586    8.181204
# 10   10 10.002716   10.002716

#### Comparing L matrices
#### L matrix estimated from myLocfit() is identical to Hat matrix estimated from linear regression.

# [,1]    [,2]   [,3]    [,4]   [,5]    [,6]   [,7]    [,8]   [,9]  [,10]
# [1,] 0.9705 -0.0691 0.1023 -0.0038 0.0000  0.0000 0.0000  0.0000 0.0000 0.0000
# [2,] 0.0000  0.5245 0.3566  0.0821 0.0369  0.0000 0.0000  0.0000 0.0000 0.0000
# [3,] 0.1433  0.3921 0.4605  0.0041 0.0000  0.0000 0.0000  0.0000 0.0000 0.0000
# [4,] 0.0000  0.0000 0.0000  0.5210 0.3964  0.1714 0.0000 -0.0888 0.0000 0.0000
# [5,] 0.0000  0.0000 0.0000  0.3935 0.3541  0.2600 0.0000 -0.0076 0.0000 0.0000
# [6,] 0.0000  0.0000 0.0000  0.2052 0.2716  0.3741 0.0000  0.1491 0.0000 0.0000
# [7,] 0.0000  0.0000 0.0000  0.0000 0.0000 -0.0228 0.3933  0.1946 0.4349 0.0000
# [8,] 0.0000  0.0000 0.0000  0.0000 0.0000  0.1464 0.2161  0.4861 0.1513 0.0000
# [9,] 0.0000  0.0000 0.0000  0.0000 0.0000 -0.0486 0.4433  0.0674 0.5378 0.0000
# [10,] 0.0000  0.0000 0.0000  0.0000 0.0000  0.0000 0.0164 -0.0476 0.0433 0.9879
# 
# [,1]    [,2]   [,3]    [,4]   [,5]    [,6]   [,7]    [,8]   [,9]  [,10]
# [1,] 0.9705 -0.0691 0.1023 -0.0038 0.0000  0.0000 0.0000  0.0000 0.0000 0.0000
# [2,] 0.0000  0.5245 0.3566  0.0821 0.0369  0.0000 0.0000  0.0000 0.0000 0.0000
# [3,] 0.1433  0.3921 0.4605  0.0041 0.0000  0.0000 0.0000  0.0000 0.0000 0.0000
# [4,] 0.0000  0.0000 0.0000  0.5210 0.3964  0.1714 0.0000 -0.0888 0.0000 0.0000
# [5,] 0.0000  0.0000 0.0000  0.3935 0.3541  0.2600 0.0000 -0.0076 0.0000 0.0000
# [6,] 0.0000  0.0000 0.0000  0.2052 0.2716  0.3741 0.0000  0.1491 0.0000 0.0000
# [7,] 0.0000  0.0000 0.0000  0.0000 0.0000 -0.0228 0.3933  0.1946 0.4349 0.0000
# [8,] 0.0000  0.0000 0.0000  0.0000 0.0000  0.1464 0.2161  0.4861 0.1513 0.0000
# [9,] 0.0000  0.0000 0.0000  0.0000 0.0000 -0.0486 0.4433  0.0674 0.5378 0.0000
# [10,] 0.0000  0.0000 0.0000  0.0000 0.0000  0.0000 0.0164 -0.0476 0.0433 0.9879
#############


#############
# # can also check the GCV if selecting the bestalpha and bestdeg...
# # In the example above, we assumed alpha=0.5 and deg=1
# # use trace(L), trace(L^T L) to get GCV compare with gcvplot
# # compute the GCV for this alpha...
# gcvalpha=(N*sum((y-yhat)^2)) / ((N-sum(diag(L)))^2)
# 
# # compute gcv from the gcvplot command
# zz=gcvplot(y ~ x, alpha= alpha, deg=porder, kern="bisq", ev=dat(), scale=TRUE)
# list(alpha=alpha, gcvmanual=gcvalpha, gcvplot=zz$values)
###############

# (iv)
# If we set alpha=1 (fit using all the data) and weights=1 (linear global fit) then we should get the same results as linear regression...
#define model parameters alpha and polynomial order....

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
  
  # and compare the L matrices from my fn and the hat matrix from linear regression...
  H1<-matrix(0,nrow=N, ncol=N)        #create an Nx1 matrix
  #H1<-hatvalues(zz)
  X<-cbind(rep(1,N),x)                 #augmented x matrix
  H1<-X %*% solve(t(X) %*% X) %*% t(X) #hat matrix
  print(round(L, digits=4))            #L matrix from modified locfit
  print(round(H1, digits=4))          #hat matrix
  identical(round(L, digits=3), round(H1, digits=3)) #compare
  
  # Plot
  par(mfrow=c(1,1))
  plot(y, yhat, pch=20, col="blue", xlab="Observed Y",ylab="Modeled Y", main="Observed vs modeled behavior\nof a generic response variable Y")
  abline(a=0, b=1)
  #now surface the lsfit estimates..
  points(y, zz$fitted.values, col="green", pch=0) #Lsfit package - green squares
  legend("topleft", legend=c("Lsfit", "First Principles"), col=c("green","blue"), pch=c(0,20), bty="n")
  
}

myLsfit(x,y,alpha=1,porder=1)
dev.copy(png, "Model Estimates from First Principles and Lsfit Package.png")
dev.off()

#### Results
#### Comparing Y observed, Yhat from myLsfit() and Yhat from lsfit().
#### Yhat estimated from myLsfit (modified myLocfit with alpha=1 and weights=diag(1)) is identical to Yhat estimated from linear regression. 
# Yobs   my_Yhat lsfit_Yhat
# 1     1 0.1250016  0.1250016
# 2     2 3.1615011  3.1615011
# 3     3 2.3727999  2.3727999
# 4     4 5.4723955  5.4723955
# 5     5 5.6459097  5.6459097
# 6     6 5.9535032  5.9535032
# 7     7 7.6176626  7.6176626
# 8     8 7.1759900  7.1759900
# 9     9 7.7280808  7.7280808
# 10   10 9.7471557  9.7471557

#### Comparing L matrices
#### L matrix estimated from myLsfit (modified myLocfit) is identical to Hat matrix estimated from linear regression.
# [,1]    [,2]    [,3]   [,4]   [,5]   [,6]    [,7]    [,8]    [,9]   [,10]
# [1,]  0.4869  0.2684  0.3251 0.1020 0.0895 0.0674 -0.0525 -0.0207 -0.0604 -0.2058
# [2,]  0.2684  0.1732  0.1979 0.1009 0.0954 0.0858  0.0337  0.0475  0.0302 -0.0330
# [3,]  0.3251  0.1979  0.2310 0.1012 0.0939 0.0810  0.0113  0.0298  0.0067 -0.0779
# [4,]  0.1020  0.1009  0.1012 0.1000 0.0999 0.0998  0.0992  0.0994  0.0992  0.0984
# [5,]  0.0895  0.0954  0.0939 0.0999 0.1003 0.1009  0.1041  0.1033  0.1044  0.1083
# [6,]  0.0674  0.0858  0.0810 0.0998 0.1009 0.1028  0.1129  0.1102  0.1135  0.1258
# [7,] -0.0525  0.0337  0.0113 0.0992 0.1041 0.1129  0.1601  0.1475  0.1632  0.2205
# [8,] -0.0207  0.0475  0.0298 0.0994 0.1033 0.1102  0.1475  0.1376  0.1500  0.1953
# [9,] -0.0604  0.0302  0.0067 0.0992 0.1044 0.1135  0.1632  0.1500  0.1665  0.2267
# [10,] -0.2058 -0.0330 -0.0779 0.0984 0.1083 0.1258  0.2205  0.1953  0.2267  0.3416

# [,1]    [,2]    [,3]   [,4]   [,5]   [,6]    [,7]    [,8]    [,9]   [,10]
# [1,]  0.4869  0.2684  0.3251 0.1020 0.0895 0.0674 -0.0525 -0.0207 -0.0604 -0.2058
# [2,]  0.2684  0.1732  0.1979 0.1009 0.0954 0.0858  0.0337  0.0475  0.0302 -0.0330
# [3,]  0.3251  0.1979  0.2310 0.1012 0.0939 0.0810  0.0113  0.0298  0.0067 -0.0779
# [4,]  0.1020  0.1009  0.1012 0.1000 0.0999 0.0998  0.0992  0.0994  0.0992  0.0984
# [5,]  0.0895  0.0954  0.0939 0.0999 0.1003 0.1009  0.1041  0.1033  0.1044  0.1083
# [6,]  0.0674  0.0858  0.0810 0.0998 0.1009 0.1028  0.1129  0.1102  0.1135  0.1258
# [7,] -0.0525  0.0337  0.0113 0.0992 0.1041 0.1129  0.1601  0.1475  0.1632  0.2205
# [8,] -0.0207  0.0475  0.0298 0.0994 0.1033 0.1102  0.1475  0.1376  0.1500  0.1953
# [9,] -0.0604  0.0302  0.0067 0.0992 0.1044 0.1135  0.1632  0.1500  0.1665  0.2267
# [10,] -0.2058 -0.0330 -0.0779 0.0984 0.1083 0.1258  0.2205  0.1953  0.2267  0.3416
###############

# (v) Replace the last value of y with an outlier.  Set alpha=0.9 and recompute the L matrix from local polynomial regression and the Hat matrix from linear regression.  Compare the weights and comment
# Now impliment my local polynomial function...
x<-1:10
ymod<-c(-1.4,2.45,1.45,5.38,5.6,5.99,8.1,7.54,8.24,10.8*10)
N<-length(ymod)

myLocfit(y=ymod, x=x, alpha=0.9, porder=1) # Prints L matrix from myLocfit and locfit package as well as yhat estimates from myLocfit and locfit package...

# Now compare with hat matrix from linear regression...
X<-cbind(rep(1,N),x) # the augmented X matrix
H1<-X %*% solve(t(X) %*% X) %*% t(X) # the Hat matrix
round(H1, digits=4)

#### Results
#### L matrix from myLocfit()
# [,1]   [,2]    [,3]    [,4]    [,5]    [,6]    [,7]    [,8]   [,9]  [,10]
# [1,] 0.5277 0.3726  0.2121  0.0725 -0.0254 -0.0699 -0.0634 -0.0261 0.0000 0.0000
# [2,] 0.3651 0.3012  0.2130  0.1204  0.0423 -0.0071 -0.0228 -0.0123 0.0000 0.0000
# [3,] 0.2155 0.2258  0.2051  0.1619  0.1086  0.0583  0.0215  0.0034 0.0000 0.0000
# [4,] 0.0862 0.1444  0.1831  0.1928  0.1722  0.1277  0.0717  0.0219 0.0000 0.0000
# [5,] 0.0000 0.0449  0.1319  0.2060  0.2344  0.2060  0.1319  0.0449 0.0000 0.0000
# [6,] 0.0000 0.0000  0.0449  0.1319  0.2060  0.2344  0.2060  0.1319 0.0449 0.0000
# [7,] 0.0000 0.0000  0.0219  0.0717  0.1277  0.1722  0.1928  0.1831 0.1444 0.0862
# [8,] 0.0000 0.0000  0.0034  0.0215  0.0583  0.1086  0.1619  0.2051 0.2258 0.2155
# [9,] 0.0000 0.0000 -0.0123 -0.0228 -0.0071  0.0423  0.1204  0.2130 0.3012 0.3651
# [10,] 0.0000 0.0000 -0.0261 -0.0634 -0.0699 -0.0254  0.0725  0.2121 0.3726 0.5277

##### Hat matrix from linear regression
# [,1]    [,2]    [,3]   [,4]   [,5]   [,6]   [,7]    [,8]    [,9]   [,10]
# [1,]  0.3455  0.2909  0.2364 0.1818 0.1273 0.0727 0.0182 -0.0364 -0.0909 -0.1455
# [2,]  0.2909  0.2485  0.2061 0.1636 0.1212 0.0788 0.0364 -0.0061 -0.0485 -0.0909
# [3,]  0.2364  0.2061  0.1758 0.1455 0.1152 0.0848 0.0545  0.0242 -0.0061 -0.0364
# [4,]  0.1818  0.1636  0.1455 0.1273 0.1091 0.0909 0.0727  0.0545  0.0364  0.0182
# [5,]  0.1273  0.1212  0.1152 0.1091 0.1030 0.0970 0.0909  0.0848  0.0788  0.0727
# [6,]  0.0727  0.0788  0.0848 0.0909 0.0970 0.1030 0.1091  0.1152  0.1212  0.1273
# [7,]  0.0182  0.0364  0.0545 0.0727 0.0909 0.1091 0.1273  0.1455  0.1636  0.1818
# [8,] -0.0364 -0.0061  0.0242 0.0545 0.0848 0.1152 0.1455  0.1758  0.2061  0.2364
# [9,] -0.0909 -0.0485 -0.0061 0.0364 0.0788 0.1212 0.1636  0.2061  0.2485  0.2909
# [10,] -0.1455 -0.0909 -0.0364 0.0182 0.0727 0.1273 0.1818  0.2364  0.2909  0.3455
#### Comments 
# The L matrix from myLocfit gives weight to only 8 of 10 observations (alpha=0.9) compared to all the observations in linear regression. # This alone helps reduce the influence of outliers
# In addition, the non-zero weights farthest from the fitted point xp are smaller in the Lmatrix than in the Hat matrix.
# Finally, if we loook at the diagonal of the L matrix, we see that y(xp) (e.g. L[i,i]) is given the most weight in estimating the function at xp, which is what you want from a local polynomial fit.  In contrast, the value of y(xp) is often not given the most weight in estimating the function at xp in linear regression, thus making the fit more susecptible to outliers.


################# End Q3
################# 
################# Start Q4
# 4. Repeat 1 with Local polynomial method.
setwd("/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013/Final")

#Import data: Mean annual precipitation at 491 locations in Colorado 1980-2002
data<-read.table(file="/Users/elliotcohen/Dropbox/Data/Climate/Rainfall/Colorado_Annual_Precip_Gridded.csv", sep=",", header=TRUE)

# i. Fit a ‘best’ local regression model
################ subset selection ##############
library(leaps)    # to provide combinations
library(MPV)  	# to help estimate PRESS and consequently, GCV
library(locfit)
#Select independent and dependent variables
X<-data[,1:3] # all the predictor set.
Y<-data[,4] # response variable Y
N = length(Y)

combs = leaps(X,Y, nbest=25)  #Get upto 25 combinations of predictor and predictants
combos = combs$which  # logical combinations of predictor variables
ncombos = length(combos[,1]) # number of combinations of predictors variables

for(i in 1:ncombos){
  xx<-as.data.frame(X[,combos[i,]])
  names(xx)<-names(X)[combos[i,]]
  bestparam(X=xx, Y=Y, family="gaussian") # find bestalpha and bestdeg for given set of predictor variables
  
  # apply bestalpha and bestdeg to fit the local polynomial
  xx<-as.matrix(xx)
  y<-as.matrix(Y)
  zz<-locfit(y~xx, alpha=bestalpha, deg=bestdeg) 
  # create vector of GCV values for each model with its own bestalpha and bestdeg
  if(i==1) {GCVs <- gcv(zz)[[4]]} else {GCVs <- rbind(GCVs,gcv(zz)[[4]])} 
}

# select the model with the overall lowest GCV and re-fit the "bestmodel"
besti<-which.min(GCVs)            #best combination of predictor variables based on GCV
names<-names(X)[combos[besti,]]     #predictors
X<-as.data.frame(X[,combos[besti,]]) #capital X = df; lowercase x = matrix
names(X)<-names
x<-as.matrix(X)
bestparam(X=X, Y=Y,family="gaussian") # alpha=0.1622, deg=1
bestmodel<-locfit(y~x, alpha=bestalpha, deg=bestdeg) #fit the best model

# Compute useful quantities for locfit...
modresid=residuals(bestmodel) # Get the residuals of the model..
nX<-dim(X)[2]
Yhat=Y-modresid  #residuals = Y - Yestimate ==> Yestimate = Y - residuals
k=dim(X)[2]      #number of regressor variables
p=k+1             #number of model parameters 
n=length(Y)       #number of observations

##### Custom Model Diagnostics (6 visual checks) #####
GLMdiagnostics(bestmodel, X, Y)
dev.copy(png,paste(bestmodel$call[1],"Custom_Model_Diagnostics", sep="_"))
dev.off()

# Comments on Model Diagnostics: 
# Q-Q plot reveals non-normality in the upper tail --> FAIL diagnostic.
# Plot of residuals vs predictor variable X1 reveals heteroscdasticity --> FAIL diagnostic.
# Plot of residuals vs model estimates of Y reveals heteroscdasticity --> FAIL diagnostic.
# Autocorrelation does not appear to be a big problem --> PASS diagnostic
# --> locfit with gaussian distribution is *not* adequate for modeling precip data.
#### END Custom MODEL DIAGNOSTICS ####

#### Model Signficance ####
# Ftest comparing locfit to lm
RSS1 = sum(residuals(bestmodel)^2)
nu1 = sum(fitted(bestmodel,what="infl"))     # trace(L)
nu2 = sum(fitted(bestmodel,what="vari"))     # trace(L^T L)
## Or
#nu1 = bestmodel$dp[6]
#nu2 = bestmodel$dp[7]
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
Ftheory = qf(0.95,(nu00-nu11), nu11)    # 95% confidence level..
Fdata>Ftheory                           # TRUE
## Ftest results:
## Fdata > Ftheory -- reject null (e.g. reject that the locfit is no different from a linear model)

# iii. Compute cross-validated and fitted estimates at each observation points and plot them against the observed values. This is to visually see how the model performs in a fitting and cross-validated mode.

# Impliment cross-validation function....
xpred<-X          #xpred<-data.frame(data$Long)
crossval(X=X, Y=Y, xpred=xpred)
dev.copy(png,paste(bestmodel$call[1],"Cross-Validation", sep="_"))
dev.off()

# Comments:  Cross-validated model estimates and estimates modeled on the full data set appear internally consistent to each other, but niether are accurate to observed precipitation values for larger values of precip (e.g. above 1000 mm/yr).


# iv. Drop 10% of observations, fit the model (i.e., the ‘best’ model from i. above) to the rest of the data and predict the dropped points. Compute RMSE and R2 and show them as boxplots. 
droptest(X=X, Y=Y, bestmodel=bestmodel, drop=0.1)
dev.copy(png,paste(bestmodel$call[1],"RMSE_Skill", sep="_"))
dev.off()

#v. Make a 3-D plot (or spatial colored/contour map) of model estimates (i.e. latitude, longitude and the model estimates) and model error. 
# (v) Spatial map of actual precipitation, model estimates, and model error
png(filename=paste(bestmodel$call[1],"ContourPlots.png" ,sep="_"), width=800, height=600, units="px", pointsize=18)
par(mfrow=c(2,2))
par(mar=c(4,4,4,6)+0.1) 

#par(mar=c(4,4,3,3)+0.1) # c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. The default is c(5, 4, 4, 2) + 0.1.

zz0<-interp(x=data$Long, y=data$Lat, z=data$Precip, duplicate="mean")
image.plot(zz0, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Observed Precipitation (mm)")
contour(zz0, add=T)
world(add=TRUE, lwd=4)

zz1<-interp(data$Long, data$Lat, fitted(bestmodel), duplicate="mean")
image.plot(zz1, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Local Polynomial (Gaussian)\nFitted Precipitation (mm)")
contour(zz1, add=T)
world(add=TRUE, lwd=4)

zz2<-interp(data$Long, data$Lat, modresid, duplicate="mean")
image.plot(zz2, col=topo.colors(20,0.5), ylab="latitude", xlab="longitude", main="Local Polynomial (Gaussian)\nResiduals (mm/yr)")
contour(zz2, add=T)
world(add=TRUE, lwd=4)

# vi. Estimate the precipitation and the standard errors on the DEM grid and show them as above.
# get DEM grid Lat-Long-Elev
#predpoints = read.table("http://cires.colorado.edu/~aslater/CVEN_6833/colo_dem.dat",sep=",")
predpoints<-read.table(file="/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013/colo_dem.txt")
names(predpoints)<-c("Lat","Long","Elev")
Xpred<-predpoints$Long #gridded lat, long and elev
ypred = predict(bestmodel,newdata=Xpred, type="response",se.fit=TRUE)

zz = interp(predpoints$Long, predpoints$Lat, ypred$fit)
surface(zz,xlab="Longitude",ylab ="Latitude", main="DEM Grid Predictions") 

#dev.copy(png,paste(bestmodel$call[1],"ContourPlots", sep="_"))
dev.off()

# Discussion:
# The local polynomial regression is better than the linear model and the generalized linear model with appropriate link function, but still not very good as can be seen in the residuals and from the model diagnostics.  We see that most of the variance is still explained by the model residuals rather than the local polynomial model, indicating that the model is insufficienct.

# Q4 complete!
# Let's repeat one more time, but now with an appropriate link function, e.g. local glm.
################# End Q4
################# 
################# Start Q5
# 5. Repeat 4 with Local Polynomial method but using the appropriate link function (i.e. ‘Local GLM’). [For the Local Polynomial approach the ‘best model’ involves fitting the best subset of predictors and the smoothing parameter, alpha. You can also compare the GCV from these four different methods.]

setwd("/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013/Final")
data<-read.table(file="/Users/elliotcohen/Dropbox/Data/Climate/Rainfall/Colorado_Annual_Precip_Gridded.csv", sep=",", header=TRUE)
X<-data[,1:3] # all the predictor set.
Y<-data[,4] # response variable Y
N = length(Y)

# choose an appropriate link function....
library(fitdistrplus)
par(mfrow=c(1,1))
plotdist(Y) # interpretation: use a gamma distribution
dev.copy(png,paste(bestmodel$call[1],"Distribution_Heuristic_1.png", sep="_"))
dev.off()

descdist(Y) # interpretation: use a gamma distribution
dev.copy(png,paste(bestmodel$call[1],"Distribution_Heuristic_2.png", sep="_"))
dev.off()

# fig1 <- fitdist(Y, "gamma")
# plot(fig1)
# summary(fig1)  
# dev.copy(png,paste(bestmodel$call[1],"Gamma_Fit", sep="_"))
# dev.off()

# i. Fit a ‘best’ local regression model
################ subset selection ##############
library(leaps)    # to provide combinations
library(MPV)  	# to help estimate PRESS and consequently, GCV
library(locfit)
#Select independent and dependent variables
X<-data[,1:3] # all the predictor set.
Y<-data[,4] # response variable Y
N = length(Y)

combs = leaps(X,Y, nbest=25)  #Get upto 25 combinations of predictor and predictants
combos = combs$which  # logical combinations of predictor variables
ncombos = length(combos[,1]) # number of combinations of predictors variables

for(i in 1:ncombos){
  xx<-as.data.frame(X[,combos[i,]])
  names(xx)<-names(X)[combos[i,]]
  bestparam(X=xx, Y=Y, family="gamma") # find bestalpha and bestdeg for given set of predictor variables
  
  # apply bestalpha and bestdeg to fit the local polynomial
  xx<-as.matrix(xx)
  y<-as.matrix(y)
  zz<-locfit(y~xx, alpha=bestalpha, deg=bestdeg, family="gamma") 
  # create vector of GCV values for each model with its own bestalpha and bestdeg
  if(i==1) {GCVs <- gcv(zz)[[4]]} else {GCVs <- rbind(GCVs,gcv(zz)[[4]])} 
}

# select the model with the overall lowest GCV and re-fit the "bestmodel"
besti<-which.min(GCVs)            #best combination of predictor variables based on GCV
names<-names(X)[combos[besti,]]     #predictors
X<-as.data.frame(X[,combos[besti,]]) #capital X = df; lowercase x = matrix
names(X)<-names
x<-as.matrix(X)
bestparam(X=X, Y=Y, family="gamma") # alpha=0.1622, deg=1
bestmodel<-locfit(y~x, alpha=bestalpha, deg=bestdeg, family="gamma") #fit the best model

# Compute useful quantities for locfit...
modresid=residuals(bestmodel) # Get the residuals of the model..
nX<-dim(X)[2]
Yhat=Y-modresid  #residuals = Y - Yestimate ==> Yestimate = Y - residuals
k=dim(X)[2]      #number of regressor variables
p=k+1             #number of model parameters 
n=length(Y)       #number of observations

##### Custom Model Diagnostics (6 visual checks) #####
GLMdiagnostics(bestmodel, X, Y)
dev.copy(png,paste(bestmodel$call[1],"Custom_Model_Diagnostics", sep="_"))
dev.off()

# Comments on Model Diagnostics: 
# Q-Q plot looks good! --> PASS diagnostic.
# Plot of residuals vs predictor variable X1 reveals heteroscdasticity --> FAIL diagnostic.
# Plot of residuals vs model estimates of Y reveals heteroscdasticity --> FAIL diagnostic.
# Autocorrelation does not appear to be a big problem --> PASS diagnostic
# --> locfit with gamma distribution performs better than locfit with gaussian distributions (residuals are normal), but is still *not* adequate for modeling precip.

#### END Custom MODEL DIAGNOSTICS ####

#### Model Signficance ####
# Ftest comparing locfit to lm
RSS1 = sum(residuals(bestmodel)^2)
nu1 = sum(fitted(bestmodel,what="infl"))     # trace(L)
nu2 = sum(fitted(bestmodel,what="vari"))     # trace(L^T L)
## Or
#nu1 = bestmodel$dp[6]
#nu2 = bestmodel$dp[7]
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
Ftheory = qf(0.95,(nu00-nu11), nu11)    # 95% confidence level..
Fdata>Ftheory                           # TRUE
## Ftest results:
## Fdata > Ftheory -- reject null (e.g. reject that the locfit is no different from a linear model)

# iii. Compute cross-validated and fitted estimates at each observation points and plot them against the observed values. This is to visually see how the model performs in a fitting and cross-validated mode.

# Impliment cross-validation function....
xpred<-X          #xpred<-data.frame(data$Long)
crossval(X=X, Y=Y, xpred=xpred)

dev.copy(png,paste(bestmodel$call[1],"Cross-Validation", sep="_"))
dev.off()

# Comments:  Cross-validated model estimates and estimates modeled on the full data set appear internally consistent to each other, but niether are accurate to observed precipitation values for larger values of precip (e.g. above 1000 mm/yr).


# iv. Drop 10% of observations, fit the model (i.e., the ‘best’ model from i. above) to the rest of the data and predict the dropped points. Compute RMSE and R2 and show them as boxplots. 
droptest(X=X, Y=Y, bestmodel=bestmodel, drop=0.1)

dev.copy(png,paste(bestmodel$call[1],"RMSE_Skill", sep="_"))
dev.off()

#v. Make a 3-D plot (or spatial colored/contour map) of model estimates (i.e. latitude, longitude and the model estimates) and model error. 
# (v) Spatial map of actual precipitation, model estimates, and model error
png(filename=paste(bestmodel$call[1],"ContourPlots.png" ,sep="_"), width=800, height=600, units="px", pointsize=18)
par(mfrow=c(2,2))
par(mar=c(4,4,4,6)+0.1) 
library(akima)
library(fields)

zz0<-interp(x=data$Long, y=data$Lat, z=data$Precip, duplicate="mean")
image.plot(zz0, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Observed Precipitation (mm)")
contour(zz0, add=T)
world(add=TRUE, lwd=4)

zz1<-interp(data$Long, data$Lat, fitted(bestmodel), duplicate="mean")
image.plot(zz1, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Local Polynomial (Gamma)\nFitted Precipitation (mm)")
contour(zz1, add=T)
world(add=TRUE, lwd=4)

zz2<-interp(data$Long, data$Lat, modresid, duplicate="mean")
image.plot(zz2, col=topo.colors(20,0.5), ylab="latitude", xlab="longitude", main="Local Polynomial (Gamma)\nResiduals (mm/yr)")
contour(zz2, add=T)
world(add=TRUE, lwd=4)

# vi. Estimate the precipitation and the standard errors on the DEM grid and show them as above.
# get DEM grid Lat-Long-Elev
#predpoints = read.table("http://cires.colorado.edu/~aslater/CVEN_6833/colo_dem.dat",sep=",")
predpoints<-read.table(file="/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013/colo_dem.txt")
names(predpoints)<-c("Lat","Long","Elev")
Xpred<-predpoints$Long #gridded lat, long and elev
ypred = predict(bestmodel,newdata=Xpred, type="response",se.fit=TRUE)

zz = interp(predpoints$Long, predpoints$Lat, ypred$fit)
surface(zz,xlab="Longitude",ylab ="Latitude", main="DEM Grid Predictions") 

#dev.copy(png,paste(bestmodel$call[1],"Contour_Plots", sep="_"))
dev.off()

# Discussion:
# The local polynomial regression with appropriate link function ("gamma) is better than the linear model, the generalized linear model with appropriate link function, and the local polynomial model with gaussian link function, but still not very good as can be seen in the residuals and from the model diagnostics.  We see that most of the variance is still explained by the model residuals rather than the local polynomial model, indicating that the model is insufficienct.

# Q5 complete!
################# End Q5
################# 
################# Start Q6

setwd("/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013/Final")

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

X<-data[,-4]  # df of Independent variable
Y<-data$Precip #vector of dependent variable (annual rainfall at lat-long position)
n<-length(Y)
nvar<-dim(X)[2]

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


######## Model 2  - Lecture 3
######## Spatial trend with a linear model + Kriging residuals #####
######## Estimated simultaneously

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
##############

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

png(filename="Observed vs Modeled Precip - 3 Kriging Models.png", width=800, height=600, units="px", pointsize=16)
par(mfrow=c(2,2))
par(mar=c(4,3,3,2)+ 0.1) 
#par(mar=c(4, 3, 3, 2) + 0.1) # numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot
zz0<-interp(x=data$Long, y=data$Lat, z=data$Precip, duplicate="mean")
image.plot(zz0, col=topo.colors(n=20,alpha=0.5), ylab="latitude", xlab="longitude", main="Observed Precipitation (mm)")
contour(zz0, add=T)
world(add=TRUE, lwd=4)

image.plot(xlon,ylat,t(zmat1a), col=topo.colors(20,0.5), ylab="latitude", xlab="longitude", main="Model 1 precip predictions (mm/yr)\nKriging")
contour(xlon,ylat,t(zmat1a),add=T)

image.plot(xlon,ylat,t(zmat2a), col=topo.colors(20,0.5), ylab="latitude", xlab="longitude", main="Model 2 precip predictions (mm/yr)\nAdditive: LM + Kriging")
contour(xlon,ylat,t(zmat2a),add=T)

image.plot(xlon,ylat,t(zmat3a), col=topo.colors(20,0.5), ylab="latitude", xlab="longitude", main="Model 3 precip predictions (mm/yr)\nHierarchical: GLM + Kriging")
contour(xlon,ylat,t(zmat3a),add=T)

#dev.copy(png,"Observed vs. Modeled Precip - 3 models.png")
dev.off()
######
###### Conclusions from Q6 and Q7...
# As expected, Kriging is effective at caputuring localized effects, and provides  much more fidelty to actual observed precipitation patterns than do linear models or even local polynomial methods (from HW questions 1, 2, 4 and 5.)
# By visual inspeciton alone, it is difficult to discern between the three spatial models... 
# Refer to the cross-validated estimates (above) for a quantitative assessment

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

dev.copy(png,"Observed vs. Modeled Precip and Residuals - Binomial Regression.png")
dev.off()

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
################# End Q10
#################
################# Skip Q11 as per instructions
################# Skip Q12 as per instructions
#################
################# Start Q13
# In this homework we explored:

# linear models: The simplest, most parsimonious regression model.  Ideal for linear datasets with a normal distribution.

# generalized linear models: Also a global, linear model akin to line regression but with the flexibility of cononical link functions.  GLM accepts any theoretical distribution, such as gaussian, gamma, poison, binomial, etc...

# local polynomial models: A non-parametric regression model that does *not* require apriori selection of a theoretical distribution.  Local polynomial methods are data-driven akin to histograms with appropriate smoothing.  Local polynomial methods use k-nearest-neighbors to estimate the function at any point xp, thus dispensig with the need for apriori bandwidth selection while allowing for local features to be captured.  

# Kriging spatial models: Kriging is based on the assumption that covariance in spatial data is a function of the distance between points. Kriging can be considered a local estimation method (e.g. nearby observations are weighted more than observations far away).  Kriging captues local features in spatial data.

# Additive spatial models and hierarchical spatial models:  Kriging can be coupled with linear regression or local polynomial methods to get "the best of both worlds".  A regression model can be fit to the observed data to capture the mean trend (e.g. increased precipitation from East to West in Colorado) and then take the residuals and apply kriging to capture non-linear spatial features.

## END HW1 !!