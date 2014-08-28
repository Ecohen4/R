#### Model Diagnostics (6 visual checks) ####
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



###### Best alpha and best polynomial order #########
# Custom function to search for the best alpha (between 0 and 1) and best polynomial order (1 or 2) for the local polynomial fit.

bestparam<-function(X, family){
  
  N = length(Y)
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
  y<<-as.matrix(Y)
  x<<-as.matrix(X)
  #dimnames(xx)<-list(rep("obs",N),names(X)[combos[besti,]])
  zz<-gcvplot(y ~ x, alpha=alpha1, deg=1, kern='bisq', ev=dat(), scale=TRUE, family=family)
  z1<-gcvplot(y ~ x, alpha=alpha2, deg=2, kern="bisq", ev=dat(), scale=TRUE, family=family)
  
  # pick the best alpha and the degree of the polynomial that gives the least GCV
  z2=order(c(zz$values,z1$values)) #order from lowest GCV to highest
  
  deg1=1
  if(z2[1] > n)deg1=2 #select degree (1 or 2) that yields lowest GCV
  bestalpha<<-alpha[z2[1]] #select alpha that yields lowest GCV and assign to global environment
  bestdeg<<-deg1 #select polynomial degree that yields lowest GCV and assign to global environment
  print(bestalpha)
  print(bestdeg)
} 
######## close function ########


####### Cross-Validated and Fitted Estimates ########
crossval<-function(X, Y, xpred){
  
  Fitted<-predict(bestmodel, newdata=xpred, type="response")
  
  par(mfrow=c(1,1))
  plot(Y, Fitted, pch=20, col="black", xlab="Observed Precip (mm)",ylab="Modeled Precip (mm)", main="Observed vs Modeled Precipitation")
  abline(a=0, b=1)
  
  # cross validated estimates (drop observations one at a time, refit the model to the remaining data (N-1) and predict at the dropped point; repeat for all observations)
  par(mfrow=c(1,1))
  n=length(Y)
  yest=1:n
  nvar=dim(X)[2]
  
  index=1:n
  for(i in 1:n){
    index1=index[index != i] #drop one observation at a time
    Xval=X[index1,] #X data less the dropped observation
    Yval=Y[index1]  #Y data less the dropped observation
    
    zz<-locfit(Yval ~ Xval, alpha=bestalpha, deg=bestdeg) #fit the model without the dropped observation
    
    #xpred=c(1,X[i,1:nvar]) #now estimate at the point that was dropped
    xpred<-X[i,1:nvar] #now estimate at the point that was dropped
    xpred<-as.numeric(xpred)
    yest[i]<-predict(zz, newdata=xpred, type="response")
  }
  #now surface the x-validated estimates..
  points(Y, yest, col="blue", pch=20)
  
} 
####### close function ########



##### Simulated RMSE Skill #####
## Drop some % of points, fit the model and predict the dropped points..
droptest<-function(X, Y, drop){
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
    newdf<-as.data.frame(cbind(yy,xx))
    zz<-locfit(yy~xx, data=newdf, alpha=bestalpha, deg=bestdeg) #fit model to remaining data
    xpred<-X[drop,] #the dropped data
    yhat<-predict(zz,newdata=xpred) #predict at dropped points using model fit w.out those points
    rmseskill[i]<-sqrt(mean(Y[drop]-yhat)^2)
    corskill[i]<-cor(Y[drop],yhat)
  }
  
  par(mfrow=c(1,2))
  boxplot(rmseskill, main="Simulated RMSE") #simple version
  boxplot(corskill, main="Simulated Cor." )  #simple version
  
  zz=myboxplot(rmseskill, main="Simulated RMSE skill",plot=FALSE)
  zz$names=rep("",length(zz$names))
  z1=bxp(zz,xlab="",ylab="RMSE",cex=1.25)
  #rmse=mean(((modresid)/sd(Y))^2) # from Balaji's code... 
  rmse<-sqrt(sum(modresid^2)/N)  #confirmed from Balaji....
  # Use this one....
  
  # Identical calculations....
  #rmse2<-sqrt(RSS1/nu11)
  #SSE = sum((modresid)^2)            
  #MSE = SSE/nu11 
  #rmse3<-sqrt(MSE)
  
  points(z1,rmse,col="red",cex=2,pch=19)  #add a point showing the true RMSE of the data.
  title(main="RMSE skill")
  
  zz=myboxplot(corskill, main="Simulated Correlation skill",plot=FALSE)
  zz$names=rep("",length(zz$names))
  z1=bxp(zz,xlab="",ylab="Cor",cex=1.25)
  title(main="Cor skill")
  
} #### close function ####




################# End Q3
################# 
################# Start Q4
# 4. Repeat 1 with Local polynomial method.
setwd("/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013")

#Import data: Mean annual precipitation at 491 locations in Colorado 1980-2002
data<-read.table(file="/Users/elliotcohen/Dropbox/Data/Climate/Rainfall/Colorado_Annual_Precip_Gridded.csv", sep=",", header=TRUE)

# i. Fit a ‘best’ local regression model
################ subset selection ##############
library(leaps)    # to provide combinations
library(MPV)		# to help estimate PRESS and consequently, GCV
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
  bestparam(X=xx, family="gaussian") # find bestalpha and bestdeg for given set of predictor variables
  
  # apply bestalpha and bestdeg to fit the local polynomial
  xx<-as.matrix(xx)
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
bestparam(X, family="gaussian") # alpha=0.1622, deg=1
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
droptest(X, Y, drop=0.1)
dev.copy(png,paste(bestmodel$call[1],"RMSE_Skill", sep="_"))
dev.off()

#v. Make a 3-D plot (or spatial colored/contour map) of model estimates (i.e. latitude, longitude and the model estimates) and model error. 
# (v) Spatial map of actual precipitation, model estimates, and model error
plot.new()
par(mfrow=c(2,2))
par(mar=c(4,3,4,6)+0.1) # c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. The default is c(5, 4, 4, 2) + 0.1.
library(akima)
library(fields)

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

dev.copy(png,paste(bestmodel$call[1],"ContourPlots", sep="_"))
dev.off()

# Discussion:
# The local polynomial regression is better than the linear model and the generalized linear model with appropriate link function, but still not very good as can be seen in the residuals and from the model diagnostics.  We see that most of the variance is still explained by the model residuals rather than the local polynomial model, indicating that the model is insufficienct.

# Q4 complete!
# Let's repeat one more time, but now with an appropriate link function, e.g. local glm.
################# End Q4
################# 
################# Start Q5
# 5. Repeat 4 with Local Polynomial method but using the appropriate link function (i.e. ‘Local GLM’). [For the Local Polynomial approach the ‘best model’ involves fitting the best subset of predictors and the smoothing parameter, alpha. You can also compare the GCV from these four different methods.]

setwd("/Users/elliotcohen/Dropbox/Advance Data Analysis/HW/HW1-2013")
data<-read.table(file="/Users/elliotcohen/Dropbox/Data/Climate/Rainfall/Colorado_Annual_Precip_Gridded.csv", sep=",", header=TRUE)
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
dev.copy(png,paste(bestmodel$call[1],"Gamma_Fit", sep="_"))
dev.off()

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
  bestparam(X=xx, family="gamma") # find bestalpha and bestdeg for given set of predictor variables
  
  # apply bestalpha and bestdeg to fit the local polynomial
  xx<-as.matrix(xx)
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
bestparam(X, family="gamma") # alpha=0.1622, deg=1
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
droptest(X, Y, drop=0.1)

dev.copy(png,paste(bestmodel$call[1],"RMSE_Skill", sep="_"))
dev.off()

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

dev.copy(png,paste(bestmodel$call[1],"Contour_Plots", sep="_"))
dev.off()

# Discussion:
# The local polynomial regression with appropriate link function ("gamma) is better than the linear model, the generalized linear model with appropriate link function, and the local polynomial model with gaussian link function, but still not very good as can be seen in the residuals and from the model diagnostics.  We see that most of the variance is still explained by the model residuals rather than the local polynomial model, indicating that the model is insufficienct.

# Q5 complete!
################# End Q5
################# 
################# Start Q6