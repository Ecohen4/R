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