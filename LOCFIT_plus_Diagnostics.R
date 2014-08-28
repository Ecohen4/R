##################
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

################## End Q3
################# 
################# 