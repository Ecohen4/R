##### Simulated RMSE Skill #####
## Drop some % of points, fit the model and predict the dropped points..
droptest<-function(X, Y, drop){
  source("/Users/elliotcohen/Dropbox/Advance Data Analysis/R source files/myboxplot.r")
  source("/Users/elliotcohen/Dropbox/Advance Data Analysis/R source files/myboxplot-stats.r")
  
  nsim = 500
  yhatcross=1:nsim
  ksecross=1:nsim  
  N = length(Y)
  N10 = round(drop*N)    #choose % of points to drop (e.g. 10%)
  index=1:N
  ########## special code for kriging....
  	nobs<-length(Y)
	nobs1<-nobs-1
	yhatcross<-1:nobs
	ksecross<-1:nobs
  
  ########### end kriging section
  
  for(i in 1:nsim){
    drop=sample(c(1:N),N10)  #sample 10% of the row indices from 1:N at random
    keep=setdiff(index,drop)  #discard values at the intersection of index and drop (e.g. drop 10% of the data, keep the rest)
    xkeep<-X[keep,]
    ykeep<-Y[keep]
    newdf<-as.data.frame(cbind(ykeep,xkeep))
    nobs<-length(ykeep)  #added this to make conformable args in computing ghat...
    nobs1<-nobs-1  #added this to make conformable args in computing ghat...
    ######### special code for kriging....
    sigma=1
    rho<-pars[1]
    theta<-pars[2]
    xp<-X[drop,] #predict at the dropped points


    # # fit of exponential by nonlinear least squares
	 # xd<- look$centers
	 # ym<- look$stats[2,]
	 # sigma<- 1.0  
	 # nls.fit<- nls( ym ~ sigma^2 + rho*( 1- exp(-xd/theta)),
	               # start= list(rho=max(look$stats[2,]), theta=quantile(xd,0.1)), control=list( maxiter=200) )
	 # pars<- coefficients(nls.fit)
	 # rho<-pars[1]
	 # theta<-pars[2]
	 # xr = round(max(xd)+sd(xd))
	 # dgrid<- seq(0,xr,,400)
	 # lines(dgrid, sigma^2 + rho*(1 - exp(-1*dgrid/theta)), col="blue", lwd=3)


	K <- sigma^2 + rho*(1 - exp(-1*rdist(xkeep,xkeep)/theta)) #dim(K)= 442 x 442
	Kstar <- sigma^2 + rho*(1 - exp(-1*rdist(as.matrix(t(xp)),xkeep)/theta)) # 3 x 442
	ghat<- Kstar %*% solve( K + diag( sigma^2, nobs)) %*% ykeep
	kse <- sqrt(var(ykeep) - diag((Kstar) %*% solve( K + diag( sigma^2, nobs)) %*% t(Kstar)))
	yhatcross[i] <- ghat
	ksecross[i] <- kse
    }
    # works up to here....
    
    
    #### or use....
	zz<- Krig(xkeep,ykeep,rho=rho,theta=theta,sigma2=sigma,m=1) #fit the model the the remaining data... 
	yest[i] = predict.Krig(zz,x=as.matrix(t(xp)),drop.Z=TRUE) #predict at the dropped points
	### error... dim(t(xp)) does not match the data used to fit the model zz...
	## try this for the standard error
	predict.se.Krig(zz,x=as.matrix(t(xp)),drop.Z=TRUE)
    
  
    ############# end special kriging section 
  
  par(mfrow=c(1,2))
  boxplot(yhatcross, main="Cross-Validated Y Estimates") #simple version
  boxplot(ksecross, main="Cross-Validated Standard Errors" )  #simple version
  
  } #### close function ####
  
  # zz=myboxplot(rmseskill, main="Simulated RMSE skill",plot=FALSE)
  # zz$names=rep("",length(zz$names))
  # z1=bxp(zz,xlab="",ylab="RMSE",cex=1.25)
  
  
  # rmse<-sqrt(sum(modresid^2)/N)
  # # alternate calculations....
  # #rmse2=mean(((modresid)/sd(Y))^2) # from Balaji's code...
  # #rmse3<-sqrt(RSS1/nu11)
  # #SSE = sum((modresid)^2)            
  # #MSE = SSE/nu11 
  # #rmse4<-sqrt(MSE)
  
  # points(z1,rmse,col="red",cex=2,pch=19)  #add a point showing the true RMSE of the data.
  # title(main="RMSE skill")
  
  # zz=myboxplot(corskill, main="Simulated Correlation skill",plot=FALSE)
  # zz$names=rep("",length(zz$names))
  # z1=bxp(zz,xlab="",ylab="Cor",cex=1.25)
  # title(main="Cor skill")
  
} #### close function ####