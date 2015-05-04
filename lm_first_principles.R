elm<-function(X,Y) {
  X <- as.matrix(X)    # predictor set X
  Y <- as.matrix(Y)    # response vector Y
  Xaug <- cbind(1,X)  # augmented matrix
  Yhat <- Xaug %*% solve(t(Xaug) %*% Xaug) %*% t(Xaug) %*% Y # least-squares solution
  residuals <- Y - Yhat
  hatm <- Xaug %*% solve(t(Xaug) %*% Xaug) %*% t(Xaug)
  beta <- solve(t(Xaug) %*% Xaug) %*% t(Xaug) %*% Y
  SST = sum((Y - mean(Y))^2)    # deviation of observed values from sample mean
  SSR = sum((Yhat - mean(Y))^2) # deviation of fitted values from sample mean
  SSE = sum((Yhat - Y)^2)       # deviation of model fit from observed response
  # SSE = sum((bestmodel$residuals)^2)     # same as above
  n<-dim(X)[1]
  k<-dim(X)[2]
  p<-length(beta)
  MSR = SSR / k                            # SSR/k predictors
  MSE = SSE/(n - p)                        # SSE/n-p
  R2 = (SST-SSE)/SST                       # R2
  PearsonR2 = (cor(Yhat,Y))^2              # equivalent to R2
  AdjustedR2 = 1-((SSE/(n-p))/(SST/(n-1))) # adjusted R2
  varExplained = SSR/SST                   # Anu suggestion Jan. 4 2014.
  Ftest=MSR/MSE                            # Ftest for multivariate regresion
  ANOVA<-data.frame(SST=SST,SSR=SSR,SSE=SSE,R2=R2,AdjustedR2=AdjustedR2,PearsonR2=PearsonR2)
  model<-data.frame(X=X, Y=Y, Yhat=Yhat, residuals=residuals)
  return(list(ANOVA=ANOVA, beta=beta, model=model))
}

## use like this
myModel<-elm(X, Y)
myModel[[1]] # to return ANOVA, same as:
myMoel$ANOVA # to return ANOVA