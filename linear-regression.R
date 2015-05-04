elm<-function(X,Y) {
  X <- as.matrix(X)    # predictor set X
  Y <- as.matrix(Y)    # response vector Y
  Xaug <- cbind(1,X)   # augmented predictor matrix
  Yhat <- Xaug %*% solve(t(Xaug) %*% Xaug) %*% t(Xaug) %*% Y # least-squares solution
  residuals <- Y - Yhat                                # error term
  hatm <- Xaug %*% solve(t(Xaug) %*% Xaug) %*% t(Xaug) # hat matrix
  beta <- solve(t(Xaug) %*% Xaug) %*% t(Xaug) %*% Y    # model coefficients
  SST = sum((Y - mean(Y))^2)    # deviation of observed values from sample mean
  SSR = sum((Yhat - mean(Y))^2) # deviation of fitted values from sample mean
  SSE = sum((Yhat - Y)^2)       # deviation of model fit from observed response
  n<-dim(X)[1]    # number of observations
  k<-dim(X)[2]    # number of predictors
  p<-length(beta) # number of model coefficients
  R2 = (SST-SSE)/SST                       # R2
  PearsonR2 = (cor(Yhat,Y))^2              # equivalent to R2
  AdjustedR2 = 1-((SSE/(n-p))/(SST/(n-1))) # adjusted R2

  # objects to return...
  ANOVA<-data.frame(SST=SST, SSR=SSR, SSE=SSE, R2=R2, AdjustedR2=AdjustedR2, PearsonR2=PearsonR2)
  model<-data.frame(X=X, Y=Y, Yhat=Yhat, residuals=residuals)
  return(list(ANOVA=ANOVA, beta=beta, model=model))
}

## example
X<-1:100
Y<-1:100 + rnorm(100, sd=10, mean=10)

myModel<-elm(X, Y)
myModel[[1]]  # return ANOVA, same as:
myModel$ANOVA # return ANOVA
myModel[[2]]  # betas
myModel$beta