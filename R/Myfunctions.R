#' @title Geometric Brownian Motion sample
#' @description Use Geometric Brownian Motion to simulate fructuation of stocks using R.
#' @param SO initial number (numeric)
#' @param mu drift(numeric vector)
#' @param sigma volatility(numeric vector)
#' @param T time interval(numeric)
#' @param numSteps length of the chain(numeric)
#' @param numRep1 number of the chain(numeric)
#' @return a Geometric Brownian Motion sample
#' @import quantmod
#' @import MASS bootstrap boot DAAG RANN energy Ball stargazer jmuOutlier xtable MAP microbenchmark
#' @importFrom stats rnorm lm predict
#' @importFrom flextable xtable_to_flextable
#' @useDynLib StatComp20093
#' @examples
#' \dontrun{
#' SO <- 50
#' mu <- 0.1
#' sigma <- 0.3
#' T <- 1
#' numSteps <- 1000
#' numRep1 <- 10
#' path = matrix(nrow=numRep1,ncol=numSteps+1)
#' path <- simGBM(SO,mu,sigma,T,numSteps,numRep1 )
#' for(i in 1:numRep1) {
#' if (i == 1)
#' plot(0:numSteps, path[i,], col = i, type='l', ylim = c(0,100))
#' else
#' lines(0:numSteps,path[i,], col = i)
#' }  
#'}
#' @export
simGBM <- function (SO,mu,sigma,T,numSteps,numRep1)
{
  dt <- T/numSteps
  nuT <- (mu-sigma^2/2)*dt
  sigmaT <- sqrt(dt)*sigma
  pathMatrix = matrix(nrow=numRep1,ncol=numSteps+1)
  pathMatrix[,1] <- SO
  for (i in 1:numRep1){
    for (j in 2:numSteps){
      
      pathMatrix[i,j] <- pathMatrix[i,j-1]*
        exp(rnorm(1,nuT,sigmaT))
      
    }
  }
  return(pathMatrix)
}

#' @title Use three inputs to predict response using R.
#' @description Calculate the cumulative abnormal return of the stocks using event study method.
#' @param r individual stock return(numeric vector)
#' @param rm market return(numeric vector)
#' @param m1 days of the start of the model estimation period before the event(numeric)
#' @param m2 days of the end of the model estimation period before the event(numeric)
#' @param c1 days of the start of the model window period before the event(numeric)
#' @param c2 days of the start of the model window period before the event(numeric)
#' @param date date of the event(date)
#' @return a random sample of size \code{n}
#' @import quantmod
#' @import MASS bootstrap boot DAAG RANN energy Ball stargazer jmuOutlier xtable MAP microbenchmark
#' @importFrom stats rnorm lm predict
#' @importFrom zoo index
#' @importFrom flextable xtable_to_flextable
#' @useDynLib StatComp20093
#' @examples
#' \dontrun{
#' library(quantmod)
#' x<-'XF'
#' y<-'HS300'
#' d1<-'2019-01-02'
#' d2<-'2020-12-12'
#' getSymbols("x",from=d1,to=d2)
#' getSymbols("y",from=d1,to=d2)
#' r <- periodReturn(X$`X.Adjusted`,period='daily',type='log')
#' rm <- periodReturn(Y$`Y.Adjusted`,period='daily',type='log')
#' do_car(r,rm,10,5,1,1,'2020-01-02')
#' }
#' @export
do_car <- function(r,rm,m1,m2,c1,c2,date) { 
  stopifnot(m1 > m2) 
  date<-as.Date(date)
  n<-max(which(index(r)<=date))
  if (n - m1 < 0) { cat("n =", n, "is too small ") } 
  else if (n + c2 > length(r)) { 
    cat("n =", n, "is too large ") } 
  else { 
    i1 <- max(1, n - m1) 
    i2 <- n - m2 
    i3 <- n - c1 
    i4 <- n + c2 
    r.model <- r[i1:i2] 
    rm.model <- rm[i1:i2] 
    r.car <- r[i3:i4] 
    rm.car <- rm[i3:i4] 
    model <- lm(r.model ~ I(r.model - rm.model)) 
    coef <- coef(model) 
    ars <- r.car - predict(model, list(r.model = r.car, rm.model = rm.car)) 
    list(date = date, coef = list(coef), ars = list(ars)) 
  }
}

