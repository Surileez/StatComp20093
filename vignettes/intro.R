## -----------------------------------------------------------------------------
lm.WH<-lm(women$weight~women$height)
summary(lm.WH)$coef

## -----------------------------------------------------------------------------
#par(mfrow=c(2,2))
plot(lm.WH)

## ----results='asis'-----------------------------------------------------------
print(xtable::xtable(head(women)),type='html')

## -----------------------------------------------------------------------------
n <- 1000 
u <- runif(n) 
x <- 2/((1-u)^0.5) # F(x) =1-(2/x)^2,x>=2
hist(x, prob = TRUE, main = expression(f(x)==8/x^3),
xlim=c(0,20),breaks=100)
y<-seq(0,20,.1)
lines(y,8/(y^3),col='blue')

## -----------------------------------------------------------------------------
n <- 10000 
u1 <- runif(n,-1,1) 
u2 <- runif(n,-1,1) 
u3 <- runif(n,-1,1) 
x<-ifelse(abs(u3)>=abs(u2)&abs(u3)>=abs(u1),u2,u3)
hist(x, prob = TRUE, main = expression(f(x)==3(1-x^2)/4),breaks=50)
y<-seq(-1,1,.1)
lines(y,3*(1-y^2)/4,col='blue')

## -----------------------------------------------------------------------------
n <- 1000 
u <- runif(n) 
x <- 2/((1-u)^0.25)-2 # F(x)=1-(2/(2+x))^4,x>=0
hist(x, prob = TRUE, main = expression(f(x)==64/(2+x)^5),
xlim=c(0,10),breaks=100)
y<-seq(0,10,.1)
lines(y,64/((2+y)^5),col='blue')

## -----------------------------------------------------------------------------
n <- 10000; 
t <- runif(n, min=0, max=pi/3) #t~U(0,pi/3)
theta_hat<- mean(pi/3*sin(t)) #Monte Carlo estimate 
print(c(theta_hat,-cos(pi/3)+cos(0)))#compare

## -----------------------------------------------------------------------------
antithetic<-function(n,anti=TRUE){
  x1 <- runif(n);
  u <- x1[1:n/2]
  v <- 1-u;
  x2<- c(u,v)
  if(anti)
  return (mean(exp(x2)))
  else
  return (mean(exp(x1)))
}

n <- 10000; 
#antithetic variate approach
theta_anti<-antithetic(n) 
#simple Monte Carlo method
theta_simple<-antithetic(n,anti=FALSE) 
m<-1000
MC_anti<-MC_simple<-numeric(m)
for(i in 1:m){
  MC_anti[i]<-antithetic(n);
  MC_simple[i]<-antithetic(n,anti=FALSE)
}
#Compare the result with the theoretical value
print(c(theta_anti,theta_simple,exp(1)-exp(0)))
#Compute an empirical estimate of the percent reduction in variance using the antithetic variate
print(c(sd(MC_anti),sd(MC_simple),sd(MC_anti)/sd(MC_simple),1-sd(MC_anti)/sd(MC_simple)))

## -----------------------------------------------------------------------------
m<- 1e4
g<-function(x){x^2/sqrt(2*pi)*exp(-x^2/2)*(x>1)}
#choose f1 as importance sampling
x<-rnorm(m)
fg1<-g(x)/dnorm(x)
theta.hat1<-mean(fg1)
#choose f2 as important sampling
x<-rexp(m,1/2)
fg2<-g(x)/dexp(x,1/2)
theta.hat2<-mean(fg2)
print(c(theta.hat1,theta.hat2))#theta results
integrate(g,1,Inf)#true integral
print(c(sd(fg1),sd(fg2)))#variance

## -----------------------------------------------------------------------------
M <- 10000
N <- 50 
k <- 5 
r <- M/k 
T5 <- numeric(k)
est <- matrix(0, N, 2)
#use reverse transform method
g<-function(x,a,b) exp(-x)/(1+x^2)*(x>a)*(x<b)
h<-function(u,a,b) -log(exp(-a)-u*(exp(-a)-exp(-b)))
fg<-function(x,a,b) g(x,a,b)/(exp(-x)/(exp(-a)-exp(-b)))
for (i in 1:N) {
u<-runif(M)
u.s<-runif(M/k)
#importance sampling
est[i, 1] <- mean(fg(h(u,0,1),0,1))
#stratified importance sampling
for(j in 1:k) T5[j]<-mean(fg(h(u.s,(j-1)/k,j/k),(j-1)/k,j/k))
est[i, 2] <- sum(T5)
}
apply(est,2,mean)
apply(est,2,sd)

## -----------------------------------------------------------------------------
n<-20
alpha<-0.05
#symmetric t-interval
t_UCL<-replicate(1000, expr = { 
  x <- rchisq(n,2)
  c(mean(x),(var(x)/n)^(0.5)*qt(1-alpha/2, df = n-1)) })
#Example 6.4, use chisq sample by assuming it as normal distribution
UCL <- replicate(1000, expr = { 
  x <- rchisq(n,2)
  (n-1) * var(x) / qchisq(alpha, df = n-1) })
#Example 6.4, use normal distribution sample
UCLN <- replicate(1000, expr = { 
  x <- rnorm(n,0,2)
  (n-1) * var(x) / qchisq(alpha, df = n-1) })
#coverage probability
mean((t_UCL[1,]+t_UCL[2,])>2&(t_UCL[1,]-t_UCL[2,])<2)
mean(UCL > 4)
mean(UCLN > 4)

## ----warning=FALSE------------------------------------------------------------
sk <- function(x) {
#computes the sample skewness coeff.
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
n <- 30
m <- 2500
a<-0.1
alpha <- c(seq(1, 100, 1))
N <- length(alpha)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qnorm(1-a/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { #for each epsilon
e <- alpha[j]
sktests <- numeric(m)
for (i in 1:m) { #for each replicate
x <- rbeta(n, e, e)
sktests[i] <- as.integer(abs(sk(x)) >= cv)
}
pwr[j] <- mean(sktests)
}
#plot power vs alpha
plot(alpha, pwr, type = "b",
xlab = bquote(alpha))
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(alpha, pwr+se, lty = 3)
lines(alpha, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
n <- 30
m <- 2500
a<-0.1
v<- c(seq(1, 100, 1))
N <- length(v)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qnorm(1-a/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { #for each epsilon
e <- v[j]
sktests <- numeric(m)
for (i in 1:m) { #for each replicate
x <- rt(n, e)
sktests[i] <- as.integer(abs(sk(x)) >= cv)
}
pwr[j] <- mean(sktests)
}
#plot power vs v
plot(v, pwr, type = "b",
xlab = bquote(v))
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(v, pwr+se, lty = 3)
lines(v, pwr-se, lty = 3)

## ----warning=FALSE------------------------------------------------------------
#Count Five test
count5test <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) > 5))
}
# generate samples under H1 to estimate power
sigma1 <- 1
sigma2 <- 1.5
m <- c(10,500,100000)
power<-length(m)
for(i in 1:3){
power[i] <- mean(replicate(m[i], expr={
x <- rnorm(20, 0, sigma1)
y <- rnorm(20, 0, sigma2)
count5test(x, y)
}))
}
print(power)

## -----------------------------------------------------------------------------
#F test
sigma1 <- 1
sigma2 <- 1.5
alpha<-0.055
m <- c(10,500,100000)
power<-length(m)
for(i in 1:3){
power[i] <- mean(replicate(m[i], expr={
x <- rnorm(20, 0, sigma1)
y <- rnorm(20, 0, sigma2)
as.integer(var.test(x,y)$p.value <= alpha)
}))
}
print(power)

## -----------------------------------------------------------------------------
#6.8
#Consider multiple normal distribution
library(MASS)
multi_sk <- function(X){
n <- nrow(X)
xbar <- colMeans(X)
sigma.hat <- cov(X) * (n - 1) / n
b <- sum(((t(t(X) - xbar))%*%solve(sigma.hat)%*%(t(X) - xbar))^3) / n^2
return (b)
}
d=2 # the dimension
n=c(10,20,30,50,100,500) # sample sizes
cv=qchisq(0.95,d*(d+1)*(d+2)/6)
p.reject=numeric(length(n)) 
m=1000 
sig<-matrix(c(10,0,0,10),nrow=2)
for (i in 1:length(n)) {
  sktests=numeric(m) 
  for (j in 1:m) {
    x=mvrnorm(n[i],rep(0,d),sig) 
    sktests[j]=as.integer(abs(n[i]*multi_sk(x)/6) >= cv )
  }
  p.reject[i]=mean(sktests) #proportion rejected
}
p.reject

## -----------------------------------------------------------------------------
#6.10
d<-2
library(MASS)
alpha <- .1
n <- 30
m <- 500
epsilon <- c(seq(0, .15, .02), seq(.15, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
cv<-qchisq(1-alpha,d*(d+1)*(d+2)/6)
for (j in 1:N) { #for each epsilon
e <- epsilon[j]
sktests <- numeric(m)
 for (i in 1:m) { #for each replicate
 sigma <- sample(c(1, 10), replace = TRUE,
 size = n, prob = c(1-e, e))
 x=matrix(0,n,d)
  for(k in 1:n){
  x[k,] <- mvrnorm(1, rep(0,2), diag(sigma[k]^2,d))}
 sktests[i] <- as.integer(n*abs(multi_sk(x))/6 >= cv)
}
pwr[j] <- mean(sktests)
}
#plot power vs epsilon
plot(epsilon, pwr, type = "b",
xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## ----waring=FALSE-------------------------------------------------------------
library(bootstrap)
n <- nrow(law) #sample size
theta_hat <- cor(law$LSAT, law$GPA)
theta_j <- numeric(n)
for (i in 1:n) {
#randomly select the indices
  x<-law[-i,]
  LSAT <- x$LSAT
  GPA <- x$GPA
  theta_j[i] <- cor(LSAT, GPA)
}
bias_jack <- (n - 1) * (mean(theta_j) - theta_hat) 
#jackknife bias
se_jack <- (n - 1) * sqrt(var(theta_j) / n)
#jackknife standard error 
cat('\n','Jackknife bias',bias_jack,'\n', 'Jackknife standard error ',se_jack)

## ----warning=FALSE------------------------------------------------------------
library(boot)
mean(aircondit[,1])#theta
N=200#boot sample size
bootstrap_result<-boot(aircondit,statistic= function(x,ind){mean(x[ind,1])},R=N)
boot.ci(bootstrap_result, conf = 0.95, type = 'all')

## ----waring=FALSE-------------------------------------------------------------
library(bootstrap)
library(boot)
lambda_hat<-eigen(cov(scor))$values
theta_hat <- lambda_hat[1] / sum(lambda_hat) 
N <- 200 #sample size
n <- nrow(scor) 
# Jackknife 
theta_j <- rep(0, n) 
for (i in 1:n) { 
   x <- scor [-i,] 
   lambda <- eigen(cov(x))$values 
   theta_j[i] <- lambda[1] / sum(lambda) 
   } 
bias_jack <- (n - 1) * (mean(theta_j) - theta_hat) 
#jackknife bias
se_jack <- (n - 1) * sqrt(var(theta_j) / n)
#jackknife standard error
cat('\n','Jackknife bias',bias_jack,'\n','Jackknife standard error',se_jack )

## ----include=FALSE------------------------------------------------------------
library(DAAG)
attach(ironslag)

## ----warning=FALSE------------------------------------------------------------
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n*(n-1)/2)
t=1
a <- seq(10, 40, .1) #sequence for plotting fits
# for n-fold cross validation
# fit models on leave-two-out samples
for (k in 1:(n-1)) 
 for(j in (k+1):n){
   y <- magnetic[c(-k,-j)]
   x <- chemical[c(-k,-j)]
   J1 <- lm(y ~ x)
   yhat1 <- J1$coef[1] + J1$coef[2]*chemical[c(k,j)]
   e1[t] <- mean(magnetic[c(k,j)]-yhat1)
   J2 <- lm(y ~ x + I(x^2))
   yhat2 <- J2$coef[1] + J2$coef[2] * chemical[c(k,j)] +J2$coef[3] * chemical[c(k,j)]^2
   e2[t] <- mean(magnetic[c(k,j)] - yhat2)
   J3 <- lm(log(y) ~ x)
   logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[c(k,j)]
   yhat3 <- exp(logyhat3)
   e3[t] <- mean(magnetic[c(k,j)] - yhat3)
   J4 <- lm(log(y) ~ log(x))
   logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[c(k,j)])
   yhat4 <- exp(logyhat4)
   e4[t] <- mean(magnetic[c(k,j)] - yhat4)
   t=t+1
}
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))


## -----------------------------------------------------------------------------
L2 <- lm(magnetic ~ chemical + I(chemical^2))
L2

## -----------------------------------------------------------------------------
detach(ironslag)

## -----------------------------------------------------------------------------
count5test <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) > 5))
}
n1 <- 20;n2 <- 30#sample sizes are not equal
mu1 <- mu2 <- 0;sigma1 <- sigma2 <- 1
m <- 1000;R <- 999;K <- 1:50
reps <- numeric(R)
tests <- replicate(m, expr = {
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
z <-c(x,y)
for (i in 1:R) {
k <- sample(K, size = n1, replace = FALSE)
x1 <- z[k]
y1 <- z[-k] #complement of x1
k <- sample(1:n2, size = n1, replace = FALSE)
y1 <- y1[k] #control the sample size
reps[i] <- count5test(x1, y1)
}
mean(reps)
} )
alphahat <- mean(tests)
alphahat

## ----warning=FALSE------------------------------------------------------------
library(RANN)
library(energy)
library(boot)
Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) 
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}

## ----warning=FALSE,include=FALSE----------------------------------------------
library(RANN)
library(energy)
library(Ball)
library(boot)
library(stargazer)

## ----warning=FALSE------------------------------------------------------------
m <- 100; k<-3; set.seed(12345)
mu1 <- mu2 <- 0#equal expectations
sigma1 <-1; sigma2 <- 2#unequal variances 
n1 <- n2 <- 20; R<-999; n <- n1+n2; N = c(n1,n2)
eqdist.nn <- function(z,sizes,k){
boot.obj <- boot(data=z,statistic=Tn,R=R,
sim = "permutation", sizes = sizes,k=k)
ts <- c(boot.obj$t0,boot.obj$t)
p.value <- mean(ts>=ts[1])
list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
x <- rnorm(n1, mu1, sigma1)#N(0,1)
y <- rnorm(n2, mu2, sigma2)#N(0,2)
z <- c(x, y) 
p.values[i,1] <- eqdist.nn(z,N,k)$p.value
p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12345)$p.value
}
alpha <- 0.1;
pow <- colMeans(p.values<alpha)
names(pow)<-c('NN','energy','Ball')

## ----results='asis'-----------------------------------------------------------
stargazer(pow,type='html')

## ----waring=FALSE-------------------------------------------------------------
mu1<-0;mu2<-1#unequal expectations
sigma1 <-1;sigma2<-2#unequal variances 
set.seed(12345)
for(i in 1:m){
x <- rnorm(n1, mu1, sigma1)#N(0,1)
y <- rnorm(n2, mu2, sigma2)#N(1,2)
z <- c(x, y) 
p.values[i,1] <- eqdist.nn(z,N,k)$p.value
p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12345)$p.value
}
pow <- colMeans(p.values<alpha)
names(pow)<-c('NN','energy','Ball')

## ----results='asis'-----------------------------------------------------------
stargazer(pow,type='html')

