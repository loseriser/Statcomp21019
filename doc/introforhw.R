## -----------------------------------------------------------------------------
#生成100以内正整数简单随机不放回抽样50例，记为y，x为1至50的自然数
set.seed(0)
x=1:50
y= sample(1:100,size=50,replace = F)


#生成y的直方图
hist(y)

#生成x,y的散点图
plot(x,y)


## -----------------------------------------------------------------------------
#生成x和y的一元线性回归模型
model<-lm(x~y)

#输出统计量表格
knitr::kable(summary(model)$coef)


## -----------------------------------------------------------------------------
set.seed(123)
n<-1000
u<-runif(n)
y<-seq(0,20,.01)

#par(mfrow=c(1,2))
#sigma=2的情形，可以看出拟合较好
sig_1<-2


x_1<-sqrt(-2*sig_1^2*log(1-u))
hist(x_1,probability = T,main=expression(sigma==2))
#添加理论密度曲线,图中为红线
lines(y,y/sig_1^2*exp(-y^2/(2*sig_1^2)),col="red")

#sigma=5的情形，可以看出拟合较好
sig_2<-5

x_2<-sqrt(-2*sig_2^2*log(1-u))
hist(x_2,probability = T,main=expression(sigma==5))
#添加理论密度曲线,图中为红线

lines(y,y/sig_2^2*exp(-y^2/(2*sig_2^2)),col="red")


## -----------------------------------------------------------------------------
bif<-function(p)
{
#定义函数bif,输入p_1=p,输出对应混合分布直方图
n<-1e3
X1<-rnorm(n,0,1)
X2<-rnorm(n,3,1)
p_1<-p;p_2<-1-p_1
r<-sample(c(0,1),n,replace = T,prob = c(p_2,p_1))
Z<-r*X1+(1-r)*X2
hist(Z,main="")
}

## -----------------------------------------------------------------------------
set.seed(123)
#输出p_1=0.75时的混合分布直方图
bif(0.75)
title(main="p=0.75",outer = T)
#输出p_1=0.1,...,1时的混合分布直方图


for(i in 1:4)
  bif(0.25*i)
title(main="p=0.1,0.2,...,1.0",outer = T)
#可以看出双峰集中在0.4~0.6
#输出p_1=0.40,...,0.60时的混合分布直方图


for(i in 1:4)
  bif(0.4+0.05*i)
title(main="p=0.40,0.42,...,0.60",outer = T)
#可以看出p_1取值在0.5附近时双峰较为明显

## -----------------------------------------------------------------------------
poi_gamma<-function(n,t,lambda,shape,rate)
{
  ##生成复合泊松分布,输入样本数n,随机过程时刻t,
  ##Poisson过程参数lambda,gamma分布参数shape和rate,
  ##输出样本的均值和方差以及理论均值和方差
  
  
  #生成n个复合泊松过程的样本
  X<-numeric(n)
  for(i in 1:n)
  {
    #生成泊松过程N(t)
    Tn<-rexp(n,lambda)
    Sn<-cumsum(Tn)
    N<-min(which(Sn>t))-1
    Y<-rgamma(N,shape=shape,rate=rate)
    X[i]<-sum(Y)
  }
  #输出样本均值和方差
print(paste0("样本均值:", round(mean(X),2),"  样本方差:", round(var(X),2)))
 #print(c(mean(X),var(X)))
 #输出理论均值和方差
print(paste0("理论均值:", round(lambda*t*shape/rate,2),"  理论方差:",round(lambda*t*shape*(shape+1)/rate^2,2)))
 #print(c(lambda*t*shape/rate,lambda*t*shape*(shape+1)/rate^2))
}

## -----------------------------------------------------------------------------
set.seed(123)
print(paste0("shape =  ", 3,"  rate =  ", 5))
poi_gamma(1000,10,2,3,5)
#模拟较好

#针对不同的shape和rate参数，
for(i in 1:3)
{
  for(j in 1:3)
  #更换不同的shape值和rate值,查看拟合情况
  {
    print(paste0("shape =  ", i,"  rate =  ", j))
    poi_gamma(1000,10,2,i,j)
  }
  
  
}
#拟合较好

## -----------------------------------------------------------------------------
beta_mc<-function(m=1e5,x)
{
  ###m为生成均匀分布样本数,默认取1e5
  ###x取值范围为0,1
  ###输出蒙特卡洛方法得到的B(x;3,3)
  if(x<0 || x>1)
    stop("x is out of range")
  t<-runif(m,min=0,max=x)
  theta.hat<-x*mean(t^2*(1-t)^2)/beta(3,3)
  return(theta.hat)
}

## -----------------------------------------------------------------------------
set.seed(0)
x<-seq(0.1,0.9,0.1)
n<-length(x)
theta.hat<-theta<-numeric(n)
for(i in 1:n)
{
  theta.hat[i]=beta_mc(x=x[i])
  theta[i]=pbeta(q=x[i],shape1 = 3,shape2 = 3)
  print(paste0("x=",x[i],"时,生成值:", round(theta.hat[i],3),",理论值:", round(theta[i],3)))
}


## -----------------------------------------------------------------------------
Rayleigh<-function(sigma, n=10000, antithetic=TRUE)
{
 u <- runif(n)#生成n/2个均匀随机数
  if(!antithetic) v <- runif(n) else
    v <- 1-u#若应用对偶变量，则令后一半均匀随机数v为1-u，否则v为n个均匀随机数且与u独立
  X1 <- sqrt(-2*sigma^2*log(u))#利用u生成一组Rayleigh分布随机数X1
  X2 <- sqrt(-2*sigma^2*log(v))#利用v生成一组Rayleigh分布随机数X2
  X <- (X1+X2)/2#取两组样本均值得到最终估计
}


## -----------------------------------------------------------------------------
set.seed(0)
sigma<-1:5
var_reduction <- matrix(nrow = 1, ncol = length(sigma))#存储结果矩阵
for (i in 1:length(sigma)) {
  simple<-Rayleigh(sigma = sigma[i],antithetic=FALSE)
  antithetic<-Rayleigh(sigma = sigma[i])
  var_reduction[1, i] <- (var(simple)-var(antithetic))/var(simple)
}
dimnames(var_reduction)[[2]] <- paste("sigma =", sigma)#重命名各列
knitr::kable(var_reduction, caption = "表2.1：对偶变量法生成Rayleigh随机数相比简单方法生成Rayleigh随机数方差减少的比例")#绘制结果表格

## -----------------------------------------------------------------------------
g<-function(x) exp(-x^2/2)*x^2/sqrt(2*pi)
f1<-function(x) rep(1,length(x))
f2<-function(x) exp(-x)/(1-exp(-1))

gs <- c(expression(g(x)==x^2*e^{-x^2/2}/sqrt(2*pi)),
            expression(f[1](x)==1),
            expression(f[2](x)==e^{-x}/(1-e^{-1})))
#计算(0,1)内理论积分值
I_g<-integrate(g,0,1)
I_g
#计算理论积分值
1/2-I_g$value   
x <- seq(0, 1, .01)
plot(x,g(x),col=1,ylim=c(0,3),lty=1,type = "l")
lines(x,f1(x),col=2,lty=2)
lines(x,f2(x),col=3,lty=3)
legend("topleft", legend = gs,
           lty = 1:3, inset = 0.02,col=1:3)

## -----------------------------------------------------------------------------
set.seed(0)
  m <- 10000
  est <- sd <- numeric(2)
  
 x <- runif(m) #using f1
  fg <- g(x)
  est[1] <- mean(fg)
  sd[1] <- sd(fg)
u <- runif(m) #f2, inverse transform method
  x <- - log(1 - u * (1 - exp(-1)))
  fg <- g(x) / (exp(-x) / (1 - exp(-1)))
  est[2] <- mean(fg)
  sd[2] <- sd(fg)

## -----------------------------------------------------------------------------
  res <- rbind(est=round(1/2-est,3), sd=round(sd,3))
  print(res)

## -----------------------------------------------------------------------------
set.seed(13)
n<-20
alpha<-0.05
m<-1e4
p<-numeric(m)

for(i in 1:m){
 x<-rchisq(n,2)
if(abs(mean(x)-2)<=sd(x)/sqrt(n)*qt(1-alpha/2,df=n-1)) p[i]=1 
}
mean(p)

## -----------------------------------------------------------------------------
set.seed(123)
n<-20
alpha<-0.05
m<-1e4
UCL<-numeric(m)
#x<-rchisq(n,2)
for(i in 1:m){
 x<-rchisq(n,2)
UCL[i]<-(n-1)*var(x)/qchisq(alpha,df=n-1)

}
mean((UCL>4))

## -----------------------------------------------------------------------------
set.seed(0)
n<-100
m<-1e4
p.val1<-p.val2<-p.val3<-numeric(m)
alpha<-.05
mu<-1
for(i in 1:m)
{
p.val1[i]<-t.test(rchisq(n,df=1),mu=mu,alternative = "two.sided")$p.value
p.val2[i]<-t.test(runif(n,min=0,max=2),mu=mu,alternative = "two.sided")$p.value
p.val3[i]<-t.test(rexp(n),mu=mu,alternative = "two.sided")$p.value
}
##可以看出一型错误率都较小，接近0.05，说明t分布足够稳健
print(round(c(mean(p.val1<=alpha),
mean(p.val2<=alpha),
mean(p.val3<=alpha)),3))


## -----------------------------------------------------------------------------
mat <-
  matrix(c(6510, 3490, 10000, 6760, 3240, 10000, 13270, 6730, 20000), 3, 3,
         dimnames = list(
           c("Rejected", "Accepted", "total"),
           c("method A", "method B", "total")
         ))
mat

## -----------------------------------------------------------------------------
n <- c(10, 20, 30, 50, 100, 500) #sample sizes
cv <- qnorm(.975, 0, sqrt(6/n)) #crit. values for each n
sk <- function(x) {
#computes the sample skewness coeff.
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}

## -----------------------------------------------------------------------------
set.seed(1234)
#n is a vector of sample sizes
#we are doing length(n) different simulations
p.reject <- numeric(length(n)) #to store sim. results
m <- 10000 #num. repl. each sim.
for (i in 1:length(n)) {
sktests <- numeric(m) #test decisions
for (j in 1:m) {
x <- rnorm(n[i])
#test decision is 1 (reject) or 0
sktests[j] <- as.integer(abs(sk(x)) >= cv[i] )
}
p.reject[i] <- mean(sktests) #proportion rejected
}
p.reject

## -----------------------------------------------------------------------------
cv <- qnorm(.975, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
round(cv, 4)


## -----------------------------------------------------------------------------
#n is a vector of sample sizes
#we are doing length(n) different simulations
p.reject <- numeric(length(n)) #to store sim. results
m <- 2500 #num. repl. each sim.
for (i in 1:length(n)) {
sktests <- numeric(m) #test decisions
for (j in 1:m) {
x <- rnorm(n[i])
#test decision is 1 (reject) or 0
sktests[j] <- as.integer(abs(sk(x)) >= cv[i] )
}
p.reject[i] <- mean(sktests) #proportion rejected
}
p.reject

## -----------------------------------------------------------------------------
alpha <- .1
n <- 30
m <- 2500
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { #for each epsilon
e <- epsilon[j]
sktests <- numeric(m)
for (i in 1:m) { #for each replicate
sigma <- sample(c(1, 10), replace = TRUE,
size = n, prob = c(1-e, e))
x <- rnorm(n, 0, sigma)
sktests[i] <- as.integer(abs(sk(x)) >= cv)
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

## -----------------------------------------------------------------------------

nn <- c(10,20,30,50)          # 样本容量
alpha <- 0.05                 # 显著性水平
d <- 2                        # 随机变量的维数
b0 <- qchisq(1-alpha,df=d*(d+1)*(d+2)/6)*6/nn  # 每种样本容量临界值向量

# 计算多元样本偏度统计量
mul.sk <- function(x){
  n <- nrow(x) # 样本个数
  xbar <- colMeans(x) 
  sigma.hat <- (n-1)/n*cov(x) # MLE估计
  
  b <- 0
  for(i in 1:nrow(x)){
    for(j in 1:nrow(x)){
      b <- b+((x[i,]-xbar)%*%solve(sigma.hat)%*%(x[j,]-xbar))^3
    }
  }
  return(b/(n^2))
}

# 计算第一类错误的经验估计
##library(mvtnorm)
set.seed(200)
p.reject <- vector(mode = "numeric",length = length(nn)) # 保存模拟结果

m <- 100

for(i in 1:length(nn)){
  mul.sktests <- vector(mode = "numeric",length = m)
  for(j in 1:m){
    data <- mvtnorm::rmvnorm(nn[i],mean = rep(0,d))
    mul.sktests[j] <- as.integer(mul.sk(data)>b0[i])
  }
  p.reject[i] <- mean(mul.sktests)
}
p.reject

summ <- rbind(nn,p.reject)
rownames(summ) <- c("n","estimate")
knitr::kable(summ)

## -----------------------------------------------------------------------------
#library(boot)
#library(bootstrap)
#library(MASS)
set.seed(0)


b.theta <- function(x,i) {
d<-x[i,]
sigma_hat<-cov(d)/nrow(x)*(nrow(x)-1)
value<-eigen(sigma_hat)$values[1]/sum(eigen(sigma_hat)$values)
return(value)
}
## theta hat为obj$t0

obj <- boot::boot(data=bootstrap::scor,statistic=b.theta,R=2000)
round(c(original=obj$t0,bias=mean(obj$t)-obj$t0,
            se=sd(obj$t)),3)

## -----------------------------------------------------------------------------
n<-nrow(bootstrap::scor)
theta.hat<-b.theta(bootstrap::scor,1:n)
theta.jack<-numeric(n)
for(i in 1:n)
{
  theta.jack[i]<-b.theta(bootstrap::scor,(1:n)[-i])
}
bias.jack<-(n-1)*(mean(theta.jack)-theta.hat)
se.jack<-sqrt((n-1)*mean((theta.jack-theta.hat)^2))
round(c(original=theta.hat,bias.jack=bias.jack,
            se.jack=se.jack),3)

## -----------------------------------------------------------------------------
#percentile confidence intervals

alpha<-c(.025,.975)
quantile(obj$t,alpha,type=6)

#BCa confidence intervals
zalpha<-qnorm(alpha)
z0<-qnorm(sum(obj$t<theta.hat)/length(obj$t))
L<-mean(theta.jack)-theta.jack
a<-sum(L^3)/(6*sum(L^2)^1.5)

adj.alpha<-pnorm(z0+(z0+zalpha)/(1-a*(z0+zalpha)))
quantile(obj$t,adj.alpha,type=6)

## -----------------------------------------------------------------------------
ske<-function(d)
{
dbar<-mean(d)
m3<-mean((d-dbar)^3)
m2<-mean((d-dbar)^2)
return(m3/m2^1.5)
}

b.skewness <- function(x,i) {
d<-x[i,]
return(ske(d))
}

## -----------------------------------------------------------------------------
#library(MASS)
set.seed(1)
nn<-1e3
x<-MASS::mvrnorm(nn,mu=rep(0,2),Sigma = diag(2))

obj <- boot::boot(data=x,statistic=b.skewness,R=2000)
round(c(original=obj$t0,bias=mean(obj$t)-obj$t0,
            se=sd(obj$t)),3)

## -----------------------------------------------------------------------------
alpha<-c(.025,.975)
#standard normal bootstrap confidence interval
sn<-obj$t0+qnorm(alpha)*sd(obj$t)
print(sn)

#basic bootstrap confidence interval
bb<-2*obj$t0-quantile(obj$t,rev(alpha),type = 1)
print(bb)

#percentile bootstrap confidence interval
pb<-quantile(obj$t,alpha,type=6)
print(pb)

## -----------------------------------------------------------------------------
set.seed(1)
m<-1e3
p.val1<-p.val2<-p.val3<-numeric(m)
p.val1l<-p.val2l<-p.val3l<-numeric(m)
p.val1r<-p.val2r<-p.val3r<-numeric(m)
for(i in 1:m)
{
  x<-rnorm(nn)
  skes<-ske(x)
  p.val1l[i]<-as.numeric(skes<sn[1])
  p.val1r[i]<-as.numeric(skes>sn[2])
  
  p.val2l[i]<-as.numeric(skes<bb[1])
  p.val2r[i]<-as.numeric(skes>bb[2])
  
  p.val3l[i]<-as.numeric(skes<pb[1])
  p.val3r[i]<-as.numeric(skes>pb[2])
  
  p.val1[i]<-1-(p.val1l[i]+p.val1r[i])
  p.val2[i]<-1-(p.val2l[i]+p.val2r[i])
  p.val3[i]<-1-(p.val3l[i]+p.val3r[i])
  
}
print(c(mean(p.val1l),mean(p.val2l),mean(p.val3l)))
print(c(mean(p.val1r),mean(p.val2r),mean(p.val3r)))

print(c(mean(p.val1),mean(p.val2),mean(p.val3)))

## -----------------------------------------------------------------------------
set.seed(1)
nn<-1e3
x<-cbind(c(rchisq(nn,df=5),rchisq(nn,df=5)))

obj <- boot::boot(data=x,statistic=b.skewness,R=2000)
round(c(original=obj$t0,bias=mean(obj$t)-obj$t0,
            se=sd(obj$t)),3)

## -----------------------------------------------------------------------------
alpha<-c(.025,.975)
#standard normal bootstrap confidence interval
sn<-obj$t0+qnorm(alpha)*sd(obj$t)
print(sn)

#basic bootstrap confidence interval
bb<-2*obj$t0-quantile(obj$t,rev(alpha),type = 1)
print(bb)

#percentile bootstrap confidence interval
pb<-quantile(obj$t,alpha,type=6)
print(pb)

## -----------------------------------------------------------------------------
set.seed(1)
m<-1e3
p.val1<-p.val2<-p.val3<-numeric(m)
p.val1l<-p.val2l<-p.val3l<-numeric(m)
p.val1r<-p.val2r<-p.val3r<-numeric(m)
for(i in 1:m)
{
  x<-rchisq(nn,df=5)
  skes<-ske(x)
  p.val1l[i]<-as.numeric(skes<sn[1])
  p.val1r[i]<-as.numeric(skes>sn[2])
  
  p.val2l[i]<-as.numeric(skes<bb[1])
  p.val2r[i]<-as.numeric(skes>bb[2])
  
  p.val3l[i]<-as.numeric(skes<pb[1])
  p.val3r[i]<-as.numeric(skes>pb[2])
  
  p.val1[i]<-1-(p.val1l[i]+p.val1r[i])
  p.val2[i]<-1-(p.val2l[i]+p.val2r[i])
  p.val3[i]<-1-(p.val3l[i]+p.val3r[i])
  
}
print(c(mean(p.val1l),mean(p.val2l),mean(p.val3l)))
print(c(mean(p.val1r),mean(p.val2r),mean(p.val3r)))
print(c(mean(p.val1),mean(p.val2),mean(p.val3)))

## -----------------------------------------------------------------------------
set.seed(0)
n<-20
R<-9999
#library(MASS)

reps<-numeric(R)
x<-rnorm(n)
y<-rnorm(n)
t0<-cor(x,y)
z<-c(x,y)
K<-1:20
for(i in 1:R){
  k<-sample(K,size=n,replace = FALSE)
  x1<-z[k]
  y1<-z[-k]
  reps[i]<-cor(x1,y1,method = "spearman")
}
p<-mean(abs(c(t0,reps))>=abs(t0))
round(c(p,cor.test(x,y)$p.value),3)

## -----------------------------------------------------------------------------
#library(RANN)
#library(energy)
#library(Ball)
#library(boot)
Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- RANN::nn2(data=z, k=k+1) # what's the first column?
  block1 <- NN$nn.idx[1:n1,-1] 
  block2 <- NN$nn.idx[(n1+1):n,-1] 
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5) 
  (i1 + i2) / (k * n)
}

## -----------------------------------------------------------------------------
m <- 50; k<-3; p<-2; mu <- 0.3; set.seed(12345)
n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot::boot(data=z,statistic=Tn,R=R,
  sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)

## -----------------------------------------------------------------------------
set.seed(1235)
for(i in 1:m){
  x <- matrix(rnorm(n1*p,0,1.37),ncol=p);
  y <- cbind(rnorm(n2),rnorm(n2));
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- energy::eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- Ball::bd.test(x=x,y=y,num.permutations=999,seed=i*12)$p.value
}
alpha <- 0.1; 
pow <- colMeans(p.values<alpha)
print(pow)

## -----------------------------------------------------------------------------
set.seed(12455)
for(i in 1:m){
  x <- matrix(rnorm(n1*p,0.1,1.37),ncol=p);
  y <- cbind(rnorm(n2),rnorm(n2));
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- energy::eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- Ball::bd.test(x=x,y=y,num.permutations=999,seed=i*123)$p.value
}
alpha <- 0.1; 
pow <- colMeans(p.values<alpha)
print(pow)

## -----------------------------------------------------------------------------
set.seed(12455)
for(i in 1:m){
  x <- matrix(rt(n1*p,1),ncol=p);
  y1 <- rnorm(n1*p,0,1);
  y2 <- rnorm(n1*p,1,2.8);
  r<-sample(c(0,1),n1*p,replace = TRUE)
  y0<-r*y1+(1-r)*y2
  y<-matrix(y0,ncol = p)
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- energy::eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- Ball::bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}
alpha <- 0.1; 
pow <- colMeans(p.values<alpha)
print(pow)

## -----------------------------------------------------------------------------
m <- 1e2; k<-3; p<-2; mu <- 0.3; set.seed(12345)
n1 <-10; n2 <- 90; R<-999; n <- n1+n2; N = c(n1,n2)
eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot::boot(data=z,statistic=Tn,R=R,
  sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- matrix(rnorm(n1*p,0,1.55),ncol=p);
  y <- cbind(rnorm(n2),rnorm(n2));
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- energy::eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- Ball::bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}
alpha <- 0.1; 
pow <- colMeans(p.values<alpha)
print(pow)

## -----------------------------------------------------------------------------
f<-function(x,theta=1,eta=0){
  stopifnot(theta>0)
  return(1/(theta*pi*(1+((x-eta)/theta)^2)))
}

## -----------------------------------------------------------------------------
set.seed(1)
m  <-  10000
sigma  <-2.7
  
x  <-  numeric(m)
x[1]  <-  rnorm(1,mean=0,sd=sigma)
#x[1]<-5
k  <-  0
u  <-  runif(m)
for  (i  in  2:m)  {
xt  <-  x[i-1]
y  <-  rnorm(1,mean=xt,sd=sigma)
num <- f(y) * dnorm(xt, mean=y,sd = sigma)
den <- f(xt) * dnorm(y,mean=xt ,sd=sigma)
if  (u[i]  <= num/den ) x[i]  <-  y  else  {

x[i]  <-  xt
k  <-  k+1         #y  is  rejected
}
}
index <- 5000:5500
y1 <- x[index]
plot(index, y1, type="l", main="", ylab="x")




## -----------------------------------------------------------------------------

###  compare the deciles of the generated observations 
###  with the deciles of the standard Cauchy distribution
b <- 1001 #discard the burnin sample
y <- x[b:m]
a <- ppoints(10)
QR <- qt(a,df=1) #quantiles of Cauchy
Q <- quantile(x, a)
qqplot(QR, Q, main="",
xlab="Cauchy Quantiles", ylab="Sample Quantiles")
qqline(Q)
hist(y, breaks="scott", main="", xlab="", freq=FALSE)
lines(QR, f(QR, theta = 1,eta=0))

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
psi <- as.matrix(psi)
n <- ncol(psi)
k <- nrow(psi)
psi.means <- rowMeans(psi) #row means
B <- n * var(psi.means) #between variance est.
psi.w <- apply(psi, 1, "var") #within variances
W <- mean(psi.w) #within est.
v.hat <- W*(n-1)/n + (B/n) #upper variance est.
r.hat <- v.hat / W #G-R statistic
return(r.hat)
}

## -----------------------------------------------------------------------------
cauchy.chain <- function(sigma, N, x1){
x  <-  numeric(N)
x[1]  <- x1
k  <-  0
u  <-  runif(N)
for  (i  in  2:N)  {
xt  <-  x[i-1]
y  <-  rnorm(1,mean=xt,sd=sigma)

num <- f(y) * dnorm(xt, mean=y,sd = sigma)
den <- f(xt) * dnorm(y,mean=xt ,sd=sigma)
if  (u[i]  <= num/den ) x[i]  <-  y  else  {
x[i]  <-  xt
k  <-  k+1         #y  is  rejected
}
}
return(x)
}

## -----------------------------------------------------------------------------
set.seed(0)
sigma <- 0.2 #parameter of proposal distribution
k <- 4 #number of chains to generate
n <- 10000 #length of chains
b <- 1000 #burn-in length
#choose overdispersed initial values
x0 <- c(-1, -0.5, 0.5, 1)
#generate the chains
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k)
X[i, ] <- cauchy.chain(sigma, n, x0[i])
#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))
#plot psi for the four chains
#par(mfrow=c(2,2))
for (i in 1:k)
plot(psi[i, (b+1):n], type="l",
xlab=i, ylab=bquote(psi))
#par(mfrow=c(1,1)) #restore default
#plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
set.seed(0)
N <- 5000 #length of chain
burn <- 1000 #burn-in length
X <- matrix(0, N, 2) #the chain, a bivariate sample
a=5
b=4
n=20
x=1;y=0.5
bichain<-function(N=5000,burn=1000,a=5,b=4,n=20,x=1,y=0.5,plot=FALSE){
  X <- matrix(0, N, 2)
  X[1, ] <- c(x,y) #initialize
for (i in 2:N) {
X[i, 1] <- rbinom(1,size=n,prob=y)
x<-X[i,1]
X[i, 2] <- rbeta(n=1,shape1 = x+a,shape2 = n-x+b)
y<-X[i,2]
}
b <- burn + 1
x <- X[b:N, ]
if(plot == TRUE)
plot(x, main="", cex=.5, xlab="x",
ylab="y", ylim=range(x[,2]))
else 
return(X)
}
bichain(plot=TRUE)

## -----------------------------------------------------------------------------
library(R2OpenBUGS)
library(MCMCglmm)
DATA1<- mcmc(data=cbind(bichain(x = 1,y=0.2,a=10,N=5000)))
DATA2<- mcmc(data=cbind(bichain(x= 10,y=0.8,a=10,N=5000)))

gelman.plot(x=c(mcmc.list(DATA1),mcmc.list(DATA2)))

## -----------------------------------------------------------------------------
ft<-function(a,n=1e7,eps=1e-9,plot="TRUE")
{ ## 输入参数：
  ## a为d维实值向量
  ## n为最大迭代次数,默认取1e7
  ## eps为精度,默认取1e-9
  
  stopifnot(is.vector(a)) #当a不为向量时报错
  d<-length(a)            #d为a的维度
  
  an<-sqrt(sum(a^2))      #当a的模长为0时,输出0
  
  if(an==0) return ((list(res=0,k=0)))
  else{
    ## 注意到k!=gamma(k-1),k>=2时成立
    ## k=0,1时,k!=1,故将k=0,1情形单独列出
     for(k in 0:1){
    temp1<-(2*k+2)*log(an)+lgamma((d+1)/2)+lgamma(k+3/2)
    temp2<-k*log(2)+log(2*k+1)+log(2*k+2)+lgamma(k+d/2+1)
    temp<-exp(temp1-temp2)
    if(temp<=eps) {
      
      return(list(res=temp,k=k))
      }
    }
    
    for(k in 2:n){
    temp1<-(2*k+2)*log(an)+lgamma((d+1)/2)+lgamma(k+3/2)
    temp2<-lgamma(k-1)+k*log(2)+log(2*k+1)+log(2*k+2)+lgamma(k+d/2+1)
    temp<-exp(temp1-temp2)
    if(temp<=eps) {
      if(plot=="TRUE")
     { cat(c("d=",d,"\n"))
      cat(c("the Euclidean norm of a=",sqrt(sum(a^2)),"\n"))
      cat(c("k=",k,"\n"))
      cat("converge\n")}
      return(list(res=temp,k=k))
      }
    
    }
    #for循环正常结束,次数说明迭代次数较少,可能需要增加n的初始值
    cat("did not converge\n")
    return(list(res=temp,k=k))
  }
  
}

## -----------------------------------------------------------------------------
ff<-function(a,n=1e5,eps=1e-9)
{
  stopifnot(is.vector(a))
  d<-length(a)
  res<-0
  an<-sqrt(sum(a^2))
  if(an==0) return (0)
  else{
     for(k in 0:1){
    temp1<-(2*k+2)*log(an)+lgamma((d+1)/2)+lgamma(k+3/2)
    temp2<-k*log(2)+log(2*k+1)+log(2*k+2)+lgamma(k+d/2+1)
    temp<-exp(temp1-temp2)
    res<-res+(-1)^k*temp
    if(temp<=eps) {
      cat(c("k=",k,"\n"))
      cat("converge\n")
      return(res)
      }
    }
    
    for(k in 2:n){
    temp1<-(2*k+2)*log(an)+lgamma((d+1)/2)+lgamma(k+3/2)
    temp2<-lgamma(k-1)+k*log(2)+log(2*k+1)+log(2*k+2)+lgamma(k+d/2+1)
    temp<-exp(temp1-temp2)
    res<-res+(-1)^k*temp
    if(temp<=eps) {
      cat(c("k=",k,"\n"))
      cat("converge\n")
      return(res)
      }
    
    }
    cat("did not converge\n")
    return(res)
  }
  
}

## -----------------------------------------------------------------------------
a<-seq(0,4,.01)
ff(a)

## -----------------------------------------------------------------------------
a<-c(1,2)
ff(a)

## -----------------------------------------------------------------------------
SS<-function(a,k){
  stopifnot(k>1)
  n<-length(a)
  res<-numeric(n)
  for(i in 1:n){
     if(abs(a[i]^2-k-1)<1e-7) res[i]<-0
      else{
          temp<-sqrt(a[i]^2*k/(k+1-a[i]^2))
          res[i]<-(1-pt(temp,df=k))
  }
  }
 return(res)
}


## -----------------------------------------------------------------------------
for(k in c(4,25,100)){
  aa<-seq(-sqrt(k),sqrt(k),1e-2)
plot(aa,SS(aa,k)-SS(aa,k-1),ylab="",type="l")
abline(h=0)
}


## -----------------------------------------------------------------------------
kk<-c(4:25,100,500,1000)
root<-matrix(0,nrow = length(kk),ncol = 3)
root[,1]<-kk
colnames(root)<-c("k","root_negative","root_positive")
for(i in  1: length(kk)){
  dsa<-function(a){
  return(SS(a,kk[i])-SS(a,kk[i]-1))
}
a<-seq(-sqrt(kk[i]),0,1e-3)
d1<-which.max(dsa(a))
d2<-which.min(dsa(a))
res<-uniroot(dsa, sort(a[c(d1,d2)]))
root[i,2]=unlist(res[1])
root[i,3]<-abs(root[i,2])
}

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(root)

## -----------------------------------------------------------------------------
da<-c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
k<-ifelse(da<1,da,0)
n<-length(da)
n_e<-n-sum(da<1)

estar<-function(lambda){
 return(n_e*(1+lambda))
}
EM<-function(k,n,max.it=1e4,eps=1e-7)
{
 
 
  i<-1
  lambda1<-5
  lambda2<-2
  
  while(abs(lambda1-lambda2)>=eps){
    lambda1<-lambda2
    lambda2<-(sum(k)+estar(lambda=lambda2))/(n)
    print(round(c(lambda2),5))
    if(i == max.it) break
    i<i+1
  }
  return(lambda2)
}
EM(k=k,n=n)

## -----------------------------------------------------------------------------
(n_e+sum(k))/(n-n_e)

## -----------------------------------------------------------------------------
set.seed(0)
trims <- c(0, 0.1, 0.2, 0.5)
x <- rcauchy(100)
lapply(trims, function(trim) mean(x, trim = trim))
lapply(trims, mean, x = x)

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared

## -----------------------------------------------------------------------------
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)

mod<-lapply(formulas, lm,data=mtcars)
r2<-lapply(mod, rsq)
round(unlist(r2),3)

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})
mod<-lapply(bootstraps, lm,formula=mpg ~ disp )
r2<-lapply(mod, rsq)
round(unlist(r2),3)

## -----------------------------------------------------------------------------
set.seed(0)
df<-data.frame(replicate(6,sample(c(1:10),10,rep=T)))
round(vapply(df, sd,FUN.VALUE = double(1)),3)

## -----------------------------------------------------------------------------
set.seed(0)
df<-data.frame(replicate(6,sample(c(1:10),10,rep=T)),y=rep("a",10))

round(vapply(df[vapply(df, is.integer, FUN.VALUE = logical(1))]
             , sd,FUN.VALUE = double(1)),3)

## -----------------------------------------------------------------------------
library(parallel)
nocore<-detectCores()-1
cl<-makeCluster(nocore)


## -----------------------------------------------------------------------------
mcsapply<-function(cl=makeCluster(detectCores()-1),X,FUN){
  parallel::parLapply(cl, X, FUN)
}

mcsapply(cl, 1:4, function(x) x+3)

## -----------------------------------------------------------------------------
mcvapply<-function(cl=makeCluster(detectCores()-1),X,FUN,USE.NAMES =TRUE){
 Y=X[sapply(X, function(x) typeof(x)==typeof(USE.NAMES))]
 if(length(Y)==0) return("ERROR")
 else(parLapply(cl, Y, FUN))
 
}
mcvapply(cl, 1:4, function(x) x+3,USE.NAMES = integer(1))
mcvapply(cl, 1:4, function(x) x+3,USE.NAMES = TRUE)

## -----------------------------------------------------------------------------
set.seed(0)
#library(Rcpp)
Rcpp::cppFunction(
'NumericMatrix gibbsC(int N) {
  NumericMatrix mat(N, 2);
  double x = 0, y = 0, a=5,b=4,n=20;
  for(int i = 0; i < N; i++) {
   
      x = rbinom(1,n,y)[0];
      y = rbeta(1,x+a,n-x+b)[0];
 
    mat(i, 0) = x;
    mat(i, 1) = y;
  }
  return(mat);
}'
)

## -----------------------------------------------------------------------------
data<-gibbsC(1000)
plot(data)

## -----------------------------------------------------------------------------
Rcpp::cppFunction('NumericVector cunif(int n,double min,double max){
            NumericVector x;
            x=runif(n,min,max);
            return(x);
            }')

## -----------------------------------------------------------------------------
set.seed(0)
cc<-cunif(n=1000,min=-1,max=1)
rr<-runif(n=1000,min=-1,max=1)
qqplot(cc,rr)

## -----------------------------------------------------------------------------
 #   library(microbenchmark)
ts <- microbenchmark::microbenchmark(cc<-cunif(n=1000,min=-1,max=1),
rr<-runif(n=1000,min=-1,max=1))
    summary(ts)[,c(1,3,5,6)]

