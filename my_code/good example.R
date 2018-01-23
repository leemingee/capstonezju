library("MSBVAR")
library("splines")
library("np")
#library("sm")
library("MASS")
library("quantreg")
library("VGAM")

n=100
nrep=1
dx=500
tau=0.5


######AR(1) covariance matrix for z
times <- 1:dx
rhoz <- 0.5
sigmaz <- 1
H <- abs(outer(times, times, "-"))
sigz <- sigmaz * rhoz^H
pp <- nrow(sigz)
#sigz[cbind(1:pp, 1:pp)] <- sigz[cbind(1:pp, 1:pp)] * sigmaz
sigz<-diag(rep(1-rhoz,dx))+matrix(rhoz,dx,dx)
sigz[4,-4]=sqrt(rhoz)
sigz[-4,4]=sqrt(rhoz)
#sigz[5,-5]=0
#sigz[-5,5]=0

muz<-rep(0,dx)
#Mstar=c(2,200,500,800,1000)
#Mstar=c(1,30,60,80,100)
Mstar=c(1,2,3,4)
d=floor(n/log(n))
jd=1

Mhat1=matrix(NA,nrow=nrep,ncol=d)

set.seed(1)
for(ii in 1:nrep)
{
  x=rmultnorm(n, muz, sigz)
  
  if(jd==1){
    z=rnorm(n)
    Qtau=qnorm(tau)
  }
  if(jd==2){
    z=rlaplace(n, location=0, scale=1)
    Qtau=qlaplace(tau)
  }
  
  #y=apply(x[,Mstar],1,sum)+z
  beta=2.5+2.5*abs(tau-0.5)
  y=beta*x[,1]+beta*x[,2]+beta*x[,3]-3*beta*sqrt(rhoz)*x[,4]+z
  #y=beta*x[,1]+beta*x[,2]+beta*x[,3]-3*beta*sqrt(rhoz)*x[,4]+0.25*beta*x[,5]+z
  
  alpha=0.05
  
  dstar=floor(2*(n/log(n))^0.5)
  cc=1
  CORM=correlation(y,x,n,dx,tau,d,Mstar,alpha,cc)
  CORM=CORM[[1]]
  
  result=PQCS1(y,x,n,dx,tau,dstar,d,Mstar,alpha,CORM)
}
result
