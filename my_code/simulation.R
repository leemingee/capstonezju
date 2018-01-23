library("MSBVAR")
library("splines")
library("np")
#library("sm")
library("MASS")
library("quantreg")
library("VGAM")

n=200
nrep=200
dx=1000
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
MhatBIC1=matrix(NA,nrow=nrep,ncol=d)
MhatBIC1C2=matrix(NA,nrow=nrep,ncol=d)

Mhat2=matrix(NA,nrow=nrep,ncol=d)
MhatBIC2=matrix(NA,nrow=nrep,ncol=d)
MhatBIC2C2=matrix(NA,nrow=nrep,ncol=d)

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
     alpha1=0.05


  dstar=floor(2*(n/log(n))^0.5)
   cc=1
  CORM=correlation(y,x,n,dx,tau,d,Mstar,alpha,cc)
  CORM=CORM[[1]]


  result=PQCS1(y,x,n,dx,tau,dstar,d,Mstar,alpha,CORM)
  RjS1[ii,]=result[[2]]
  MS1[ii]=result[[3]]
  pjS1[ii,]=result[[4]]
  paS1[ii]=result[[5]]
  A=result[[1]]
  Mhat1[ii,]=A
  BIC1=BIC(y,x,n,dx,tau,dstar,d,Mstar,alpha,A)
  MhatBIC=BIC1[[1]]
  MhatBIC1[ii,1:length(MhatBIC)]=MhatBIC
  BIC1=BIC2(y,x,n,dx,tau,dstar,d,Mstar,alpha,A)
  MhatBIC=BIC1[[1]]
  MhatBIC1C2[ii,1:length(MhatBIC)]=MhatBIC
  
  result=PQCS2(y,x,n,dx,tau,dstar,d,Mstar,alpha,CORM)
  #result=PQCS2new(y,x,n,dx,tau,dstar,d,Mstar,alpha,cc)
  RjS2[ii,]=result[[2]]
  MS2[ii]=result[[3]]
  pjS2[ii,]=result[[4]]
  paS2[ii]=result[[5]]
  A=result[[1]]
  Mhat2[ii,]=A
  BIC22=BIC(y,x,n,dx,tau,dstar,d,Mstar,alpha,A)
  MhatBIC=BIC22[[1]]
  MhatBIC2[ii,1:length(MhatBIC)]=MhatBIC
  BIC22=BIC2(y,x,n,dx,tau,dstar,d,Mstar,alpha,A)
  MhatBIC=BIC22[[1]]
  MhatBIC2C2[ii,1:length(MhatBIC)]=MhatBIC
}


CBIC=seq(0,0,length=nrep)
OBIC=seq(0,0,length=nrep)

CBIC2=seq(0,0,length=nrep)
OBIC2=seq(0,0,length=nrep)

TPBIC=seq(0,0,length=nrep)
FPBIC=seq(0,0,length=nrep)
TPBIC2=seq(0,0,length=nrep)
FPBIC2=seq(0,0,length=nrep)

for(ii in 1:nrep){
BIC1=MhatBIC2[ii,]
BIC1=na.omit(BIC1)
INM=intersect(BIC1,Mstar)
CBIC[ii]=(length(INM)==length(Mstar))*(length(INM)==length(BIC1))
OBIC[ii]=(length(INM)==length(Mstar))*(length(INM)<length(BIC1))
TPBIC[ii]=length(INM)
FPBIC[ii]=length(BIC1)-length(INM)

BIC1=MhatBIC2C2[ii,]
BIC1=na.omit(BIC1)
INM=intersect(BIC1,Mstar)
CBIC2[ii]=(length(INM)==length(Mstar))*(length(INM)==length(BIC1))
OBIC2[ii]=(length(INM)==length(Mstar))*(length(INM)<length(BIC1))
TPBIC2[ii]=length(INM)
FPBIC2[ii]=length(BIC1)-length(INM)

}

mean(CBIC)
mean(OBIC)
mean(CBIC2)
mean(OBIC2)
mean(TPBIC)
mean(FPBIC)
mean(TPBIC2)
mean(FPBIC2)
