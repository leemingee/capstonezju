PQCSIN2BIC2<-function(n,p,ro,tau,alpha)
{
  #sigma
  sigma<-matrix(0,p,p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma[i,j]=ro^abs(i-j)
    }
  }
  library(MASS)
  library(quantreg)
  d0=floor(2*(n/log(n))^0.5)
  d1=floor(n/log(n)) #选出predictors的最终维数,least
  am=matrix(nrow=re,ncol=d1)
  place=matrix(nrow=re,ncol=3)
  z=c(1:p) #原始变量集
  
    xi=mvrnorm(n, rep(0,p), sigma)
    x=xi^2
    e=rnorm(n)
    yi=yi=3*xi[,1]*xi[,2]+3*xi[,1]*xi[,3]+e  #2*xi[,1]+3*xi[,1]*xi[,2]+3*xi[,1]*xi[,3]+e  
    y=yi^2
    CORM=correlation(y,x,n,p,alpha)
    ##screening procedure
    a=vector(length=0) #起始的选出变量集
    for (i in 1:d1) {
      b=setdiff(z,a) #该轮要与y做qpcor的变量集
      qpcor=vector(length=length(b))
      if(i<=d0){
        for (j in b) {
          sj=CORM[j,]
          sj=na.omit(sj)
          sj=union(sj,a) #计算qpcor时的条件集合
          sj=unique(sj)
          xsj=cbind(1,x[,sj])
          thetaj=ginv(t(xsj)%*%xsj)%*%t(xsj)%*%x[,j] #coef(lm(x[,j]~x[,sj])) 
          piej=coef(rq(y~x[,sj],tau=tau))
          thegmaj=sum((x[,j]-xsj%*%thetaj)^2)/n
          qpcorj=sum((tau-funci(y-xsj%*%piej))*(x[,j]-xsj%*%thetaj))/(n*sqrt((tau-tau^2)*thegmaj))
          qpcor[which(b==j)]=qpcorj
        }
        a[i]=b[which.max(abs(qpcor))]}
      else{
        for (j in b) {
          sj=CORM[j,]
          sj=na.omit(sj)
          sj=union(sj,a[1:d0]) #计算qpcor时的条件集合
          sj=unique(sj)
          thetaj=coef(lm(x[,j]~x[,sj]))
          piej=coef(rq(y~x[,sj],tau=tau))
          xsj=cbind(1,x[,sj])
          thegmaj=sum((x[,j]-xsj%*%thetaj)^2)/n
          qpcorj=sum((tau-funci(y-xsj%*%piej))*(x[,j]-xsj%*%thetaj))/(n*sqrt((tau-tau^2)*thegmaj))
          qpcor[which(b==j)]=qpcorj
        }
        a[i]=b[which.max(abs(qpcor))]}
    }
   ABIC=BIC2(y,x,n,tau,a)
   return(list(ABIC,a))
}



