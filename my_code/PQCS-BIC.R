PQCSBIC<-function(n,p,ro,tau,alpha)
{
  sigma<-matrix(ro,p,p)
  for (i in 1:p) {
    sigma[i,4]=ro^0.5
    sigma[4,i]=ro^0.5
    sigma[i,i]=1
  }
  library(MASS)
  library(quantreg)
  d0=floor(2*(n/log(n))^0.5)
  d1=floor(n/log(n)) #选出predictors的最终维数,least
  belta=2.5*(1+abs(tau-0.5))
  z=c(1:p) #原始变量集
  
    x=mvrnorm(n, rep(0,p), sigma)
    e=rnorm(n)
    y=belta*x[,1]+belta*x[,2]+belta*x[,3]-3*belta*ro^0.5*x[,4]+e
    CORM=correlation(y,x,n,p,alpha)
    ##screening procedure
    A=vector(length=0) #起始的选出变量集
    for (i in 1:d1) {
      b=setdiff(z,A) #该轮要与y做qpcor的变量集
      qpcor=vector(length=length(b))
      if(i<=d0){
        for (j in b) {
          sj=CORM[j,]
          sj=na.omit(sj)
          sj=union(sj,A) #计算qpcor时的条件集合
          sj=unique(sj)
          xsj=cbind(1,x[,sj])
          thetaj=coef(lm(x[,j]~x[,sj])) #ginv(t(xsj)%*%xsj)%*%t(xsj)%*%x[,j] 
          piej=coef(rq(y~x[,sj],tau=tau))
          thegmaj=sum((x[,j]-xsj%*%thetaj)^2)/n
          qpcorj=sum((tau-funci(y-xsj%*%piej))*(x[,j]-xsj%*%thetaj))/(n*sqrt((tau-tau^2)*thegmaj))
          qpcor[which(b==j)]=qpcorj
        }
        A[i]=b[which.max(abs(qpcor))]}
      else{
        for (j in b) {
          sj=CORM[j,]
          sj=na.omit(sj)
          sj=union(sj,A[1:d0]) #计算qpcor时的条件集合
          sj=unique(sj)
          thetaj=coef(lm(x[,j]~x[,sj]))
          piej=coef(rq(y~x[,sj],tau=tau))
          xsj=cbind(1,x[,sj])
          thegmaj=sum((x[,j]-xsj%*%thetaj)^2)/n
          qpcorj=sum((tau-funci(y-xsj%*%piej))*(x[,j]-xsj%*%thetaj))/(n*sqrt((tau-tau^2)*thegmaj))
          qpcor[which(b==j)]=qpcorj
        }
        A[i]=b[which.max(abs(qpcor))]}
    }
    ABIC=BIC(y,x,n,tau,A)
    return(list(ABIC,A))
}










