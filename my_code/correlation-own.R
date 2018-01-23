correlation<-function(y,x,n,p,alpha)
{
  i=diag(n) 
  rou=h=matrix(NA,p,p)
  bd=floor(cp*sqrt(n/log(n)))
  sj=matrix(NA,p,ncol=bd)
  for (t in 1:p) {
    rou[t,]=cor(x[,t],x[,]) 
    h[t,]=order(abs(rou[t,]),decreasing=TRUE) #相关系数由大到小排序
    h=h[,-1]
    mn=min(col(h),bd)
    for (j in 1:p) {
      for (m in 1:ncol(h)) {
        d=h[j,1:m] #选出的cor前m个的变量集
        f=setdiff(h[j,],d) #因cor未被选出的变量集
        jyo=vector(length=0)
        xs=x[,d]
        q=i-xs%*%solve(t(xs)%*%xs)%*%t(xs)
        for (k in f) {
          cjk=sqrt(sum((q%*%x[,k])^2))*sqrt(sum((q%*%x[,j])^2))/n
          prou=t(x[,j])%*%q%*%x[,k]/(cjk*n) #变量j与变量k的偏相关系数
          fjk=0.5*log((1+prou)/(1-prou))
          jy=(n-length(d)-3)^0.5*abs(fjk) #检验统计量
          jyo[which(f==k)]=jy #算出每一个jk检验，并排序
        }
        jyo=na.omit(jyo)
        if(max(jyo)<qnorm(1-alpha/2)){
          break
        }
      }
      mj=m
      cp=1 #参数
      if(mj<cp*sqrt(n/log(n))){
        mj=mj
      }else{
        mj=floor(cp*sqrt(n/log(n)))
      }
      sj[j,1:mj]=h[j,1:mj]#计算qpcor时的条件集合
    }
    return(sj)
  }
  
}
  