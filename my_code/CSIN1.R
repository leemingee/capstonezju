CSIN1<-function(n,p,ro,re)
{
  #sigma
  sigma<-matrix(0,p,p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma[i,j]=ro^abs(i-j)
    }
  }
  d1=floor(n/log(n)) #选出predictors的最终维数
  am=matrix(nrow=re,ncol=d1)
  place=matrix(nrow=re,ncol=3)
  
  for (g in 1:re) {
    x=mvrnorm(n, rep(0,p), sigma)
    e=rnorm(n)
    y=y=2*x[,1]+3*x[,1]*x[,2]+3*x[,1]*x[,3]+e #3*x[,1]*x[,2]+3*x[,1]*x[,3]+e   
    CORM=cor(y,x)
    a=order(abs(CORM),decreasing=TRUE) #相关系数由大到小排序
    am[g,1:d1]=a[1:d1]
    place[g,1:3]=c(which(a==1),which(a==2),which(a==3))
  }
  rankave=apply(place, 2, mean)
  return(list(rankave,place,am))
}



