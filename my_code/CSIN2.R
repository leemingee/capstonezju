CSIN2<-function(n,p,ro,re)
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
    xi=mvrnorm(n, rep(0,p), sigma)
    x=xi^2
    e=rnorm(n)
    yi=3*xi[,1]*xi[,2]+3*xi[,1]*xi[,3]+e   #yi=2*xi[,1]+3*xi[,1]*xi[,2]+3*xi[,1]*xi[,3]+e
    y=yi^2
    CORM=cor(y,x)
    h=order(abs(CORM),decreasing=TRUE) #相关系数由大到小排序
    }
    am[g,1:d1]=h[1:d1]
  #place[g,1:3]=c(which(a==1),which(a==2),which(a==3))
  #rankave=apply(place, 2, mean)
  return(list(am))
}



