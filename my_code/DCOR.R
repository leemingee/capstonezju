DCOR<-function(y,x,n,p)
{
  dcor=vector(length=p)
  #for A
  a=matrix(NA,n,n)
  for (j in 1:n) {
    for (k in 1:n) {
      a[j,k]=sqrt((y[j]-y[k])^2)
    }
  }
  armean=apply(a,1,mean)
  acmean=apply(a,2,mean)
  amean=mean(armean)
  A=matrix(NA,n,n)
  for (j in 1:n) {
    for (k in 1:n) {
      A[j,k]=a[j,k]-armean[j]-acmean[k]-amean
    }
  }
  dcov2yy=sum(A^2)/(n^2)
  #for B(p个)
  for (i in 1:p) {
    b=matrix(NA,n,n)
    for (j in 1:n) {
      for (k in 1:n) {
        b[j,k]=sqrt((x[j,i]-x[k,i])^2)
      }
    }
    brmean=apply(b,1,mean)
    bcmean=apply(b,2,mean)
    bmean=mean(brmean)
    B=matrix(NA,n,n)
    for (j in 1:n) {
      for (k in 1:n) {
        B[j,k]=b[j,k]-brmean[j]-bcmean[k]-bmean
      }
    }
    dcov2yi=sum(A*B)/(n^2)
    dcov2ii=sum(A^2)/(n^2)
    
  }
  
}





