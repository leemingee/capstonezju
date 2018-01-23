BIC<-function(y,x,n,dx,tau,dstar,d,Mstar,alpha,A)
{
  da=length(A)
  BIC=seq(0,0,length=da)
  #Cn=log(da)
  Cn=max(log(log(da)),1)
  #Cn=1
  for(i in 1:da){
    xs=x[,A[1:i]]
    xss=cbind(1,xs)
    qlm=rq(y~xs,tau)
    phhat=qlm$coef
    yres=y-xss%*%phhat
    yres=as.vector(yres)
    rho=(tau-(yres<0))*yres
    BIC[i]=log(mean(rho))+i*log(n)*Cn/(2*n)
    #BIC[i]=log(mean(rho))+i*(log(n)+2*log(d))/(2*n)
  }
  ii=which(BIC==min(BIC))
  MhatBIC=A[1:ii]
   return(list(MhatBIC))

}