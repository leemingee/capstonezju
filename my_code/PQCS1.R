PQCS1<-function(y,x,n,dx,tau,dstar,d,Mstar,alpha,CORM)
{

 A=NULL
  indx=seq(1:dx)

   pj=seq(0,0,length=length(Mstar))
   step=0
   while(step<=(dstar-1))
   	{
   step=step+1

   if(length(A)==0){
   xr=x
   ind=indx
   }else{
   xr=x[,-A]
   ind=indx[-A]
   }
   dxr=ncol(xr)  
  
  QPCOR=seq(0,0,length=dxr)
  for(j in 1:dxr){
   xj=xr[,j]
   ins=CORM[ind[j],]
   ins=na.omit(ins)
   Sj=c(ins,A)
   Sj=unique(Sj)
   if(length(Sj)>0)
   {
   xs=x[,Sj]
   xss=cbind(1,xs)
   vartheta=ginv(t(xss)%*%xss)%*%t(xss)%*%xj
   z=xj-xss%*%vartheta
   qlm=rq(y~xs,tau)
   phhat=qlm$coef
   yres=y-xss%*%phhat
   yres=as.vector(yres)
   }else{
   yres=y-quantile(y, tau)
   z=xj-mean(xj)
   }
   sig2=mean(z^2)
   ph=(tau-(yres<0))*z
   ph=as.matrix(ph)
   ph=apply(ph,2,mean)
   QPCOR[j]=abs(ph/sqrt((tau-tau^2)*sig2))
   }
   jstar=ind[which(QPCOR==max(QPCOR))]
   A=c(A,jstar)
}

Astar=A
step=0
   while(step<=(d-1-dstar)||sum(pj)<length(Mstar))
   	{
   step=step+1

   xr=x[,-A]
   ind=indx[-A]
   dxr=ncol(xr)  
  
  QPCOR=seq(0,0,length=dxr)
  for(j in 1:dxr){
   xj=xr[,j]
   ins=CORM[ind[j],]
   ins=na.omit(ins)
   Sj=c(ins,Astar)
   Sj=unique(Sj)
   if(length(Sj)>0)
   {
   xs=x[,Sj]
   xss=cbind(1,xs)
   vartheta=ginv(t(xss)%*%xss)%*%t(xss)%*%xj
   z=xj-xss%*%vartheta
   qlm=rq(y~xs,tau)
   phhat=qlm$coef
   yres=y-xss%*%phhat
   yres=as.vector(yres)
   }else{
   yres=y-quantile(y, tau)
   z=xj-mean(xj)
   }
   sig2=mean(z^2)
   ph=(tau-(yres<0))*z
   ph=as.matrix(ph)
   ph=apply(ph,2,mean)
   QPCOR[j]=abs(ph/sqrt((tau-tau^2)*sig2))
   }
   jstar=ind[which(QPCOR==max(QPCOR))]
   A=c(A,jstar)
  pj=seq(0,0,length=length(Mstar))
      for(jj in 1:length(Mstar))
    {
   pj[jj]=sum(A==Mstar[jj])
   } 
}
     Mhat=A[1:d]
   Rj=seq(0,0,length=length(Mstar))  ##average rank of Xj
   pj=seq(0,0,length=length(Mstar))  ##proportion of Xj being selected
   for(jj in 1:length(Mstar)){
   Rj[jj]=which(A==Mstar[jj])
   pj[jj]=sum(Mhat==Mstar[jj])
   }
   M=max(Rj)  ##largest rank of true predictors
   pa=1-as.numeric(mean(pj)<1)  ##proportion of all true predictors being selected
   return(list(Mhat,Rj,M,pj,pa,A))

}
