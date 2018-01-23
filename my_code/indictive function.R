##indicative function
funci=function(w){
  l=length(w)
  I=vector(length=l)
  for (u in 1:l) {
    if(w[u]<0){
      I[u]=1
    }
    else{
      I[u]=0
    }
  }
  return(I)
}