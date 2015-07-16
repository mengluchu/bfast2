filter.st.median<-function(a) #if evi is more than 1 or less than -1, replace with median of its surrounding pixels
{
  
  index1<-which(abs(a)>1,arr.ind=TRUE)
  
  
  for (i1 in (1:dim(index1)[1]))
  {
    
    i<-as.integer(index1[i1,1]) 
    
    j<-as.integer(index1[i1,2]) 
    t<-as.integer(index1[i1,3])
    
    ih<-i+1
    il<-i-1
    
    jh<-j+1
    jl<-j-1
    
    th<-t+1
    tl<-t-1
    
    if(ih>max(index1[,1]))
      ih=max(index1[,1])
    if(il<min(index1[,1]))
      il=min(index1[,1])
    
    if(jh>max(index1[,2]))
      jh=max(index1[,2])
    if(jl<min(index1[,2]))
      jl=min(index1[,2])
    
    if(th>max(index1[,3]))
      th=max(index1[,3])
    if(tl<min(index1[,3]))
      tl=min(index1[,3])
    
    
    a[i,j,t]<-median(c(
      a[ih,j,t],a[il,j,t],a[i,jh,t],a[i,jh,t],
      a[ih,j,tl],a[il,j,tl],a[i,jh,tl],a[i,jh,tl],
      a[ih,j,th],a[il,j,th],a[i,jh,th],a[i,jh,th]))
    
  }
  return(a)   
}
