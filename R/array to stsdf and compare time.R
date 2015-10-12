#1) STFDF:  covert array to stsdf: 
#2) STSDF:  array store the change point
#           deter spatial polygon list; prodes spatial polygon list
#3) way 1:  manually build spatial and temporal index
#4) way 2:  convert from array to STFDF, then to STIDF, then to STSDF (from STIDF to STSDF the endtime is by default)
#5) polygontopoint:           rasterize polygons to spatialpoint

# arraytoSTFDF2 :           1)
# deterpolytopointstoSTSDF  2) 3)
# changearraytoSTSDF:       2) 3)
# waystfdf2stsdf:           4)

#hexspt comparetime plottimediff

arraytoSTFDF2<-function(array,attr.name="value", alltime=bt, x=c(58930:59079),y=c(48210:48359),crs,months=0.3)
{
  
  t<-c(1:dim(array)[3])
  tt2<-unique(as.POSIXct(alltime[t]))  # array index to time
  y1<-rep(y,each=length(x))
  x1<-rep(x,length(y))
  xyd<-as.data.frame(cbind(x1,y1))
  xyc<-getxyMatrix(xyd,231.6564)
  xyc<-as.data.frame(xyc)
  
  names(xyc)<-c("x","y")
  coordinates(xyc)<-~x+y
  proj4string(xyc)<-crs
  
  data2<-as.vector(array) 
  
  data2<-  data.frame(data2)
  if(months!=0){
    ltime<-  tt2+months*3600*24*30
    etime<-  tt2-months*3600*24*30
    
    stfdf1<-STFDF(xyc, etime, data2,ltime) }else{ 
      stfdf1<-STFDF(xyc, tt2, data2 )  }
  return(stfdf1)
}

polygontopoint<-function(x=c(58930:59079),y=c(48210:48359),sppolygon,crs=MODISCRS) {
  #sppolygon<-prodes67808
  
  x1<-rep(x,each=length(y))
  y1<-rep(y,length(x))
  
  xyd<-as.data.frame(cbind(x1,y1))
  xyc<-getxyMatrix(xyd,231.6564)
  xyc<-as.data.frame(xyc)
  coordinates(xyc)<-~x+y
  SpatialPoints(xyc)
  #as.character(round(unique(xyc@coords[,2])))
  
  proj4string(xyc)<-crs
  
  proj4string(sppolygon)<-CRS("+proj=utm +zone=21 +south  +ellps=WGS84 ")
  sppolygon=spTransform(sppolygon,  CRS(crs) )
  points<-xyc[sppolygon,]
  
  return( points)
}
#p1<-polygontopoint(c(58930:59079),c(48210:48359),spatialPolygons)

deterpolytopointstoSTSDF<-function(tdall2=tdall,pp01=pp0,rpp2=rpp, crs=MODISCRS,attr.name="value",months=0.3)
{
  
  detertime<-tdall2
  dtime<-as.POSIXct(detertime )
  
  uniquetime<-unique(dtime)
  rp1<-table(dtime)
  tt2<-uniquetime[order(uniquetime)]
  
  lt<-length(uniquetime)
  
  tn<-unlist(lapply(lt:1, function(i) rep(i,rp1[i])))
  
  tn2<- unlist(lapply(1:length(tn), function(i) rep(tn[i],rpp2[[i]])))  #rpp2: the number of points
  
  sn2<-c(1:length(pp01))
  index1<-as.matrix(cbind( sn2 , tn2  ))
  
  data2<-data.frame(rep(1,length(pp01)))
  names( data2)<-attr.name
  proj4string(pp01)<-CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs')
  
  ltime<- tt2+months*3600*24*30 
  etime<- tt2-months*3600*24*30
  stsdf2<-STSDF(pp01,etime,index=index1,data2,ltime ) 
  
  #stpolygon<-spTransform(stsdf2, CRS(crs))
  
  return(stsdf2)
}

changearraytoSTSDF<-function(array,crs,attr.name="value", alltime=bt, x=c(58930:59079),y=c(48210:48359),months=6)
{
  
 
  #woud be easier to get the unordered index of time than get the right index of space
  change7<-which(!is.na(array),arr.ind=TRUE) #0.05
  change8<-which(!is.na(array))
  x11<-change7[,1]+x[1]-1 # for the second 150 by 150 array   
  y11<-change7[,2]+y[1]-1
  t<-change7[,3]
  
  xyd<-data.frame(cbind(x11,y11))
  xyc<-getxyMatrix(xyd,231.6564)
  xyc<-data.frame(xyc)
  names(xyc)<-c('x','y')
  xyc$ID<-c(1:length(x11))
  coordinates(xyc)<-~x+y 
  
  tt2<-unique(as.POSIXct(alltime[t]))  # array index to time
  
  zdxyc<-zerodist(xyc)     # non-unique spatial points   
  
  sel = lapply(2:length(zdxyc[,2]),  function(i) !identical( zdxyc[i-1,2], zdxyc[i,2]))
  sel<- unlist(sel)
  sel2<-c(TRUE,sel)
  sel2<-sel2[-length(sel2)]
  zuni = zdxyc[sel2,]
  length(zdxyc[,2])
  
  spl<-length(xyc)        
  
  uniquexyc<- remove.duplicates(xyc)
  dunixyc<-data.frame(uniquexyc)
  dunixyc$IDnew<-c(1:length(uniquexyc))
  
  noID<-data.frame(cbind(dunixyc$ID,dunixyc$IDnew))
  
  remo<- match(zuni[,1],dunixyc$ID)  # return the index of ID
  
  oldindex<-c(dunixyc$ID,zuni[,2])
  newindex<-c(dunixyc$IDnew,noID[remo,2])
  
  oldandnewindex<-data.frame(cbind(newindex,oldindex))
  names(oldandnewindex)<-c("new","old")
  oldnew<-oldandnewindex[order(oldindex),] 
  lt<-length(tt2)
  #spl<-length(newone2)
  #runi<-rank(unique())
  tn<-lapply(1:lt, function(i) rep(i,table(t)[i]))
  #sn<-lapply(1:spl, function(i) rep(runi[i],newone2[i]))
  index1<-as.matrix(cbind(oldnew$new,unlist(tn) ))
  
  data2<-na.omit(as.vector(array)) 
  names( data2)<-attr.name
  data2<-as.data.frame(data2)
  
  
  ltime<-  tt2+months*3600*24*30
  etime<-  tt2-months*3600*24*30
  spxyc<-as( uniquexyc,"SpatialPoints")
  proj4string( spxyc)<-crs
  stsdf1<-STSDF(spxyc, etime,index=index1,data2,ltime) 
  return(stsdf1)
} 

#not so computationally efficient, but easiest way, pay attention to the x and y sequence
waystfdf2stsdf<-function(data=t3darrbfamul2, alltime=bt,x=c(58930:59079),y=c(48210:48359),crs=MODISCRS,months=0.3)  
{
  changestfdf<-arraytoSTFDF2(array=data, alltime=alltime,x=x,y=y,crs=crs,months=months)   #change array, points
  #which(!is.na(edivisive1[,,1:160]),arr.ind=TRUE)
  stsdfchangepoints1<-as(changestfdf,"STIDF")
  stsdfchangepoints2<-as(stsdfchangepoints1,"STSDF")
  return(stsdfchangepoints2)
}
 
hexspt<-function(stdf,name)
{
  
  hedb<-hexbin(data.frame(stdf)$x,data.frame(stdf)$y)
  
  plot(hedb,main=name,xlab='x',ylab='y')
 
}

comparetime<-function(sts1=changests,sts2=deterpoinf  )
{
  x<-geometry(sts1) #bfast
  y<-geometry(sts2) #2125 deter points
  
  xspin<-na.omit(over(y@sp,x@sp))  
  bfastindeter<-x[xspin,]
  
  yspin<-na.omit(over(x@sp,y@sp)) # deter in bfast  (more replicated points?) 1370 # 873?
  deterinbfast<-y[yspin,]
 
  deterinbfastspid<-over(bfastindeter@sp, deterinbfast@sp)
  bfastindeterspid<-over(deterinbfast@sp,bfastindeter@sp)
  
  
  lt1<-c()
  lt2<-c()
  ini<-array(c(0,0),c(1,2))
  for (i in 1:length(deterinbfastspid))
  {
    lt1[i]<-length(time(deterinbfast[deterinbfastspid[i], ] ))
    lt2[i]<-length(time(bfastindeter[i, ] ))
    ini<-  rbind(ini, cbind(as.integer( time(deterinbfast[deterinbfastspid[i], ] ))
                            ,as.integer( time(bfastindeter[i, ]))))
  }
  ini<-ini[-1,]
  return(ini)
}


plottimediff<-function(sts1=changests,sts2=deterpoinf , timedf=ini, xlab="BFAST",ylab="DETER",title, bufferdays1=150,bufferdays2=1 )
{ 
  x<-geometry(sts1) #bfast
  y<-geometry(sts2) #2125 deter points
  
  xspin<-na.omit(over(y@sp,x@sp)) #bfast in deter 1040 y@sp[x@sp,] more than 873 because deter points are replicated
  
  #xspin1<-x@index[,1][xspin] # index in x
  #length(unique(xspin)) is the same as the length of yspin 
  bfastindeter<-x[xspin,]
  
  yspin<-na.omit(over(x@sp,y@sp)) # deter in bfast  (more replicated points?) 1370 # 873?
  deterinbfast<-y[yspin,]
  #yspin1<-y@index[,1][yspin]
  #plot(x@sp[y@sp,],col='skyblue') # mind they are not the same since there are duplicated spatial points # 1370
  #points(y@sp[x@sp,],col='pink') # 1040
  #table(over(y@sp[x@sp,], x@sp[y@sp,]))
  
  deterinbfastspid<-over(bfastindeter@sp, deterinbfast@sp)
  bfastindeterspid<-over(deterinbfast@sp,bfastindeter@sp)
  
  
  lt1<-c()
  lt2<-c()
  ini<-array(c(0,0),c(1,2))
  for (i in 1:length(deterinbfastspid))
  {
    lt1[i]<-length(time(deterinbfast[deterinbfastspid[i], ] ))
    lt2[i]<-length(time(bfastindeter[i, ] ))
    ini<-  rbind(ini, cbind(as.integer( time(deterinbfast[deterinbfastspid[i], ] ))
                            ,as.integer( time(bfastindeter[i, ]))))
  }
  ini<-ini[-1,]
  t1<-as.POSIXct(ini[,1],origin='1970-01-01') #deter 
  t2<-as.POSIXct(ini[,2],origin='1970-01-01') #bfast        
  t22<-t2-3600*24*bufferdays1
  t23<-t2+3600*24*bufferdays1
  t111<-as.integer(t1)
  t222<-as.integer(t22)
  t223<-as.integer(t23)
  hedb<-hexbin(cbind(ini[,1],ini[,2]))
  hedb2<-hexbin( t111,t222)
  hedb3<-hexbin( t111,t223)
  hvp <- hexViewport(hedb)
  hvp2 <- hexViewport(hedb2)
  hvp3 <- hexViewport(hedb3)
  maxini<-max(ini)
  minini<-min(ini)
  
  #jpeg(paste( i,"bfasttime vs dtertime.jpg "), height=4, width=7, res=400,unit="in")
  
  a<- hexbinplot(ini[,1]~ini[,2],xlab=xlab, ylab=ylab,aspect = 1
                 ,xlim=c(minini,maxini+30000000) ,style='nested.centroids',
                 ylim =c(minini,maxini+30000000),main=title,
  )
  
  # hexVP.abline(hvp,a=0,b=1,col='orange')
  #  hexVP.abline(hvp2,a=0,b=1,col='red')
  # hexVP.abline(hvp3,a=0,b=1,col='skyblue')
  #legend("bottomright",c("bfast a year earlier", 
  #                 "bfast same as deter","bfast a year later"),
  #       col=c("red","orange","skyblue"),pch=1)
  
  #dev.off()
  
  
  hexspt(bfastindeter, name="bfast in deter")
  hexspt(deterinbfast, name="deter in bfast")
  ##return the number of points that are overlapped
  tend1<-t1+3600*24*bufferdays1 # half year late and early buffer
  tb1<-t1-3600*24*bufferdays1
  tend2<-t2+3600*24*bufferdays2
  tb2<-t2-3600*24*bufferdays2
  interval1<-Intervals(cbind(tb1,tend1), closed = c(TRUE, FALSE))
  interval2<-Intervals(cbind(tb2,tend2), closed = c(TRUE, FALSE))
  # 
  rett<-c()
  for (i in 1:length(t1))
  {  
    ret = interval_overlap(interval2[  i] ,interval1[i] )
    if(length(unlist(ret))!=0)
      rett[i]<-ret[[1]]
  }  
  lenofover<-length(which(!is.na(rett))) # 76/ (1370) match for half year buffer
  return(list(a,lenofover))
  
}
