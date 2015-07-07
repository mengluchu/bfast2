
load("t3darrbfamul3.Rdata") 
load("t3darrbfamul2.Rdata") 
library("raster")
library("rasterVis")
valueoffirstorlastbreak<- function(changearray,timeofbreak="earliest", x=c(58930:59079),y=c(48210:48359)) # map the changes from the array that stores the changes. multiple changes are mapped as one change point
{  
  
  change7<-which(!is.na(changearray ),arr.ind=TRUE) #0.05
  change<-which(!is.na(changearray ) ) #0.05
  
  #xct1<-change7[,1]+x[1]-1 # for the second 150 by 150 array
  xct2<-change7[,1] 
  
  #yct1<-change7[,2]+y[1]-1
  yct2<-change7[,2] 
  
  dfallxyt<-as.data.frame(cbind(xct2,yct2))
  names(dfallxyt)<-c('x','y')
  coordinates(dfallxyt)<-~x+y #make the time value for searching
  bdt<-dfallxyt 
  
  zbdt<-zerodist(bdt) 
  
  change8<-change7[-zbdt[,1],]
  xct1<-change8[,1]+x[1]-1 # for the second 150 by 150 array   
  yct1<-change8[,2]+y[1]-1
  
  xychange<-cbind(xct1,yct1)
  changein<-getxyMatrix(xychange,231.6564)
  if (timeofbreak=="earliest")
  { changedata<-changearray[change][-zbdt[,2]] }
  else {changedata<-changearray[change][-zbdt[,1]]} #"latest"
  
  spmodist51<-SpatialPointsDataFrame(coordinates(changein),data=data.frame(changedata))
  proj4string(spmodist51)<-'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
  
  return(spmodist51)
  
} 

classofbl<-valueoffirstorlastbreak(t3darrbfamul2,timeofbreak='lastest') #bfast classes
classofbe<-valueoffirstorlastbreak(t3darrbfamul2)                       # times

timeofbl<-valueoffirstorlastbreak(t3darrbfamul3,timeofbreak='lastest')
timeofbe<-valueoffirstorlastbreak(t3darrbfamul3)


gridded(timeofbe)<-TRUE
rastert1<-raster(timeofbe)
levelplot(rastert1,margin = FALSE,col.regions =bpy.colors(),main="Time of Earliest Breakpoints")

gridded(timeofbl)<-TRUE
rastert2<-raster(timeofbl)
levelplot(rastert2,margin = FALSE,col.regions =bpy.colors(),main="Time of Latest Breakpoints")


gridded(classofbe)<-TRUE
rastert3<-raster(classofbe)
levelplot(rastert3,margin = FALSE,col.regions =bpy.colors(),main="Types of of Earliest Breakpoints")

gridded(classofbl)<-TRUE
rastert4<-raster(classofbl)
levelplot(rastert4,margin = FALSE,col.regions =bpy.colors(),main="Types of of Latest Breakpoints")

# plot together
rasterst2<-stack(rastert1,rastert2)
timeob<-levelplot(rasterst2,col.regions =bpy.colors(),names.attr=c("Time of Earliest Breakpoints","Time of Latest Breakpoints"))

classes<-c("NA","NA","monotonic decrease","monotonic increase", "setback (interupted gradual change)","boost(interupted gradual change)","greening (reversal)","browning (reversal)")
rastert3r = ratify(rastert3)
rastert4r = ratify(rastert4)

rat <- levels(rastert3r)[[1]]
rat$class <-classes[rat$ID]
levels(rastert3)<-rat
levels(rastert4)<-rat
rasterst<-stack(rastert3,rastert4)
ty<-levelplot(rasterst,col.regions =bpy.colors(),layout=c(1, 2),names.attr=c("Types of Earliest Breakpoints","Types of Latest Breakpoints"))

jpeg('typesofbreakpoints.jpg', height=12, width=12, res=400,unit="in")
plot(ty)
dev.off()
jpeg('timeofbreakpoints.jpg', height=12, width=12, res=400,unit="in")
plot(timeob)
dev.off()
plot(rasterst)
