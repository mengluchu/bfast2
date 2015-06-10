install_github("bfast2","mengluchu")

install_github("strucchange","mengluchu",build_vignettes = FALSE)
library(devtools)
library(strucchange)
library(bfast)


setwd("C:/Users/m_lu0002/Dropbox/mengluchu/bfast")
list.files()
load("a1.saved")
load("fevi8.Rdata")

fevi20b20<-fevi8[20:50,20:50,] # input example array
fevi150b150<-fevi8[1:150,1:150,]
fevi9b9<-fevi8[1:9,1:9,]
fevi3b3<-fevi8[1:3,1:3,]
#fevi3b4<-fevi8[1:3,1:4,]

#output1<-bfmarray(fevi20b20,dates=a1,aggre='month',start=9) # bfast monitor
########## try spatial CAR bfast

#get STF

array<-fevi3b3
itrydf<-as.data.frame.table(array) #y x t interate y ->x ->t: (x1,y1,t1)(x1,y2,t1)(x2,y1,t1)(2,2,t1)(1,1,t2)(1,2,t2)

# deseasonality
fevi3b312t<-apply(fevi3b3,c(1,2),function(x) (ts(x,start=c(2000,1),frequency=46)-stl(ts(x,start=c(2000,1),frequency=46),'per')$time.series[,"trend"]))
dim(fevi3b312t)
f2<-aperm(fevi3b312t,c(2,3,1))
#
itrydf<-as.data.frame.table(f2)
aa2<-itrydf$Freq

#
aa3<-as.data.frame(aa2)
eday <- as.Date("2000-01-30")           # date 
e8day <- seq(eday, length.out=636, by="8 days")

x<-c(1:3)
y<-c(1:3)
x1<-rep(x,length(y))
y1<-rep(y,each=length(x))
xyd<-as.data.frame(cbind(x1,y1))
coordinates(xyd)<-~x1+y1

stfdf3b3<-STFDF(xyd,e8day, aa3) 
length()
# get neighbor 

cn<-cell2nb(3,3, type ="queen",torus =FALSE)
neigh1<-nbMult(cn, stfdf3b3, addT = FALSE, addST = FALSE) # only spatial neighbours are added for eath time replicate

listcn636<-nb2listw(neigh1)
fevi3b312<-ts(fevi3b3[i,j,],start=c(2000,1),frequency=46)
try2<-spautolm(aa2~c(rep(c(1:636),each=9)),family="SAR",listw=listcn636)
#try1<-spautolm(aa2~c(rep(c(1:636),each=12)),family="CAR",listw=listcn636)
rn<-lapply(1:9,function(i) {residuals(try2)[seq(i,636*9-(9-i),9)]})


p1<-array(,c(3,3))
p2<-array(,c(3,3))
p3<-array(,c(3,3))
dim(f2) #deseasonalized 
for (j in 1:3) 
{
  for(i in 1:3)  
  { 
    ii<-i+(j-1)*3
    
    fevi3b312t<-ts(f2[i,j,],start=c(2000,1),frequency=46)
    p.Vt2 <- sctest(efp(fevi3b312t ~ ti, h = 0.15, type = "OLS-CUSUM" ))
    
    p.Vt2$p.value
    p.Vt <- sctest(efp(fevi3b312t ~ ti, h = 0.15, type = "OLS-CUSUM", spatial1=as.numeric(rn[[ii]]))  )
    p.Vt$p.value
    p1[i,j]<- p.Vt$p.value # spautolm residuals
    p2[i,j]<-p.Vt2$p.value # linear regression residuals
    
    #p3[i,j]<-bfast(fevi3b312t,max.iter=1,level=0.05,h=0.15, type = "OLS-CUSUM")$output[[1]]$Vt.bp[1]
  }
}

}

p1
#p2
#p1: two more breakpoints: i=1:2 j=1
i=2; j=2
ii<-i+(j-1)*2
fevi3b312<-ts(f2[i,j,],start=c(2000,1),frequency=46)
plot(stl(fevi3b312,"per"))
plot(efp(fevi3b312 ~ ti, h = 0.15, type = "OLS-CUSUM", spatial1=as.numeric(rn[[ii]])))
plot(efp(fevi3b312 ~ ti, h = 0.15, type = "OLS-CUSUM"))
# rec-cusum cusum mosum most breakpoints ->least breakpoints
#plot(bfast(fevi3b312,max.iter=1,level=0.05,h=0.15))
#p1 spatial
#p2 original
#tsa1<-ts(aa2[(seq(2,636*9-7,9))],start=c(2000,8),frequency=46,end=c(2013,45))
#remainder<- stl(tsa1,"period")$time.series[,"remainder"]
#plot(residuals(lm(remainder~c(rep(1:636,1)))),typ="l") # seasonality is still here ?
#plot(aa2[(seq(1,636*9-8,9))],typ="l")
#summary(residuals(tryll))
#plot(aa2[(seq(1,636*9-8,9))],typ="l")

#a<-seq(1, 636*9-8,9)
#b<-seq(2, 636*9-7,9)
#c<-seq(3, 636*9-6,9)
#d<-seq(4, 636*9-5,9)

#plot(fitted(try1)[a],typ="l")
#plot(fitted(try1)[b],typ="l")
#plot(fitted(try1)[c],typ="l")
#plot(fitted(try1)[d],typ="l")


#r1<-residuals(try1)[a]
#r2<-residuals(try1)[b]
#r3<-residuals(try1)[c]
#r4<-residuals(try1)[d]
#r5<-residuals(try1)[seq(5,636*9-4,9)]
#r6<-residuals(try1)[seq(6,636*9-3,9)]
#r7<-residuals(try1)[seq(7,636*9-2,9)]
#r8<-residuals(try1)[seq(8,636*9-1,9)]
#r9<-residuals(try1)[seq(9,636*9,9)]
#rr<-list(r1,r2,r3,r4,r5,r6,r7,r8,r9)



#fevi3b312<-ts(fevi3b3[1,2,],start=c(2000,1),frequency=46)

#p3<-bfast(fevi3b312,max.iter=1,level=0.05,h=0.15)$output[[1]]$bp.Vt

#edit(bfast)
#fevi3b3<-fevi8[40:42,40:42,]
#ti<-1:636

 


str(pvalue.efp(fevi3b312 ~ ti, h = 0.15,lim.process="Brownian motion",alt.boundary=FALSE)
    
    
    help(sctest)
    edit(efp)
    efp()
    ## convert time to time dimension
    
    mtime<-monthofyears(output1[1,,]) #timearr
    
    dimx<-dim(output1)[2]
    dimy<-dim(output1)[3]
    dimt<-max(mtime[!is.na(mtime)])
    
    newarray<-array(,c(dimx,dimy,dimt)) #newarr
    
    mag<- output1[2,,]   # variablearr (magnitude)
    
    #load("C:/Users/m_lu0002/Dropbox/mengluchu/bfast/app1.R")
    
    t3darr<-timeasdim(newarray,mtime,mag) # 3 dimensional array with magnitude 
    
    which(!is.na(t3darr),arr.ind=TRUE)
    image(t3darr[,,176]) #test: bfastmonitor magnitude
    
    #test bfast: save only time
    t1<-fevi20b20
    efp()
    dimt=23*12
    dimx<-dim(output3)[1]
    dimy<-dim(output3)[2]
    
    newarray<-array(,c(dimx,dimy,dimt)) #newarr
    
    mag<- array(100,c(dimx,dimy))
    
    t3dbfaarr<-timeasdimbfa(newarray,output3,mag)
    which(!is.na(t3dbfaarr),arr.ind=TRUE)
    image(t3dbfaarr[,,126])
    ########################################### bfast with save other values
    p<-proc.time()
    
    dim(t1)
    t1<-fevi20b20[1:20,1:20,]
    fevi8
    p1<-proc.time()
    output1 <-bfaarray2(fevi150b150,dates=a1,aggre="month",season="harmonic",max.iter=1,level=0.05)
    tp1<-proc.time()-p
    
    save(tp,file='tp.Rdata')
    output2 <-bfaarray5(t1,dates=a1,aggre='month',season="harmonic",max.iter=1,level=0.05)
    outputsarmosum <-bfaarraypvalue(fevi8,dates=a1,aggre='month',season="harmonic",max.iter=1,level=0.05,pvat=tssarar4,ijk=1,pvas=1)
    outputmosum <-bfaarraypvalue(fevi8,dates=a1,aggre='month',season="harmonic",max.iter=1,level=0.05,pvat=tssarar2,ijk=1,pvas=1)
    output1 <-bfaarray2(fevi8,dates=a1,aggre="month",season="harmonic",max.iter=1,level=0.05)
    
    #str(output2)
    allbreakpoints<-output1[1:6,,] # breakpoint
    allmagnitudes<-output1[7:12,,] # magnitude
    breakpointclass<-output1[13:18,,] # class
    
    #which(allbreakpoints!=0)
    mtime<-allbreakpoints #timearr
    
    dimx<-dim(output2)[2]
    dimy<-dim(output2)[3]
    dimt<-max(mtime[ mtime!=0])
    
    newarray<-array(,c(dimx,dimy,dimt)) #newarr
    
    mag<- allmagnitudes  # variablearr (magnitude)
    
    t3darrbfamul<-tasd.bfa.mul(newarray,mtime,mag) # 3 dimensional array with magnitude 
    
    t3darrbfamul2<-tasd.bfa.mul(newarray,mtime, breakpointclass) # 3 dimensional array with classes 
    t3darrbfamul3<-tasd.bfa.mul(newarray,mtime, allbreakpoints) # 3 dimensional array with time 
    
    proc.time()-p
    #t3darrbfamul[which(is.na(t3darrbfamul))]<-0
    
    which(!is.na(t3darrbfamul2),arr.ind=TRUE)
    which(!is.na(t3darrbfamul2),arr.ind=TRUE)
    plot(raster(na.omit(t3darrbfamul2[,,126]))) #test: bfastmonitor magnitude
   
    
    
    raster126<-raster(t3darrbfamul2[,,]) 
    #plot classies for each time step
    library(rasterVis)
    library(raster) 

    raster126<-raster(t3darrbfamul2[,,126]) 
    classes<-c("NA","NA","monotonic decrease","monotonic increase", "setback (interupted gradual change)","boost(interupted gradual change)","greening (reversal)","browning (reversal)")
    raster126r = ratify(raster126)
    rat <- levels(raster126r)[[1]]
    
    rat$class <-classes[rat$ID]
    levels(raster126)<-rat
    cols = c("#ffff99", "#7fc97f", "#1f78b4", "#fdc086","#e0f3db"
             ,"#fde0dd","#fa9fb5","#c51b8a", "#a8ddb5", "#43a2ca")
    
    l<- levelplot(raster126,col.regions=cols, pretty=TRUE, margin=F)
    plot(l)    
)
    
    
    levelplot(raster126)
    
    a<-bfastchangepoint(t3darrbfamul2)
    
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
classofbl<-valueoffirstorlastbreak(t3darrbfamul2,timeofbreak='lastest')
classofbe<-valueoffirstorlastbreak(t3darrbfamul2)
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



save(t3darrbfamul,file="t3darrbfamul.Rdata")
   which.max( table(t))
    t3darrbfamul[1,17,140]
    