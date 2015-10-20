# validation tools: confusion matrix, hexbin plot, comparison in space and time
# using space time object, raster, array

# compare in space, get a list of SptialPoints: True Positive, fn, tn, fn
#input: spatial points (with the same reference system) : change points, reference points, totalpoints
require(spacetime)
require(raster)
tfpn<-function(changepoint, reference.sppoint=pdd,totalp=19167){
  #changepoint<-pvpoint
  
  changepointf<-changepoint[prodespoints00]
  deterpoints<-reference.sppoint
  result<-data.frame()
  deterpoints2<-c()
  TP<-changepointf[reference.sppoint,]
  FN<-reference.sppoint[is.na(over(reference.sppoint,changepointf))]
  FP<-changepoint[is.na(over(changepointf,reference.sppoint))]
  PRODESNE<-prodespoints00[ which(is.na(over(prodespoints00,reference.sppoint)))]
  BFASTNE<-prodespoints00[ which(is.na(over(prodespoints00,changepoint)))]
  TN<-PRODESNE[BFASTNE]
  return(list(TP,FN,FP,TN))
}
#example: TP of MOSUM and SAR MOSUM 
#cmsarm<-tfpn(psarmosum,reference.sppoint = pmosum) #TP, FN, FP, TN
#TPts1<-stfdfevi8c[cmsarm[[1]],]# sar mosum and mosum agreed TP
# example to run:
#load("/Users/lumeng/Dropbox/mengluchu/bfast/tssarar2.Rdata")
#load("/Users/lumeng/Dropbox/mengluchu/bfast/pdd.Rdata")
#sptssarar2<-pvaluepoint2sp(tssarar2,x=c(58930:59079),y=c(48210:48359),xoff=1,yoff=1,pvalue=0.05,crs=CRS("+proj=utm +zone=21 +south")) #map changes from an array of p-values
#tfpn1<-tfpn(sptssarar2,pdd)
#TPsp<-tfpn1[[1]]
#FNsp<-tfpn1[[2]]
#FPsp<-tfpn1[[3]]
#TNsp<-tfpn1[[4]]

# tools for space-time objects

hexspt <- function(stdf, name) {
  
  hedb <- hexbin(data.frame(stdf)$x, data.frame(stdf)$y)
  
  plot(hedb, main = name, xlab = "x", ylab = "y")
  
}

comparetime <- function(sts1 = changests, sts2 = deterpoinf) {
  x <- geometry(sts1)  #bfast
  y <- geometry(sts2)  #2125 deter points
  
  xspin <- na.omit(over(y@sp, x@sp))
  bfastindeter <- x[xspin, ]
  
  yspin <- na.omit(over(x@sp, y@sp))  # deter in bfast  (more replicated points?) 1370 # 873?
  deterinbfast <- y[yspin, ]
  
  deterinbfastspid <- over(bfastindeter@sp, deterinbfast@sp)
  bfastindeterspid <- over(deterinbfast@sp, bfastindeter@sp)
  
  
  lt1 <- c()
  lt2 <- c()
  ini <- array(c(0, 0), c(1, 2))
  for (i in 1:length(deterinbfastspid)) {
    lt1[i] <- length(time(deterinbfast[deterinbfastspid[i], ]))
    lt2[i] <- length(time(bfastindeter[i, ]))
    ini <- rbind(ini, cbind(as.integer(time(deterinbfast[deterinbfastspid[i], ])), as.integer(time(bfastindeter[i, 
                                                                                                                ]))))
  }
  ini <- ini[-1, ]
  return(ini)
}


plottimediff <- function(sts1 = changests, sts2 = deterpoinf, timedf = ini,
                         xlab = "BFAST", ylab = "DETER", 
                         title, bufferdays1 = 150, bufferdays2 = 1) {
  x <- geometry(sts1)  #bfast
  y <- geometry(sts2)  #2125 deter points
  
  xspin <- na.omit(over(y@sp, x@sp))  #bfast in deter 1040 y@sp[x@sp,] more than 873 because deter points are replicated
  
  # xspin1<-x@index[,1][xspin] # index in x length(unique(xspin)) is the same as the length of yspin
  bfastindeter <- x[xspin, ]
  
  yspin <- na.omit(over(x@sp, y@sp))  # deter in bfast  (more replicated points?) 1370 # 873?
  deterinbfast <- y[yspin, ]
  # yspin1<-y@index[,1][yspin] plot(x@sp[y@sp,],col='skyblue') # mind they are not the same since there are
  # duplicated spatial points # 1370 points(y@sp[x@sp,],col='pink') # 1040 table(over(y@sp[x@sp,],
  # x@sp[y@sp,]))
  
  deterinbfastspid <- over(bfastindeter@sp, deterinbfast@sp)
  bfastindeterspid <- over(deterinbfast@sp, bfastindeter@sp)
  
  
  lt1 <- c()
  lt2 <- c()
  ini <- array(c(0, 0), c(1, 2))
  for (i in 1:length(deterinbfastspid)) {
    lt1[i] <- length(time(deterinbfast[deterinbfastspid[i], ]))
    lt2[i] <- length(time(bfastindeter[i, ]))
    ini <- rbind(ini, cbind(as.integer(time(deterinbfast[deterinbfastspid[i], ])), as.integer(time(bfastindeter[i, 
                                                                                                                ]))))
  }
  ini <- ini[-1, ]
  t1 <- as.POSIXct(ini[, 1], origin = "1970-01-01")  #deter 
  t2 <- as.POSIXct(ini[, 2], origin = "1970-01-01")  #bfast        
  t22 <- t2 - 3600 * 24 * bufferdays1
  t23 <- t2 + 3600 * 24 * bufferdays1
  t111 <- as.integer(t1)
  t222 <- as.integer(t22)
  t223 <- as.integer(t23)
  hedb <- hexbin(cbind(ini[, 1], ini[, 2]))
  hedb2 <- hexbin(t111, t222)
  hedb3 <- hexbin(t111, t223)
  hvp <- hexViewport(hedb)
  hvp2 <- hexViewport(hedb2)
  hvp3 <- hexViewport(hedb3)
  maxini <- max(ini)
  minini <- min(ini)
  
  # jpeg(paste( i,'bfasttime vs dtertime.jpg '), height=4, width=7, res=400,unit='in')
  
  a <- hexbinplot(ini[, 1] ~ ini[, 2], xlab = xlab, ylab = ylab, aspect = 1, xlim = c(minini, maxini + 3e+07), 
                  style = "nested.centroids", ylim = c(minini, maxini + 3e+07), main = title, )
  
  # hexVP.abline(hvp,a=0,b=1,col='orange') hexVP.abline(hvp2,a=0,b=1,col='red')
  # hexVP.abline(hvp3,a=0,b=1,col='skyblue') legend('bottomright',c('bfast a year earlier', 'bfast same as
  # deter','bfast a year later'), col=c('red','orange','skyblue'),pch=1)
  
  # dev.off()
  
  
  hexspt(bfastindeter, name = "bfast in deter")
  hexspt(deterinbfast, name = "deter in bfast")
  ## return the number of points that are overlapped
  tend1 <- t1 + 3600 * 24 * bufferdays1  # half year late and early buffer
  tb1 <- t1 - 3600 * 24 * bufferdays1
  tend2 <- t2 + 3600 * 24 * bufferdays2
  tb2 <- t2 - 3600 * 24 * bufferdays2
  interval1 <- Intervals(cbind(tb1, tend1), closed = c(TRUE, FALSE))
  interval2 <- Intervals(cbind(tb2, tend2), closed = c(TRUE, FALSE))
  # 
  rett <- c()
  for (i in 1:length(t1)) {
    ret = interval_overlap(interval2[i], interval1[i])
    if (length(unlist(ret)) != 0) 
      rett[i] <- ret[[1]]
  }
  lenofover <- length(which(!is.na(rett)))  # 76/ (1370) match for half year buffer
  return(list(a, lenofover))
  
} 

timedif<-function(sts1=changests,sts2 =deterpoinf) # time differences 
{
  ini<- comparetime(sts1=sts1,sts2=sts2)
  
  time1<-as.POSIXlt(ini[,1],origin="1970-01-01")
  
  time2<-as.POSIXlt(ini[,2],origin="1970-01-01")
  
  diffd<- time1- time2
  
  return(diffd)
}

#replot histogram of time by having frequency as percentage
hist.time<-function(object,title){
  
  v.hist<-hist(object,breaks=20, freq=TRUE,plot=FALSE) 
  v.hist$counts <- v.hist$counts/sum(v.hist$counts)
  aaa<-as.Date(strftime(as.POSIXlt(v.hist$breaks,origin="1970-01-01"), format="%Y-%m-%d"))
  v.hist$breaks<- aaa 
  plot(v.hist,xaxt="n",xlab="",ylab="frequency",main=title)
  axis.Date(1,at=aaa,labels=format(aaa,"%Y-%b "),las=2)
}

#example:
#get the time 
#timeroccu<-time(mochangepointcu@time[mochangepointcu@index[,2]] )
#timestablecu<-time(mochangepoint2cu@time[mochangepoint2cu@index[,2]] )
#plotting the time histogram 
#jpeg("time cusum roc stable.jpeg ",width = 480,height = 480)
#par(mfrow=c(1,2))
#par(mar=c(5,2,2,2))
#hist.time(timeroccu,title=" time of cusum roc detected change")
#hist.time(timestablecu,title=" time of cusum stable detected change")
#dev.off()
#
#time differences
#bfa.de<-timedif(sts1=changests,sts2=deterpoinf)
#edi.de<-timedif(sts1=stfdfedivi,sts2=deterpoinf)
#bfa.edi<-timedif(sts1=changests,sts2=stfdfedivi)
#bfa.de1<-bfa.de/365
#edi.de1<-edi.de/365
#bfa.edi1<-bfa.edi/3600/24/365
#hist(as.numeric(bfa.de1 ),,breaks=10,main="Differences in time: Bfast Vs.Deter",xlab="year")
#v.hist$counts <- v.hist$counts/sum(v.hist$counts)
#plot(v.hist,ylab="probability",main="Differences in time: Bfast Vs.DETER",xlab="year")

# from generate confusion matrix 
#
#cts<-rbind(ttssarar4,ttssarar6,ttssarar2,ttssarar22，ttssarar3，ttssarar5，ttssarar1，ttssarar11)
#names1<-c("OLS-MOSUM p-value: 0.05","AR(1) OLS-MOSUM pvalue: 0.05","SAR OLS-MOSUM p-value: 0.05", "SAR OLS-MOSUM p-value: 0.1","OLS-CUSUM p-value: 0.05","AR(1) OLS-CUSUM p-value: 0.05","SAR OLS-CUSUM p-value: 0.05","SAR OLS-CUSUM p-value: 0.2")
#barplotcm(cts, n=8,names1=names1)
#generatemapRAS1(ptssarar2,groundtruth,prodespoints00)
