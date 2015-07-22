library("devtools")
install_github("strucchange","mengluchu",build_vignettes = TRUE) # how to build on modified package
install_github("bfast2","mengluchu",build_vignettes = TRUE)
library("bfast")
library('rgdal')
library('raster')
library("spacetime")
library("scidb")
library("raster")
#scidbconnect("server_name", port, "user_name", "password")
data(prodespoints00) # mask forest
data(pdd)            # validation

#data(fevi8) # array saved in Rdata

#output from SCIDB
sciref<-scidb("resarefpscidb2")
tssarar1s<-array(as.double(sciref[][]$spcu),c(150,150))[2:149,2:149] 
tssarar2s<-array(as.double(sciref[][]$spmo),c(150,150))[2:149,2:149] 
tssarar3s<-array(as.double(sciref[][]$cu),c(150,150))[2:149,2:149]
tssarar4s<-array(as.double(sciref[][]$mo),c(150,150))[2:149,2:149]
tssarar5s<-array(as.double(sciref[][]$arcu),c(150,150))[2:149,2:149]
tssarar6s<-array(as.double(sciref[][]$armo),c(150,150))[2:149,2:149]
aaa<-tssarar1ss[2:149,2:149]-tssarar1
tssarar1ss<-aperm(tssarar1s,c(2,1))
length(
str(  tssarar1s[which(tssarar1s)])
  )

summary(aaa)
#stored results

data(tssarar1) # corrected st model with original array sar cusum
data(tssarar2)#sar mosum
data(tssarar3)# cusum
data(tssarar4)# mosum
data(tssarar5)#ar cusum
data(tssarar6)#ar mosum

length(which(tssarar2s<0.05,arr.ind=TRUE))
#takes a week to run,checked  
#tsall<-SARefpdents(inputarray=fevi8, le=636)
tssarar1<-tsall[[1]]
tssarar2<-tsall[[2]]
tssarar3<-tsall[[3]]
tssarar4<-tsall[[4]]
tssarar5<-tsall[[5]]
tssarar6<-tsall[[6]]
 
groundtruth<-pdd    
ttssarar1<-generatecmpvalue(tssarar1,pdd ,pv=0.05) #sarcu
ttssarar11<-generatecmpvalue(result.array=tssarar1 ,reference.sppoints=pdd ,pv=0.2)
ttssarar2<-generatecmpvalue(tssarar2, pdd,pv=0.05) #sarmo
ttssarar22<-generatecmpvalue(tssarar2, pdd,pv=0.1) 
ttssarar3<-generatecmpvalue(tssarar3, pdd,pv=0.05) #cu
ttssarar4<-generatecmpvalue(tssarar4, pdd,pv=0.05) # mo
ttssarar5<-generatecmpvalue(tssarar5, pdd,pv=0.05) #arcu
ttssarar6<-generatecmpvalue(tssarar6, pdd,pv=0.05) #armo



cts<-rbind(ttssarar4,ttssarar6,ttssarar2,ttssarar22，ttssarar3，ttssarar5，ttssarar1，ttssarar11)
names1<-c("OLS-MOSUM p-value: 0.05","AR(1) OLS-MOSUM pvalue: 0.05","SAR OLS-MOSUM p-value: 0.05", "SAR OLS-MOSUM p-value: 0.1","OLS-CUSUM p-value: 0.05","AR(1) OLS-CUSUM p-value: 0.05","SAR OLS-CUSUM p-value: 0.05","SAR OLS-CUSUM p-value: 0.2")

barplotcm(cts, n=8,names1=names1)


ptssarar1<-generateppvalue(tssarar1s,pv=0.05)
ptssarar11<-generateppvalue(tssarar1s,pv=0.005)
ptssarar2<-generateppvalue(tssarar2s, pv=0.05)
ptssarar22<-generateppvalue(tssarar2,pv=0.025)
ptssarar3<-generateppvalue(tssarar3,pv=0.05)
ptssarar4<-generateppvalue(tssarar2,pv=0.05)
ptssarar5<-generateppvalue(tssarar5s,pv=0.05)
ptssarar6<-generateppvalue(tssarar6s,pv=0.05)
generatemapRAS1(ptssarar1,groundtruth,prodespoints00)
