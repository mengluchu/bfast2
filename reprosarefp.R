install_github("strucchange","mengluchu",build_vignettes = TRUE) # how to build on modified package
install_github("bfast2","mengluchu",build_vignettes = TRUE)
library("bfast")
library('rgdal')
library('raster')
library("spacetime")
data(prodespoints00) # mask forest
data(pdd)            # validation

#data(fevi8) # array saved in Rdata

#output from SCIDB
scidbarray<-resarefpscidb 
sciref<-scidb("scidbarray")

tssarar1<-array(sciref[,,][]$cu,c(150,150))
tssarar2<-array(sciref[,,][]$mo,c(150,150))
tssarar3<-array(sciref[,,][]$spcu,c(150,150))
tssarar4<-array(sciref[,,][]$spmo,c(150,150))
tssarar5<-array(sciref[,,][]$arcu,c(150,150))
tssarar6<-array(sciref[,,][]$armo,c(150,150))

#stored results
#data(tssarar1) # corrected st model with original array cusum
#data(tssarar2)#mosum
#data(tssarar3)#sar cusum
#data(tssarar4)#sar mosum
#data(tssarar5)#ar cusum
#data(tssarar6)#ar mosum


#takes a week to run
#tsall<-SARefpdents(inputarray=fevi8, le=636)
#tssarar1<-tsall[[1]]
#tssarar2<-tsall[[2]]
#tssarar3<-tsall[[3]]
#tssarar4<-tsall[[4]]
#tssarar5<-tsall[[5]]
#tssarar6<-tsall[[6]]

groundtruth<-pdd    
ttssarar1<-generatecmpvalue(tssarar1,pdd,pv=0.05)
ttssarar11<-generatecmpvalue(result.array=tssarar1,reference.sppoints=pdd ,pv=0.2)
ttssarar2<-generatecmpvalue(tssarar2,pdd,pv=0.05)
ttssarar22<-generatecmpvalue(tssarar2,pdd,pv=0.1)
ttssarar3<-generatecmpvalue(tssarar3,pdd,pv=0.05)
ttssarar4<-generatecmpvalue(tssarar4,pdd,pv=0.05)
ttssarar5<-generatecmpvalue(tssarar5,pdd,pv=0.05)
ttssarar6<-generatecmpvalue(tssarar6,pdd,pv=0.05)

cts<-rbind(ttssarar4,ttssarar6,ttssarar2,ttssarar22，ttssarar3，ttssarar5，ttssarar1，ttssarar11)
names1<-c("OLS-MOSUM p-value: 0.05","AR(1) OLS-MOSUM pvalue: 0.05","SAR OLS-MOSUM p-value: 0.05", "SAR OLS-MOSUM p-value: 0.1","OLS-CUSUM p-value: 0.05","AR(1) OLS-CUSUM p-value: 0.05","SAR OLS-CUSUM p-value: 0.05","SAR OLS-CUSUM p-value: 0.2")

barplotcm(cts, n=8,names1=names1)
