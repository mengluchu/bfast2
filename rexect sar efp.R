iquery("
store(
  project(
       apply(      
       repart(
       apply(
       subarray(
       
       MOD09Q1_JUARA,58930,48210,6,59079,48359,643),
       evi2,2.5*((nir*0.0001-red*0.0001)/(nir*0.0001+2.4*red*0.0001+1.0))),
       <red:int16,nir:int16,quality:uint16,evi2:double>[col_id=0:149,1,1,row_id=0:149,1,1,time_id=0:637,638,0]),
       the_i,double(col_id),the_j,double(row_id),the_t,double(time_id)),
       the_i, the_j,the_t, evi2),
       repro1)")

# sar epf 
iquery(" r_exec(repro1, 'output_attrs=8',
       
       'expr= 
       
       library(strucchange)
       library(spdep)
       library(nlme)
       
       load(\"/home/menglu/X.Rdata\") # coefficient matrix 
       load(\"/home/menglu/listcn636.Rdata\") # neighbor
       
       
       dim1<-length(unique(the_i)) 
       dim2<-length(unique(the_j)) 
       dim3<-length(unique(the_t))
       tl=1:636
       w=1/46       
       co <- cos(2*pi*tl*w)
       si <- sin(2*pi*tl*w)
       co2 <- cos(2*pi*tl*w*2)
       si2 <- sin(2*pi*tl*w*2)  
       co3 <- cos(2*pi*tl*w*3)
       si3 <- sin(2*pi*tl*w*3) 
       newarray<-array(evi2,c(dim1,dim2,636))          # scidb array to r array                
       
       
       fevi3b3<-newarray    
       fevi3b3[is.na(fevi3b3)] <- 0
       
       #if(length(which(fevi3b3==0))<100)
       # {
       aa2<-as.vector(fevi3b3)
       aa2[aa2==0]<-NA
       fevi3b3[fevi3b3==0]<-median(fevi3b3)    
       
       # }  
       if(dim1<3 || dim2<3)
{       
       if(the_i<2|| the_j <2) # fisrt row or first collumn
{
       fevi3b312t1<-ts(fevi3b3[1,1,],start=c(2000,1),frequency=46) # reconstruct the time series with its own
       rcol<-min(the_i) 
       rrow<-min(the_j)                 
} else {
       fevi3b312t1<-ts(fevi3b3[1,2,],start=c(2000,1),frequency=46) # reconstruct the time series with its own
       rcol<-max(the_i) 
       rrow<-max(the_j)        
}
       
       
       resar1<-coredata(residuals(gls(fevi3b312t1 ~ tl+co+co2+co3+si+si2+si3,correlation=corAR1())))
       
       p.Vt3 <- sctest(efp(fevi3b312t1 ~  tl+co+co2+co3+si+si2+si3,   h = 0.15, type = \"OLS-CUSUM\" )) 
       p.Vt4 <- sctest(efp(fevi3b312t1 ~  tl+co+co2+co3+si+si2+si3,   h = 0.15, type = \"OLS-MOSUM\" ))
       p.Vt5 <- sctest(efp(fevi3b312t1 ~  tl+co+co2+co3+si+si2+si3,   h = 0.15, type = \"OLS-CUSUM\" ,spatial1=as.numeric(resar1)) ) 
       p.Vt6 <- sctest(efp(fevi3b312t1 ~  tl+co+co2+co3+si+si2+si3,   h = 0.15, type = \"OLS-MOSUM\" ,spatial1=as.numeric(resar1)) )
       
       spcusum1<-p.Vt3$p.value # spautolm residuals CUSUM
       spmosum1 <-p.Vt4$p.value # spautolm residuals  MOSUM  
       cusum1 <-p.Vt3$p.value # CUSUM
       mosum1<-p.Vt4$p.value # MOSUM
       cusumar1 <-p.Vt5$p.value # CUSUM ar 1
       mosumar1 <-p.Vt6$p.value # MOSUM ar1
       rcol<-min(the_i)*100
       rrow<-min(the_j)*100
} else {  
       
       
       fevi3b312t1<-ts(fevi3b3[2,2,],start=c(2000,1),frequency=46) # reconstruct the time series
       
       try2<-spautolm(aa2~. , data.frame(aa2,X),family=\"SAR\",method= \"Matrix\", listw=listcn636,na.action=na.exclude,zero.policy=TRUE)
       
       rn<-lapply(1:9,function(i) {residuals(try2)[seq(i,636*9-(9-i),9)]})
       #get residuals for each time series
       
       ii<-5   # get the middle pixel (5 for 3*3 matrix)
       
       resar1<-coredata(residuals(gls(fevi3b312t1 ~ tl+co+co2+co3+si+si2+si3,correlation=corAR1())))
       
       p.Vt1  <- sctest(efp(fevi3b312t1 ~ tl+co+co2+co3+si+si2+si3,    h = 0.15, type = \"OLS-CUSUM\", spatial1=as.numeric(rn[[ii]]))  )
       p.Vt2  <- sctest(efp(fevi3b312t1 ~ tl+co+co2+co3+si+si2+si3,    h = 0.15, type = \"OLS-MOSUM\", spatial1=as.numeric(rn[[ii]]))  )
       p.Vt3 <-  sctest(efp(fevi3b312t1 ~  tl+co+co2+co3+si+si2+si3,   h = 0.15, type = \"OLS-CUSUM\" )) 
       p.Vt4 <-  sctest(efp(fevi3b312t1 ~  tl+co+co2+co3+si+si2+si3,   h = 0.15, type = \"OLS-MOSUM\" ))
       p.Vt5 <-  sctest(efp(fevi3b312t1 ~  tl+co+co2+co3+si+si2+si3,   h = 0.15, type = \"OLS-CUSUM\" ,spatial1=as.numeric(resar1)) ) 
       p.Vt6 <-  sctest(efp(fevi3b312t1 ~  tl+co+co2+co3+si+si2+si3,   h = 0.15, type = \"OLS-MOSUM\" ,spatial1=as.numeric(resar1)) )
       
       spcusum1  <-p.Vt1$p.value # spautolm residuals CUSUM
       spmosum1  <-p.Vt2$p.value # spautolm residuals  MOSUM  
       cusum1    <-p.Vt3$p.value # CUSUM
       mosum1    <-p.Vt4$p.value # MOSUM
       cusumar1 <-p.Vt5$p.value # CUSUM ar 1
       mosumar1 <-p.Vt6$p.value # MOSUM ar1
        rcol<-(min(the_i)+1)
        rrow <-(min(the_j)+1)  
              }
 
     
list(spcusum1,spmosum1,cusum1,mosum1,cusumar1,mosumar1,rcol,rrow)
   
 ') ",
       return=TRUE
       )


