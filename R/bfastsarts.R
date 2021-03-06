
bfastsar <- function(starr, LE=636, h=0.15 , season =c("dummy","harmonic","none"), max.iter = 3, breaks = NULL, hpc = "none", level = 0.05, type= "OLS-MOSUM")
{
  require('spdep')
  le=LE
  tl=1:le
  w=1/46
  #harmonic
  
  Yt<-ts(starr[2,2,],start=c(2000,1),frequency=46)
  
  co <- cos(2*pi*tl*w); si <- sin(2*pi*tl*w)
  co2 <- cos(2*pi*tl*w*2);si2 <- sin(2*pi*tl*w*2)  
  co3 <- cos(2*pi*tl*w*3);si3 <- sin(2*pi*tl*w*3) 
  
  X = matrix(0, le * 9, 9*8)
  
  for( i in 1:9)
  {
    
    X[seq(i,by=9,length.out=le),1+(i-1)*8] = 1 
    X [seq(i,by=9,length.out=le),2+(i-1)*8] = tl 
    X[seq(i,by=9,length.out=le),3+(i-1)*8] =co
    X[seq(i,by=9,length.out=le),4+(i-1)*8] =co2
    X[seq(i,by=9,length.out=le),5+(i-1)*8] =co3
    X[seq(i,by=9,length.out=le),6+(i-1)*8] =si
    X[seq(i,by=9,length.out=le),7+(i-1)*8] =si2
    X[seq(i,by=9,length.out=le),8+(i-1)*8] =si3
  }
  
  
  colnames(X) = paste0("v", 1:(9*8))
  
  
  
  #sar
  season <- match.arg(season)
  level = rep(level, length.out = 2)
  ti <- time(Yt)
  
  f <- frequency(Yt)      # on cycle every f time points (seasonal cycle)
  if(class(Yt)!="ts")
    stop ("Not a time series object")
  ## return value
  output <- list()
  Tt <- 0
  
  
  
  # seasonal model setup
  if (season=="harmonic") {
    w <- 1/f # f = 23 when freq=23 :-)
    tl <- 1:length(Yt)
    co <- cos(2*pi*tl*w); si <- sin(2*pi*tl*w)
    co2 <- cos(2*pi*tl*w*2);si2 <- sin(2*pi*tl*w*2)
    co3 <- cos(2*pi*tl*w*3);si3 <- sin(2*pi*tl*w*3)
    smod <- Wt ~ co+si+co2+si2+co3+si3
    # Start the iterative procedure and for first iteration St=decompose result
    St <- stl(Yt, "periodic")$time.series[, "seasonal"]
    
  } else if (season=="dummy") {
    # Start the iterative procedure and for first iteration St=decompose result
    St <- stl(Yt, "periodic")$time.series[, "seasonal"]
    D <- seasonaldummy(Yt)
    D[rowSums(D) == 0,] <- -1
    smod <- Wt ~ -1 + D
  } else if (season == "none") {
    print("No seasonal model will be fitted!")
    St <- 0
  } else stop("Not a correct seasonal model is selected ('harmonic' or 'dummy') ")
  
  # number/timing of structural breaks in the trend/seasonal component
    fevi3b312t<-apply(starr,c(1,2),function(x) (ts(x,start=c(2000,1),frequency=46)-St))
    f2<-aperm(fevi3b312t,c(2,3,1))  
    aa2<-as.vector(f2) 
    try2<-spautolm(aa2~. , data.frame(aa2,X),family="SAR",method= "Matrix", listw=listcn636)   
    rn<-lapply(1:9,function(i) {residuals(try2)[seq(i,le*9-(9-i),9)]})
    Vt<-ts(f2[2,2,],start=c(2000,1),frequency=46) # reconstruct the time series
    
    
    p.Vt <- sctest(efp(Vt ~ ti, h=h, type=type,spatial1=rn[[5]]))
  
  
  Vt.bp <- 0
  Wt.bp <- 0 
  CheckTimeTt <- 1
  CheckTimeSt <- 1
  i <- 0
  while ( (!identical(CheckTimeTt,Vt.bp) | !identical(CheckTimeSt,Wt.bp)) & i < max.iter)
  {
    CheckTimeTt <- Vt.bp
    CheckTimeSt <- Wt.bp
    # TREND
    # Vt <- Yt-St
  
    
    if (p.Vt$p.value <= level[1]) 
    {
      bp.Vt <- breakpoints(Vt ~ ti, h=h,breaks=breaks, hpc = hpc)
      nobp.Vt <- is.na(breakpoints (bp.Vt)[1])
    } 
    else 
    {
      nobp.Vt <- TRUE
      bp.Vt <- NA       
    }
    if (nobp.Vt)
    {
      fm0 <- lm(Vt ~  ti)
      Vt.bp <- 0      # no breaks times
      Tt <- ts(fitted(fm0))     # Data minus trend
      tsp(Tt) <- tsp(Yt)
      ci.Vt <- NA
    } 
    else
    {
      fm1 <- lm(Vt ~ breakfactor(bp.Vt)/ti)
      ci.Vt <- confint(bp.Vt, het.err = FALSE)
      Vt.bp <- ci.Vt$confint[,2]
      Tt <- ts(fitted(fm1))     # Data minus trend
      tsp(Tt) <- tsp(Yt)
    }
    
    # SEASONAL COMPONENT
    if (season=="none") {
      Wt <- 0
      St <- 0
      bp.Wt <- NA; ci.Wt <- NA; nobp.Wt<- TRUE
    } else
    {
      Wt <- Yt - Tt
      if (p.Vt$p.value <= level[2]) # OR statement 
      {
        bp.Wt <- breakpoints(smod, h=h,breaks=breaks, hpc = hpc) # Breakpoints in the seasonal component
        nobp.Wt <- is.na(breakpoints (bp.Wt)[1])
      } 
      else 
      {
        nobp.Wt <- TRUE
        bp.Wt <- NA       
      }
      if (nobp.Wt)
      {
        sm0 <- lm(smod)
        St <- ts(fitted(sm0))  #  The fitted seasonal component
        tsp(St) <- tsp(Yt)
        Wt.bp <- 0             # no seasonal breaks
        ci.Wt <- NA
      } 
      else
      {
        if(season=="dummy") sm1 <-lm(Wt ~ -1+D %in% breakfactor(bp.Wt))
        if(season=="harmonic") sm1 <- lm(Wt ~ (co+si+co2+si2+co3+si3) %in% breakfactor(bp.Wt)) 
        St <- ts(fitted(sm1))  #  The fitted seasonal component
        tsp(St) <- tsp(Yt)
        ci.Wt <- confint(bp.Wt, het.err = FALSE)
        Wt.bp <- ci.Wt$confint[,2] 
      }
    }
    i <- i+1
    output[[i]] <- list(Tt=Tt,St=St,Nt=Yt-Tt-St,
                        Vt=Vt, bp.Vt=bp.Vt, Vt.bp=Vt.bp, ci.Vt=ci.Vt,
                        Wt=Wt, bp.Wt=bp.Wt, Wt.bp=Wt.bp, ci.Wt=ci.Wt)
  }
  if (!nobp.Vt) # probably only works well for dummy model!
  {
    Vt.nrbp <- length(bp.Vt$breakpoints)
    co <- coef(fm1) # final fitted trend model
    Mag <- matrix(NA,Vt.nrbp,3)
    for (r in 1:Vt.nrbp) 
    {
      if (r==1) 
        y1 <- co[1]+co[r+Vt.nrbp+1]*ti[Vt.bp[r]]
      else 
        y1 <- co[1]+co[r]+co[r+Vt.nrbp+1]*ti[Vt.bp[r]]
      y2 <- (co[1]+co[r+1])+co[r+Vt.nrbp+2]*ti[Vt.bp[r]+1]
      Mag[r,1] <- y1
      Mag[r,2] <- y2
      Mag[r,3] <- y2-y1   
    }
    index <- which.max(abs(Mag[,3]))
    m.x <- rep(Vt.bp[index],2)
    m.y <- c(Mag[index,1],Mag[index,2]) #Magnitude position
    Magnitude <- Mag[index,3] # Magnitude of biggest change
    Time <- Vt.bp[index]
  } 
  else 
  {
    m.x <- NA; m.y <- NA
    Magnitude <- 0  # if we do not detect a break then the magnitude is zero
    Time <- NA # if we do not detect a break then we have no timing of the break
    Mag <- 0
  }
  return(structure(list(Yt=Yt,output=output,nobp=list(Vt=nobp.Vt,Wt=nobp.Wt),Magnitude=Magnitude,Mags=Mag,
                        Time=Time,jump=list(x=ti[m.x],y=m.y)),class="bfast"))  
}
#i=1
#starr<-fevi8[i:(i+2),j:(j+2),]
#
#bfastsar(starr,season='harmonic', max.iter = 3)
