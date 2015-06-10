bfastmore<-function (Yt, h = 0.15, season = c("dummy", "harmonic", "none"), 
                     max.iter = NULL, breaks = NULL, hpc = "none", level = 0.05, 
                     type = "OLS-MOSUM") 
{
  require(strucchange)
  season <- match.arg(season)
  level = rep(level, length.out = 2)
  ti <- time(Yt)
  f <- frequency(Yt)
  if (class(Yt) != "ts") 
    stop("Not a time series object")
  output <- list()
  Tt <- 0
  if (season == "harmonic") {
    w <- 1/f
    tl <- 1:length(Yt)
    co <- cos(2 * pi * tl * w)
    si <- sin(2 * pi * tl * w)
    co2 <- cos(2 * pi * tl * w * 2)
    si2 <- sin(2 * pi * tl * w * 2)
    co3 <- cos(2 * pi * tl * w * 3)
    si3 <- sin(2 * pi * tl * w * 3)
    smod <- Wt ~ co + si + co2 + si2 + co3 + si3
    St <- stl(Yt, "periodic")$time.series[, "seasonal"]
  }
  else if (season == "dummy") {
    St <- stl(Yt, "periodic")$time.series[, "seasonal"]
    D <- seasonaldummy(Yt)
    D[rowSums(D) == 0, ] <- -1
    smod <- Wt ~ -1 + D
  }
  else if (season == "none") {
    print("No seasonal model will be fitted!")
    St <- 0
  }
  else stop("Not a correct seasonal model is selected ('harmonic' or 'dummy') ")
  Vt.bp <- 0
  Wt.bp <- 0
  CheckTimeTt <- 1
  CheckTimeSt <- 1
  i <- 0
  while ((!identical(CheckTimeTt, Vt.bp) | !identical(CheckTimeSt, 
                                                      Wt.bp)) & i < max.iter) {
    CheckTimeTt <- Vt.bp
    CheckTimeSt <- Wt.bp
    Vt <- Yt - St
    p.Vt <- sctest(efp(Vt ~ ti, h = h, type = type))
    if (p.Vt$p.value <= level[1]) {
      bp.Vt <- breakpoints(Vt ~ ti, h = h, breaks = breaks, 
                           hpc = hpc)
      nobp.Vt <- is.na(breakpoints(bp.Vt)[1])
    }
    else {
      nobp.Vt <- TRUE
      bp.Vt <- NA
    }
    if (nobp.Vt) {
      fm0 <- lm(Vt ~ ti)
      Vt.bp <- 0
      Tt <- ts(fitted(fm0))
      tsp(Tt) <- tsp(Yt)
      ci.Vt <- NA
    }
    else {
      fm1 <- lm(Vt ~ breakfactor(bp.Vt)/ti)
      ci.Vt <- confint(bp.Vt, het.err = FALSE)
      Vt.bp <- ci.Vt$confint[, 2]
      Tt <- ts(fitted(fm1))
      tsp(Tt) <- tsp(Yt)
    }
    if (season == "none") {
      Wt <- 0
      St <- 0
      bp.Wt <- NA
      ci.Wt <- NA
      nobp.Wt <- TRUE
    }
    else {
      Wt <- Yt - Tt
      p.Wt <- sctest(efp(smod, h = h, type = type))
      if (p.Wt$p.value <= level[2]) {
        bp.Wt <- breakpoints(smod, h = h, breaks = breaks, 
                             hpc = hpc)
        nobp.Wt <- is.na(breakpoints(bp.Wt)[1])
      }
      else {
        nobp.Wt <- TRUE
        bp.Wt <- NA
      }
      if (nobp.Wt) {
        sm0 <- lm(smod)
        St <- ts(fitted(sm0))
        tsp(St) <- tsp(Yt)
        Wt.bp <- 0
        ci.Wt <- NA
      }
      else {
        if (season == "dummy") 
          sm1 <- lm(Wt ~ -1 + D %in% breakfactor(bp.Wt))
        if (season == "harmonic") 
          sm1 <- lm(Wt ~ (co + si + co2 + si2 + co3 + 
                            si3) %in% breakfactor(bp.Wt))
        St <- ts(fitted(sm1))
        tsp(St) <- tsp(Yt)
        ci.Wt <- confint(bp.Wt, het.err = FALSE)
        Wt.bp <- ci.Wt$confint[, 2]
      }
    }
    i <- i + 1
    output[[i]] <- list(Tt = Tt, St = St, Nt = Yt - Tt - 
                          St, Vt = Vt, bp.Vt = bp.Vt, Vt.bp = Vt.bp, ci.Vt = ci.Vt, 
                        Wt = Wt, bp.Wt = bp.Wt, Wt.bp = Wt.bp, ci.Wt = ci.Wt)
  }
  if (!nobp.Vt) {
    Vt.nrbp <- length(bp.Vt$breakpoints)
    co <- coef(fm1)
    Mag <- matrix(NA, Vt.nrbp, 3)
    for (r in 1:Vt.nrbp) {
      if (r == 1) 
        y1 <- co[1] + co[r + Vt.nrbp + 1] * ti[Vt.bp[r]]
      else y1 <- co[1] + co[r] + co[r + Vt.nrbp + 1] * 
        ti[Vt.bp[r]]
      y2 <- (co[1] + co[r + 1]) + co[r + Vt.nrbp + 2] * 
        ti[Vt.bp[r] + 1]
      Mag[r, 1] <- y1
      Mag[r, 2] <- y2
      Mag[r, 3] <- y2 - y1
    }
    index <- which.max(abs(Mag[, 3]))
    m.x <- rep(Vt.bp[index], 2)
    m.y <- c(Mag[index, 1], Mag[index, 2])
    Magnitude <- Mag[index, 3]
    Time <- Vt.bp[index]
  }
  else {
    m.x <- NA
    m.y <- NA
    Magnitude <- 0
    Time <- NA
    Mag <- 0
  }
  return(structure(list(Yt = Yt, output = output, nobp = list(Vt = nobp.Vt, 
                                                              Wt = nobp.Wt), Magnitude = Magnitude, Mags = Mag, co = co, Time = Time, 
                        jump = list(x = ti[m.x], y = m.y)), class = "bfast"))
}


bfaclass<-function(fit)
{
  co<-fit$co
  mag<-fit$Mags
  ttrend<-fit$output[[1]]$Vt.bp  # time of trend
  # tsea<-fit$output[[1]]$Wt.bp  #time of seasonality
  
  nbp<-length(fit$output[[1]]$Vt.bp)
  
  out <- rep(NA, nbp)
  

  if(ttrend[1] == 0) {     
    for(r in 1:nbp)
    {
    slope <- co[2]# slope
    if(slope > 0) out[r] <- 1
    if(slope < 0) out[r] <- 2
    }
  } else {
    ## if break, list segment and break point parameters (p$..)
    ToB <- as.numeric(ttrend)  # time of break
    s<-c()
    m<-c()
    for(r in 1: (nbp+1))
    {
      
      s[r] <- co[nbp+1+r] # slope segment 1
      
    }
    for(r in 1:nbp)
    {
      m  <- mag[r,3] # magnitude of abrupt change
      if(s[r]> 0 && s[r+1] > 0 && m > 0) out[r] <- 3
      if(s[r] < 0 && s[r+1] < 0 && m < 0) out[r] <- 4
      # interrupted gradual change (setback or boost)
      if(s[r] > 0 &&  s[r+1] > 0 && m < 0) out[r] <- 5
      if(s[r] < 0 &&  s[r+1] < 0 && m > 0) out[r] <- 6
      # trend reversal (greening to browning v.v.)
      if(s[r] > 0 &&  s[r+1] < 0) out[r] <- 7
      if(s[r] < 0 &&  s[r+1] > 0) out[r] <- 8
    }
  }
  return(out)}