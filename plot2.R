 x<-bfa
plot2.bfast <- function (x) 
{
   
  out <- x$output 
  TSt <- out$Yt
  noise <- x$out$Nt
  
  ft <- cbind(Tseasonal = out$TSt, remainder = out$Nt)
    tsp(ft) <- tsp(x$Yt)
    ft <- list(time.series = ft)
    # fit = x passes the BFAST object to seasonal() for ANOVA
par(mfrow=c(3,1))
  plot(out$Yt,ylab='time series' ) 
  plot(out$TSt,ylab='trend and seasonality')
  if (!x$nobp) {
          lines(out$bp.Vt)
          lines(out$ci.Vt)
          legend("topright", paste("Time of BP(s)", paste(out$Vt.bp, 
                                                          collapse = ",")), col = 2)
                 }
plot(out$Nt,ylab='noise')
par(mfrow=c(1,1))
}
plot2.bfast(tr2)
plot(tr1)
par(mfrow=c(1,1))
var(tr2$output$Nt) #var 0.01525689
var(tr1$output[[1]]$Nt)
plot(tr2$output$Nt[1:150],typ='l')
lines(tr1$output[[1]]$Nt[1:150],col='red') # bfast residual  var 0.01507316
