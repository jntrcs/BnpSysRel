##bspMoments.R

E1E2 <- function(bsp) {

  base<-bsp$centeringMeasure
  prec<-bsp$precisionAt
  n<-length(base)
  ## Calculate the 2nd Moment
  temp <- ((1-base[-1])*(prec[-1]*(1-base[-1])+1))/
    ((1-base[-n])*(prec[-1]*(1-base[-n])+1))
  E2 <- cumprod(temp)-1+2*base[-1]
  #if (any(is.nan(E2)))print(paste("E2 had a NAN", E2))
  #E2[is.nan(E2)] <- 0

  ## Return list E1, E2, and times
  return( data.frame(E1=base, E2=c(0, E2)))

}

b=E1E2(bsp)
varx = b$E2 - b$E1^2
