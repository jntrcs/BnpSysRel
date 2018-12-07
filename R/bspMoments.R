##bspMoments.R


#' Calculates the first and second moments of a bsp object at each spot on the support
#'
#' @param bsp The bsp object
#'
#' @return what should this be??
#' @export
#'
#' @examples
#' bsp=bsp(c(1:3), centeringMeasure = c(.1,.9, .98), precision = 2)
#' E1E2(bsp)
#'
#'
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

#b=E1E2(bsp)
#varx = b$E2 - b$E1^2


################################################################################
##
## Find First and Second Moment of Merged BSPs for Components in Series.
##
################################################################################
E1E2_series <- function(bsp1, bsp2) {

  times = unique(sort(c(bsp1$support,bsp2$support)))
  n <- length(times)

  new_C1_E1 <- rep(NA,n)
  new_C1_E2 <- rep(NA,n)
  for (i in 1:n) {
    if (length(which(times[i]>=C1_times))>0) {
      ind <- max(which(times[i]>=C1_times))
      new_C1_E1[i] <- C1_E1[ind]
      new_C1_E2[i] <- C1_E2[ind]
    } else {
      new_C1_E1[i] <- 0
      new_C1_E2[i] <- 0
    }
  }

  new_C2_E1 <- rep(NA,n)
  new_C2_E2 <- rep(NA,n)
  for (i in 1:n) {
    if (length(which(times[i]>=C2_times))>0) {
      ind <- max(which(times[i]>=C2_times))
      new_C2_E1[i] <- C2_E1[ind]
      new_C2_E2[i] <- C2_E2[ind]
    } else {
      new_C2_E1[i] <- 0
      new_C2_E2[i] <- 0
    }
  }

  E1 <- 1-(1-new_C1_E1)*(1-new_C2_E1)
  #Why are we multiplying by 1-new_c1_E1, page nine seems to say times by the base
  E2 <- 1-2*(1-new_C1_E1)*(1-new_C2_E1)+(1-2*new_C1_E1+new_C1_E2)*(1-2*new_C2_E1+new_C2_E2)

  m <- which.max(E1)

  ## Return list E1, E2, and times
  return( list("E1" = E1[1:m], "E2" = E2[1:m], "times"= times[1:m]) )

}

################################################################################
##
## Find First and Second Moment of Merged BSPs for Components in Parallel.
##
################################################################################
E1E2_parallel <- function(C1_E1, C1_E2, C1_times, C2_E1, C2_E2, C2_times) {

  times = unique(sort(c(C1_times,C2_times)))
  n <- length(times)

  new_C1_E1 <- rep(NA,n)
  new_C1_E2 <- rep(NA,n)
  for (i in 1:n) {
    if (length(which(times[i]>=C1_times))>0) {
      ind <- max(which(times[i]>=C1_times))
      new_C1_E1[i] <- C1_E1[ind]
      new_C1_E2[i] <- C1_E2[ind]
    } else {
      new_C1_E1[i] <- 0
      new_C1_E2[i] <- 0
    }
  }

  new_C2_E1 <- rep(NA,n)
  new_C2_E2 <- rep(NA,n)
  for (i in 1:n) {
    if (length(which(times[i]>=C2_times))>0) {
      ind <- max(which(times[i]>=C2_times))
      new_C2_E1[i] <- C2_E1[ind]
      new_C2_E2[i] <- C2_E2[ind]
    } else {
      new_C2_E1[i] <- 0
      new_C2_E2[i] <- 0
    }
  }

  E1 <- new_C1_E1*new_C2_E1

  E2 <- new_C1_E2*new_C2_E2

  m <- n-which.min(rev(E1))+1

  if (E1[m]==0 & E2[m]==0)m=m+1

  ## Return list E1, E2, and times
  return( list("E1" = E1[m:n], "E2" = E2[m:n], "times"= times[m:n]) )

}
