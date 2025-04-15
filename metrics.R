library(philentropy)
wasserstein_dist <- function(a, b, p=1, wa=NULL, wb=NULL) {
  m <- length(a)
  n <- length(b)
  stopifnot(m > 0 && n > 0)
  if (m == n && is.null(wa) && is.null(wb)) {
    return(mean(abs(sort(b)-sort(a))^p)^(1/p))
  }
  stopifnot(is.null(wa) || length(wa) == m)
  stopifnot(is.null(wb) || length(wb) == n)
  if (is.null(wa)) {
    wa <- rep(1,m)
  } else { # remove points with zero weight
    wha <- wa > 0
    wa <- wa[wha]
    a <- a[wha]
    m <- length(a)
  }
  if (is.null(wb)) {
    wb <- rep(1,n)
  } else { # remove points with zero weight
    whb <- wb > 0
    wb <- wb[whb]
    b <- b[whb]
    n <- length(b)
  }
  
  orda <- order(a)
  ordb <- order(b)
  a <- a[orda]
  b <- b[ordb]
  wa <- wa[orda]
  wb <- wb[ordb]
  ua <- (wa/sum(wa))[-m]
  ub <- (wb/sum(wb))[-n]
  cua <- c(cumsum(ua))  
  cub <- c(cumsum(ub))  
  arep <- hist(cub, breaks = c(-Inf, cua, Inf), plot = FALSE)$counts + 1
  brep <- hist(cua, breaks = c(-Inf, cub, Inf), plot = FALSE)$counts + 1
  # we sum over rectangles with cuts on the vertical axis each time one of the two ecdfs makes a jump
  # arep and brep tell us how many times each of the a and b data have to be repeated in order to get the points on the horizontal axis
  # note that sum(arep)+sum(brep) = m+n-1 (we do not count the height-zero final rectangle where both ecdfs jump to 1)
  
  aa <- rep(a, times=arep)
  bb <- rep(b, times=brep)
  
  uu <- sort(c(cua,cub))
  uu0 <- c(0,uu)
  uu1 <- c(uu,1)
  areap <- sum((uu1-uu0)*abs(bb-aa)^p)^(1/p)
  #  print(rbind(uu1-uu0, pmax(aa,bb)-pmin(aa,bb)))
  return(areap)
}

JSD_dist <- function(a, b) {
  nbins=20
  a_b_breaks <- seq(min(c(a,b), na.rm = T), max(c(a,b), na.rm = T), length.out=nbins+1)
  a_counts <-  hist(a, breaks=a_b_breaks, plot=F)
  b_counts <-  hist(b, breaks=a_b_breaks, plot=F)
  
  prob_counts_a <- a_counts$counts / sum(a_counts$counts)
  prob_counts_b <- b_counts$counts / sum(b_counts$counts)
  
  prob_mat <- rbind(prob_counts_a,
                    prob_counts_b)
  
  d <- JSD(prob_mat)
  names(d) <- c()
  return(d)  
}

KL_dist <- function(a, b) {
  nbins=20
  a_b_breaks <- seq(min(c(a,b), na.rm = T), max(c(a,b), na.rm=T), length.out=nbins+1)
  a_counts <-  hist(a, breaks=a_b_breaks, plot=F)
  b_counts <-  hist(b, breaks=a_b_breaks, plot=F)
  
  prob_counts_a <- a_counts$counts / sum(a_counts$counts)
  prob_counts_b <- b_counts$counts / sum(b_counts$counts)
  
  prob_mat <- rbind(prob_counts_a,
                    prob_counts_b)
  
  d <- KL(prob_mat)
  
  names(d) <- c()
  return(d)  
}

dprime_dist <- function(a , b) {
  return(abs(mean(a, na.rm=T) - mean(b, na.rm=T)) / (0.5 * sd(a, na.rm=T) + 0.5 * sd(b, na.rm=T)))
}

mse <- function(a, b) {
  return(mean((a - b) ** 2))
}
