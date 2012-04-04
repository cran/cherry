pickFisher <- function(p, select = seq_along(p), alpha=0.05, silent=FALSE) {

  rej <- p[select]
  if (any(is.na(rej))) stop("invalid selection or NA in p-values")
  nr <- setdiff(p, rej)
  nr <- nr[nr > min(rej)]
  lrej <- -2*log(sort(rej, decreasing=TRUE))
  lnr <- -2*log(sort(nr, decreasing=TRUE))
  cum.r <- cumsum(lrej)
  cum.nr <- c(0, cumsum(lnr))
  crit.val <- sapply(1:length(cum.r), function(st) {
    max(qchisq(1-alpha, df=2*(0:(length(cum.nr)-1) + st), lower.tail=TRUE) - cum.nr)
  })
  out <- max(c(0,which(cum.r <= crit.val)))
  if (!silent) {
    cat("Rejected ", length(rej), " hypotheses at confidence level ", 1-alpha, ".\n", sep="")
    cat("Correct rejections >= ", length(rej)-out, "; ", sep="")
    cat("False rejections <= ", out, ".\n", sep="")
    invisible(length(rej)-out)
  } else
    length(rej)-out
}


curveFisher <- function(p, select = seq_along(p), order, alpha=0.05, plot = TRUE) {

  selected <- !missing(select) || missing(order)
  ordered <- !missing(order)
  if (ordered & selected) 
    stop("please provide either select or order, but not both")
  if (selected)
    ranks <- sort(rank(p, ties.method="first")[select])
  else
    ranks <- rank(p, ties.method="first")[order]
  if (length(ranks)==0 || any(is.na(ranks)))
    stop("invalid selection or NA in p-values")
  others <- setdiff(length(p):1, ranks)
  lpv <- -2*log(sort(p))
  res <- numeric(length(ranks))
  chisqs <- qchisq(1-alpha, df=2*1:length(lpv), lower.tail=TRUE)
  st <- 1
  for (ed in 1:length(ranks)) {
    if (selected)
      ins <- ranks[seq(st,ed)]
    else
      ins <- sort(ranks[1:ed])[seq(st,ed)]
    outs <- setdiff(length(lpv):min(ins), ins)
    cr.v <- max(chisqs[ed-st+1+0:length(outs)] - cumsum(c(0,lpv[outs])))
    rej <- (sum(lpv[ins]) >= cr.v)
    if (rej)
      st <- st+1
    res[ed] <- st-1
  }
  names(res) <- names(ranks)
  if (plot) {
    false <- c(0, res)
    xs <- 1:length(false)-.5
    tots <- 0:length(res)
    plot(xs, tots, type="S", xlab="number of rejections", ylab="number of rejections", lty=2)
    lines(xs, false, type="S")
    legend("topleft", c(paste("correct rejections (", 100*(1-alpha), "% conf.)", sep=""),"others"), lty=1:2)
    invisible(res)
  } else
    res

}


pickSimes <- function(p, select = seq_along(p), alpha=0.05, hommel=FALSE, silent=FALSE) {

  ranks <- sort(rank(p, ties.method="first")[select])
  p <- sort(p)
  others <- setdiff(length(p):1, ranks)
  st <- 1
  ed <- length(ranks)
  ready <- FALSE
  while (!ready) {
    ins <- seq_along(p) %in% ranks[seq(st,ed)]
    outs <- (!logical(length(p))) & (cummax(ins)==1) & (!ins)
    participate <- numeric(length(p))
    participate[ins] <- 1+sum(outs)
    participate[outs] <- seq_len(sum(outs))
    maxlag <- cumsum(outs)
    rej <- TRUE
    i <- 0
    while (rej && (i <= sum(outs))) {
      bottom.is <- (participate > i)
      K <- sum(bottom.is)
      if (hommel)
        lag <- floor(1:K - p[bottom.is]/(alpha/(K*sum(1/1:K))))
      else
        lag <- floor(1:K - p[bottom.is]/(alpha/K))
      if (any(lag >= 0 & lag >= maxlag[bottom.is] - i & ins[bottom.is]))
        i <- Inf
      else if (any(lag >= 0))
         i <- i + 1 + max(pmin(lag, maxlag[bottom.is] - i))
      else
        rej <- FALSE
    }
    if (rej) {
      st <- st+1
      ready <- st > ed
    } else
      ready <- TRUE
  }
  out <- ed-st+1
  if (!silent) {
    cat("Rejected ", length(ranks), " hypotheses. At confidence level ", 1-alpha, ":\n", sep="")
    cat("Correct rejections >= ", length(ranks)-out, "; ", sep="")
    cat("False rejections <= ", out, ".\n", sep="")
    invisible(length(ranks)-out)
  } else
    length(ranks)-out
}


curveSimes <- function(p, select = seq_along(p), order, alpha=0.05, hommel=FALSE, plot = TRUE) {

  selected <- !missing(select) || missing(order)
  ordered <- !missing(order)
  if (ordered & selected) 
    stop("please provide either select or order, but not both")
  if (selected) {
    ranks <- sort(rank(p, ties.method="first")[select])
    endpoint <- length(ranks) - pickSimes(p, select, alpha, hommel, silent=TRUE)
  } else {
    ranks <- rank(p, ties.method="first")[order]
    endpoint <- length(ranks) - pickSimes(p, order, alpha, hommel, silent=TRUE)
  }
  if (length(ranks)==0 || any(is.na(ranks)))
    stop("invalid selection or NA in p-values")
  p <- sort(p)
  others <- setdiff(length(p):1, ranks)
  res <- numeric(length(ranks))
  st <- 1
  for (ed in 1:length(ranks)) {
    if (selected)
      ins <- seq_along(p) %in% ranks[seq(st,ed)]
    else
      ins <- seq_along(p) %in% sort(ranks[1:ed])[seq(st,ed)]
    outs <- (!logical(length(p))) & (cummax(ins)==1) & (!ins)
    participate <- numeric(length(p))
    participate[ins] <- 1+sum(outs)
    participate[outs] <- seq_len(sum(outs))
    maxlag <- cumsum(outs)
    rej <- TRUE
    i <- 0
    while (rej && (i <= sum(outs))) {
      bottom.is <- (participate > i)
      K <- sum(bottom.is)
      if (hommel)
        lag <- floor(1:K - p[bottom.is]/(alpha/(K*sum(1/1:K))))
      else
        lag <- floor(1:K - p[bottom.is]/(alpha/K))
      if (any(lag >= 0 & lag >= maxlag[bottom.is] - i & ins[bottom.is]))
        i <- Inf
      else if (any(lag >= 0))
         i <- i + 1 + max(pmin(lag, maxlag[bottom.is] - i))
      else
        rej <- FALSE
    }
    if (rej)
      st <- st+1
    res[ed] <- st-1
    if (st > endpoint) {
      res[ed:length(ranks)] <- endpoint
      break
    }
  }
  names(res) <- names(p[ranks])
  if (plot) {
    false <- c(0, res)
    xs <- 1:length(false)-.5
    tots <- 0:length(res)
    plot(xs, tots, type="S", xlab="number of rejections", ylab="number of rejections", lty=2)
    lines(xs, false, type="S")
    legend("topleft", c(paste("correct rejections (", 100*(1-alpha), "% conf.)", sep=""),"others"), lty=1:2)
    invisible(res)
  } else
    res
}









