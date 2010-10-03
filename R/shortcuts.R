pickFisher <- function(all, select = seq_along(all), alpha=0.05, silent=FALSE) {
  rej <- all[select]
  nr <- setdiff(all, rej)
  nr <- nr[nr > min(rej)]
  lrej <- -2*log(sort(rej, decreasing=TRUE))
  lnr <- -2*log(sort(nr, decreasing=TRUE))
  cum.r <- cumsum(lrej)
  cum.nr <- c(0, cumsum(lnr))
  crit.val <- sapply(1:length(cum.r), function(st) {
    max(qchisq(1-alpha, df=2*(0:(length(cum.nr)-1) + st), lower=TRUE) - cum.nr)
  })
  out <- max(c(0,which(cum.r <= crit.val)))
  if (!silent) {
    cat("Rejected ", length(rej), " hypotheses at confidence level ", 1-alpha, ".\n", sep="")
    cat("Correct rejections >= ", length(rej)-out, "; ", sep="")
    cat("False rejections <= ", out, ".\n", sep="")
    invisible(out)
  } else
    out
}

pickSimes <- function(all, select = seq_along(all), alpha=0.05, hommel = FALSE, silent=FALSE) {
  rej <- sort(all[select], decreasing=TRUE)
  nr <- sort(setdiff(all, rej), decreasing=TRUE)
  nr <- nr[nr > min(rej)]
  if (!hommel)
    res <- which(!sapply(1:length(rej), function(s)
      all(sapply(0:length(nr), function(i)
        any(sort(c(rej[1:s], nr[seq_len(i)])) <= alpha*(1:(s+i))/(s+i) )
      ))
    ))
  else
    res <- which(!sapply(1:length(rej), function(s)
      all(sapply(0:length(nr), function(i)
        any(sort(c(rej[1:s], nr[seq_len(i)]))*sum(1/(1:(s+i))) < alpha*(1:(s+i))/(s+i) )
      ))
    ))
  out <- max(c(0,res))
  if (!silent) {
    cat("Rejected ", length(rej), " hypotheses at confidence level ", 1-alpha, ".\n", sep="")
    cat("Correct rejections >= ", length(rej)-out, "; ", sep="")
    cat("False rejections <= ", out, ".\n", sep="")
    invisible(out)
  } else
    out
}

