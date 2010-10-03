require(bitops)

setClass("closure",
  representation(
    defining = "integer",
    hypotheses = "character",
    alpha = "numeric"
  )
)

closed <- function(test, hypotheses, alpha = 0.05) {

  N <- length(hypotheses)
  Nmax <- log2(.Machine$integer.max+1)
  if (N > Nmax)
    stop("no more than", Nmax, "hypotheses supported in full closed testing.\n Use a shortcut-based test.")
  closure <- 1:(2^N-1)
  base <- 2^(1:N-1)

  # finds offspring hypotheses of a hypothesis (NB including self)
  offspring <- function(x) {
    res <- bitAnd(x, closure)
    res[res != 0]
  }

  # sort the closure to decreasing number of participating hypotheses
  lengths <- rowSums(sapply(base, function(bs) bitAnd(closure, bs) != 0))
  closure <- closure[sort.list(lengths, decreasing = TRUE)]
  sort.closure <- sort.list(closure)
  
  # perform closed testing
  futile <- rep(FALSE, length(closure))
  for (i in closure) {
    if (!futile[i] && test(hypotheses[.bit2boolean(i, N)]) > alpha) {
      futile[offspring(i)] <- TRUE
    }
  }
  rejected <- which(!futile)
  
  # reduce: find defining rejections
  def <- .defining(rejected, N)
  
  # return
  out <- new("closure", defining = def, hypotheses = hypotheses, alpha = alpha)
  out
}

setMethod("show", "closure", function(object) {
  cat("Closed testing result on ", length(object@hypotheses), " hypotheses at significance level ", object@alpha, ".\n", sep="")
  res <- pick(object, object@hypotheses, silent=TRUE)
  cat("False hypotheses >= ", length(object@hypotheses) - res, "; ", sep="")
  cat("True hypotheses <= ", res, ".\n", sep="")
  object
})

setMethod("summary", "closure", function(object) {
  cat("Closed testing result on ", length(object@hypotheses), " hypotheses at significance level ", object@alpha, ".\n", sep="")
  res <- pick(object, object@hypotheses, silent=TRUE)
  cat("False hypotheses >= ", length(object@hypotheses) - res, "; ", sep="")
  cat("True hypotheses <= ", res, ".\n\n", sep="")
  cat("Defining rejections:\n")
  print(defining(object))
  invisible(object)
})

setGeneric("defining", function(object, ...) standardGeneric("defining"))
setMethod("defining", "closure", function(object, ...) {
  .num2names(object@defining, object@hypotheses)
})

setGeneric("hypotheses", function(object, ...) standardGeneric("hypotheses"))
setMethod("hypotheses", "closure", function(object, ...) {
  object@hypotheses
})

# gets the defining rejections from all rejections
.defining <- function(rejected, N) {
  closure <- 1:(2^N-1)
  ancestors <- function(x) {
    bitOr(x, closure)
  }
  isdone <- integer(0)
  todo <- rejected
  while (length(todo) > 0) {
    isdone <- c(setdiff(isdone, ancestors(todo[1])), todo[1])
    todo <- setdiff(todo, ancestors(todo[1]))
  }
  isdone
}

# converts from integer to boolean (as binary)
.bit2boolean <- function(x, N) {
  base <- 2^(1:N-1)
  bitAnd(x, base) != 0
}

.num2names <- function(rejected, vars) {
  N <- length(vars)
  bools <- lapply(rejected, .bit2boolean, N=N)
  lapply(bools, function(b) vars[b])
}

pick <- function(closure, reject, silent=FALSE) {
  N <- length(closure@hypotheses)
  reject <- which(closure@hypotheses %in% reject)
  clos <- 1:(2^N-1)
  interest <- unique(bitAnd(sum(2^(reject-1)), clos))
  interest <- interest[interest>0]
  isAncestor <- function(x,y) { # is x an ancestor of y?
    bitOr(x,y) == x
  }
  interest <- interest[!apply(outer(interest, closure@defining, isAncestor), 1, any)]
  base <- 2^(1:N-1)
  lengths <- lapply(base, function(bs) bitAnd(interest, bs) != 0)
  lengths <- rowSums(do.call(cbind, lengths))
  out <- max(lengths)
  if (!silent) {
    cat("Rejected ", length(reject), " hypotheses at confidence level ", 1-closure@alpha, ".\n", sep="")
    cat("Correct rejections >= ", length(reject)-out, "; ", sep="")
    cat("False rejections <= ", out, ".\n", sep="")
    invisible(out)
  } else
    out
}

