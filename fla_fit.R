## =============================================================================
## fla_fit: full likelihood approach (FLA) for logistic regression with
## misclassified binary responses and repeated classifications (Section 2.2
## and Supplementary Section S4 of the paper). Base R only.
##
## Usage (run this file once, or source("fla_fit.R"), then):
##   fit <- fla_fit(X, S, m)
##     X : model matrix, with a column of ones for the intercept,
##         e.g. X <- model.matrix(~ x1 + x2, data)
##     S : number of positive classifications of each unit (0, ..., m)
##     m : number of classifications per unit (a number, or a vector with one
##         value per unit when only part of the units is classified repeatedly)
##     w : optional frequency weights
##   fit$coef   estimates of (beta, e1, e2)
##   fit$se     standard errors from the observed information
##   fit$vcov   covariance matrix of the estimates
##   fit$ok     TRUE if the maximisation converged, the observed
##              information is positive definite and the estimates of e1
##              and e2 are not at the boundary (below 1e-4)
##   fit$loglik log-likelihood at the maximum
## Wald tests and confidence intervals follow from fit$coef and fit$se; see
## Supp_S4_fla_example.R and Section5_lens.R.
## =============================================================================

expit <- function(z) 1 / (1 + exp(-z))
`%||%` <- function(a, b) if (is.null(a)) b else a

## ----------------------------------------------------------------------------
##    m may be a scalar or a vector (number of classifications of each unit),
##    so that samples in which only part of the units is classified repeatedly
##    (m_i = 1 for the others) are allowed; w = optional frequency weights.
## -----------------------------------------------------------------------------
fla_nll <- function(theta, X, S, m, w = NULL) {
  w  <- w %||% rep(1, nrow(X))
  p  <- ncol(X); b <- theta[1:p]; e1 <- theta[p + 1]; e2 <- theta[p + 2]
  if (e1 <= 0 || e2 <= 0 || e1 + e2 >= 1) return(1e10)
  F  <- expit(drop(X %*% b))
  -sum(w * log(F * dbinom(S, m, 1 - e1) + (1 - F) * dbinom(S, m, e2)))
}

## score of fla_nll with respect to (beta, e1, e2)
.fla_grad <- function(theta, X, S, m, w) {
  p <- ncol(X); b <- theta[1:p]; e1 <- theta[p + 1]; e2 <- theta[p + 2]
  F  <- expit(drop(X %*% b))
  b1 <- dbinom(S, m, 1 - e1); b0 <- dbinom(S, m, e2)
  P  <- F * b1 + (1 - F) * b0
  c(-colSums(X * (w * F * (1 - F) * (b1 - b0) / P)),
    -sum(w * F * b1 * ((m - S) / e1 - S / (1 - e1)) / P),
    -sum(w * (1 - F) * b0 * (S / e2 - (m - S) / (1 - e2)) / P))
}

## starting values: beta = 0 and the (e1, e2) pairs below, independent of the SM
fla_starts <- rbind(c(0.05, 0.05), c(0.10, 0.10), c(0.20, 0.20),
                    c(0.05, 0.20), c(0.20, 0.05))

fla_fit <- function(X, S, m, w = NULL, start_beta = NULL, starts = fla_starts,
                    warn = TRUE) {
  w <- w %||% rep(1, nrow(X))
  p <- ncol(X)
  if (warn && all(m == 1))
    warning("m = 1 for all units: e1 and e2 are identified only through the ",
            "logistic form of the model; estimates may be severely biased and ",
            "highly variable. Classify at least part of the units repeatedly.")
  start_beta <- start_beta %||% rep(0, p)
  ## optimise on an unconstrained scale: e_k = 0.5 * expit(eta_k), e_k in (0, 0.5)
  tr <- function(t) c(t[1:p], 0.5 * expit(t[p + 1]), 0.5 * expit(t[p + 2]))
  f  <- function(t) fla_nll(tr(t), X, S, m, w)
  gr <- function(t) {
    th <- tr(t); g <- .fla_grad(th, X, S, m, w)
    c(g[1:p], g[p + 1] * th[p + 1] * (1 - 2 * th[p + 1]),
      g[p + 2] * th[p + 2] * (1 - 2 * th[p + 2]))
  }
  best <- NULL
  for (k in seq_len(nrow(starts))) {
    st <- c(start_beta, qlogis(2 * starts[k, 1]), qlogis(2 * starts[k, 2]))
    r  <- try(optim(st, f, gr, method = "BFGS",
                    control = list(maxit = 2000, reltol = 1e-12)), silent = TRUE)
    if (!inherits(r, "try-error") && is.finite(r$value) &&
        (is.null(best) || r$value < best$value)) best <- r
  }
  if (is.null(best)) return(list(coef = rep(NA, p + 2), se = rep(NA, p + 2), ok = FALSE))
  theta <- tr(best$par)
  names(theta) <- c(colnames(X) %||% paste0("b", 1:p - 1), "e1", "e2")
  ## observed information on the original scale
  H <- try(suppressWarnings(optimHess(theta, function(th) fla_nll(th, X, S, m, w),
                     function(th) .fla_grad(th, X, S, m, w))), silent = TRUE)
  V <- if (inherits(H, "try-error") || any(!is.finite(H))) NULL else
       tryCatch(solve(H), error = function(e) NULL)
  ok <- best$convergence == 0 && !is.null(V) && all(is.finite(diag(V))) &&
        all(diag(V) > 0) && all(theta[p + 1:2] > 1e-4)
  list(coef = theta, se = if (ok) sqrt(diag(V)) else rep(NA, p + 2),
       vcov = if (ok) V else NULL, ok = ok, loglik = -best$value)
}
