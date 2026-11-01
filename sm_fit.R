## =============================================================================
## sm_fit: simple majority (SM) method for logistic regression with
## misclassified binary responses and repeated classifications (Sections 2.3,
## 2.4 and 4.5 of the paper). Base R, plus the package lpSolve for the
## detection of separation.
##
## Usage (run this file once, or source("sm_fit.R"), then):
##   fit <- sm_fit(X, S, m)
##     X : model matrix, with a column of ones for the intercept,
##         e.g. X <- model.matrix(~ x1 + x2, data)
##     S : number of positive classifications of each unit (0, ..., m)
##     m : number of classifications per unit (odd, m >= 3)
##     w : optional frequency weights
##   fit$coef          estimates of beta (logistic regression on the majority)
##   fit$se_usual      usual standard errors, as in glm()
##   fit$se_sandwich   sandwich standard errors
##   fit$ok            FALSE if the fit did not converge or there is separation
## The same estimates and usual standard errors are given by
##   glm(as.integer(S > m / 2) ~ x1 + x2, family = binomial, data = data)
##
##   unan_fit(X, S, m)   the same fit restricted to the units with unanimous
##                       classifications (Section 4.5), for confidence
##                       intervals in large samples
##   moment_e(S, m)      estimates of e1 and e2 with standard errors
##                       (Section 2.4)
##   a_probs(m, e1, e2)  probabilities that the majority is wrong, a1 and a2;
##                       with the estimates of e1 and e2, the diagnostic of
##                       the adequacy of the SM (Section 2.4)
## See Section5_lens.R for an example with Wald tests and odds ratios.
## =============================================================================

expit <- function(z) 1 / (1 + exp(-z))
`%||%` <- function(a, b) if (is.null(a)) b else a

## -----------------------------------------------------------------------------
## 1. Detection of complete or quasi-complete separation by linear programming
##    (Konis, 2007): separation exists if some beta != 0 satisfies
##    (2y_i - 1) x_i'beta >= 0 for all i. We maximise sum_i (2y_i - 1) x_i'beta
##    under these constraints and |beta_j| <= 1; a positive optimum indicates
##    separation.
## -----------------------------------------------------------------------------
separation <- function(X, y) {
  s  <- 2 * y - 1
  A  <- X * s
  p  <- ncol(X)
  ## beta = bp - bm, bp, bm >= 0
  obj <- c(colSums(A), -colSums(A))
  con <- rbind(cbind(A, -A), cbind(diag(p), diag(p)))
  dir <- c(rep(">=", nrow(A)), rep("<=", p))
  rhs <- c(rep(0, nrow(A)), rep(1, p))
  r <- lpSolve::lp("max", obj, con, dir, rhs)
  r$status == 0 && r$objval > 1e-7
}

## -----------------------------------------------------------------------------
## 2. SM estimator with usual and sandwich standard errors (Section 2.3)
##    S = number of positive classifications (0..m); X = model matrix;
##    w = optional frequency weights
## -----------------------------------------------------------------------------
sm_fit <- function(X, S, m, w = NULL, check_separation = TRUE) {
  w   <- w %||% rep(1, nrow(X))
  y   <- as.integer(S > m / 2)                          # majority response
  sep <- if (check_separation) separation(X[w > 0, , drop = FALSE], y[w > 0]) else FALSE
  fit <- suppressWarnings(glm.fit(X, y, weights = w, family = binomial()))
  b   <- fit$coefficients
  p   <- expit(drop(X %*% b))
  A   <- crossprod(X * (w * p * (1 - p)), X)
  B   <- crossprod(X * (w * (y - p)^2), X)
  Ai  <- tryCatch(solve(A), error = function(e) NULL)
  ok  <- fit$converged && !sep && !is.null(Ai) && all(is.finite(b))
  if (is.null(Ai)) Ai <- matrix(NA, ncol(X), ncol(X))
  Vs  <- Ai %*% B %*% Ai
  list(coef = b, se_usual = sqrt(diag(Ai)), se_sandwich = sqrt(diag(Vs)),
       vcov_usual = Ai, vcov_sandwich = Vs, ok = ok, separation = sep)
}

## Same fit restricted to the units with unanimous classifications (Section 4.5)
unan_fit <- function(X, S, m, w = NULL) {
  keep <- S == 0 | S == m
  r <- sm_fit(X[keep, , drop = FALSE], S[keep], m, w = if (is.null(w)) NULL else w[keep])
  r$n_unan <- sum(if (is.null(w)) keep else w[keep])
  r
}

## -----------------------------------------------------------------------------
## 3. Moment estimators of e1, e2 (Section 2.4) and delta-method standard errors
## -----------------------------------------------------------------------------
.moment_map <- function(mu) {
  d  <- mu[2] - mu[1]^2
  c1 <- (mu[3] - mu[1] * mu[2]) / d
  c0 <- (mu[1] * mu[3] - mu[2]^2) / d
  r  <- sqrt(c1^2 - 4 * c0)
  list(e = c(e1 = 1 - (c1 + r) / 2, e2 = (c1 - r) / 2), c1 = c1, c0 = c0, d = d, r = r)
}

## Maximum likelihood of the two-component binomial mixture by the
## EM algorithm; used when the moment estimators are not available
mixture_em <- function(S, m, start = c(Fbar = 0.5, e1 = 0.1, e2 = 0.1),
                       tol = 1e-10, maxit = 5000) {
  Fb <- start[1]; e1 <- start[2]; e2 <- start[3]
  for (it in 1:maxit) {
    p1 <- Fb * dbinom(S, m, 1 - e1); p0 <- (1 - Fb) * dbinom(S, m, e2)
    tau <- p1 / (p1 + p0)                          # P(true = 1 | S)
    new <- c(mean(tau), 1 - sum(tau * S) / (m * sum(tau)),
             sum((1 - tau) * S) / (m * sum(1 - tau)))
    if (max(abs(new - c(Fb, e1, e2))) < tol) break
    Fb <- new[1]; e1 <- new[2]; e2 <- new[3]
  }
  th  <- c(Fbar = Fb, e1 = e1, e2 = e2)
  nll <- function(th) -sum(log(th[1] * dbinom(S, m, 1 - th[2]) +
                               (1 - th[1]) * dbinom(S, m, th[3])))
  H   <- try(suppressWarnings(optimHess(th, nll)), silent = TRUE)
  V   <- if (inherits(H, "try-error")) NULL else tryCatch(solve(H), error = function(e) NULL)
  se  <- if (is.null(V) || any(!is.finite(diag(V)) | diag(V) <= 0)) rep(NA, 3) else sqrt(diag(V))
  list(e = th[2:3], se = se[2:3], Fbar = unname(th[1]), iterations = it)
}

## Moment estimators (Section 2.4); if the discriminant is negative or the
## estimates fall outside (0, 0.5), the EM estimator above is used instead
## (method = "EM" in the output).
moment_e <- function(S, m, fallback = TRUE) {
  stopifnot(m >= 3)
  n  <- length(S)
  w  <- cbind(S / m,
              S * (S - 1) / (m * (m - 1)),
              S * (S - 1) * (S - 2) / (m * (m - 1) * (m - 2)))
  mu <- colMeans(w)
  mm <- suppressWarnings(.moment_map(mu))
  ok <- is.finite(mm$r) && all(mm$e > 0 & mm$e < 0.5)
  if (!ok) {
    if (!fallback) return(list(e = c(e1 = NA, e2 = NA), se = c(NA, NA), Fbar = NA,
                               ok = FALSE, method = "moments"))
    r <- mixture_em(S, m)
    return(c(r, list(ok = all(r$e > 0 & r$e < 0.5), method = "EM")))
  }
  ## gradients of c1 and c0 (Section 2.4) and of e1, e2
  gc1 <- c(2 * mu[1] * mm$c1 - mu[2], -mu[1] - mm$c1, 1) / mm$d
  gc0 <- c(mu[3] + 2 * mu[1] * mm$c0, -2 * mu[2] - mm$c0, mu[1]) / mm$d
  e   <- mm$e
  J   <- rbind((gc0 - (1 - e[1]) * gc1) / mm$r,
               (gc0 - e[2] * gc1) / mm$r)
  V   <- J %*% cov(w) %*% t(J) / n
  list(e = e, se = sqrt(diag(V)), vcov = V,
       Fbar = unname((mu[1] - e[2]) / (1 - e[1] - e[2])),
       ok = TRUE, method = "moments")
}

## -----------------------------------------------------------------------------
## 4. Probabilities that the majority of m classifications is wrong (Section 2.3)
## -----------------------------------------------------------------------------
a_probs <- function(m, e1, e2) {
  k <- (m - 1) / 2
  c(a1 = pbinom(k, m, 1 - e1),          # P(majority = 0 | true = 1)
    a2 = 1 - pbinom(k, m, e2))          # P(majority = 1 | true = 0)
}
