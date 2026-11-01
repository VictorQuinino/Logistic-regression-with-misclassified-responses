## =============================================================================
## Functions used by the programs of the paper
## "Logistic regression with misclassified binary responses: a simple majority
## approach based on repeated classifications".
## The SM is in sm_fit.R and the FLA in fla_fit.R, both loaded below; this file
## adds the binomial regression of S, the asymptotic quantities and the
## summaries of the finite-sample check. Base R, plus lpSolve for the
## detection of separation.
## This file is loaded by the other programs; it does not need to be run.
## =============================================================================
source("sm_fit.R")
source("fla_fit.R")

## -----------------------------------------------------------------------------
## 1. Binomial logistic regression of S on X with m trials (Section S2)
## -----------------------------------------------------------------------------
binom_fit <- function(X, S, m) {
  fit <- suppressWarnings(glm.fit(X, cbind(S, m - S), family = binomial()))
  b   <- fit$coefficients
  p   <- expit(drop(X %*% b))
  V   <- tryCatch(solve(crossprod(X * (m * p * (1 - p)), X)), error = function(e) NULL)
  list(coef = b, se = if (is.null(V)) rep(NA, ncol(X)) else sqrt(diag(V)),
       ok = fit$converged && !is.null(V))
}

## -----------------------------------------------------------------------------
## 2. Covariate distribution G as a quadrature rule (no simulation)
## -----------------------------------------------------------------------------
.gauss_hermite <- function(K) {          # weight: standard normal density
  J <- matrix(0, K, K)
  for (k in 1:(K - 1)) J[k, k + 1] <- J[k + 1, k] <- sqrt(k)
  ev <- eigen(J, symmetric = TRUE)
  list(x = ev$values, w = ev$vectors[1, ]^2)
}
.gauss_legendre01 <- function(K) {       # weight: uniform on (0,1)
  J <- matrix(0, K, K)
  for (k in 1:(K - 1)) J[k, k + 1] <- J[k + 1, k] <- k / sqrt(4 * k^2 - 1)
  ev <- eigen(J, symmetric = TRUE)
  list(x = (ev$values + 1) / 2, w = ev$vectors[1, ]^2)
}

## Design of Hausman et al. (1998): x1 ~ lognormal(0,1), x2 ~ Bernoulli(1/3),
## x3 ~ U(0,1), independent ("indep"), or log(x1) and qnorm(x3) bivariate normal
## with correlation rho ("corr").
covariate_grid <- function(type = c("indep", "corr"), rho = 0.5, K = 60) {
  type <- match.arg(type)
  gh <- .gauss_hermite(K)
  x2 <- c(0, 1); w2 <- c(2 / 3, 1 / 3)
  if (type == "indep") {
    gl <- .gauss_legendre01(30)
    g  <- expand.grid(i1 = seq_along(gh$x), i2 = 1:2, i3 = seq_along(gl$x))
    X  <- cbind(1, exp(gh$x[g$i1]), x2[g$i2], gl$x[g$i3])
    w  <- gh$w[g$i1] * w2[g$i2] * gl$w[g$i3]
  } else {
    g  <- expand.grid(i1 = seq_along(gh$x), i2 = 1:2, iu = seq_along(gh$x))
    z1 <- gh$x[g$i1]; z3 <- rho * z1 + sqrt(1 - rho^2) * gh$x[g$iu]
    X  <- cbind(1, exp(z1), x2[g$i2], pnorm(z3))
    w  <- gh$w[g$i1] * w2[g$i2] * gh$w[g$iu]
  }
  colnames(X) <- c("(Intercept)", "x1", "x2", "x3")
  list(X = X, w = w / sum(w))
}

## -----------------------------------------------------------------------------
## 3. Asymptotic quantities (Section 3)
## -----------------------------------------------------------------------------
## Pseudo-true value of a logistic fit to a binary response with success
## probability pi(x), and its sandwich covariance (Section 3.1)
.pseudo_true <- function(X, w, pi, b) {
  for (it in 1:200) {                     # Newton-Raphson for E[x (pi - F*)] = 0
    Fs <- expit(drop(X %*% b))
    U  <- colSums(X * (w * (pi - Fs)))
    A  <- crossprod(X * (w * Fs * (1 - Fs)), X)
    db <- solve(A, U); b <- b + db
    if (max(abs(db)) < 1e-12) break
  }
  Fs <- expit(drop(X %*% b))
  A  <- crossprod(X * (w * Fs * (1 - Fs)), X)
  B  <- crossprod(X * (w * (pi * (1 - pi) + (pi - Fs)^2)), X)
  Ai <- solve(A)
  list(bstar = b, V = Ai %*% B %*% Ai, Ainv = Ai)
}

## SM: misclassification of the majority (a1, a2)
asym_sm <- function(beta, m, e1, e2, G) {
  a  <- a_probs(m, e1, e2)
  F  <- expit(drop(G$X %*% beta))
  r  <- .pseudo_true(G$X, G$w, F * (1 - a[1]) + (1 - F) * a[2], beta)
  r$a <- a; r
}

## Binomial regression of S/m (Section S2): pseudo-true value with (e1, e2)
asym_binom <- function(beta, e1, e2, G) {
  F <- expit(drop(G$X %*% beta))
  .pseudo_true(G$X, G$w, F * (1 - e1) + (1 - F) * e2, beta)
}

## SM fitted to the unanimous units (Section 4.5). Returns the pseudo-true value
## (intercept shifted by log(q1/q0)), the covariance per selected unit and the
## expected proportion of unanimous units.
asym_unan <- function(beta, m, e1, e2, G) {
  q1 <- (1 - e1)^m + e1^m; q0 <- (1 - e2)^m + e2^m
  u1 <- e1^m / q1;         u2 <- e2^m / q0
  F  <- expit(drop(G$X %*% beta))
  sel <- F * q1 + (1 - F) * q0            # P(unanimous | x)
  wu  <- G$w * sel; pu <- sum(wu); wu <- wu / pu
  Fu  <- F * q1 / sel                     # P(true = 1 | x, unanimous)
  b0  <- beta; b0[1] <- b0[1] + log(q1 / q0)
  r <- .pseudo_true(G$X, wu, Fu * (1 - u1) + (1 - Fu) * u2, b0)
  r$target <- b0; r$prop <- pu; r$u <- c(u1 = u1, u2 = u2); r
}

## FLA: Fisher information per unit (Section S1)
asym_fla <- function(beta, m, e1, e2, G) {
  X <- G$X; w <- G$w; p <- ncol(X)
  F <- expit(drop(X %*% beta))
  I <- matrix(0, p + 2, p + 2)
  for (s in 0:m) {
    b1 <- dbinom(s, m, 1 - e1); b0 <- dbinom(s, m, e2)
    db1 <- b1 * ((m - s) / e1 - s / (1 - e1))
    db0 <- b0 * (s / e2 - (m - s) / (1 - e2))
    P  <- F * b1 + (1 - F) * b0
    g  <- cbind(X * (F * (1 - F) * (b1 - b0)), F * db1, (1 - F) * db0)
    I  <- I + crossprod(g * (w / P), g)
  }
  Vall <- solve(I)
  list(V = Vall[1:p, 1:p], Vall = Vall, I = I)
}

## Reference without misclassification (WCE)
asym_wce <- function(beta, G) {
  F <- expit(drop(G$X %*% beta))
  solve(crossprod(G$X * (G$w * F * (1 - F)), G$X))
}

## Coverage of a nominal 95% Wald interval with bias b and standard error se
cover <- function(b, se, z = qnorm(0.975)) pnorm(z - b / se) - pnorm(-z - b / se)

## All comparison measures for one scenario (Section 3.3)
compare_scenario <- function(beta, m, e1, e2, n, G, z = qnorm(0.975), size = TRUE) {
  sm <- asym_sm(beta, m, e1, e2, G); fl <- asym_fla(beta, m, e1, e2, G)
  vs <- diag(sm$V); vt <- diag(fl$V); b <- sm$bstar - beta
  se <- sqrt(vs / n); delta <- b / se; rho <- sqrt(diag(sm$Ainv) / vs)
  rej <- function(d) pnorm(-z + d) + pnorm(-z - d)
  sz <- if (size) sapply(2:length(beta), function(j) {   # size of test of beta_j = 0
    b0 <- beta; b0[j] <- 0
    s0 <- asym_sm(b0, m, e1, e2, G)
    rej(s0$bstar[j] / sqrt(s0$V[j, j] / n))
  }) else rep(NA, length(beta) - 1)
  data.frame(coef = names(beta) %||% paste0("b", seq_along(beta) - 1),
             m = m, e1 = e1, e2 = e2, n = n,
             a1 = sm$a[1], a2 = sm$a[2],
             bstar = sm$bstar, bias = b, rel_bias = b / beta,
             var_SM = vs, var_FLA = vt, var_WCE = diag(asym_wce(beta, G)),
             mse_SM = b^2 + vs / n, mse_FLA = vt / n,
             mse_ratio = (b^2 + vs / n) / (vt / n),
             nstar = ifelse(vt > vs, (vt - vs) / b^2, 0),
             rho = rho,
             cov_sandwich = pnorm(z - delta) - pnorm(-z - delta),
             cov_usual = pnorm(z * rho - delta) - pnorm(-z * rho - delta),
             power_SM = rej(sm$bstar / se),
             power_FLA = rej(beta / sqrt(vt / n)),
             size_SM = c(NA, sz),
             row.names = NULL)
}

## -----------------------------------------------------------------------------
## 4. Sequential classification (Section 4.4): expected number of
##    classifications until (m + 1)/2 concordant ones, for a unit whose
##    classifications are correct with probability 1 - e
## -----------------------------------------------------------------------------
expected_sequential <- function(m, e) {
  k <- (m + 1) / 2; q <- 1 - e
  ## the process stops at classification N = k + j, j = 0..k-1, when one of the
  ## two outcomes reaches k for the first time
  N <- k + 0:(k - 1)
  pN <- sapply(0:(k - 1), function(j)
    choose(k - 1 + j, j) * (q^k * e^j + e^k * q^j))
  sum(N * pN)
}

## -----------------------------------------------------------------------------
## 5. Summary measures of the finite-sample check for one scenario (Section 4.3).
## A fit is counted as a failure if flagged as such or if the estimate or its
## standard error is not finite. FLA summaries use the successful FLA fits;
## SM summaries use the fits without separation. The "_all" summaries keep the
## fits flagged as failures, using every finite estimate, as an analyst would.
## -----------------------------------------------------------------------------
summarize_finite <- function(s, beta, j, G, z = qnorm(0.975)) {
  m <- s$m[1]; e1 <- s$e1[1]; e2 <- s$e2[1]; n <- s$n[1]; b <- beta[j]
  fin <- function(...) Reduce(`&`, lapply(list(...), is.finite))
  ok  <- s$sm_ok == 1 & fin(s$sm_est, s$sm_seu, s$sm_ses)
  okf <- s$fl_ok %in% 1 & fin(s$fl_est, s$fl_se)
  oku <- s$un_ok == 1 & fin(s$un_est, s$un_se)
  okb <- s$bi_ok == 1 & fin(s$bi_est)
  ok0 <- s$sm0_ok == 1 & fin(s$sm0_tu, s$sm0_ts)
  do_fla <- any(!is.na(s$fl_ok))
  oka  <- is.finite(s$sm_est)
  okfa <- is.finite(s$fl_est)
  asy <- compare_scenario(beta, m, e1, e2, n, G, size = FALSE)[j, ]
  mse <- function(x) mean((x - b)^2)
  cv  <- function(est, se) mean(abs(est - b) <= z * se)
  data.frame(
    m = m, e1 = e1, e2 = e2, n = n, R = nrow(s),
    mse_ratio_asy = asy$mse_ratio,
    mse_ratio_fs  = if (do_fla) mse(s$sm_est[ok]) / mse(s$fl_est[okf]) else NA,
    ## the same, keeping the fits flagged as failures (all finite estimates)
    mse_ratio_fs_all = if (do_fla) mse(s$sm_est[oka]) / mse(s$fl_est[okfa]) else NA,
    bias_fla_all = if (do_fla) mean(s$fl_est[okfa]) - b else NA,
    mae_sm_all  = median(abs(s$sm_est[oka] - b)),
    mae_fla_all = if (do_fla) median(abs(s$fl_est[okfa] - b)) else NA,
    max_abs_fla_failed = if (do_fla && any(okfa & !okf)) max(abs(s$fl_est[okfa & !okf])) else NA,
    bias_sm = mean(s$sm_est[ok]) - b,
    bias_fla = if (do_fla) mean(s$fl_est[okf]) - b else NA,
    mae_sm  = median(abs(s$sm_est[ok] - b)),
    mae_fla = if (do_fla) median(abs(s$fl_est[okf] - b)) else NA,
    fail_sm = 100 * mean(!ok), fail_fla = if (do_fla) 100 * mean(!okf) else NA,
    cov_sm_asy = asy$cov_sandwich,
    cov_sm  = cv(s$sm_est[ok], s$sm_seu[ok]), cov_sm_sand = cv(s$sm_est[ok], s$sm_ses[ok]),
    cov_fla = if (do_fla) cv(s$fl_est[okf], s$fl_se[okf]) else NA,
    size_sm = mean(abs(s$sm0_tu[ok0]) > z),
    size_sm_sand = mean(abs(s$sm0_ts[ok0]) > z),
    size_fla = if (do_fla) mean(abs(s$fl0_t[s$fl0_ok %in% 1 & is.finite(s$fl0_t)]) > z) else NA,
    pow_sm  = mean(abs(s$sm_est[ok] / s$sm_seu[ok]) > z),
    pow_fla = if (do_fla) mean(abs(s$fl_est[okf] / s$fl_se[okf]) > z) else NA,
    cov_unan = cv(s$un_est[oku], s$un_se[oku]),
    width_ratio_unan = median(s$un_se[oku]) / median(s$sm_seu[ok]),
    fail_unan = 100 * mean(!oku),
    mse_ratio_binom_sm = mse(s$bi_est[okb]) / mse(s$sm_est[ok]),
    size_binom = mean(abs(s$bi0_t[s$bi0_ok == 1 & is.finite(s$bi0_t)]) > z))
}

## -----------------------------------------------------------------------------
## 6. Replicates of the finite-sample check with reproducible random streams.
##    Replicate r uses stream k = 1 + (r - 1) %% 2 of the L'Ecuyer-CMRG
##    generator, derived from the current seed, and the replicates of each
##    stream are run in order. This reproduces the streams that
##    parallel::mclapply(..., mc.cores = 2, mc.set.seed = TRUE) assigns, so the
##    results are the same on any operating system and number of processors.
##    On Linux and macOS the two streams run in parallel.
## -----------------------------------------------------------------------------
run_replicates <- function(R, f, streams = 2) {
  seed0 <- get(".Random.seed", envir = .GlobalEnv)
  seeds <- vector("list", streams); s <- seed0
  for (k in seq_len(streams)) { s <- parallel::nextRNGStream(s); seeds[[k]] <- s }
  one_stream <- function(k) {
    assign(".Random.seed", seeds[[k]], envir = .GlobalEnv)
    lapply(seq(k, R, by = streams), f)
  }
  res <- if (.Platform$OS.type == "unix")
    parallel::mclapply(seq_len(streams), one_stream, mc.cores = streams,
                       mc.set.seed = FALSE)
  else lapply(seq_len(streams), one_stream)
  assign(".Random.seed", seed0, envir = .GlobalEnv)
  out <- vector("list", R)
  for (k in seq_len(streams)) out[seq(k, R, by = streams)] <- res[[k]]
  out
}
