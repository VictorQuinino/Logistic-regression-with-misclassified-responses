## =============================================================================
## Supplementary Section S4 (Fitting the FLA): example of fla_fit. With a
## single classification per unit the FLA estimates may be severely biased and
## highly variable; classifying part of the units repeatedly removes most of
## the problem. Running time: about 2 minutes.
## Output: printed on the screen
## =============================================================================
## ---- Start: use the folder of this program, which must also contain
##      functions.R and fla_fit.R. Works with Rscript and with "Source" in
##      RStudio; otherwise set the working directory to this folder. ---------
here <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  f <- if (length(a)) sub("^--file=", "", a[1]) else
    tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(f) && requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) f <- rstudioapi::getSourceEditorContext()$path
  if (is.null(f) || !nzchar(f)) getwd() else dirname(normalizePath(f))
})
if (!file.exists(file.path(here, "functions.R")))
  stop("functions.R not found. Keep all programs in the same folder and set it ",
       "as the working directory (RStudio: Session > Set Working Directory > ",
       "To Source File Location).")
setwd(here)
source("functions.R")
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)
eps <- function(f) {                       # EPS copy of a PDF figure (Poppler)
  st <- suppressWarnings(tryCatch(system2("pdftops", c("-eps", file.path("figures",
          paste0(f, ".pdf")), file.path("figures", paste0(f, ".eps")))),
          error = function(e) 1))
  if (st != 0) message("pdftops not found: EPS copy of ", f, " not created")
}
## ----------------------------------------------------------------------------
set.seed(2026)
beta <- c(-1, 0.2, 1.5, -0.6); e1 <- 0.10; e2 <- 0.10; n <- 500; R <- 200

sim <- function(frac_rep, m_rep = 3) {
  X  <- cbind(1, rlnorm(n), rbinom(n, 1, 1/3), runif(n))
  Yt <- rbinom(n, 1, expit(drop(X %*% beta)))
  mi <- ifelse(seq_len(n) <= round(frac_rep * n), m_rep, 1)   # classifications per unit
  S  <- ifelse(Yt == 1, rbinom(n, mi, 1 - e1), rbinom(n, mi, e2))
  f  <- suppressWarnings(fla_fit(X, S, mi))
  if (isTRUE(f$ok)) f$coef[3] else NA
}
res <- sapply(c(0, 0.3), function(fr) replicate(R, sim(fr)))
colnames(res) <- c("m = 1 for all units", "m = 3 for 30% of units")
summ <- apply(res, 2, function(x) c(failures = mean(is.na(x)),
                                    bias = mean(x, na.rm = TRUE) - beta[3],
                                    sd = sd(x, na.rm = TRUE),
                                    median_abs_error = median(abs(x - beta[3]), na.rm = TRUE)))
print(round(summ, 3))

## Single data set: a warning is issued when m = 1 for all units
X <- cbind(1, rlnorm(n), rbinom(n, 1, 1/3), runif(n))
S <- rbinom(n, 1, ifelse(rbinom(n, 1, expit(drop(X %*% beta))) == 1, 1 - e1, e2))
fit <- fla_fit(X, S, 1)

## ---- Use of fla_fit on one data set with m = 3: estimates, standard errors,
##      Wald tests and 95% confidence intervals ---------------------------------
X   <- cbind("(Intercept)" = 1, x1 = rlnorm(n), x2 = rbinom(n, 1, 1/3), x3 = runif(n))
Yt  <- rbinom(n, 1, expit(drop(X %*% beta)))
S   <- ifelse(Yt == 1, rbinom(n, 3, 1 - e1), rbinom(n, 3, e2))
fit <- fla_fit(X, S, 3)
z   <- fit$coef / fit$se
print(round(cbind(estimate = fit$coef, se = fit$se, z = z, p = 2 * pnorm(-abs(z)),
                  lower95 = fit$coef - 1.96 * fit$se, upper95 = fit$coef + 1.96 * fit$se), 3))
