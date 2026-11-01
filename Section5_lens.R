## =============================================================================
## Section 5 (Application to lens inspection data): Tables 6 and 7, and
## Table S7 of the Supplementary Material, with the numbers quoted in
## Section 5. The data are those of Table 5 of the paper.
## Running time: about 2 minutes (parametric bootstrap of the goodness of fit).
## Output: results/lens_results.txt
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
## Number of lenses with s = 0..5 'perfect' classifications, by (L, R)
counts <- rbind("L0R0" = c(13, 11, 5,  3, 20, 48),
                "L0R1" = c( 8,  9, 4,  9, 27, 43),
                "L1R0" = c( 5,  2, 2,  5, 31, 55),
                "L1R1" = c( 1,  1, 2, 11, 34, 51))
m   <- 5
grp <- data.frame(L = c(0, 0, 1, 1), R = c(0, 1, 0, 1), row.names = rownames(counts))
lens <- do.call(rbind, lapply(rownames(counts), function(g)
  data.frame(L = grp[g, "L"], R = grp[g, "R"], S = rep(0:m, counts[g, ]))))
X <- model.matrix(~ L + R, lens); S <- lens$S

## ---- SM: the analyst needs only these two lines ------------------------------
lens$Y <- as.integer(lens$S >= 3)                       # majority of 5
fit    <- glm(Y ~ L + R, family = binomial, data = lens)
print(summary(fit))
sm <- sm_fit(X, S, m)                                   # adds sandwich SEs

## ---- SM restricted to the units with unanimous classifications (Section 4.5) --
lens$U <- as.integer(lens$S %in% c(0, m))
fitU   <- glm(Y ~ L + R, family = binomial, data = lens, subset = U == 1)
un     <- unan_fit(X, S, m)

## ---- moment estimators of e1, e2 and probabilities that the majority errs ----
me   <- moment_e(S, m)
ahat <- a_probs(m, me$e[1], me$e[2])

## ---- FLA ----------------------------------------------------------------------
fl <- fla_fit(X, S, m)
pub <- c(0.853, 1.584, 0.531, fl$coef[4:5])             # Valadares et al. (2021)

## ---- Table 6: estimates (standard errors) --------------------------------------
f3 <- function(est, se) sprintf("%.3f (%.3f)", est, se)
tab6 <- data.frame(
  parameter   = c("beta0", "beta1 (L)", "beta2 (R)", "e1", "e2"),
  SM          = c(f3(sm$coef, sm$se_usual),    f3(me$e, me$se)),
  SM_sandwich = c(f3(sm$coef, sm$se_sandwich), f3(me$e, me$se)),
  FLA         = f3(fl$coef, fl$se))

## ---- Table 7: odds ratios, 95% Wald intervals and p-values -------------------
orow <- function(method, b, se, k) {
  data.frame(effect = c("Laboratory", "Region")[k - 1], method = method,
             OR = sprintf("%.2f", exp(b[k])),
             CI = sprintf("(%.2f; %.2f)", exp(b[k] - 1.96 * se[k]), exp(b[k] + 1.96 * se[k])),
             p  = sprintf("%.3f", 2 * pnorm(-abs(b[k] / se[k]))))
}
tab7 <- do.call(rbind, lapply(2:3, function(k) rbind(
  orow("SM", sm$coef, sm$se_usual, k), orow("SM (unan.)", un$coef, un$se_usual, k),
  orow("FLA", fl$coef[1:3], fl$se[1:3], k))))

## ---- Table S7: P(perfect lens) with 95% CI (logit scale) ----------------------
Xg <- cbind(1, grp$L, grp$R)
ci <- function(b, V) {
  eta <- drop(Xg %*% b); se <- sqrt(rowSums((Xg %*% V) * Xg))
  sprintf("%.3f (%.3f; %.3f)", expit(eta), expit(eta - 1.96 * se), expit(eta + 1.96 * se))
}
tabS7 <- data.frame(group = rownames(grp), SM = ci(sm$coef, sm$vcov_usual),
                    FLA = ci(fl$coef[1:3], fl$vcov[1:3, 1:3]))

## ---- goodness of fit of the binomial mixture: Pearson X^2, 15 df -------------
pearson <- function(theta, cnt) {
  F <- expit(drop(Xg %*% theta[1:3]))
  E <- t(sapply(1:4, function(g) sum(cnt[g, ]) *
         (F[g] * dbinom(0:m, m, 1 - theta[4]) + (1 - F[g]) * dbinom(0:m, m, theta[5]))))
  sum((cnt - E)^2 / E)
}
X2 <- pearson(fl$coef, counts)
## parametric bootstrap p-value
set.seed(20261018); B <- 1000
X2b <- replicate(B, {
  F <- expit(drop(Xg %*% fl$coef[1:3]))
  cb <- t(sapply(1:4, function(g) {
    Yt <- rbinom(sum(counts[g, ]), 1, F[g])
    Sb <- ifelse(Yt == 1, rbinom(length(Yt), m, 1 - fl$coef[4]), rbinom(length(Yt), m, fl$coef[5]))
    tabulate(Sb + 1, nbins = m + 1) }))
  Sb <- rep(rep(0:m, 4), as.vector(t(cb)))
  Xb <- Xg[rep(rep(1:4, each = m + 1), as.vector(t(cb))), , drop = FALSE]
  fb <- fla_fit(Xb, Sb, m)
  if (isTRUE(fb$ok)) pearson(fb$coef, cb) else NA })
p_boot <- mean(X2b >= X2, na.rm = TRUE)

## ---- asymptotic comparison for this design (Section 5) -----------------------
Gl <- list(X = Xg, w = rep(1 / 4, 4))
bl <- setNames(fl$coef[1:3], c("b0", "b1", "b2"))
asl <- compare_scenario(bl, m, fl$coef[4], fl$coef[5], 400, Gl, size = FALSE)

## ---- output ----------------------------------------------------------------------
sink("results/lens_results.txt")
cat("Table 6\n"); print(tab6)
cat("\nTable 7\n"); print(tab7)
cat("\nTable S7\n"); print(tabS7)
cat(sprintf("\nUnanimous lenses: %d\n", un$n_unan))
cat(sprintf("a1hat = %.3f, a2hat = %.3f\n", ahat[1], ahat[2]))
cat(sprintf("SM/FLA estimate of the laboratory effect: %.3f\n", sm$coef[2] / fl$coef[2]))
cat(sprintf("Interval widths, SM/FLA (log OR): L %.3f, R %.3f\n",
            sm$se_usual[2] / fl$se[2], sm$se_usual[3] / fl$se[3]))
cat(sprintf("Interval widths, unanimous/SM (log OR): L %.3f, R %.3f\n",
            un$se_usual[2] / sm$se_usual[2], un$se_usual[3] / sm$se_usual[3]))
cat(sprintf("FLA log-likelihood at ML: %.2f; at published beta: %.2f\n",
            fl$loglik, -fla_nll(pub, X, S, m)))
cat(sprintf("Pearson X2 = %.1f (15 df); parametric bootstrap p = %.2f (B = %d, %d valid)\n",
            X2, p_boot, B, sum(!is.na(X2b))))
cat("Max. difference between estimated probabilities (SM vs FLA):",
    sprintf("%.3f", max(abs(expit(drop(Xg %*% sm$coef)) - expit(drop(Xg %*% fl$coef[1:3]))))), "\n")
cat("\nAsymptotic comparison for this design (n = 400), FLA estimates as true values\n")
print(asl[, c("coef", "mse_ratio", "nstar", "cov_sandwich")], digits = 3)
sink()
cat(readLines("results/lens_results.txt"), sep = "\n")
