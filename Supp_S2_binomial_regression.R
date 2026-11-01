## =============================================================================
## Supplementary Section S2 (Binomial regression of the number of positive
## classifications): asymptotic relative bias and finite-sample comparison
## with the SM.
## The finite-sample part uses the replicates saved by
## Section4_3_finite_sample.R; run that program first. Running time: under
## 1 minute.
## Output: results/numbers_supp_S2.txt
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

rds <- "results/finite_sample_reps.rds"
if (!file.exists(rds))
  stop("Run Section4_3_finite_sample.R first: it saves the replicates used here.")

beta <- c(b0 = -1, b1 = 0.2, b2 = 1.5, b3 = -0.6)
G    <- covariate_grid("indep")
ev   <- c(0.05, 0.10, 0.20)

## ---- asymptotic relative bias of beta2: binomial regression and SM (m = 3) -----
binom_rb <- sapply(ev, function(e) (asym_binom(beta, e, e, G)$bstar[3] - beta[3]) / beta[3])
sm_rb    <- sapply(ev, function(e)
  subset(compare_scenario(beta, 3, e, e, 50, G), coef == "b2")$rel_bias)

## ---- finite sample, from the saved replicates ------------------------------
reps <- readRDS(rds)
d <- do.call(rbind, lapply(split(reps, list(reps$m, reps$e1, reps$e2, reps$n), drop = TRUE),
  summarize_finite, beta = beta, j = 3, G = G))
d5 <- subset(d[order(d$m, d$e1, d$e2, d$n), ], n <= 500)
s  <- function(x) paste(range(x, na.rm = TRUE), collapse = " to ")
e2 <- function(a, b) sprintf("(%.2f,%.2f)", a, b)

summ <- list(
  "rel. bias of beta2, binomial regression, e = 0.05, 0.10, 0.20" = binom_rb,
  "rel. bias of beta2, SM with m = 3, e = 0.05, 0.10, 0.20" = sm_rb,
  "scenarios with MSE(binomial) < MSE(SM)" =
    with(subset(d5, mse_ratio_binom_sm < 1), paste0("m=", m, ",", e2(e1, e2), ",n=", n)),
  "size of the test, binomial regression (%)" = s(100 * d5$size_binom),
  "size of the test, SM (%)" = s(100 * d5$size_sm))
sink("results/numbers_supp_S2.txt"); print(summ, digits = 4); sink()
print(summ, digits = 4)
cat("Section S2 done: see results/\n")
