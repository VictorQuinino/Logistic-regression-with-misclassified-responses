## =============================================================================
## Section 4.5 (Confidence intervals from the unanimous classifications):
## asymptotic results quoted in the text and Table S5 of the Supplementary
## Material (finite sample).
## Table S5 uses the replicates saved by Section4_3_finite_sample.R; run that
## program first. Running time: about 1 minute.
## Output: results/analytic_unanimous.csv, results/numbers_section4_5.txt,
##         results/tables_section4_5.txt (LaTeX rows of Table S5)
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

beta  <- c(b0 = -1, b1 = 0.2, b2 = 1.5, b3 = -0.6)
egrid <- c(0.05, 0.10, 0.15, 0.20)
slope <- c("b1", "b2", "b3")
G     <- covariate_grid("indep")

## ---- asymptotic results for the unanimous units ------------------------------
unan <- do.call(rbind, lapply(c(3, 5), function(m)
  do.call(rbind, lapply(egrid, function(e1) do.call(rbind, lapply(egrid, function(e2) {
    u <- asym_unan(beta, m, e1, e2, G)
    do.call(rbind, lapply(c(500, 1000), function(n) {
      j <- 2:4; nu <- n * u$prop
      data.frame(m = m, e1 = e1, e2 = e2, n = n, coef = slope, prop = u$prop,
                 u1 = unname(u$u[1]), u2 = unname(u$u[2]),
                 shift = log(((1 - e1)^m + e1^m) / ((1 - e2)^m + e2^m)),
                 cov = cover(u$bstar[j] - u$target[j], sqrt(diag(u$V)[j] / nu)))
    }))
  }))))))
write.csv(unan, "results/analytic_unanimous.csv", row.names = FALSE)

## ---- finite sample (Table S5), from the saved replicates -----------------------
reps <- readRDS(rds)
d <- do.call(rbind, lapply(split(reps, list(reps$m, reps$e1, reps$e2, reps$n), drop = TRUE),
  summarize_finite, beta = beta, j = 3, G = G))
d <- d[order(d$m, d$e1, d$e2, d$n), ]; rownames(d) <- NULL

pc <- function(x) formatC(100 * x, format = "f", digits = 1)
p1 <- function(x) formatC(x, format = "f", digits = 1)
f2 <- function(x) ifelse(abs(x) < 0.1, formatC(x, format = "f", digits = 3),
                         formatC(x, format = "f", digits = 2))
e2 <- function(a, b) sprintf("(%.2f,%.2f)", a, b)
rows <- with(d, paste0(m, " & ", e2(e1, e2), " & ", n, " & ", pc(cov_sm), " & ",
  pc(cov_unan), " & ", f2(width_ratio_unan), " & ", p1(fail_unan), "\\\\"))
k <- seq(4, length(rows) - 1, by = 4); rows[k] <- paste0(rows[k], "\n\\addlinespace")
writeLines(c("% Table S5", rows), "results/tables_section4_5.txt")

## ---- numbers quoted in Section 4.5 ----------------------------------------
s <- function(x) paste(range(x, na.rm = TRUE), collapse = " to ")
summ <- list(
  "proportion of unanimous units, m = 3, e1 = e2 = 0.05 and 0.20" =
    unique(subset(unan, m == 3 & e1 == e2 & e1 %in% c(.05, .2))[, c("e1", "prop", "u1")]),
  "intercept shift, (0.05,0.20), m = 3" =
    unique(subset(unan, m == 3 & e1 == .05 & e2 == .2)$shift),
  "min asymptotic coverage of slopes, unanimous units, n = 500 and 1000" = min(unan$cov),
  "finite sample: unanimous coverage (%), n>=500" = s(100 * subset(d, n >= 500)$cov_unan),
  "finite sample: SM coverage (%), minimum for n>=500" = min(100 * subset(d, n >= 500)$cov_sm),
  "finite sample: width ratio unanimous/SM, n>=500" = s(subset(d, n >= 500)$width_ratio_unan),
  "finite sample: width ratio unanimous/SM, n=50" = s(subset(d, n == 50)$width_ratio_unan),
  "finite sample: unanimous failures (%), n=50" = s(subset(d, n == 50)$fail_unan))
sink("results/numbers_section4_5.txt"); print(summ, digits = 4); sink()
print(summ, digits = 4)
cat("Section 4.5 done: see results/\n")
