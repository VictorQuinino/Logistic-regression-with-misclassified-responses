## =============================================================================
## Section 4.4 (Number of classifications and sequential classification):
## Table 3, and Table S4 of the Supplementary Material.
## No simulation. Running time: under 1 minute.
## Output: results/table3.csv (Table 3, minimum asymptotic coverage of the SM,
##         in %), results/tableS3.csv (Table S4, expected number of
##         classifications per unit under sequential classification),
##         results/tables_section4_4.txt (LaTeX rows)
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

beta  <- c(b0 = -1, b1 = 0.2, b2 = 1.5, b3 = -0.6)
egrid <- c(0.05, 0.10, 0.15, 0.20)
slope <- c("b1", "b2", "b3")
G     <- covariate_grid("indep")

## ---- Table 3: minimum coverage over (e1, e2) and slopes, by m and n ---------
cov_m <- do.call(rbind, lapply(c(3, 5, 7, 9, 11), function(m)
  do.call(rbind, lapply(egrid, function(e1) do.call(rbind, lapply(egrid, function(e2)
    do.call(rbind, lapply(c(100, 500, 1000), function(n)
      compare_scenario(beta, m, e1, e2, n, G, size = FALSE)))))))))
tab3 <- with(subset(cov_m, coef %in% slope), tapply(cov_sandwich, list(m = m, n = n), min))
write.csv(round(100 * tab3, 1), "results/table3.csv")
print(round(100 * tab3, 1))

## ---- Table S4: expected number of classifications, sequential --------------
tabS3 <- t(sapply(c(0.10, 0.20), function(e)
  sapply(c(3, 5, 7, 9, 11), expected_sequential, e = e)))
dimnames(tabS3) <- list(e = c(0.10, 0.20), m = c(3, 5, 7, 9, 11))
write.csv(round(tabS3, 2), "results/tableS3.csv")
print(round(tabS3, 2))

## ---- LaTeX rows ----------------------------------------------------------------
t3 <- round(100 * tab3, 1)
out <- c("% Table 3",
  paste0(rownames(t3), " & ", apply(t3, 1, function(r)
    paste(formatC(r, format = "f", digits = 1), collapse = " & ")), "\\\\"), "",
  "% Table S4",
  paste0(formatC(as.numeric(rownames(tabS3)), format = "f", digits = 2), " & ",
    apply(round(tabS3, 2), 1, function(r) paste(formatC(r, format = "f", digits = 2),
                                                collapse = " & ")), "\\\\"))
writeLines(out, "results/tables_section4_4.txt")
cat("Section 4.4 done: see results/\n")
