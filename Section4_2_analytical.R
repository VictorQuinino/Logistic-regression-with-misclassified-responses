## =============================================================================
## Section 4.2 (Analytical results): Figure 1, Table 1, and Figures S1 and S2
## of the Supplementary Material, with the numbers quoted in Section 4.2.
## No simulation: expectations over the covariate distribution by Gauss
## quadrature. Running time: about 1 minute.
## Output: figures/fig_mse_ratio.pdf (Figure 1), results/table1.csv (Table 1),
##         figures/fig_summary.pdf (Figure S1), figures/fig_coverage.pdf
##         (Figure S2), results/analytic_all.csv (all scenarios),
##         results/summary_by_n_m.csv, results/numbers_section4_2.txt
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
mgrid <- c(3, 5, 7)
ngrid <- c(50, 100, 500)
slope <- c("b1", "b2", "b3")

grid_run <- function(G, mg = mgrid, ng = ngrid, size = TRUE)
  do.call(rbind, lapply(mg, function(m)
    do.call(rbind, lapply(egrid, function(e1)
      do.call(rbind, lapply(egrid, function(e2)
        do.call(rbind, lapply(ng, function(n)
          compare_scenario(beta, m, e1, e2, n, G, size = size)))))))))

## ---- all scenarios: independent and correlated covariates ------------------
G    <- covariate_grid("indep")
res  <- grid_run(G);  res$design <- "indep"
Gc   <- covariate_grid("corr", rho = 0.5)
resc <- grid_run(Gc); resc$design <- "corr"
write.csv(rbind(res, resc), "results/analytic_all.csv", row.names = FALSE)

## ---- Table 1: a1 = a2 and n* for beta2, e1 = e2 = e -------------------------
tabA <- sapply(mgrid, function(m) sapply(egrid, function(e) a_probs(m, e, e)[1]))
sub  <- subset(res, coef == "b2" & e1 == e2 & n == 500)
tabN <- sapply(mgrid, function(m) sub$nstar[sub$m == m][order(sub$e1[sub$m == m])])
tab1 <- data.frame(e = egrid, a_m3 = tabA[, 1], a_m5 = tabA[, 2], a_m7 = tabA[, 3],
                   nstar_m3 = tabN[, 1], nstar_m5 = tabN[, 2], nstar_m7 = tabN[, 3])
write.csv(tab1, "results/table1.csv", row.names = FALSE)
print(tab1, digits = 3)

## ---- numbers quoted in Section 4.2 ----------------------------------------
sl   <- subset(res, coef %in% slope)
d500 <- subset(res, n == 500)
amax <- pmax(d500$a1, d500$a2)
summ <- list(
  "max |rel. bias|, e1,e2<=0.10, m>=5" =
    max(abs(subset(res, e1 <= .1 & e2 <= .1 & m >= 5)$rel_bias)),
  "max attenuation of slopes (rel. bias), by m" =
    aggregate(rel_bias ~ m, sl, min),
  "max |rel. bias| of intercept, (0.05,0.20), m = 3" =
    max(abs(subset(res, coef == "b0" & e1 == .05 & e2 == .2 & m == 3)$rel_bias)),
  "max MSE ratio by n and m" = aggregate(mse_ratio ~ n + m, res, max),
  "max MSE ratio, n = 500, max(a1,a2) <= 0.03" = max(d500$mse_ratio[amax <= 0.03]),
  "min max(a1,a2) among n = 500 scenarios with MSE ratio > 1.5" =
    min(amax[d500$mse_ratio > 1.5]),
  "n* for beta2, e1 = e2 (Table 1)" = tab1,
  "range of rho by m" = aggregate(rho ~ m, res, range),
  "max |1 - rho_j| under beta_j = 0 (independent covariate)" = max(sapply(
    c(2, 4), function(j) { b0 <- beta; b0[j] <- 0
      max(sapply(mgrid, function(m) max(sapply(egrid, function(e1) max(sapply(egrid,
        function(e2) { s <- asym_sm(b0, m, e1, e2, G)
          abs(1 - sqrt(s$Ainv[j, j] / s$V[j, j])) }))))))})),
  "min coverage (sandwich), slopes, e1,e2<=0.10" =
    min(subset(sl, e1 <= .1 & e2 <= .1)$cov_sandwich),
  "min coverage, slopes, n = 50" = min(subset(sl, n == 50)$cov_sandwich),
  "min coverage, slopes, by n and m" = aggregate(cov_sandwich ~ n + m, sl, min),
  "size (sandwich), independent design" = range(sl$size_SM),
  "size (sandwich), correlated design" =
    range(subset(resc, coef %in% slope)$size_SM),
  "max loss of power SM vs FLA (percentage points), by m" =
    aggregate(I(100 * (power_FLA - power_SM)) ~ m, sl, max))
sink("results/numbers_section4_2.txt"); print(summ, digits = 4); sink()

## ---- Figures (black and white; legends below the panels) --------------------
sty <- list(lty = c(b1 = 1, b2 = 2, b3 = 3), pch = c(b1 = 19, b2 = 2, b3 = 15))
diag_res <- res[res$e1 == res$e2, ]

bottom_legend <- function(labels, lty, pch, ncol = length(labels)) {
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", labels, lty = lty, pch = pch, lwd = 1.5, bty = "n",
         ncol = ncol, cex = 1.3, seg.len = 3, x.intersp = 0.8, text.width = 0.16,
         xpd = TRUE)
}

panel_plot <- function(file, ylab, yfun, ref, ylim) {
  pdf(file, width = 9, height = 8)
  par(mfrow = c(3, 3), mar = c(4, 4, 2, 1), oma = c(4, 0, 0, 0))
  for (mm in mgrid) for (nn in ngrid) {
    d <- diag_res[diag_res$m == mm & diag_res$n == nn, ]
    plot(NA, xlim = range(egrid), ylim = ylim, xlab = expression(paste(e[1], " = ", e[2])),
         ylab = ylab, main = paste0("m = ", mm, ", n = ", nn), las = 1)
    abline(h = ref, col = "grey60")
    for (cf in slope) {
      dd <- d[d$coef == cf, ]
      lines(dd$e1, yfun(dd), lty = sty$lty[cf], pch = sty$pch[cf], type = "b", lwd = 1.5)
    }
  }
  bottom_legend(c(expression(beta[1]), expression(beta[2]), expression(beta[3])),
                sty$lty, sty$pch)
  invisible(dev.off())
}

## Figure 1: MSE_SM / MSE_FLA
panel_plot("figures/fig_mse_ratio.pdf", "MSE ratio SM/FLA", function(d) d$mse_ratio, 1,
           range(c(0, diag_res$mse_ratio[diag_res$coef %in% slope])))
## Figure S2: coverage of the SM 95% interval (sandwich SE)
panel_plot("figures/fig_coverage.pdf", "Coverage (SM)", function(d) d$cov_sandwich, 0.95,
           c(min(diag_res$cov_sandwich[diag_res$coef %in% slope]), 1))

## Figure S1: SD and MSE of SM, FLA and WCE, averaged over the 16 (e1, e2)
res$sd_SM  <- sqrt(res$var_SM / res$n);  res$sd_FLA <- sqrt(res$var_FLA / res$n)
res$sd_WCE <- sqrt(res$var_WCE / res$n); res$mse_WCE <- res$var_WCE / res$n
agg <- aggregate(cbind(sd_SM, sd_FLA, sd_WCE, mse_SM, mse_FLA, mse_WCE) ~ coef + n + m,
                 data = res, FUN = mean)
write.csv(agg, "results/summary_by_n_m.csv", row.names = FALSE)

pdf("figures/fig_summary.pdf", width = 10, height = 5.8)
par(mfrow = c(2, 4), mar = c(4.2, 4.2, 2, 0.8), mgp = c(2.6, 0.7, 0), oma = c(4, 0, 0, 0))
cn   <- c(b0 = expression(beta[0]), b1 = expression(beta[1]),
          b2 = expression(beta[2]), b3 = expression(beta[3]))
mlty <- c(SM = 1, FLA = 2, WCE = 3)
mpch <- c(SM = 19, FLA = 2, WCE = 0)
draw <- function(cf, vars, ylab) {
  d  <- agg[agg$coef == cf, ]; d <- d[order(d$n, d$m), ]
  plot(NA, xlim = c(0.7, 9.3), ylim = range(unlist(d[vars])), xaxt = "n", xlab = "",
       ylab = ylab, main = cn[[cf]], log = "y", las = 1)
  axis(1, at = 1:9, labels = rep(mgrid, 3), cex.axis = 0.8)
  mtext(paste("n =", ngrid), side = 1, line = 2.3, at = c(2, 5, 8), cex = 0.6)
  abline(v = c(3.5, 6.5), col = "grey85")
  for (k in seq_along(vars)) for (blk in 0:2) {
    idx <- blk * 3 + 1:3
    lines(idx, d[[vars[k]]][idx], lty = mlty[k], pch = mpch[k], type = "b", lwd = 1.5)
  }
}
for (cf in names(cn)) draw(cf, c("sd_SM", "sd_FLA", "sd_WCE"), "SD")
for (cf in names(cn)) draw(cf, c("mse_SM", "mse_FLA", "mse_WCE"), "MSE")
bottom_legend(names(mlty), mlty, mpch)
invisible(dev.off())

for (f in c("fig_mse_ratio", "fig_coverage", "fig_summary")) eps(f)
cat("Section 4.2 done: see results/ and figures/\n")
