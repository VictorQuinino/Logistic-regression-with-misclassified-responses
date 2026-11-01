## =============================================================================
## Section 4.3 (Finite-sample check): Table 2, and Tables S1 and S2 and
## Figure S3 of the Supplementary Material, with the numbers quoted in
## Section 4.3.
## The simulation also stores the fits used in Section 4.5 (Table S5) and in
## Supplementary Section S2; run this program before Section4_5_unanimous.R
## and Supp_S2_binomial_regression.R.
## Running time: about 40 minutes (about 20 minutes on Linux or macOS, where
## the two random streams run in parallel). If the replicates were already
## saved, they are read instead; set rerun <- TRUE to simulate again.
## Output: results/finite_sample_reps.rds (all replicates),
##         results/finite_sample.csv (summaries by scenario),
##         results/tables_section4_3.txt (LaTeX rows of Tables 2, S1, S2 and
##         numbers quoted in the text), figures/fig_finite.pdf (Figure S3)
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

R     <- 1000        # replicates per scenario
rerun <- FALSE       # TRUE: simulate again even if the replicates were saved
rds   <- "results/finite_sample_reps.rds"

beta  <- c(b0 = -1, b1 = 0.2, b2 = 1.5, b3 = -0.6)
j     <- 3           # coefficient reported: beta2 (binary covariate x2)
jnull <- 2           # size: test of beta1 = 0, data generated under beta1 = 0
G     <- covariate_grid("indep")

## ---- simulation -------------------------------------------------------------
## FLA summaries use the replicates in which the FLA converged with a positive
## definite observed information; SM summaries use the replicates without
## separation. Failure rates are reported.
gen_X <- function(n) cbind("(Intercept)" = 1, x1 = rlnorm(n), x2 = rbinom(n, 1, 1/3),
                           x3 = runif(n))
gen_S <- function(X, b, m, e1, e2) {
  Yt <- rbinom(nrow(X), 1, expit(drop(X %*% b)))
  ifelse(Yt == 1, rbinom(nrow(X), m, 1 - e1), rbinom(nrow(X), m, e2))
}

one_rep <- function(m, e1, e2, n, do_fla) {
  X  <- gen_X(n); S <- gen_S(X, beta, m, e1, e2)
  sm <- sm_fit(X, S, m); un <- unan_fit(X, S, m); bi <- binom_fit(X, S, m)
  fl <- if (do_fla) fla_fit(X, S, m) else list(coef = rep(NA, 6), se = rep(NA, 6), ok = NA)
  b0 <- beta; b0[jnull] <- 0
  X0 <- gen_X(n); S0 <- gen_S(X0, b0, m, e1, e2)
  sm0 <- sm_fit(X0, S0, m); bi0 <- binom_fit(X0, S0, m)
  fl0 <- if (do_fla) fla_fit(X0, S0, m) else list(coef = rep(NA, 6), se = rep(NA, 6), ok = NA)
  c(sm_est = sm$coef[[j]], sm_seu = sm$se_usual[[j]], sm_ses = sm$se_sandwich[[j]],
    sm_ok = sm$ok,
    un_est = un$coef[[j]], un_se = un$se_usual[[j]], un_ok = un$ok,
    bi_est = bi$coef[[j]], bi_ok = bi$ok,
    fl_est = fl$coef[[j]], fl_se = fl$se[[j]], fl_ok = fl$ok,
    sm0_tu = sm0$coef[[jnull]] / sm0$se_usual[[jnull]],
    sm0_ts = sm0$coef[[jnull]] / sm0$se_sandwich[[jnull]], sm0_ok = sm0$ok,
    bi0_t  = bi0$coef[[jnull]] / bi0$se[[jnull]], bi0_ok = bi0$ok,
    fl0_t  = fl0$coef[[jnull]] / fl0$se[[jnull]], fl0_ok = fl0$ok,
    ## intercept (Table S3 of the Supplementary Material)
    sm_est_int = sm$coef[[1]], sm_seu_int = sm$se_usual[[1]],
    sm_ses_int = sm$se_sandwich[[1]], fl_est_int = fl$coef[[1]], fl_se_int = fl$se[[1]])
}

if (rerun || !file.exists(rds)) {
  scen <- expand.grid(n = c(50, 100, 500, 1000), e = c("0.05,0.05", "0.05,0.20", "0.20,0.20"),
                      m = c(3, 5), stringsAsFactors = FALSE)
  RNGkind("L'Ecuyer-CMRG"); set.seed(20261017)
  reps <- vector("list", nrow(scen))
  for (i in seq_len(nrow(scen))) {
    m <- scen$m[i]; n <- scen$n[i]; e <- as.numeric(strsplit(scen$e[i], ",")[[1]])
    cat(sprintf("Scenario %d of %d: m = %d, (e1, e2) = (%.2f, %.2f), n = %d\n",
                i, nrow(scen), m, e[1], e[2], n))
    s <- as.data.frame(do.call(rbind, run_replicates(R, function(r)
           one_rep(m, e[1], e[2], n, do_fla = n <= 500))))
    reps[[i]] <- cbind(m = m, e1 = e[1], e2 = e[2], n = n, s)
  }
  saveRDS(do.call(rbind, reps), rds)
  RNGkind("default")
}

## ---- summaries by scenario -------------------------------------------------
reps <- readRDS(rds)
d <- do.call(rbind, lapply(split(reps, list(reps$m, reps$e1, reps$e2, reps$n), drop = TRUE),
  summarize_finite, beta = beta, j = j, G = G))
d <- d[order(d$m, d$e1, d$e2, d$n), ]; rownames(d) <- NULL
write.csv(d, "results/finite_sample.csv", row.names = FALSE)

## intercept: the same summaries with the intercept estimates (Table S3)
int_cols <- function(s) {
  s$sm_est <- s$sm_est_int; s$sm_seu <- s$sm_seu_int; s$sm_ses <- s$sm_ses_int
  s$fl_est <- s$fl_est_int; s$fl_se <- s$fl_se_int; s
}
di <- do.call(rbind, lapply(split(reps, list(reps$m, reps$e1, reps$e2, reps$n), drop = TRUE),
  function(s) summarize_finite(int_cols(s), beta = beta, j = 1, G = G)))
di <- di[order(di$m, di$e1, di$e2, di$n), ]; rownames(di) <- NULL
write.csv(di, "results/finite_sample_intercept.csv", row.names = FALSE)

## ---- LaTeX rows of Tables 2, S1 and S2 -----------------------------------------
f2 <- function(x) ifelse(x > 0 & x < 0.001, "$<$0.001",
                   ifelse(abs(x) < 0.1, formatC(x, format = "f", digits = 3),
                          formatC(x, format = "f", digits = 2)))
pc <- function(x) formatC(100 * x, format = "f", digits = 1)
p1 <- function(x) formatC(x, format = "f", digits = 1)
e2 <- function(a, b) sprintf("(%.2f,%.2f)", a, b)
brk <- function(r, every) { k <- seq(every, length(r) - 1, by = every)
  r[k] <- paste0(r[k], "\n\\addlinespace"); r }
out <- character(0); add <- function(title, rows) out <<- c(out, paste("%", title), rows, "")

d5  <- subset(d, n <= 500)
lab <- with(d5, paste0(m, " & ", e2(e1, e2), " & ", n))

add("Table 2 (m = 3)", brk(with(subset(d5, m == 3), paste0(e2(e1, e2), " & ", n, " & ",
  f2(mse_ratio_asy), " & ", f2(mse_ratio_fs), " & ", f2(mse_ratio_fs_all), " & ",
  f2(mae_sm), " & ", f2(mae_fla), " & ", f2(mae_fla_all),
  " & ", p1(fail_fla), " & ", pc(cov_sm), " & ", pc(cov_fla), " & ", pc(size_sm),
  " & ", pc(size_fla), "\\\\")), 3))

add("Table S1", brk(with(d5, paste0(lab, " & ", f2(mse_ratio_asy), " & ", f2(mse_ratio_fs),
  " & ", f2(mse_ratio_fs_all), " & ", f2(bias_sm), " & ", f2(bias_fla), " & ", f2(bias_fla_all),
  " & ", f2(mae_sm), " & ", f2(mae_fla), " & ", f2(mae_fla_all),
  " & ", p1(fail_sm), " & ", p1(fail_fla), "\\\\")), 3))

add("Table S2", brk(with(d5, paste0(lab, " & ", pc(cov_sm_asy), " & ", pc(cov_sm),
  " & ", pc(cov_fla), " & ", pc(size_sm), " & ", pc(size_sm_sand), " & ", pc(size_fla),
  " & ", pc(pow_sm), " & ", pc(pow_fla), "\\\\")), 3))

di5 <- subset(di, n <= 500)
add("Table S3 (intercept)", brk(with(di5, paste0(m, " & ", e2(e1, e2), " & ", n,
  " & ", f2(mse_ratio_asy), " & ", f2(mse_ratio_fs), " & ", f2(mse_ratio_fs_all),
  " & ", f2(bias_sm), " & ", f2(bias_fla), " & ", f2(bias_fla_all),
  " & ", pc(cov_sm), " & ", pc(cov_fla), "\\\\")), 3))

## ---- numbers quoted in Section 4.3 ----------------------------------------
s <- function(x) paste(range(x, na.rm = TRUE), collapse = " to ")
small <- subset(d5, n <= 100)
txt <- c(
  paste("FLA failures (%), n<=100, m=3:", s(subset(small, m == 3)$fail_fla),
        "; m=5:", s(subset(small, m == 5)$fail_fla)),
  paste("SM failures (%), all:", s(d$fail_sm)),
  paste("max FLA bias, n<=100:", max(small$bias_fla), "; SM bias there:",
        small$bias_sm[which.max(small$bias_fla)]),
  paste("FS MSE ratio <= asymptotic for all n<=100:", all(small$mse_ratio_fs <= small$mse_ratio_asy)),
  paste("MAE SM <= MAE FLA (+0.01) for all n<=100:", all(small$mae_sm <= small$mae_fla + 0.01)),
  paste("size SM usual (%):", s(100 * d5$size_sm), "; FLA:", s(100 * d5$size_fla),
        "; SM sandwich:", s(100 * d5$size_sm_sand)),
  paste("max |coverage FS - asy| SM (pp), n<=100:",
        max(abs(100 * (small$cov_sm - small$cov_sm_asy)))),
  paste("coverage FLA (%):", s(100 * d5$cov_fla)),
  paste("power SM - FLA (pp), n<=100:", s(100 * (small$pow_sm - small$pow_fla))),
  paste("keeping FLA failures: MSE ratio, n<=100:", s(signif(small$mse_ratio_fs_all, 2)),
        "; max |FLA estimate| among failed fits:", round(max(d5$max_abs_fla_failed, na.rm = TRUE), 1),
        "; MAE FLA all minus successful, max:", round(max(d5$mae_fla_all - d5$mae_fla), 3)),
  paste("intercept, n = 500: asy MSE ratio", s(round(subset(di5, n == 500)$mse_ratio_asy, 2)),
        "; FS", s(round(subset(di5, n == 500)$mse_ratio_fs, 2)),
        "; coverage SM (%)", s(round(100 * subset(di5, n == 500)$cov_sm, 1))),
  paste("intercept, n <= 100: FS MSE ratio (successful)", s(round(subset(di5, n <= 100)$mse_ratio_fs, 2)),
        "; keeping failures", s(signif(subset(di5, n <= 100)$mse_ratio_fs_all, 2))))
add("Numbers quoted in Section 4.3", paste("%", txt))
writeLines(out, "results/tables_section4_3.txt")
cat(out, sep = "\n")

## ---- Figure S3: bias and MSE ratio SM/FLA for beta2 -----------------------------
ec <- unique(d5[, c("e1", "e2")])
nn <- sort(unique(d5$n))
pch_m <- c("3" = 19, "5" = 2)

pdf("figures/fig_finite.pdf", width = 10, height = 7)
layout(matrix(c(1, 2, 3, 4, 4, 4, 5, 6, 7, 8, 8, 8), ncol = 3, byrow = TRUE),
       heights = c(1, 0.13, 1, 0.13))
par(mar = c(4, 4.3, 2.2, 0.8), mgp = c(2.5, 0.7, 0))
leg <- function(labels, lty, pch) {
  par(mar = c(0, 0, 0, 0)); plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("center", labels, lty = lty, pch = pch, lwd = 1.5, bty = "n",
         ncol = length(labels), seg.len = 3, x.intersp = 0.8, text.width = 0.25, cex = 1.3)
  par(mar = c(4, 4.3, 2.2, 0.8))
}
yl <- range(c(0, d5$bias_sm, d5$bias_fla))
for (k in seq_len(nrow(ec))) {
  sk <- d5[d5$e1 == ec$e1[k] & d5$e2 == ec$e2[k], ]
  plot(NA, xlim = c(1, 3), ylim = yl, xaxt = "n", xlab = "n", ylab = "Bias",
       main = sprintf("(e1, e2) = (%.2f, %.2f)", ec$e1[k], ec$e2[k]), las = 1)
  axis(1, at = 1:3, labels = nn); abline(h = 0, col = "grey60")
  for (m in c(3, 5)) {
    sm <- sk[sk$m == m, ]; x <- match(sm$n, nn); mm <- as.character(m)
    lines(x, sm$bias_sm,  lty = 1, pch = pch_m[mm], type = "b", lwd = 1.5)
    lines(x, sm$bias_fla, lty = 2, pch = pch_m[mm], type = "b", lwd = 1.5)
  }
}
leg(c("SM", "FLA", "m = 3", "m = 5"), c(1, 2, NA, NA), c(NA, NA, 19, 2))
yl <- range(c(d5$mse_ratio_asy, d5$mse_ratio_fs))
for (k in seq_len(nrow(ec))) {
  sk <- d5[d5$e1 == ec$e1[k] & d5$e2 == ec$e2[k], ]
  plot(NA, xlim = c(1, 3), ylim = yl, log = "y", xaxt = "n", yaxt = "n", xlab = "n",
       ylab = "MSE ratio SM/FLA (log scale)")
  at <- c(0.001, 0.01, 0.1, 0.5, 1, 2, 5)
  axis(2, at = at, labels = format(at, scientific = FALSE, drop0trailing = TRUE), las = 1)
  axis(1, at = 1:3, labels = nn); abline(h = 1, col = "grey60")
  for (m in c(3, 5)) {
    sm <- sk[sk$m == m, ]; x <- match(sm$n, nn); mm <- as.character(m)
    lines(x, sm$mse_ratio_asy, lty = 3, pch = pch_m[mm], type = "b", lwd = 1.5)
    lines(x, sm$mse_ratio_fs,  lty = 1, pch = pch_m[mm], type = "b", lwd = 1.5)
  }
}
leg(c("Finite sample", "Asymptotic", "m = 3", "m = 5"), c(1, 3, NA, NA), c(NA, NA, 19, 2))
invisible(dev.off())
eps("fig_finite")
cat("Section 4.3 done: see results/ and figures/\n")
