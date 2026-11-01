## =============================================================================
## Supplementary Section S3 (Estimators of e1 and e2): finite-sample comparison
## of the closed-form moment estimators of Section 2.4 with the maximum
## likelihood estimators of the FLA, in the scenarios of Section 4.3.
## Running time: about 15 minutes (about 8 minutes on Linux or macOS, where
## the two random streams run in parallel). If the replicates were already
## saved, they are read instead; set rerun <- TRUE to simulate again.
## Output: results/e_estimators_reps.rds (all replicates),
##         results/tables_supp_S3_e.txt (LaTeX rows of the table and numbers
##         quoted in the text)
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
dir.create("results", showWarnings = FALSE)
## ----------------------------------------------------------------------------
rerun <- FALSE
rds   <- "results/e_estimators_reps.rds"
R     <- 1000
beta  <- c(b0 = -1, b1 = 0.2, b2 = 1.5, b3 = -0.6)
z     <- qnorm(0.975)

gen_X <- function(n) cbind("(Intercept)" = 1, x1 = rlnorm(n), x2 = rbinom(n, 1, 1/3),
                           x3 = runif(n))
gen_S <- function(X, b, m, e1, e2) {
  Yt <- rbinom(nrow(X), 1, expit(drop(X %*% b)))
  ifelse(Yt == 1, rbinom(nrow(X), m, 1 - e1), rbinom(nrow(X), m, e2))
}

## moment estimators (with the EM fallback of Section 2.4) and FLA estimators
one_rep <- function(m, e1, e2, n) {
  X  <- gen_X(n); S <- gen_S(X, beta, m, e1, e2)
  mo <- tryCatch(moment_e(S, m), error = function(err)
          list(e = c(NA, NA), se = c(NA, NA), ok = FALSE, method = "error"))
  if (is.null(mo$se)) mo$se <- c(NA, NA)
  fl <- tryCatch(fla_fit(X, S, m), error = function(err)
          list(coef = rep(NA, 6), se = rep(NA, 6), ok = FALSE))
  c(mo_e1 = mo$e[[1]], mo_e2 = mo$e[[2]], mo_se1 = mo$se[[1]], mo_se2 = mo$se[[2]],
    mo_ok = isTRUE(mo$ok), mo_em = identical(mo$method, "EM"), mo_err = identical(mo$method, "error"),
    fl_e1 = fl$coef[[5]], fl_e2 = fl$coef[[6]], fl_se1 = fl$se[[5]], fl_se2 = fl$se[[6]],
    fl_ok = isTRUE(fl$ok))
}

if (rerun || !file.exists(rds)) {
  scen <- expand.grid(n = c(50, 100, 500), e = c("0.05,0.05", "0.05,0.20", "0.20,0.20"),
                      m = c(3, 5), stringsAsFactors = FALSE)
  RNGkind("L'Ecuyer-CMRG"); set.seed(20261031)
  reps <- vector("list", nrow(scen))
  for (i in seq_len(nrow(scen))) {
    m <- scen$m[i]; n <- scen$n[i]; e <- as.numeric(strsplit(scen$e[i], ",")[[1]])
    cat(sprintf("Scenario %d of %d: m = %d, (e1, e2) = (%.2f, %.2f), n = %d\n",
                i, nrow(scen), m, e[1], e[2], n))
    s <- as.data.frame(do.call(rbind, run_replicates(R, function(r)
           suppressWarnings(one_rep(m, e[1], e[2], n)))))
    reps[[i]] <- cbind(m = m, e1 = e[1], e2 = e[2], n = n, s)
  }
  saveRDS(do.call(rbind, reps), rds)
  RNGkind("default")
}
reps <- readRDS(rds)

## ---- summaries: moment estimators over all replicates in which they are
##      available; FLA over its successful fits ----------------------------------
summ <- function(est, se, ok, true) {
  k <- ok & is.finite(est)
  ks <- k & is.finite(se)
  c(bias = mean(est[k]) - true, rmse = sqrt(mean((est[k] - true)^2)),
    cover = mean(abs(est[ks] - true) <= z * se[ks]))
}
rmse <- function(x, t) sqrt(mean((x - t)^2))
A <- function(s) s$mo_ok == 1 & is.finite(s$fl_e1) & is.finite(s$fl_e2)
d <- do.call(rbind, lapply(split(reps, list(reps$m, reps$e1, reps$e2, reps$n), drop = TRUE),
  function(s) {
    t1 <- s$e1[1]; t2 <- s$e2[1]
    data.frame(m = s$m[1], e1 = t1, e2 = t2, n = s$n[1],
      t(summ(s$mo_e1, s$mo_se1, s$mo_ok == 1, t1)) |> `colnames<-`(paste0("mo1_", c("bias", "rmse", "cov"))),
      t(summ(s$fl_e1, s$fl_se1, s$fl_ok == 1, t1)) |> `colnames<-`(paste0("fl1_", c("bias", "rmse", "cov"))),
      t(summ(s$mo_e2, s$mo_se2, s$mo_ok == 1, t2)) |> `colnames<-`(paste0("mo2_", c("bias", "rmse", "cov"))),
      t(summ(s$fl_e2, s$fl_se2, s$fl_ok == 1, t2)) |> `colnames<-`(paste0("fl2_", c("bias", "rmse", "cov"))),
      mo_fail = 100 * mean(s$mo_ok != 1), mo_em = 100 * mean(s$mo_em == 1), mo_err = 100 * mean(s$mo_err == 1),
      fl_fail = 100 * mean(s$fl_ok != 1),
      ## same replicates for both: moment estimates available and FLA estimates
      ## returned, including those of the FLA fits flagged as failures
      mo1a = rmse(s$mo_e1[A(s)], t1), fl1a = rmse(s$fl_e1[A(s)], t1),
      mo2a = rmse(s$mo_e2[A(s)], t2), fl2a = rmse(s$fl_e2[A(s)], t2),
      bias_mo_a = max(abs(c(mean(s$mo_e1[A(s)]) - t1, mean(s$mo_e2[A(s)]) - t2))),
      bias_fl_a = max(abs(c(mean(s$fl_e1[A(s)]) - t1, mean(s$fl_e2[A(s)]) - t2))))
  }))
d <- d[order(d$m, d$e1, d$e2, d$n), ]
rownames(d) <- NULL
write.csv(d, "results/e_estimators.csv", row.names = FALSE)

f3 <- function(x) sprintf("%.3f", x); p1 <- function(x) sprintf("%.1f", 100 * x)
rows <- with(d, paste0(m, " & (", sprintf("%.2f", e1), ",", sprintf("%.2f", e2), ") & ", n,
  " & ", f3(mo1a), " & ", f3(fl1a), " & ", sprintf("%.2f", mo1a / fl1a),
  " & ", f3(mo2a), " & ", f3(fl2a), " & ", sprintf("%.2f", mo2a / fl2a),
  " & ", p1(mo1_cov), " & ", p1(mo2_cov),
  " & ", sprintf("%.1f", mo_em), " & ", sprintf("%.1f", fl_fail), "\\\\"))
out <- c("% Table: RMSE of moment (Mo.) and FLA estimators of e1 and e2, ratio Mo./FLA,",
         "% coverage (%) of 95% Wald intervals of the moment estimators, % EM fallback, % FLA failures",
         "% (RMSE on the same replicates, FLA including its estimates at the point of failure)",
         rows, "",
         paste("ratio RMSE Mo./FLA, same replicates incl. FLA failures: range",
               paste(range(round(c(d$mo1a / d$fl1a, d$mo2a / d$fl2a), 2)), collapse = "-"),
               "; m = 3, n = 50:", paste(range(round(c(d$mo1a / d$fl1a, d$mo2a / d$fl2a)[c(d$m, d$m) == 3 & c(d$n, d$n) == 50], 2)), collapse = "-")),
         paste("max |bias| same replicates, moment:", f3(max(d$bias_mo_a)), "; FLA:", f3(max(d$bias_fl_a))),
         "Successful FLA fits only:",
         paste("ratio RMSE Mo./FLA, e1: range", paste(range(round(d$mo1_rmse / d$fl1_rmse, 2)), collapse = "-"),
               "; e2: range", paste(range(round(d$mo2_rmse / d$fl2_rmse, 2)), collapse = "-")),
         paste("same, m = 3, n = 50:", paste(range(round(c(d$mo1_rmse / d$fl1_rmse, d$mo2_rmse / d$fl2_rmse)[c(d$m, d$m) == 3 & c(d$n, d$n) == 50], 2)), collapse = "-")),
         paste("same, n = 500:", paste(range(round(c(d$mo1_rmse / d$fl1_rmse, d$mo2_rmse / d$fl2_rmse)[c(d$n, d$n) == 500], 2)), collapse = "-")),
         paste("coverage moment intervals (%), n = 500:", paste(range(p1(c(d$mo1_cov, d$mo2_cov)[c(d$n, d$n) == 500])), collapse = "-"),
               "; n <= 100:", paste(range(p1(c(d$mo1_cov, d$mo2_cov)[c(d$n, d$n) <= 100])), collapse = "-")),
         paste("coverage FLA intervals (%):", paste(range(p1(c(d$fl1_cov, d$fl2_cov))), collapse = "-")),
         paste("max |bias| moment:", f3(max(abs(c(d$mo1_bias, d$mo2_bias)))),
               "; FLA:", f3(max(abs(c(d$fl1_bias, d$fl2_bias))))),
         paste("EM fallback (%), max:", max(d$mo_em), "; moment failures (%), max:", max(d$mo_fail),
               "; second largest:", sort(d$mo_fail, decreasing = TRUE)[2]),
         paste("FLA failures (%), m = 3, n <= 100, max:", max(d$fl_fail[d$m == 3 & d$n <= 100]),
               "; m = 5, max:", max(d$fl_fail[d$m == 5])))
writeLines(out, "results/tables_supp_S3_e.txt")
cat(out, sep = "\n")
print(d, digits = 3)

## ---- availability and comparison on the replicates where both estimators
##      were available ---------------------------------------------------------
both <- reps$mo_ok == 1 & reps$fl_ok == 1
avail <- c(only_moments = sum(reps$mo_ok == 1 & reps$fl_ok != 1),
           only_FLA = sum(reps$mo_ok != 1 & reps$fl_ok == 1),
           neither = sum(reps$mo_ok != 1 & reps$fl_ok != 1))
rc <- do.call(rbind, lapply(split(reps, list(reps$m, reps$e1, reps$e2, reps$n), drop = TRUE),
  function(s) {
    k <- s$mo_ok == 1 & s$fl_ok == 1
    rm <- function(x, t) sqrt(mean((x - t)^2))
    c(rm(s$mo_e1[k], s$e1[1]) / rm(s$fl_e1[k], s$e1[1]),
      rm(s$mo_e2[k], s$e2[1]) / rm(s$fl_e2[k], s$e2[1]))
  }))
out2 <- c(paste("replicates with only moment / only FLA / neither available:",
                paste(avail, collapse = " / ")),
          paste("moments available in", round(100 * avail[1] / (avail[1] + avail[3]), 1),
                "% of the replicates in which the FLA failed"),
          paste("ratio RMSE Mo./FLA on replicates where both available: range",
                paste(range(round(rc, 2)), collapse = "-")))
write(out2, "results/tables_supp_S3_e.txt", append = TRUE)
cat(out2, sep = "\n")

## ---- replicates in which the FLA failed and the moment estimates were
##      available: RMSE of the FLA estimates at the point of failure relative
##      to the moment estimates (scenarios with at least 20 such replicates) ----
ff <- subset(reps, fl_ok != 1 & mo_ok == 1 & is.finite(fl_e1) & is.finite(fl_e2))
rf <- do.call(rbind, lapply(split(ff, list(ff$m, ff$e1, ff$e2, ff$n), drop = TRUE),
  function(s) if (nrow(s) >= 20) {
    rm <- function(x, t) sqrt(mean((x - t)^2))
    c(rm(s$fl_e1, s$e1[1]) / rm(s$mo_e1, s$e1[1]), rm(s$fl_e2, s$e2[1]) / rm(s$mo_e2, s$e2[1]))
  }))
out3 <- c(paste("FLA failed, moments available: ratio RMSE FLA/Mo., range",
                paste(range(round(rf, 2)), collapse = "-"),
                "; FLA estimate of e1 or e2 below 1e-4 in",
                round(100 * mean(ff$fl_e1 < 1e-4 | ff$fl_e2 < 1e-4), 1), "% of them"))
write(out3, "results/tables_supp_S3_e.txt", append = TRUE)
cat(out3, sep = "\n")
