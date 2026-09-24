################################################################################
## diag_cb.R -- why does Connected Breeding saturate in the current
## implementation, when it retains a large additive variance?
##
## Instruments the CB cycle to separate three candidate explanations:
##   (a) the donor reservoir falls progressively further behind the elite pool,
##       so late-cycle bridge material is increasingly destructive;
##   (b) bridge lines admitted on predicted merit are not elite-equivalent in
##       TRUE genetic value, because the prediction model is elite-trained
##       (winner's curse on out-of-reference material);
##   (c) the "retained additive variance" is largely between-family structure
##       created by admixture, not usable allelic diversity -- in which case
##       genic variance (sum 2pq a^2) will NOT be elevated the way var(gv) is.
################################################################################

Sys.setenv(OMP_NUM_THREADS = "1")

## ------------------------------------------------------------------ paths --
## Locate the bundle this script belongs to, so it works whether it is run with
## Rscript, sourced, or run from the RStudio editor, and whatever the working
## directory happens to be.  The script's OWN folder is tried first, so if two
## copies of the bundle exist (for example a newer one unzipped inside an older
## one) each copy uses its own data.  Set CB_SIM_ROOT to override.
.cb_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f  <- grep("^--file=", ca, value = TRUE)
  if (length(f)) return(dirname(sub("^--file=", "", f[1])))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- tryCatch(get0("ofile", envir = sys.frame(i), inherits = FALSE),
                   error = function(e) NULL)
    if (is.character(of) && length(of) == 1L) return(dirname(of))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      isTRUE(tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE))) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                  error = function(e) "")
    if (is.character(p) && nzchar(p)) return(dirname(p))
  }
  NA_character_
}

.cb_is_root <- function(p) file.exists(file.path(p, "R", "01_functions.R"))

.cb_candidates <- function() {
  cand <- character(0)
  env <- Sys.getenv("CB_SIM_ROOT")
  if (nzchar(env)) cand <- c(cand, env)
  sd <- .cb_script_dir()
  if (!is.na(sd)) cand <- c(cand, file.path(sd, ".."), sd)
  cand <- c(cand, getwd(), file.path(getwd(), ".."),
            file.path(getwd(), "..", ".."))
  cand <- unique(normalizePath(cand, mustWork = FALSE))
  cand[vapply(cand, .cb_is_root, logical(1))]
}

.cb_roots <- .cb_candidates()
if (!length(.cb_roots)) {
  stop("Could not locate the simulation bundle.\n",
       "  Expected a folder containing R/01_functions.R at or above: ", getwd(),
       "\n  Run this script from the bundle folder, e.g.\n",
       "      setwd('/path/to/connected-breeding-salinity-simulation')\n",
       "      source('R/08_figures.R')\n",
       "  or set it explicitly:\n",
       "      Sys.setenv(CB_SIM_ROOT = '/path/to/connected-breeding-salinity-simulation')",
       call. = FALSE)
}
ROOT   <- .cb_roots[1]
OUTDIR <- file.path(ROOT, "out")
FIGDIR <- file.path(ROOT, "figures")
RDIR   <- file.path(ROOT, "R")
message("bundle root: ", ROOT)
## ---------------------------------------------------------------------------

suppressMessages(source(file.path(RDIR, "01_functions.R")))

args  <- commandArgs(trailingOnly = TRUE)
nRep  <- if (length(args) >= 1) as.integer(args[1]) else 2
seed0 <- if (length(args) >= 2) as.integer(args[2]) else 5000

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
p <- defaultParams()
rows <- list()

for (r in seq_len(nRep)) {
  seed <- seed0 + r
  f <- makeFounders(p, seed = seed)
  SP <- f$SP
  donors <- drawDonors(f, p, 0)

  for (strategy in c("GSRS", "CB")) {
    set.seed(seed)
    elite <- setPheno(f$elite0, h2 = p$h2, simParam = SP)
    trainPop <- elite
    parents <- selectInd(elite, nInd = p$nParent, use = "pheno", simParam = SP)
    donorPool <- if (strategy == "CB") donors$pop else NULL

    for (cy in seq_len(20)) {
      cand <- randCross(parents, nCrosses = p$nCross, nProgeny = p$nProgeny,
                        simParam = SP)
      trIdx <- sample.int(nInd(cand), p$nTrain)
      trNew <- setPheno(cand[trIdx], h2 = p$h2, simParam = SP)
      trainPop <- c(trainPop, trNew)
      if (nInd(trainPop) > p$trainCap)
        trainPop <- trainPop[(nInd(trainPop) - p$trainCap + 1):nInd(trainPop)]

      mod  <- RRBLUP(trainPop, use = "pheno", simParam = SP)
      cand <- setEBV(cand, mod, simParam = SP)
      e    <- as.numeric(ebv(cand)); g <- as.numeric(gv(cand))
      accuracy <- stats::cor(e, g)
      eliteBar <- mean(sort(e, decreasing = TRUE)[seq_len(p$nParent)])
      top20gv  <- mean(as.numeric(gv(selectInd(cand, p$nParent, use = "ebv",
                                               simParam = SP))))

      bridgeGvMean <- NA_real_; admittedGv <- NA_real_
      nElig <- NA_integer_; donorGvMean <- NA_real_; bridgeAcc <- NA_real_

      if (strategy == "CB") {
        topElite <- selectInd(cand, nInd = 10, use = "ebv", simParam = SP)
        b1 <- randCross2(donorPool, topElite, nCrosses = p$nBridgeCross1,
                         nProgeny = p$nBridgeProg1, simParam = SP)
        b1 <- setEBV(b1, mod, simParam = SP)
        b1t <- selectInd(b1, nInd = 20, use = "ebv", simParam = SP)
        b2 <- randCross2(b1t, topElite, nCrosses = p$nBridgeCross2,
                         nProgeny = p$nBridgeProg2, simParam = SP)
        bPh <- setPheno(b2[sample.int(nInd(b2), p$nBridgePheno)], h2 = p$h2,
                        simParam = SP)
        trainPop <- c(trainPop, bPh)
        if (nInd(trainPop) > p$trainCap)
          trainPop <- trainPop[(nInd(trainPop) - p$trainCap + 1):nInd(trainPop)]
        b2 <- setEBV(b2, mod, simParam = SP)
        eb <- as.numeric(ebv(b2)); gb <- as.numeric(gv(b2))
        bridgeAcc    <- stats::cor(eb, gb)
        bridgeGvMean <- mean(gb)
        donorGvMean  <- mean(as.numeric(gv(donorPool)))
        elig <- b2[eb >= eliteBar]
        nElig <- nInd(elig)
        if (nElig > 0) {
          adm <- selectInd(elig, nInd = min(p$nBridge, nElig), use = "ebv",
                           simParam = SP)
          admittedGv <- mean(as.numeric(gv(adm)))
          nB <- nInd(adm)
          parents <- c(selectInd(cand, nInd = p$nParent - nB, use = "ebv",
                                 simParam = SP), adm)
        } else {
          parents <- selectInd(cand, nInd = p$nParent, use = "ebv", simParam = SP)
        }
        donorPool <- randCross(donorPool, nCrosses = nInd(donorPool),
                               nProgeny = 1, simParam = SP)
      } else {
        parents <- selectInd(cand, nInd = p$nParent, use = "ebv", simParam = SP)
      }

      rows[[length(rows) + 1]] <- data.frame(
        rep = r, strategy = strategy, cycle = cy, year = cy * p$cycleYears,
        accuracy = accuracy,
        bridgeAcc = bridgeAcc,
        candMeanY   = toYield(mean(g), f),
        top20MeanY  = toYield(top20gv, f),
        parentMeanY = toYield(mean(as.numeric(gv(parents))), f),
        bridgeMeanY = if (is.na(bridgeGvMean)) NA_real_ else toYield(bridgeGvMean, f),
        admittedY   = if (is.na(admittedGv)) NA_real_ else toYield(admittedGv, f),
        donorMeanY  = if (is.na(donorGvMean)) NA_real_ else toYield(donorGvMean, f),
        nElig = nElig,
        ## total genetic variance of the candidate population (includes any
        ## between-family structure created by admixture)
        varG_total = toYieldVar(stats::var(g), f),
        ## genic variance: sum over QTL of 2pq a^2 -- allelic diversity only,
        ## unaffected by linkage disequilibrium or family structure
        varG_genic = toYieldVar(as.numeric(genicVarA(cand, simParam = SP)), f)
      )
    }
  }
  message("diag rep ", r, "/", nRep, " done")
}

d <- do.call(rbind, rows)
write.csv(d, file.path(OUTDIR, "diag_cb.csv"), row.names = FALSE)

ag <- aggregate(cbind(accuracy, bridgeAcc, candMeanY, top20MeanY, parentMeanY,
                      bridgeMeanY, admittedY, donorMeanY, nElig,
                      varG_total, varG_genic) ~ strategy + year,
                d, mean, na.action = na.pass, na.rm = TRUE)
ag <- ag[order(ag$strategy, ag$year), ]

fmt <- function(x, k = 2) formatC(x, format = "f", digits = k, width = 7)
cat("\n===== GS-RS =====\n")
cat(sprintf("%5s %8s %9s %9s %9s %10s %10s\n", "year", "acc",
            "candMean", "top20", "parents", "varTotal", "varGenic"))
for (i in which(ag$strategy == "GSRS")) with(ag[i, ],
  cat(sprintf("%5d %8s %9s %9s %9s %10s %10s\n", year, fmt(accuracy),
              fmt(candMeanY), fmt(top20MeanY), fmt(parentMeanY),
              fmt(varG_total, 3), fmt(varG_genic, 3))))

cat("\n===== Connected Breeding (current implementation) =====\n")
cat(sprintf("%5s %6s %6s %8s %8s %8s %8s %8s %8s %6s %9s %9s\n", "year",
            "acc", "brAcc", "candMean", "top20", "parents", "bridge",
            "admitted", "donor", "nElig", "varTotal", "varGenic"))
for (i in which(ag$strategy == "CB")) with(ag[i, ],
  cat(sprintf("%5d %6s %6s %8s %8s %8s %8s %8s %8s %6.0f %9s %9s\n", year,
              fmt(accuracy, 2), fmt(bridgeAcc, 2), fmt(candMeanY),
              fmt(top20MeanY), fmt(parentMeanY), fmt(bridgeMeanY),
              fmt(admittedY), fmt(donorMeanY), nElig,
              fmt(varG_total, 3), fmt(varG_genic, 3))))

cat("\n--- diagnosis summary ---\n")
cb <- ag[ag$strategy == "CB", ]; gs <- ag[ag$strategy == "GSRS", ]
cat(sprintf("Elite-donor gap: year 3 %.2f t/ha -> year 60 %.2f t/ha\n",
            cb$candMeanY[1] - cb$donorMeanY[1],
            cb$candMeanY[nrow(cb)] - cb$donorMeanY[nrow(cb)]))
cat(sprintf("Admitted bridge TRUE gv minus top-20 elite TRUE gv: year 3 %+.2f -> year 60 %+.2f t/ha\n",
            cb$admittedY[1] - cb$top20MeanY[1],
            cb$admittedY[nrow(cb)] - cb$top20MeanY[nrow(cb)]))
cat(sprintf("Prediction accuracy, elite %.2f vs bridge %.2f (mean over cycles)\n",
            mean(cb$accuracy), mean(cb$bridgeAcc, na.rm = TRUE)))
cat(sprintf("GENIC variance at year 60: CB %.3f vs GS-RS %.3f (ratio %.2f)\n",
            cb$varG_genic[nrow(cb)], gs$varG_genic[nrow(gs)],
            cb$varG_genic[nrow(cb)] / gs$varG_genic[nrow(gs)]))
cat(sprintf("TOTAL variance at year 60: CB %.3f vs GS-RS %.3f (ratio %.2f)\n",
            cb$varG_total[nrow(cb)], gs$varG_total[nrow(gs)],
            cb$varG_total[nrow(cb)] / gs$varG_total[nrow(gs)]))
cat(sprintf("Late-horizon slope (yr 45-60), CB %+.4f vs GS-RS %+.4f t/ha/yr\n",
            coef(lm(candMeanY ~ year, cb[cb$year >= 45, ]))[2],
            coef(lm(candMeanY ~ year, gs[gs$year >= 45, ]))[2]))
cat("\nwrote", file.path(OUTDIR, "diag_cb.csv"), "\n")
