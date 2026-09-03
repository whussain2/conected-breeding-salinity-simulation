################################################################################
## 07_run_baseline.R -- extra replicates of the baseline pair only.
##
## The headline contrast (Connected Breeding versus closed-pool genomic recurrent
## selection under baseline assumptions) is the quantity most worth replicating,
## and it costs only two programme runs per replicate.  These replicates are
## pooled with the grid runs by 05_analyse.R.
##
## Usage: Rscript R/07_run_baseline.R <nRep> <firstSeed> <outfile>
################################################################################

Sys.setenv(OMP_NUM_THREADS = "1")
suppressMessages(source("R/01_functions.R"))

args    <- commandArgs(trailingOnly = TRUE)
nRep    <- if (length(args) >= 1) as.integer(args[1]) else 20
seed0   <- if (length(args) >= 2) as.integer(args[2]) else 7000
outfile <- if (length(args) >= 3) args[3] else "out/raw_Cbase.csv"

base  <- defaultParams()
first <- TRUE

for (r in seq_len(nRep)) {
  seed <- seed0 + r
  f <- makeFounders(base, seed = seed)
  donors <- drawDonors(f, base, 0)
  for (sc in c("GSRS_base", "CB_base")) {
    st <- if (sc == "CB_base") "CB" else "GSRS"
    res <- try(runProgram(st, f, base,
                          donors = if (st == "CB") donors else NULL, seed = seed),
               silent = TRUE)
    if (inherits(res, "try-error")) {
      message("FAILED ", sc, " rep ", r, ": ", as.character(res)); next
    }
    res$scenario     <- sc
    res$axis         <- "baseline"
    res$level        <- "baseline"
    res$rep          <- r
    res$seed         <- seed
    res$donorImprove <- 0L
    res$donorMeanObs <- donors$realisedMean
    res$donorGapObs  <- donors$realisedGap
    res$burnInUsed   <- f$burnInUsed
    res$eliteVarObs  <- f$eliteVarObs
    res$h2           <- base$h2
    res$cycleYears   <- base$cycleYears
    res$nBridge      <- base$nBridge
    res$bridgeBar    <- base$bridgeBar
    res$bridgeLag    <- base$bridgeLag
    res$varGxERatio  <- 0
    write.table(res, outfile, sep = ",", row.names = FALSE,
                col.names = first, append = !first)
    first <- FALSE
  }
  message(sprintf("baseline rep %d/%d done", r, nRep))
}
message("finished -> ", outfile)
