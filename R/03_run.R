################################################################################
## 03_run.R -- driver: runs every scenario across replicates and writes
##             out/raw_results.csv
##
## Usage:  Rscript R/03_run.R <nRep> <firstSeed> <outfile>
################################################################################

Sys.setenv(OMP_NUM_THREADS = "1")
suppressMessages({
  source("R/01_functions.R")
  source("R/02_scenarios.R")
})

args    <- commandArgs(trailingOnly = TRUE)
nRep    <- if (length(args) >= 1) as.integer(args[1]) else 20
seed0   <- if (length(args) >= 2) as.integer(args[2]) else 1000
outfile <- if (length(args) >= 3) args[3] else "out/raw_results.csv"

S    <- buildScenarios()
base <- defaultParams()

## G x E levels that need their own founder object (trait must be created with
## the appropriate varGxE at SimParam time)
gxeLevels <- sort(unique(c(0, vapply(S, function(s) s$needsGxE, numeric(1)))))

message(sprintf("scenarios = %d | replicates = %d | GxE founder builds/rep = %d",
                length(S), nRep, length(gxeLevels)))

first <- TRUE
for (r in seq_len(nRep)) {
  seed <- seed0 + r
  repT0 <- Sys.time()

  ## founder objects for this replicate, one per required G x E level
  founders <- list()
  for (gg in gxeLevels) {
    pf <- base; pf$varGxERatio <- gg
    founders[[as.character(gg)]] <- makeFounders(pf, seed = seed)
  }

  for (s in S) {
    f <- founders[[as.character(s$needsGxE)]]
    p <- s$mod(base)
    p$varGxERatio <- s$needsGxE

    donors <- if (s$strategy == "CB") drawDonors(f, p, s$donorImprove) else NULL

    res <- try(runProgram(s$strategy, f, p, donors = donors, seed = seed),
               silent = TRUE)
    if (inherits(res, "try-error")) {
      message("FAILED: ", s$name, " rep ", r, " : ", as.character(res))
      next
    }

    res$scenario     <- s$name
    res$axis         <- s$axis
    res$level        <- s$level
    res$rep          <- r
    res$seed         <- seed
    res$donorImprove <- s$donorImprove
    res$donorMeanObs <- if (is.null(donors)) NA_real_ else donors$realisedMean
    res$donorGapObs  <- if (is.null(donors)) NA_real_ else donors$realisedGap
    res$burnInUsed   <- f$burnInUsed
    res$eliteVarObs  <- f$eliteVarObs
    res$h2           <- p$h2
    res$cycleYears   <- p$cycleYears
    res$nBridge      <- p$nBridge
    res$bridgeBar    <- p$bridgeBar
    res$bridgeLag    <- p$bridgeLag
    res$varGxERatio  <- p$varGxERatio

    write.table(res, outfile, sep = ",", row.names = FALSE,
                col.names = first, append = !first)
    first <- FALSE
  }

  message(sprintf("rep %d/%d done in %.1f min", r, nRep,
                  as.numeric(difftime(Sys.time(), repT0, units = "mins"))))
}

message("finished -> ", outfile)
