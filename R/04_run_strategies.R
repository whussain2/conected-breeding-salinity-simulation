################################################################################
## 04_run_strategies.R -- additional established strategies simulated under the
## SAME parameterisation, so that the comparative entries in Supplemental
## Table S5 are outputs of this simulation rather than unattributed estimates.
##
##   PSRS  : closed-pool recurrent phenotypic selection
##   MAGIC : one-time broad-base multi-parent population, then closed genomic
##           recurrent selection (no recurrent donor inflow)
##
## Usage: Rscript R/04_run_strategies.R <nRep> <firstSeed> <outfile>
################################################################################

Sys.setenv(OMP_NUM_THREADS = "1")
suppressMessages(source("R/01_functions.R"))

args    <- commandArgs(trailingOnly = TRUE)
nRep    <- if (length(args) >= 1) as.integer(args[1]) else 10
seed0   <- if (length(args) >= 2) as.integer(args[2]) else 9000
outfile <- if (length(args) >= 3) args[3] else "out/raw_strategies.csv"

base <- defaultParams()
STRATS <- c("PSRS", "MAGIC")

first <- TRUE
for (r in seq_len(nRep)) {
  seed <- seed0 + r
  f <- makeFounders(base, seed = seed)
  donors <- drawDonors(f, base, 0)
  for (st in STRATS) {
    res <- try(runProgram(st, f, base, donors = donors, seed = seed),
               silent = TRUE)
    if (inherits(res, "try-error")) {
      message("FAILED ", st, " rep ", r, ": ", as.character(res)); next
    }
    res$scenario <- st
    res$axis     <- "strategy_comparison"
    res$level    <- st
    res$rep      <- r
    res$seed     <- seed
    write.table(res, outfile, sep = ",", row.names = FALSE,
                col.names = first, append = !first)
    first <- FALSE
  }
  message(sprintf("rep %d/%d done", r, nRep))
}
message("finished -> ", outfile)
