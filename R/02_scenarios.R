################################################################################
## 02_scenarios.R -- scenario definitions
##
## Each scenario is a named list:
##   strategy   "GSRS" or "CB"
##   axis       the sensitivity axis it belongs to
##   level      the label of the level on that axis
##   mod        a function(p) returning a modified parameter list
##   donorImprove  generations of donor pre-improvement (donor quality)
##   needsGxE   varGxE/varA ratio required at founder-simulation time (0 = none)
################################################################################

sc <- function(name, strategy, axis, level, mod = identity,
               donorImprove = 0, needsGxE = 0) {
  list(name = name, strategy = strategy, axis = axis, level = level,
       mod = mod, donorImprove = donorImprove, needsGxE = needsGxE)
}

buildScenarios <- function() {
  S <- list()

  ## ---- baseline -----------------------------------------------------------
  S[[length(S) + 1]] <- sc("GSRS_base", "GSRS", "baseline", "baseline")
  S[[length(S) + 1]] <- sc("CB_base",   "CB",   "baseline", "baseline")

  ## ---- (1) donor quality --------------------------------------------------
  ## GS-RS is a closed pool and is unaffected by donor quality, so only CB is
  ## re-run; the baseline GS-RS trajectory is the common reference.
  for (di in c(1L, 2L, 3L)) {
    S[[length(S) + 1]] <- sc(sprintf("CB_donorImp%d", di), "CB",
                             "donor_quality",
                             sprintf("%d cycle(s) of donor pre-improvement", di),
                             donorImprove = di)
  }

  ## ---- (2) number of elite-equivalent bridge lines per cycle --------------
  for (nb in c(1, 4, 8)) {
    S[[length(S) + 1]] <- sc(paste0("CB_bridge", nb), "CB", "n_bridge",
                             paste0(nb, " per cycle"),
                             mod = local({
                               nb_ <- nb
                               function(p) { p$nBridge <- nb_; p }
                             }))
  }

  ## ---- (3) bridge-quality bar (how "elite-equivalent" is defined) ---------
  for (bb in c("poolMean", "q75")) {
    lab <- c(poolMean = "elite pool mean", q75 = "elite top quartile")[bb]
    S[[length(S) + 1]] <- sc(paste0("CB_bar_", bb), "CB", "bridge_bar",
                             lab,
                             mod = local({
                               bb_ <- bb
                               function(p) { p$bridgeBar <- bb_; p }
                             }))
  }

  ## ---- (4) heritability ---------------------------------------------------
  for (h in c(0.15, 0.50)) {
    for (st in c("GSRS", "CB")) {
      S[[length(S) + 1]] <- sc(sprintf("%s_h2_%.2f", st, h), st, "heritability",
                               sprintf("h2 = %.2f", h),
                               mod = local({
                                 h_ <- h
                                 function(p) { p$h2 <- h_; p }
                               }))
    }
  }

  ## ---- (5) cycle length ---------------------------------------------------
  for (cl in c(2, 4)) {
    for (st in c("GSRS", "CB")) {
      S[[length(S) + 1]] <- sc(sprintf("%s_cyc%d", st, cl), st, "cycle_length",
                               paste0(cl, " yr"),
                               mod = local({
                                 cl_ <- cl
                                 function(p) { p$cycleYears <- cl_; p }
                               }))
    }
  }

  ## ---- (6) genotype x environment interaction ----------------------------
  for (gg in c(0.5, 1.0)) {
    for (st in c("GSRS", "CB")) {
      S[[length(S) + 1]] <- sc(sprintf("%s_gxe%.1f", st, gg), st, "gxe",
                               sprintf("varGxE/varA = %.1f", gg),
                               mod = local({
                                 gg_ <- gg
                                 function(p) { p$varGxERatio <- gg_; p }
                               }),
                               needsGxE = gg)
    }
  }

  ## ---- (7) resource-matched and time-penalised comparisons ---------------
  ## Reviewer point: the baseline comparison gives Connected Breeding an extra
  ## bridging programme (crossing, genotyping, phenotyping) and assumes bridge
  ## lines are available immediately.  These scenarios remove both advantages.

  ## (7a) equal total budget: the bridging programme is paid for out of the
  ##      elite programme, so CB genotypes 1,200 lines and phenotypes 200
  ##      plots per cycle -- exactly the GS-RS budget.
  eqBudget <- function(p) {
    p$nCross       <- 15   # 15 x 40 = 600 elite candidates
    p$nBridgeProg1 <- 30   # 10 x 30 = 300 bridge lines, step 1
    p$nBridgeProg2 <- 30   # 10 x 30 = 300 bridge lines, step 2
    p$nTrain       <- 140
    p$nBridgePheno <- 60
    p
  }
  S[[length(S) + 1]] <- sc("CB_equalBudget", "CB", "resource_matched",
                           "equal genotyping + phenotyping budget",
                           mod = eqBudget)

  ## (7b) development lag: bridge lines take two extra cycles (6 years) to
  ##      reach elite-equivalent status before they can be recycled.
  S[[length(S) + 1]] <- sc("CB_lag2", "CB", "resource_matched",
                           "2-cycle bridge development lag",
                           mod = function(p) { p$bridgeLag <- 2L; p })

  ## (7c) both constraints together (strictest test)
  S[[length(S) + 1]] <- sc("CB_equalBudget_lag2", "CB", "resource_matched",
                           "equal budget + 2-cycle lag",
                           mod = function(p) {
                             p <- eqBudget(p); p$bridgeLag <- 2L; p
                           })

  ## (7d) time cost: running the bridging programme lengthens the CB cycle to
  ##      4 years while GS-RS stays at 3 years.
  S[[length(S) + 1]] <- sc("CB_cyc4_vs_GSRS3", "CB", "resource_matched",
                           "CB 4-yr cycle vs GS-RS 3-yr cycle",
                           mod = function(p) { p$cycleYears <- 4; p })

  S
}
