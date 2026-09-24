################################################################################
## 02_scenarios.R -- scenario definitions
##
## Each scenario is a named list:
##   strategy      "GSRS" or "CB"
##   axis          the sensitivity axis it belongs to
##   level         the label of the level on that axis
##   group         "core" (heavily replicated) or "peripheral"
##   mod           function(p) returning a modified parameter list
##   donorImprove  generations of donor pre-improvement (donor quality)
##   needsGxE      varGxE/varA required at founder-simulation time (0 = none)
################################################################################

sc <- function(name, strategy, axis, level, mod = identity,
               donorImprove = 0, needsGxE = 0, group = "peripheral") {
  list(name = name, strategy = strategy, axis = axis, level = level,
       mod = mod, donorImprove = donorImprove, needsGxE = needsGxE,
       group = group)
}

setp <- function(...) {
  args <- list(...)
  function(p) { for (n in names(args)) p[[n]] <- args[[n]]; p }
}

buildScenarios <- function(group = NULL) {
  S <- list()

  ## ---- baseline -----------------------------------------------------------
  S[[length(S) + 1]] <- sc("GSRS_base", "GSRS", "baseline", "baseline",
                           group = "core")
  S[[length(S) + 1]] <- sc("CB_base", "CB", "baseline",
                           "2 bridge lines per cycle (as specified)",
                           group = "core")

  ## ---- (1) donor inflow rate: the immigration ladder ----------------------
  ## This is the axis that governs whether the response saturates.  Two
  ## elite-equivalent bridge lines per cycle out of 20 recycled parents is a
  ## 10% immigration rate; the ladder spans 5% to 40%.
  for (nb in c(1L, 4L, 8L)) {
    S[[length(S) + 1]] <- sc(sprintf("CB_nB%d", nb), "CB", "inflow_rate",
                             sprintf("%d bridge line%s per cycle", nb,
                                     if (nb == 1) "" else "s"),
                             mod = setp(nBridge = nb), group = "core")
  }

  ## ---- (2) mechanism controls --------------------------------------------
  ## Each removes one component of the bridging design, to show which parts of
  ## it are load-bearing.
  S[[length(S) + 1]] <- sc("CB_noRecurrent", "CB", "mechanism",
                           "bridge pool rebuilt each cycle (not recurrent)",
                           mod = setp(bridgeRecurrent = FALSE), group = "core")
  S[[length(S) + 1]] <- sc("CB_noWithinFam", "CB", "mechanism",
                           "bridge pool selected across families",
                           mod = setp(bridgeWithinFam = FALSE), group = "core")
  S[[length(S) + 1]] <- sc("CB_noDistinctFam", "CB", "mechanism",
                           "highest-ranking eligible lines admitted",
                           mod = setp(admitDistinctFam = FALSE), group = "core")

  ## The failure mode identified during development: a strict elite-equivalence
  ## threshold combined with a bridge population that is rebuilt from raw donors
  ## each cycle.  The bridge pool then falls progressively further behind the
  ## rising elite mean until no line can clear the threshold, and the programme
  ## reverts in practice to a closed pool.
  S[[length(S) + 1]] <- sc("CB_strictBar_noRecurrent", "CB", "mechanism",
                           "strict bar + non-recurrent bridge pool",
                           mod = setp(bridgeBar = "parentMean",
                                      bridgeRecurrent = FALSE),
                           group = "collapse")

  ## ---- (3) stringency of the elite-equivalence bar -----------------------
  for (bb in c("q75", "parentMean")) {
    lab <- c(q75 = "elite top quartile",
             parentMean = "mean of the selected elite parents")[bb]
    S[[length(S) + 1]] <- sc(paste0("CB_bar_", bb), "CB", "bridge_bar", lab,
                             mod = setp(bridgeBar = bb))
  }

  ## ---- (4) donor quality --------------------------------------------------
  for (di in c(1L, 2L, 3L)) {
    S[[length(S) + 1]] <- sc(sprintf("CB_donorImp%d", di), "CB",
                             "donor_quality",
                             sprintf("donors pre-improved %s cycle%s",
                                     c("one", "two", "three")[di],
                                     if (di == 1) "" else "s"),
                             donorImprove = di)
  }

  ## ---- (5) heritability ---------------------------------------------------
  for (h in c(0.15, 0.50)) {
    for (st in c("GSRS", "CB")) {
      S[[length(S) + 1]] <- sc(sprintf("%s_h2_%.2f", st, h), st, "heritability",
                               sprintf("h2 = %.2f", h), mod = setp(h2 = h))
    }
  }

  ## ---- (6) cycle length ---------------------------------------------------
  for (cl in c(2, 4)) {
    for (st in c("GSRS", "CB")) {
      S[[length(S) + 1]] <- sc(sprintf("%s_cyc%d", st, cl), st, "cycle_length",
                               paste0(cl, " yr"), mod = setp(cycleYears = cl))
    }
  }

  ## ---- (7) genotype x environment interaction ----------------------------
  for (gg in c(0.5, 1.0)) {
    for (st in c("GSRS", "CB")) {
      S[[length(S) + 1]] <- sc(sprintf("%s_gxe%.1f", st, gg), st, "gxe",
                               sprintf("varGxE/varA = %.1f", gg),
                               mod = setp(varGxERatio = gg), needsGxE = gg)
    }
  }

  ## ---- (8) resource-matched and time-penalised comparisons ---------------
  eqBudget <- function(p) {
    p$nCross       <- 15   # 15 x 40 = 600 elite candidates
    p$nDonorProg   <- 15   # 10 x 15 = 150 fresh donor lines
    p$nBridgeProg2 <- 15   # 10 x 15 = 150 improvement lines
    p$nTrain       <- 140
    p$nBridgePheno <- 60
    p
  }
  S[[length(S) + 1]] <- sc("CB_equalBudget", "CB", "resource_matched",
                           "equal genotyping + phenotyping budget",
                           mod = eqBudget)
  S[[length(S) + 1]] <- sc("CB_lag2", "CB", "resource_matched",
                           "2-cycle bridge development lag",
                           mod = setp(bridgeLag = 2L))
  S[[length(S) + 1]] <- sc("CB_equalBudget_lag2", "CB", "resource_matched",
                           "equal budget + 2-cycle lag",
                           mod = function(p) { p <- eqBudget(p)
                                               p$bridgeLag <- 2L; p })
  S[[length(S) + 1]] <- sc("CB_cyc4_vs_GSRS3", "CB", "resource_matched",
                           "CB 4-yr cycle vs GS-RS 3-yr cycle",
                           mod = setp(cycleYears = 4))

  ## `group` selects either a replication tier ("core"/"peripheral") or a single
  ## sensitivity axis by name, so a subset can be re-run on its own.
  if (!is.null(group)) {
    groups <- unique(vapply(S, function(x) x$group, character(1)))
    S <- if (group %in% groups)
      S[vapply(S, function(x) x$group == group, logical(1))]
    else
      S[vapply(S, function(x) x$axis == group, logical(1))]
  }
  S
}
