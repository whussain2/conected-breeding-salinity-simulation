################################################################################
## Connected Breeding vs closed-pool genomic recurrent selection
## Stochastic simulation, sensitivity analyses and resource-matched comparisons
##
## Companion code for:
##   "An Integrated Framework to Develop and Deliver Salt-Tolerant Rice
##    Varieties for Coastal Ecologies"
##
## 01_functions.R -- founder simulation, breeding-programme engines, helpers
##
## Requires: AlphaSimR (>= 1.5)
################################################################################

suppressMessages(library(AlphaSimR))

## -----------------------------------------------------------------------------
## Default parameters (baseline scenario; matches Supplemental Table S6)
## -----------------------------------------------------------------------------
defaultParams <- function() {
  list(
    ## Genome / trait architecture
    nChr            = 12,      # chromosomes
    genLen          = 1.2,     # Morgans per chromosome (~1,440 cM total)
    segSites        = 400,     # segregating sites simulated per chromosome
    nQtlPerChr      = 83,      # ~1,000 additive QTL genome-wide
    nSnpPerChr      = 100,     # 1,200-SNP panel (~1K-RiCA)
    splitGen        = 100,     # generations since elite/donor pool split
    nFounderPerPool = 200,     # coalescent founders simulated per pool

    ## Trait scaling (elite base population)
    eliteMean       = 6.0,     # t ha-1
    eliteVarA       = 0.5,     # (t ha-1)^2


    ## G x E
    varGxERatio     = 0,       # varGxE / varA ; 0 = no GxE (baseline)

    ## Programme size
    nEliteFounder   = 60,
    nDonorFounder   = 60,
    maxBurnIn       = 25,      # cap on elite burn-in generations
    nDonorDrift     = 5,       # generations of random mating in the donor pool
    eliteDonorGap   = 2.5,     # t ha-1 by which the elite pool must exceed the
                               # donor reservoir at the start of the comparison
    nDonorImprove   = 0,       # generations of upward selection applied to the
                               # donor reservoir before use ("donor quality":
                               # 0 = raw exotic donors, >0 = pre-improved)
    nParent         = 20,      # parents recycled per cycle
    nCross          = 30,      # elite crosses per cycle
    nProgeny        = 40,      # progeny per elite cross (30 x 40 = 1,200)
    nTrain          = 200,     # phenotyped training lines per cycle
    trainCap        = 1000,    # sliding-window cap on training records
    h2              = 0.30,    # narrow-sense heritability of training phenotypes

    ## Connected Breeding: donor-elite (DE) bridging
    nBridge         = 2,       # elite-equivalent bridge lines entering per cycle
    nBridgeCross1   = 10,      # donor x elite crosses
    nBridgeProg1    = 40,      # progeny per donor x elite cross  (400 lines)
    nBridgeCross2   = 10,      # (best bridge) x elite crosses
    nBridgeProg2    = 40,      # progeny per second-step cross    (400 lines)
    nBridgePheno    = 100,     # bridge lines phenotyped per cycle
    bridgeLag       = 0,       # cycles of delay before bridges enter elite pool
    bridgeBar       = "parentMean",  # eligibility bar for an "elite-equivalent"
                               # bridge line: "parentMean" = mean GEBV of the
                               # selected elite parents (the manuscript's
                               # elite-equivalent assumption); "poolMean" =
                               # mean GEBV of the whole elite candidate pool;
                               # "q75" = mean GEBV of the top quartile

    ## Horizon
    cycleYears      = 3,       # years per cycle
    horizonYears    = 60       # total horizon
  )
}

## -----------------------------------------------------------------------------
## Founder simulation.  One call per replicate; the resulting object is reused
## by every scenario within that replicate so that scenarios are compared on
## identical founder genomes (paired design).
## -----------------------------------------------------------------------------
makeFounders <- function(p, seed) {
  set.seed(seed)

  founderPop <- runMacs2(
    nInd     = 2 * p$nFounderPerPool,
    nChr     = p$nChr,
    segSites = p$segSites,
    genLen   = p$genLen,
    split    = p$splitGen,
    Ne       = 100
  )

  SP <- SimParam$new(founderPop)
  SP$setTrackPed(FALSE)
  ## single-threaded execution: AlphaSimR parallelises recombination with
  ## per-thread random-number streams, so nThreads must be 1 for bit-for-bit
  ## reproducibility from a given seed
  SP$nThreads <- 1L
  SP$restrSegSites(minQtlPerChr = p$nQtlPerChr, minSnpPerChr = p$nSnpPerChr,
                   overlap = FALSE)
  ## Trait is created on the standardised scale (mean 0, var 1) across all
  ## founders; results are rescaled to t ha-1 afterwards so that the elite base
  ## population has mean 6.0 and additive variance 0.5.
  if (p$varGxERatio > 0) {
    SP$addTraitAG(nQtlPerChr = p$nQtlPerChr, mean = 0, var = 1,
                  varGxE = p$varGxERatio)
  } else {
    SP$addTraitA(nQtlPerChr = p$nQtlPerChr, mean = 0, var = 1)
  }
  SP$addSnpChip(nSnpPerChr = p$nSnpPerChr)
  SP$setVarE(h2 = p$h2)

  pool1 <- newPop(founderPop[1:p$nFounderPerPool], simParam = SP)
  pool2 <- newPop(founderPop[(p$nFounderPerPool + 1):(2 * p$nFounderPerPool)],
                  simParam = SP)

  ## ---- Trait scale --------------------------------------------------------
  ## The trait is simulated on a standardised scale (mean 0, variance 1 in the
  ## unselected founder population).  One base-population additive standard
  ## deviation is defined as sqrt(eliteVarA) = 0.707 t ha-1, so that the base
  ## additive variance is 0.5 (t ha-1)^2 as specified.  Because this scale is a
  ## fixed constant rather than a per-replicate quantity, the elite-donor
  ## performance differential is expressed in identical units in every
  ## replicate.
  kappa <- sqrt(p$eliteVarA)

  ## ---- Donor reservoir: diverged and unimproved for yield -----------------
  donorPool <- pool2
  for (g in seq_len(p$nDonorDrift)) {
    donorPool <- randCross(donorPool, nCrosses = 100, nProgeny = 2,
                           simParam = SP)
  }
  donorRawMean <- mean(as.numeric(gv(donorPool)))

  ## ---- Elite pool: burn-in of recurrent selection -------------------------
  ## Elite germplasm is by definition the product of prior selection.  Burn-in
  ## continues until the elite pool exceeds the donor reservoir by the target
  ## differential (default 2.5 t ha-1, i.e. a 6.0 t ha-1 elite pool against a
  ## 3.5 t ha-1 donor pool), so the differential is set by construction rather
  ## than left to drift between the coalescent sub-populations.
  pop <- pool1
  gUsed <- 0L
  for (g in seq_len(p$maxBurnIn)) {
    if (kappa * (mean(as.numeric(gv(pop))) - donorRawMean) >= p$eliteDonorGap) break
    pop <- setPheno(pop, h2 = p$h2, simParam = SP)
    par <- selectInd(pop, nInd = 60, use = "pheno", simParam = SP)
    pop <- randCross(par, nCrosses = 100, nProgeny = 2, simParam = SP)
    gUsed <- g
  }
  elite0 <- pop[sample.int(nInd(pop), p$nEliteFounder)]

  gvE  <- as.numeric(gv(elite0))
  beta <- p$eliteMean - kappa * mean(gvE)

  list(SP = SP, elite0 = elite0, donorPool0 = donorPool,
       kappa = kappa, beta = beta, seed = seed,
       burnInUsed  = gUsed,
       eliteVarObs = (kappa^2) * stats::var(gvE),
       donorRawMean = donorRawMean)
}

## Convert simulated genetic values to the t ha-1 scale
toYield    <- function(g, f) f$beta + f$kappa * g
toYieldVar <- function(v, f) (f$kappa^2) * v

## -----------------------------------------------------------------------------
## Draw the donor founders at a given performance level.
## Donor quality is varied by sampling the donor reservoir at different
## quantiles of genetic value, which changes the donor mean while leaving the
## donor pool's allelic content and its divergence from the elite pool intact.
## Donor quality is set by the number of generations of pre-improvement applied
## to the donor reservoir (0 = raw exotic donors).  The realised donor mean and
## its deficit relative to the elite pool are recorded for every run.
## -----------------------------------------------------------------------------
drawDonors <- function(f, p, nImprove = 0) {
  SP   <- f$SP
  pool <- f$donorPool0
  ## nImprove generations of truncation selection inside the donor reservoir:
  ## level 0 is a raw, unimproved exotic/landrace pool; higher levels represent
  ## donors that have already been through pre-breeding.  The ordering is
  ## therefore monotone within every replicate, and the realised donor mean is
  ## reported for each level rather than being forced to a fixed value.
  if (nImprove > 0) {
    for (g in seq_len(nImprove)) {
      top  <- selectInd(pool, nInd = 60, use = "gv", simParam = SP)
      pool <- randCross(top, nCrosses = 100, nProgeny = 2, simParam = SP)
    }
  }
  idx <- sample.int(nInd(pool), p$nDonorFounder)
  y   <- toYield(as.numeric(gv(pool[idx])), f)
  list(pop = pool[idx],
       nImprove     = nImprove,
       realisedMean = mean(y),
       realisedGap  = p$eliteMean - mean(y))
}

## -----------------------------------------------------------------------------
## Fit RR-BLUP on the current training population and return GEBVs
## -----------------------------------------------------------------------------
predictGEBV <- function(trainPop, targetPop, SP) {
  mod <- RRBLUP(trainPop, use = "pheno", simParam = SP)
  as.numeric(setEBV(targetPop, mod, simParam = SP)@ebv)
}

## -----------------------------------------------------------------------------
## Breeding-programme engine
##
## strategy = "GSRS"  : closed-pool recurrent genomic selection
## strategy = "CB"    : Connected Breeding (same GS engine + DE bridging)
## strategy = "PSRS"  : closed-pool recurrent PHENOTYPIC selection (all
##                      candidates phenotyped; no genomic prediction)
## strategy = "MAGIC" : a one-time broad-base multi-parent population founded
##                      from the elite AND donor founders together, thereafter
##                      improved by closed-pool genomic recurrent selection with
##                      no further inflow of donor diversity
##
## Returns a data.frame with one row per cycle:
##   year, meanG (t ha-1), gain (t ha-1 above elite base), varG (t ha-1)^2,
##   varGpct (% of base additive variance), plus the realised resource use
##   (genotyped lines, phenotyped plots, crosses) accumulated to that cycle.
## -----------------------------------------------------------------------------
runProgram <- function(strategy, f, p, donors = NULL, seed = 1) {
  set.seed(seed)
  SP <- f$SP

  nCycle <- floor(p$horizonYears / p$cycleYears)

  ## Base population statistics (cycle 0).  These are always taken from the
  ## ELITE base population so that gain and variance retention are expressed on
  ## a common reference scale across all strategies.
  baseMeanG <- mean(as.numeric(gv(f$elite0)))
  baseVarG  <- stats::var(as.numeric(gv(f$elite0)))

  ## ---- Initial training set and parents ----------------------------------
  ## MAGIC founds a single broad-base population from elite AND donor founders
  ## and is then closed: the donor diversity is supplied once, not recurrently.
  startPop <- if (strategy == "MAGIC") c(f$elite0, donors$pop) else f$elite0
  elite <- setPheno(startPop, h2 = p$h2, simParam = SP)
  trainPop <- elite
  parents  <- selectInd(elite, nInd = p$nParent, use = "pheno", simParam = SP)

  donorPool <- if (strategy == "CB") donors$pop else NULL
  bridgeQueue <- list()   # bridges awaiting entry when bridgeLag > 0

  ## resource counters
  nGeno <- nInd(elite); nPheno <- nInd(elite); nCrossTot <- 0
  nBridgeIn <- 0

  out <- vector("list", nCycle)

  for (cy in seq_len(nCycle)) {

    ## Target environment for this cycle.  When G x E is simulated, selection
    ## acts on phenotypes expressed in one realised environment, whereas the
    ## recorded genetic value is the main effect across the target population
    ## of environments -- i.e. G x E degrades selection accuracy exactly as it
    ## does in a real multi-environment programme.
    pEnv <- if (p$varGxERatio > 0) stats::runif(1) else NULL

    ## ---- 1. elite crossing block ----------------------------------------
    cand <- randCross(parents, nCrosses = p$nCross, nProgeny = p$nProgeny,
                      simParam = SP)
    nCrossTot <- nCrossTot + p$nCross
    nGeno     <- nGeno + nInd(cand)

    ## ---- 2a. phenotypic recurrent selection ------------------------------
    ## All candidates are phenotyped (h2 = p$h2) and parents chosen directly on
    ## phenotype.  This costs far more plots per cycle than genomic selection
    ## and is less accurate at this heritability, but requires no genotyping.
    if (strategy == "PSRS") {
      cand   <- setPheno(cand, h2 = p$h2, p = pEnv, simParam = SP)
      nPheno <- nPheno + nInd(cand)
      parents <- selectInd(cand, nInd = p$nParent, use = "pheno", simParam = SP)
      gvc <- as.numeric(gv(cand))
      out[[cy]] <- data.frame(
        cycle = cy, year = cy * p$cycleYears,
        meanG = toYield(mean(gvc), f),
        gain  = toYield(mean(gvc), f) - toYield(baseMeanG, f),
        varG  = toYieldVar(stats::var(gvc), f),
        varGpct = 100 * stats::var(gvc) / baseVarG,
        nGeno = 0, nPheno = nPheno, nCross = nCrossTot, nBridgeIn = 0)
      next
    }

    ## ---- 2b. phenotype the training slice --------------------------------
    trIdx <- sample.int(nInd(cand), min(p$nTrain, nInd(cand)))
    trNew <- setPheno(cand[trIdx], h2 = p$h2, p = pEnv, simParam = SP)
    nPheno <- nPheno + nInd(trNew)
    trainPop <- c(trainPop, trNew)
    if (nInd(trainPop) > p$trainCap) {
      keep <- (nInd(trainPop) - p$trainCap + 1):nInd(trainPop)
      trainPop <- trainPop[keep]
    }

    ## ---- 3. genomic prediction and elite selection ----------------------
    mod       <- RRBLUP(trainPop, use = "pheno", simParam = SP)
    cand      <- setEBV(cand, mod, simParam = SP)
    eliteEbv  <- as.numeric(ebv(cand))
    eliteBar  <- switch(p$bridgeBar,
      poolMean   = mean(eliteEbv),
      q75        = mean(eliteEbv[eliteEbv >= stats::quantile(eliteEbv, 0.75)]),
      parentMean = mean(sort(eliteEbv, decreasing = TRUE)[seq_len(p$nParent)]),
      stop("unknown bridgeBar")
    )

    ## ---- 4. Connected Breeding: donor-elite bridging --------------------
    newBridges <- NULL
    if (strategy == "CB") {
      topElite <- selectInd(cand, nInd = 10, use = "ebv", simParam = SP)

      ## step 1: donor x elite
      b1 <- randCross2(females = donorPool, males = topElite,
                       nCrosses = p$nBridgeCross1, nProgeny = p$nBridgeProg1,
                       simParam = SP)
      nCrossTot <- nCrossTot + p$nBridgeCross1
      nGeno     <- nGeno + nInd(b1)
      b1  <- setEBV(b1, mod, simParam = SP)
      b1t <- selectInd(b1, nInd = min(20, nInd(b1)), use = "ebv", simParam = SP)

      ## step 2: (best bridge) x elite -- raises the bridge mean toward elite
      b2 <- randCross2(females = b1t, males = topElite,
                       nCrosses = p$nBridgeCross2, nProgeny = p$nBridgeProg2,
                       simParam = SP)
      nCrossTot <- nCrossTot + p$nBridgeCross2
      nGeno     <- nGeno + nInd(b2)

      ## bridge lines are phenotyped too (a real cost of the strategy)
      if (p$nBridgePheno > 0) {
        bIdx <- sample.int(nInd(b2), min(p$nBridgePheno, nInd(b2)))
        bPh  <- setPheno(b2[bIdx], h2 = p$h2, p = pEnv, simParam = SP)
        nPheno   <- nPheno + nInd(bPh)
        trainPop <- c(trainPop, bPh)
        if (nInd(trainPop) > p$trainCap) {
          keep <- (nInd(trainPop) - p$trainCap + 1):nInd(trainPop)
          trainPop <- trainPop[keep]
        }
      }

      b2 <- setEBV(b2, mod, simParam = SP)
      ## "elite-equivalent" pre-selection: only bridges reaching the mean
      ## genomic breeding value of the elite candidate pool are eligible
      elig <- b2[as.numeric(ebv(b2)) >= eliteBar]
      if (nInd(elig) > 0) {
        newBridges <- selectInd(elig, nInd = min(p$nBridge, nInd(elig)),
                                use = "ebv", simParam = SP)
      }

      ## maintain the donor reservoir (random mating, no selection)
      donorPool <- randCross(donorPool, nCrosses = nInd(donorPool),
                             nProgeny = 1, simParam = SP)
    }

    ## ---- 5. assemble parents for the next cycle -------------------------
    if (strategy == "CB") {
      ## honour the bridging lag, if any
      bridgeQueue[[length(bridgeQueue) + 1]] <- newBridges
      dueIdx <- length(bridgeQueue) - p$bridgeLag
      entering <- if (dueIdx >= 1) bridgeQueue[[dueIdx]] else NULL

      nB <- if (is.null(entering)) 0 else nInd(entering)
      nBridgeIn <- nBridgeIn + nB
      topE <- selectInd(cand, nInd = p$nParent - nB, use = "ebv", simParam = SP)
      parents <- if (nB > 0) c(topE, entering) else topE
    } else {
      parents <- selectInd(cand, nInd = p$nParent, use = "ebv", simParam = SP)
    }

    ## ---- 6. record ------------------------------------------------------
    gvc <- as.numeric(gv(cand))
    out[[cy]] <- data.frame(
      cycle   = cy,
      year    = cy * p$cycleYears,
      meanG   = toYield(mean(gvc), f),
      gain    = toYield(mean(gvc), f) - toYield(baseMeanG, f),
      varG    = toYieldVar(stats::var(gvc), f),
      varGpct = 100 * stats::var(gvc) / baseVarG,
      nGeno   = nGeno,
      nPheno  = nPheno,
      nCross  = nCrossTot,
      nBridgeIn = nBridgeIn
    )
  }

  res <- do.call(rbind, out)
  res$strategy <- strategy
  res
}
