simpact4ABC_demographics <- function(inputvector){
  library(RSimpactCyan)
  source("/user/data/gent/vsc400/vsc40070/simpact-test/code/postsimscript.R")
  source("/user/data/gent/vsc400/vsc40070/simpact-test/code/datainput_demographics.R")
  source("/user/data/gent/vsc400/vsc40070/simpact-test/code/modeloutput_demographics.R")
  
  ABC_DestDir <- "/user/data/gent/vsc400/vsc40070/simpact-test/simoutput"
  
  cfgABC["person.agegap.man.dist.normal.sigma"] <- inputvector[2]
  cfgABC["person.agegap.woman.dist.normal.sigma"] <- inputvector[2]
  cfgABC["formation.hazard.agegapry.numrel_man"] <- inputvector[3]
  cfgABC["formation.hazard.agegapry.numrel_woman"] <- inputvector[3]
  cfgABC["formation.hazard.agegapry.numrel_diff"] <- inputvector[4]
  cfgABC["person.eagerness.dist.gamma.a"] <- inputvector[5]
  cfgABC["person.eagerness.dist.gamma.b"] <- inputvector[6]
  cfgABC["formation.hazard.agegapry.eagerness_diff"] <- inputvector[7]
  cfgABC["formation.hazard.agegapry.gap_factor_man_exp"] <- inputvector[8]
  cfgABC["formation.hazard.agegapry.gap_factor_woman_exp"] <- inputvector[8]
  cfgABC["formation.hazard.agegapry.gap_factor_man_age"] <- inputvector[8] * inputvector[9]
  cfgABC["formation.hazard.agegapry.gap_factor_woman_age"] <- inputvector[8] * inputvector[9]
  cfgABC["conception.alpha_base"] <- inputvector[10]
  cfgABC["conception.alpha_agewoman"] <- inputvector[11]  
  cfgABC["dissolution.alpha_4"] <- inputvector[12]

  seedid <- inputvector[1]  

  simTime <- cfgABC[["population.simtime"]]
  if(is.null(simTime))
    stop("population.simtime was not present in the settings")

  maxEvents = 1000000
  cfgABC["population.maxevents"] <- maxEvents

  results <- simpact.run.wrapper(cfgABC,
                         ABC_DestDir,
                         agedist = agedist.data.frame,
                         intervention = NULL,
                         release=TRUE,
                         slowalg=FALSE,
                         parallel=FALSE,
                         seed=seedid)

if (results["simulationtime"] < simTime)
  {
    # Ik kan op dit moment twee redenen bedenken waarom de simulatie te vroeg zou stoppen
    #  - maximaal aantal events is bereikt, kan op gecheckt worden dmv ret["eventsexecuted"]
    #  - geen events meer, gebeurt bvb als populatie uitsterft. In dat geval zal het maximaal
    #    aan
    if (results["eventsexecuted"] >= maxEvents-1)
    {
      # Ik doe hier een -1 omdat in R getallen standaard voorgesteld worden als floating
      # point getallen, en echte gelijkheden daarmee nogal gevaarlijk zijn.
      unlink(paste0(file.path(ABC_DestDir, results$id), "personlog.csv"))
      unlink(paste0(file.path(ABC_DestDir, results$id), "relationlog.csv"))
      unlink(paste0(file.path(ABC_DestDir, results$id), "eventlog.csv"))
      unlink(paste0(file.path(ABC_DestDir, results$id), "periodiclog.csv"))
      unlink(paste0(file.path(ABC_DestDir, results$id), "treatmentlog.csv"))
      unlink(paste0(file.path(ABC_DestDir, results$id), "agedist.csv"))
      stop("MAXEVENTS: Simulation stopped prematurely, max events reached")
    }
    else
    {
      # Misschien moet er in dit geval niet gestopt worden
      unlink(paste0(file.path(ABC_DestDir, results$id), "personlog.csv"))
      unlink(paste0(file.path(ABC_DestDir, results$id), "relationlog.csv"))
      unlink(paste0(file.path(ABC_DestDir, results$id), "eventlog.csv"))
      unlink(paste0(file.path(ABC_DestDir, results$id), "periodiclog.csv"))
      unlink(paste0(file.path(ABC_DestDir, results$id), "treatmentlog.csv"))
      unlink(paste0(file.path(ABC_DestDir, results$id), "agedist.csv"))
      stop("Simulation stopped prematurely, probably ran out of events")
    }
  }
  
  summaryStats <- modeloutput.all(results)
  unlink(paste0(file.path(ABC_DestDir, results$id), "personlog.csv"))
  unlink(paste0(file.path(ABC_DestDir, results$id), "relationlog.csv"))
  unlink(paste0(file.path(ABC_DestDir, results$id), "eventlog.csv"))
  unlink(paste0(file.path(ABC_DestDir, results$id), "periodiclog.csv"))
  unlink(paste0(file.path(ABC_DestDir, results$id), "treatmentlog.csv"))
  (paste0(file.path(ABC_DestDir, results$id), "agedist.csv"))
  return(summaryStats)
}


