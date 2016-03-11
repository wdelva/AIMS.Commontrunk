simpact4ABC_demographics <- function(inputvector){
  library(RSimpactCyan)
  source("/Users/wimdelva/Documents/AIMS Essays/2016/AIMS.Commontrunk/code/postsimscript.R")
  source("/Users/wimdelva/Documents/AIMS Essays/2016/AIMS.Commontrunk/code/datainput_demographics.R")
  source("/Users/wimdelva/Documents/AIMS Essays/2016/AIMS.Commontrunk/code/modeloutput_demographics.R")
  
  ABC_DestDir <- "/Users/wimdelva/Documents/AIMS Essays/2016/AIMS.Commontrunk/simoutput"
  
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
  cfgABC["formation.hazard.agegapry.meanage"] <- inputvector[10]
  cfgABC["dissolution.alpha_0"] <- inputvector[11]
  cfgABC["dissolution.alpha_4"] <- inputvector[12]

  seedid <- inputvector[1]  

  simTime <- cfgABC[["population.simtime"]]
  if(is.null(simTime))
    stop("population.simtime was not present in the settings")

  maxEvents = 250000
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
      stop("MAXEVENTS: Simulation stopped prematurely, max events reached")
    }
    else
    {
      # Misschien moet er in dit geval niet gestopt worden
      stop("Simulation stopped prematurely, probably ran out of events")
    }
  }
  
  summaryStats <- modeloutput.all(results)
  unlink(paste0(file.path(ABC_DestDir, results$id), "*"))
  return(summaryStats)
}


