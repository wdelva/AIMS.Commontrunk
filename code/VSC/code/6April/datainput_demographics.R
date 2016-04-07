agedist.data.frame <- agedist.create(shape = 5, scale = 65)

cfg <- list()
cfg["population.eyecap.fraction"] <- 0.5 #0.25# 0.25 #0.1

cfg["population.simtime"] <- 25# 25 #10 #110#110#110 # 20 for ABC model fitting
cfg["population.numwomen"] <- 2500#2500#2000
cfg["population.nummen"] <- 2500#2500#2000

cfg["periodiclogging.interval"] <- 1
cfg["periodiclogging.starttime"] <- 0
cfg["syncrefyear.interval"] <- 1
cfg["syncpopstats.interval"] <- 1

cfg["hivseed.time"] <- 10#10#1000        # So there is no HIV in this simulation for ABC model fitting to demographics
cfg["hivseed.type"] <- "amount"
cfg["hivseed.amount"] <- 25#20     # 50 out of 10000 is 0.5% of the population
cfg["hivseed.age.min"] <- 20
cfg["hivseed.age.max"] <- 30 
cfg["conception.alpha_base"] <- -2 # Default is -3 which gives an intra-couple conception rate of 0.05. -2 gives a rate of 0.14
cfg["formation.hazard.type"] <- "agegapry"
cfg["formation.hazard.agegapry.gap_factor_man_const"] <- 0
cfg["formation.hazard.agegapry.gap_factor_woman_const"] <- 0
cfg["person.agegap.man.dist.type"] <- "normal" #person.agegap.man.dist.type[model]
cfg["person.agegap.woman.dist.type"] <- "normal" #person.agegap.woman.dist.type[model]
cfg["person.agegap.man.dist.normal.mu"] <- -4 #-1 # Not to be varied anymore in ABC fitting of reference model
cfg["person.agegap.woman.dist.normal.mu"] <- -4 #-1 # Not to be varied anymore in ABC fitting of reference model
cfg["formation.hazard.agegapry.baseline"] <- 0 # 5.5 # Not to be varied anymore in ABC fitting of reference model
cfg["formation.hazard.agegapry.numrel_scale_man"] <- 0.0 # Not to be varied anymore in ABC fitting of reference model
cfg["formation.hazard.agegapry.numrel_scale_woman"] <- 0.0 # Not to be varied anymore in ABC fitting of reference model
cfg["formation.hazard.agegapry.gap_agescale_man"] <- 0.28 # Not to be varied anymore in ABC fitting of reference model
cfg["formation.hazard.agegapry.gap_agescale_woman"] <- 0.28 # Not to be varied anymore in ABC fitting of reference model
cfg["person.eagerness.dist.type"] <- "gamma"
cfg["formation.hazard.agegapry.eagerness_sum"] <- 1 # Not to be varied anymore in ABC fitting of reference model
cfg["formation.hazard.agegapry.beta"] <- 0 # Not to be varied anymore in ABC fitting of reference model
cfg["dissolution.beta"] <- 0 # duration of relationship effect # Not to be varied anymore in ABC fitting of reference model
cfg["person.survtime.logoffset.dist.type"] = "normal"
cfg["person.survtime.logoffset.dist.normal.mu"] = 0
cfg["person.survtime.logoffset.dist.normal.sigma"] = 0.1
cfg["person.art.accept.threshold.dist.type"] = "fixed"
cfg["person.art.accept.threshold.dist.fixed.value"] = 0 # Make sure a person never accepts treatment
cfg["diagnosis.baseline"] <- -1000 # This will result in timing of HIV diagnosis way beyond the simulation period.
cfg["transmission.param.a"] <- -1.0352239#-0.3#-1.3997
cfg["transmission.param.b"] <- -89.3399940#-8#-12.0220 
cfg["transmission.param.c"] <- 0.4948478#0.1649
cfg["transmission.param.d1"] <- 0
cfg["transmission.param.d2"] <- 0
cfg["transmission.param.f1"] <- log(5) # ~1.6 such that the hazard is x 5 in 15 yo
cfg["transmission.param.f2"] <- log(log(2.5) / log(5)) / 5 #~ -0.11 and x 2.5 in 20 yo,compared to the reference (>>25)


# The following parameters will be overridden in the ABC loop

cfg["person.agegap.man.dist.normal.sigma"] <- 4.5
cfg["person.agegap.woman.dist.normal.sigma"] <- 4.5
cfg["formation.hazard.agegapry.numrel_man"] <- -0.2
cfg["formation.hazard.agegapry.numrel_woman"] <- -0.2
cfg["formation.hazard.agegapry.numrel_diff"] <- -0.1
cfg["person.eagerness.dist.gamma.a"] <- 0.125 # k parameter on wiki page
cfg["person.eagerness.dist.gamma.b"] <- 8 # theta parameter on wiki page
cfg["formation.hazard.agegapry.eagerness_diff"] <- -1 #0
cfg["formation.hazard.agegapry.gap_factor_man_exp"] <- -0.2
cfg["formation.hazard.agegapry.gap_factor_woman_exp"] <- -0.2
cfg["formation.hazard.agegapry.gap_factor_man_age"] <- 0.1 # New meaning is a multiplication factor, relative to gap_factor_gender_exp
cfg["formation.hazard.agegapry.gap_factor_woman_age"] <- 0.1 # New meaning is a multiplication factor, relative to gap_factor_gender_exp 
cfg["formation.hazard.agegapry.meanage"] <- -0.1
cfg["dissolution.alpha_0"] <- 0.1 # baseline
cfg["dissolution.alpha_4"] <- -0.05 # average age effect

targetwindow.end <- as.numeric(cfg["population.simtime"])
targetwindow.start <- targetwindow.end - 1   # We use the data over 1 year to calculate the PTR
targetwindow.duration <- targetwindow.end - targetwindow.start
targetwindow.midpoint <- targetwindow.start + targetwindow.duration/2


cfgABC <- cfg
cfgABC.times <- list(population.simtime = cfgABC$population.simtime, hivseed.time = cfgABC$hivseed.time)

save(cfgABC.times, file=paste0("/user/data/gent/vsc400/vsc40070/simpact-test/data/", "cfgABC.times.RData"))

