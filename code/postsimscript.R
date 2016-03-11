library(data.table)
library(shape)
library(reshape2)
library(nlme)
library(phylobase)
#library(sna)
library(ggplot2)
library(stringr)
library(plyr)
library(RColorBrewer)
library(survival) # for Cox proportional hazards regression models
library(AlgDesign)
library(lhs)
library(stringi)
#install.packages("readcsvcolumns", repos="http://research.edm.uhasselt.be/jori/")
library(readcsvcolumns)
#library(emulator)

#library(MIICD) # for Cox proportional hazards regression models
#library(interval)

#library(sqldf)

# We want to replace the default age distribution with one that is consistent with the non-HIV-related mortality hazard.
# mortality.normal.weibull.shape                               = 4
# mortality.normal.weibull.scale                               = 70
# mortality.normal.weibull.genderdiff                          = 0 (default is 5 and will need to be changed to 0 further down when the cfg object is created)


agedist.create <- function(shape = 5, scale = 65){
  agebins <- seq(0.5, 99.5, 1)
  probofstillalive <- 1 - pweibull(agebins, shape = shape, scale = scale)
  #plot(agebins,probofstillalive, type="l")
  #meanweib <- 70 * gamma(1+ (1/4))
  fractionsinagebins <- 100 * probofstillalive/sum(probofstillalive)
  simple.age.data.frame <- data.frame(Age = c(agebins, 100.5), Percent.Male = c(fractionsinagebins, 0), Percent.Female = c(fractionsinagebins, 0))
  
  return(simple.age.data.frame)
}

create.cfg <- function(inputvector = inputvector){
  cfg <- list()
  cfg["population.eyecap.fraction"] <- inputvector[1]
  cfg["population.simtime"] <- inputvector[2]
  cfg["population.numwomen"] <- inputvector[3]
  cfg["population.nummen"] <- inputvector[4]
  
  cfg["periodiclogging.interval"] <- inputvector[5]
  cfg["periodiclogging.starttime"] <- inputvector[6]
  cfg["syncrefyear.interval"] <- inputvector[7]
  cfg["syncpopstats.interval"] <- inputvector[8]
  
  cfg["hivseed.time"] <- inputvector[9]
  cfg["hivseed.type"] <- inputvector[10]
  cfg["hivseed.amount"] <- inputvector[11]
  cfg["hivseed.age.min"] <- inputvector[12]
  cfg["hivseed.age.max"] <- inputvector[13]
  cfg["formation.hazard.type"] <- inputvector[14]
  cfg["formation.hazard.agegapry.gap_factor_man_const"] <- inputvector[15]
  cfg["formation.hazard.agegapry.gap_factor_woman_const"] <- inputvector[16]
  cfg["person.agegap.man.dist.type"] <- inputvector[17]
  cfg["person.agegap.woman.dist.type"] <- inputvector[18]
  cfg["person.agegap.man.dist.normal.mu"] <- inputvector[19]
  cfg["person.agegap.woman.dist.normal.mu"] <- inputvector[20]
  cfg["formation.hazard.agegapry.baseline"] <- inputvector[21]
  cfg["formation.hazard.agegapry.numrel_scale_man"] <- inputvector[22]
  cfg["formation.hazard.agegapry.numrel_scale_woman"] <- inputvector[23]
  cfg["formation.hazard.agegapry.gap_agescale_man"] <- inputvector[24]
  cfg["formation.hazard.agegapry.gap_agescale_woman"] <- inputvector[25]
  
  cfg["person.eagerness.dist.type"] <- inputvector[26]
  cfg["formation.hazard.agegapry.eagerness_sum"] <- inputvector[27]
  
  
  cfg["person.survtime.logoffset.dist.type"] <- inputvector[28]
  cfg["person.survtime.logoffset.dist.normal.mu"] <- inputvector[29]
  cfg["person.survtime.logoffset.dist.normal.sigma"] <- inputvector[30]
  cfg["person.art.accept.threshold.dist.type"] <- inputvector[31]
  cfg["person.art.accept.threshold.dist.fixed.value"] <- inputvector[32]
  
  # Nobody will start ART during the simulations
  cfg["diagnosis.baseline"] <- inputvector[33]
  # We are using the transmission parameters as estimated by fitting non-linear model to Fraser model output
  cfg["transmission.param.a"] <- inputvector[34]
  cfg["transmission.param.b"] <- inputvector[35]
  cfg["transmission.param.c"] <- inputvector[36]
  
  cfg["transmission.param.d1"] <- inputvector[37]
  cfg["transmission.param.d2"] <- inputvector[38]
  cfg["transmission.param.f1"] <- inputvector[39]
  cfg["transmission.param.f2"] <- inputvector[40]
  
  cfg["dissolution.alpha_0"] <- inputvector[41]
  cfg["dissolution.alpha_4"] <- inputvector[42]
  cfg["dissolution.beta"] <- inputvector[43]
  
  
  
  cfg["person.agegap.man.dist.normal.sigma"] <- inputvector[44]
  cfg["person.agegap.woman.dist.normal.sigma"] <- inputvector[44] #!
  cfg["formation.hazard.agegapry.numrel_man"] <- inputvector[45]
  cfg["formation.hazard.agegapry.numrel_woman"] <- inputvector[45] #!
  cfg["formation.hazard.agegapry.numrel_diff"] <- inputvector[46]
  
  cfg["person.eagerness.dist.gamma.a"] <- inputvector[47]
  cfg["person.eagerness.dist.gamma.b"] <- inputvector[48]
  cfg["formation.hazard.agegapry.eagerness_diff"] <- inputvector[49]
  
  cfg["formation.hazard.agegapry.gap_factor_man_exp"] <- inputvector[50]
  cfg["formation.hazard.agegapry.gap_factor_woman_exp"] <- inputvector[50] #!
  cfg["formation.hazard.agegapry.gap_factor_man_age"] <- inputvector[50] * inputvector[51]
  cfg["formation.hazard.agegapry.gap_factor_woman_age"] <- inputvector[50] * inputvector[51] #!
  cfg["formation.hazard.agegapry.meanage"] <- inputvector[52]
  cfg["formation.hazard.agegapry.beta"] <- inputvector[53]
  

  
 
  
  
  
  
  
  return(cfg)
}




is.positive <- function(number){
  is.positive <- number > 0
  return(is.positive)
}

errFunction <- function(e)
{
  length.targetSS <- 11+16 # 36 HIV prevalence measurements + 11 other summary statistics
  largeNumber <- NA#0.5#1e1
  if (length(grep("NaN",e$message)) != 0)
    return(summaryStats = rep(largeNumber, length.targetSS))
  
  
  largeNumber <- NA#0.5#2e1
  if (length(grep("MAXEVENTS",e$message)) != 0)
    return(summaryStats = rep(largeNumber, length.targetSS))
  
  largeNumber <- NA#0.5#3e1
  if (length(grep("OUTPUTVECTOR", e$message)) != 0)
    return(summaryStats = rep(largeNumber, length.targetSS))
  
  
  # Een andere foutmelding dan NaN of MAXEVENTS of OUTPUTVECTOR, zorg dat we toch stoppen
  # (tenzij er nog een andere tryCatch is om ook dit op te vangen)
  stop(e)
}

agemixing.errorFunction <- function(e)
{
  if (length(grep("system is computationally singular: reciprocal condition number", e$message)) != 0)
    return(list())
  
  if (length(grep("non-numeric argument to mathematical function", e$message)) !=0)
    return(list())
  
  if (length(grep("Lapack routine dgesv: system is exactly singular", e$message)) !=0)
    return(list())
  
  
  stop(e)
}

agemixing.lme.errFunction <- function(e)
{
  return(list())
}



simpact.run.wrapper <- function(configParams, destDir, agedist = "${SIMPACT_DATA_DIR}sa_2003.csv", intervention = NULL, release = TRUE, slowalg = FALSE, parallel=FALSE, seed=-1, dryrun = FALSE, identifierFormat = "%T-%y-%m-%d-%H-%M-%S_%p_%r%r%r%r%r%r%r%r-" )
{
  cfg <- configParams
  if (is.null(cfg))
    cfg <- list()
  
  intv <- intervention
  if (is.null(intv))
    intv <- list()
  
  ret <- simpact.run(cfg, destDir, agedist, intv, release, slowalg, parallel, seed, dryrun, identifierFormat)
  
  inputParams <- list()
  inputParams["config"] <- list(cfg)
  inputParams["destination"] <- destDir
  inputParams[["agedist"]] <- agedist
  inputParams["intervention"] <- list(intv)
  inputParams["alg.release"] <- release
  inputParams["alg.slow"] <- slowalg
  inputParams["seed"] <- seed
  inputParams["dryrun"] <- dryrun
  inputParams["identifierFormat"] <- identifierFormat
  
  ret["input.params"] <- list(inputParams)
  return(ret)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

readthedata <- function(modeloutput){
  path <- as.character(modeloutput["outputfile"])
  outputID <- as.character(modeloutput["id"])
  DestDir <- sub(pattern = paste0(outputID, "output.txt"), replacement = "", x = path, fixed=T)
  personlogfilename <- paste0(DestDir, outputID, "personlog.csv")
  relationlogfilename <- paste0(DestDir, outputID, "relationlog.csv")
  eventlogfilename <- paste0(DestDir, outputID, "eventlog.csv")
  treatmentlogfilename <- paste0(DestDir, outputID, "treatmentlog.csv")
  periodiclogfilename <- paste0(DestDir, outputID, "periodiclog.csv")
    
  ptable <- fread(personlogfilename, sep = ",", skip = 0)
  rtable <- fread(relationlogfilename, sep = ",", skip = 0)
  etable <- setDT(read.csv.columns(eventlogfilename, has.header = FALSE, column.types = "rssiirsiir"))
  #etable <- fread(eventlogfilename, sep = ",", skip = 0)
  # setnames(etable, c("col_001", "col_002", "col_003", "col_004", "col_005", "col_006", "col_007", "col_008", "col_009", "col_010"),
  #          c("eventtime", "eventname", "p1name", "p1ID", "p1gender", "p1age", "p2name", "p2ID", "p2gender", "p2age"))
  ttable <- fread(treatmentlogfilename, sep  = ",", skip = 0)

  if (file.exists(periodiclogfilename)){
    ltable <- fread(periodiclogfilename, sep = ",", skip = 0)
    outputtables <- list(ptable = ptable, rtable = rtable, etable = etable, ttable = ttable, ltable = ltable)
  } else {
    outputtables <- list(ptable = ptable, rtable = rtable, etable = etable, ttable = ttable)
  }
  return(outputtables)

}

####
# To create heatplot, we need to loop through time steps, and calculate age- and gender-specific prevalence and incidence rates.
prevalenceheatmapdata <- function(modeloutput){ #(DT = datalist$ptable, cfg = configfile){
  datalist <- readthedata(modeloutput)
  DT <- datalist$ptable
  # inp <- modeloutput[["input.params"]][["config"]]
  # cfg <- inp[["config"]] # Should be 'cfg' again
  cfg <- modeloutput[["input.params"]][["config"]]
  output <- data.table()
  timestep.width <- 1
  decplaces <- decimalplaces(timestep.width)
  timesteps <- seq(as.numeric(cfg["hivseed.time"]), as.numeric(cfg["population.simtime"]), by = timestep.width)
  for (time_i in timesteps){
    DTalive.infected <- alive.infected(DT = datalist$ptable, time = time_i, dec.places = decplaces) # First we only take the data of people who were alive at time_i
    Prevalencetimestepdata <- DTalive.infected[ ,Prevalence := sum(Infected) / sum(!is.na(Infected)),by = "Gender,Age"] # And we calculate HIV prevalence, by gender and age group
    #setnames(Prevalencetimestepdata, "V1", "Prevalence")
    Prevalencetimestepdata <- cbind(time_i, Prevalencetimestepdata)
    output <- rbind(output, Prevalencetimestepdata)
  }
  return(output)
}

prevalenceheatmapdata.list <- function(datalist){ #(DT = datalist$ptable, cfg = configfile){
  DT <- datalist$ptable
  DTL <- datalist$ltable
  output <- data.table()
  timestep.width <- 1
  decplaces <- decimalplaces(timestep.width)
  timesteps <- DTL$Time
  for (time_i in timesteps){
    DTalive.infected <- alive.infected(DT = datalist$ptable, time = time_i, dec.places = decplaces) # First we only take the data of people who were alive at time_i
    DTalive.infected$Age5 <- cut_interval(DTalive.infected$Age, length = 5, right = FALSE)
    Prevalencetimestepdata <- DTalive.infected[ ,Prevalence := sum(Infected) / sum(!is.na(Infected)),by = "Gender,Age"] # And we calculate HIV prevalence, by gender and age group
    Prevalencetimestepdata <- DTalive.infected[ ,Prevalence5 := sum(Infected) / sum(!is.na(Infected)),by = "Gender,Age5"]
    #setnames(Prevalencetimestepdata, "V1", "Prevalence")
    Prevalencetimestepdata <- cbind(time_i, Prevalencetimestepdata)
    output <- rbind(output, Prevalencetimestepdata)
  }
  return(output)
}

incidenceheatmapdata <- function(modeloutput){ #DTP = datalist$ptable, DTE = datalist$etable, cfg = configfile){
  datalist <- readthedata(modeloutput)
  DTP <- datalist$ptable
  DTE <- datalist$etable
  cfg <- modeloutput[["input.params"]][["config"]]
  output <- data.table()
  timestep.width <- 1
  decplaces <- decimalplaces(timestep.width)
  timesteps <- seq(as.numeric(cfg["hivseed.time"]), as.numeric(cfg["population.simtime"]), by = timestep.width)
  for (time_i in head(timesteps, -1)){
    DTcases <- incidentcase(DT = DTE, interval = c(time_i, (time_i + timestep.width)))  # New infections that happened at time x with: time_i < x <= time_(i+1)
    DTnonHIVdeaths <- nonHIVdeath(DT = DTE, interval = c(time_i, (time_i + timestep.width))) # New non-HIV deaths that happened at time x with: time_i < x <= time_(i+1)
    DTcases$ID <- DTcases$p2ID
    DTcases$Gender <- DTcases$p2gender
    DTnonHIVdeaths$ID <- DTnonHIVdeaths$p1ID
    DTnonHIVdeaths$Gender <- DTnonHIVdeaths$p1gender
    
    # Personyears at risk during the interval
    DTalive.infected <- alive.infected(DT = DTP, time = time_i, dec.places = decplaces) # First we only take the data of people who were alive at time_i
    DTalive.uninfected <- DTalive.infected[Infected == "FALSE"] # We only keep those that were not infected at the start of the interval
    # A merge with etable will find those that died and/or got infected during the interval
    # To do this, we need to add a key to DTalive.uninfected
    DTPkeycols <- c("ID", "Gender")
    DTcaseskeycols <- c("ID", "Gender") # The ID and gender of the newly infected person  
    DTnonHIVdeathskeycols <- c("ID", "Gender") # The ID and gender of the newly died person (non-HIV death)
    setkeyv(DTalive.uninfected, DTPkeycols)
    setkeyv(DTcases, DTcaseskeycols)
    setkeyv(DTnonHIVdeaths, DTnonHIVdeathskeycols)
    DTforincidencecalc_a <- merge(DTalive.uninfected, DTcases, all.x=TRUE)
    DTforincidencecalc <- merge(DTforincidencecalc_a, DTnonHIVdeaths, all.x=TRUE, suffixes = c(".incident", ".dead"))
    DTforincidencecalc$eventtime.incident[is.na(DTforincidencecalc$eventtime.incident)] <- time_i + timestep.width
    DTforincidencecalc$eventtime.dead[is.na(DTforincidencecalc$eventtime.dead)] <- time_i + timestep.width
    
    DTforincidencecalc$PYend <- pmin(DTforincidencecalc$eventtime.incident, DTforincidencecalc$eventtime.dead)
    DTforincidencecalc$PY <- DTforincidencecalc$PYend - time_i
    
    Incidencetimestepdata <- DTforincidencecalc[ ,Incidence_G_A := sum(eventname.incident=="transmission", na.rm = TRUE) / sum(PY),by = "Gender,Age"] # And we calculate HIV Incidence, by gender and age group
    Incidencetimestepdata <- DTforincidencecalc[ ,PY_G_A := sum(PY),by = "Gender,Age"] # And we calculate PY, by gender and age group
    Incidencetimestepdata <- DTforincidencecalc[ ,Incidence_G := sum(eventname.incident=="transmission", na.rm = TRUE) / sum(PY),by = "Gender"] # And we calculate HIV Incidence, by gender and age group
    Incidencetimestepdata <- DTforincidencecalc[ ,PY_G := sum(PY),by = "Gender"] # And we calculate PY, by gender and age group
    Incidencetimestepdata <- DTforincidencecalc[ ,Incidence_All := sum(eventname.incident=="transmission", na.rm = TRUE) / sum(PY)] # And we calculate HIV Incidence, by gender and age group
    Incidencetimestepdata <- DTforincidencecalc[ ,PY_All := sum(PY)] # And we calculate PY, by gender and age group
    # (so we can be selective when plotting, and not plot if PY per age-gender group is too low)
    # setnames(Incidencetimestepdata, "V1", "Incidence") # Cases per personyear at risk
    Incidencetimestepdata <- cbind(time_i, Incidencetimestepdata)
    output <- rbind(output, Incidencetimestepdata)
  }
  incid.table <- unique(subset(output, select = c(time_i, Incidence_All, PY_All)))
  incid.table$cases <- incid.table$Incidence_All * incid.table$PY_All
  mean.incidence <- sum(incid.table$cases) / sum(incid.table$PY_All)
  cumul.incidence <- 1 - exp(-mean.incidence*(length(timesteps)-1))
  outputlist <- list(output=output,
                     incid.table=incid.table,
                     mean.incidence=mean.incidence,
                     cumul.incidence=cumul.incidence)
  return(outputlist)
}

AgeGap.last <- function(datatable){
  tail(datatable$AgeGap, 1)
}

# Let's analyse the HIV incidence data as in the Harling 2014 JAIDS paper, using Cox proportional hazard models, but let's also do a non-parametric estimation.
# If we are going to use the tmerge and neardate functions from the survival package, we need to several create datasets, which we will merge.

# For each gender, we prepare a set of dataframes.
# The set mimics the data obtained from a scenario in which all individuals test annually for HIV.
# They also complete a survey at each testing round and provide the age of the most recently acquired new partner (MRP),
# their own age, and whether they had multiple partners in the last year.
# We can also record the age of the new partner that was acquired most recently before HIV acquisition,
# which may be after the self-reported MRP in the survey they took at their last negative HIV test.
# This last variable can be used to measure how much of the signal is lost by using MRP at the annual surveys
# as a proxy for who is likely to be the person that caused the HIV infection.

# m.baseline.df: one row per person at risk of HIV infection, with baseline, non-time-varying covariates.
# These are the unique identifier (IDm) and the birth date (TOB), from which we can derive Calendar year.
# m.timeline.df: one row per period of exposure time. We can use 1-year intervals as periods of exposure
# except when there is a non-AIDS death or a transmission event, the interval would be shorter.
# m.timeline.df should also include an indicator variable for the (Calendar) time of observation.
# IDm, age0, age1, time0, time1.
# m.outcome.hiv.postest.df: one row with the time of first HIV+ test ##newname = event(y,x) Mark an event at time y, with value x. We can choose x=3
# InfectTime, status (==3)
# m.outcome.hiv.pos.df: one row with the exact time of HIV acquisition ##newname = event(y,x) Mark an event at time y, with value x. We can choose x=1 
# InfectTime, status (==1)
# m.outcome.survey.ageMRP.df: one row with the age of the MRP at the time of the survey ##newname = tdc(y, x)
# The argument y is assumed to be on the scale of the start and end time, if x is present the count is set to the value of x (e.g. age of MRP). 
# SurveyTime, survey.ageMRP
# m.outcome.infect.ageMRP.df: one row with the age of the MRP at the time of the infection ##newname = tdc(y, x)
# The argument y is assumed to be on the scale of the start and end time, if x is present the count is set to the value of x (e.g. age of MRP). 
# InfectTime, infect.ageMRP

# 1. baseline.df: Everybody who ever existed in the population after HIV was introduced, minus the hiv seeds themselves.
baseline.datamaker <- function(datalist = datalist, hivseed.time = as.numeric(cfgABC["hivseed.time"])){
  DTP <- datalist$ptable
  baseline.df <- subset(DTP, TOD > hivseed.time & InfectType != 0) # This has more variables than we strictly need, but we'll leave this for now.
  baseline.df$ID.Gender <- paste(baseline.df$ID, baseline.df$Gender, sep=".")
  return(baseline.df)
}
# 2. timeline.df: one row per period of exposure time. We can use 1-year intervals as periods of exposure
# except when there is a non-AIDS death or a transmission event, the interval would be shorter.
# m.timeline.df should also include an indicator variable for the (Calendar) time of observation.
# IDm, age0, age1, time0, time1.
timeline.datamaker <- function(datalist = datalist,
                               endfollowup = as.numeric(cfgABC["population.simtime"]),
                               hivseed.time = as.numeric(cfgABC["hivseed.time"])){
  baseline.df <- baseline.datamaker(datalist, hivseed.time)
  exposure.end <- pmin(pmin(baseline.df$TOD, endfollowup), ceiling(baseline.df$InfectTime)) # End of exposure time: given by death, end of simulation, or (when) infection (may be picked up by an annual HIV test).
  exposure.start <- pmin(baseline.df$TODebut, exposure.end)
  timeline.df <- data.table(ID.Gender=baseline.df$ID.Gender, TOB=baseline.df$TOB, exposure.start, exposure.end)
  someexposure.indicator <- exposure.start != exposure.end # No exposure time is not allowed
  timeline.df <- subset(timeline.df, someexposure.indicator) # No exposure time is not allowed
  baseline.df <- subset(baseline.df, someexposure.indicator) # No exposure time is not allowed
  baseline.timeline.df <- tmerge(data1 = baseline.df, data2 = timeline.df, id = ID.Gender, tstart = timeline.df$exposure.start, tstop = timeline.df$exposure.end)
  return(baseline.timeline.df)
}
# 3. outcome.df
# m.outcome.hiv.postest.df: one row with the time of first HIV+ test ##newname = event(y,x) Mark an event at time y, with value x. We can choose x=3
# InfectTime, status (==3)
# m.outcome.hiv.infect.df: one row with the exact time of HIV acquisition ##newname = event(y,x) Mark an event at time y, with value x. We can choose x=1 
# InfectTime, status (==1)
# m.outcome.survey.ageMRP.df: one row with the age of the MRP at the time of the survey ##newname = tdc(y, x)
# The argument y is assumed to be on the scale of the start and end time, if x is present the count is set to the value of x (e.g. age of MRP). 
# SurveyTime, survey.ageMRP
# m.outcome.infect.ageMRP.df: one row with the age of the MRP at the time of the infection ##newname = tdc(y, x)
# The argument y is assumed to be on the scale of the start and end time, if x is present the count is set to the value of x (e.g. age of MRP). 
# InfectTime, infect.ageMRP
outcome.datamaker <- function(datalist = datalist,
                              endfollowup = as.numeric(cfgABC["population.simtime"]),
                              hivseed.time = as.numeric(cfgABC["hivseed.time"])){
  baseline.timeline.df <- timeline.datamaker(datalist, endfollowup, hivseed.time)
  DTL <- datalist$ltable
  DTR <- datalist$rtable
  DTE <- datalist$etable
  oneyear.intervals.end <- tail(datalist$ltable$Time, -1) # c(1, 2, 3, ..., population.simtime) It does not make sense to have a survey at time=0 because nothing happened before that time.
  outcome.survey <- expand.grid(ID.Gender = baseline.timeline.df$ID.Gender, oneyear.intervals.end = oneyear.intervals.end)
  
  
  Surv.df <- tmerge(data1 = baseline.timeline.df, data2 = outcome.survey, id = ID.Gender, survey.number = event(y = oneyear.intervals.end, x = oneyear.intervals.end)) # adding survey rounds
  Surv.df <- tmerge(data1 = Surv.df, data2 = Surv.df, id = ID.Gender, hiv.infect = event(y = InfectTime)) # adding HIV infection events
  Surv.df$HIVposTime <- ceiling(Surv.df$InfectTime)
  Surv.df <- tmerge(data1 = Surv.df, data2 = Surv.df, id = ID.Gender, hiv.postest = event(y = HIVposTime, x = rep(3, nrow(Surv.df)))) # adding first HIV positive test events
  
  #outcome.survey.ageMRP.df: one row with the age of the MRP at the time of the survey ##newname = tdc(y, x)
  # The argument y is assumed to be on the scale of the start and end time, if x is present the count is set to the value of x (e.g. age of MRP). 
  DTR.m <- DTR.w <- DTR
  DTR.m$ID.Gender <- paste(DTR$IDm, 0, sep = ".")
  DTR.w$ID.Gender <- paste(DTR$IDw, 1, sep = ".")
  # index of MRP
  indx.m <- neardate(id1 = Surv.df$ID.Gender, id2 = DTR.m$ID.Gender, y1 = Surv.df$survey.number, y2 = DTR.m$FormTime, best = "prior")
  # WE MUST NOT FORGET TO CHECK IF THE PARTICIPANT WAS IN A RELATIONSHIP WITH THE MRP AT THE TIME OF THE SURVEY.
  # Harling et al. did not do this check.
  indx.m <- ifelse((DTR.m$DisTime[indx.m] - Surv.df$survey.number) > 0, indx.m, NA) # Now we know the relationship with the MRP was ongoing.
  indx.w <- neardate(id1 = Surv.df$ID.Gender, id2 = DTR.w$ID.Gender, y1 = Surv.df$survey.number, y2 = DTR.w$FormTime, best = "prior")
  # WE MUST NOT FORGET TO CHECK IF THE PARTICIPANT WAS IN A RELATIONSHIP WITH THE MRP AT THE TIME OF THE SURVEY.
  # Harling et al. did not do this check.
  indx.w <- ifelse((DTR.w$DisTime[indx.w] - Surv.df$survey.number) > 0, indx.w, NA) # Now we know the relationship with the MRP was ongoing.
  AgeGaps.m <- DTR.m[indx.m, AgeGap]
  AgeGaps.w <- DTR.w[indx.w, AgeGap]
  AgeGaps.m.and.w <- AgeGaps.m
  AgeGaps.m.and.w[is.na(AgeGaps.m.and.w) & !is.na(AgeGaps.w)] <- AgeGaps.w[is.na(AgeGaps.m.and.w) & !is.na(AgeGaps.w)]
  outcome.survey.agegapMRP.df <- data.table(ID.Gender= Surv.df$ID.Gender, survey.number=Surv.df$survey.number, AgeGapMRP=AgeGaps.m.and.w)
  TestSurv.df <- tmerge(data1 = Surv.df, data2 = outcome.survey.ageMRP.df, id = ID.Gender,
                    AgeGapMRP = tdc(y = outcome.survey.agegapMRP.df$survey.number, x = outcome.survey.agegapMRP.df$AgeGapMRP))
  
  # Preparing data for "MIICD" coxph model with interval censored data
  #   The data must contain at last two columns: left and right. For interval censored data, the left
  #   and the right columns indicates lower and upper bounds of intervals respectively. Inf in the right
  #   column stands for right censored observations.
  
  # First, we forget about the hiv.infect event
  #TestSurv.df <- subset(TestSurv.df, hiv.infect != 1)
  
  # we treat each observation as a right censored one, except the ones that led to hiv.postest===3
  #TestSurv.df$left <- TestSurv.df$tstop - TestSurv.df$tstart
  #TestSurv.df$right <- Inf
  #TestSurv.df$right[TestSurv.df$hiv.postest==3] <- 1
  
  
  # tmerge.dfs <- list(baseline.timeline.df=baseline.timeline.df,
  #                    outcome.survey=outcome.survey)
  return(TestSurv.df)
}















incidence.Harling <- function(datalist = datalist, endtime = as.numeric(cfgABC["population.simtime"]), window = 10){
  DTP <- datalist$ptable
  DTR <- datalist$rtable
  DTE <- datalist$etable
  timepoint <- endtime
  time_i <- timepoint - window # We look at the incidence in the window before the timepoint. The window can be much wider than 1 year.
  timestep.width <- window
  
  DTcases <- incidentcase(DT = DTE, interval = c(time_i, (time_i + timestep.width)))  # New infections that happened at time x with: time_i < x <= time_(i+1)
  DTnonHIVdeaths <- nonHIVdeath(DT = DTE, interval = c(time_i, (time_i + timestep.width))) # New non-HIV deaths that happened at time x with: time_i < x <= time_(i+1)
  DTcases$ID <- DTcases$p2ID
  DTcases$Gender <- DTcases$p2gender
  DTnonHIVdeaths$ID <- DTnonHIVdeaths$p1ID
  DTnonHIVdeaths$Gender <- DTnonHIVdeaths$p1gender
  
  # Personyears at risk during the interval
  DTalive.infected <- alive.infected(DT = DTP, time = time_i, dec.places = 0) # First we only take the data of people who were alive at time_i
  DTalive.uninfected <- DTalive.infected[Infected == "FALSE"] # We only keep those that were not infected at the start of the interval
  # A merge with etable will find those that died and/or got infected during the interval
  # To do this, we need to add a key to DTalive.uninfected
  DTPkeycols <- c("ID", "Gender")
  DTcaseskeycols <- c("ID", "Gender") # The ID and gender of the newly infected person  
  DTnonHIVdeathskeycols <- c("ID", "Gender") # The ID and gender of the newly died person (non-HIV death)
  setkeyv(DTalive.uninfected, DTPkeycols)
  setkeyv(DTcases, DTcaseskeycols)
  setkeyv(DTnonHIVdeaths, DTnonHIVdeathskeycols)
  DTforincidencecalc_a <- merge(DTalive.uninfected, DTcases, all.x=TRUE)
  DTforincidencecalc <- merge(DTforincidencecalc_a, DTnonHIVdeaths, all.x=TRUE, suffixes = c(".incident", ".dead"))
  DTforincidencecalc$eventtime.incident[is.na(DTforincidencecalc$eventtime.incident)] <- time_i + timestep.width
  DTforincidencecalc$eventtime.dead[is.na(DTforincidencecalc$eventtime.dead)] <- time_i + timestep.width
  
  DTforincidencecalc$PYend <- pmin(DTforincidencecalc$eventtime.incident, DTforincidencecalc$eventtime.dead)
  DTforincidencecalc$PY <- DTforincidencecalc$PYend - time_i
  DTforincidencecalc.Male <- subset(DTforincidencecalc, Gender==0)
  DTforincidencecalc.Female <- subset(DTforincidencecalc, Gender==1)
  
  DTforincidencecalc.Male$IDm <- DTforincidencecalc.Male$ID
  DTforincidencecalc.Female$IDw <- DTforincidencecalc.Female$ID
  setkey(DTforincidencecalc.Male, IDm)
  setkey(DTforincidencecalc.Female, IDw)
  
  # The data in DTforincidencecalc can give us the age differences in the relationships that led to transmission.
  # But we also need the age differences of all other relationships that were ongoing in the year before the timepoint.
  
  OngoingRels.lastyear <- DTR[FormTime <= timepoint & DisTime > time_i]
  OngoingRels.lastyear.Male.maxAgeGap <- data.table(ddply(OngoingRels.lastyear, .(IDm), summarise, max.AgeGap = max(AgeGap, na.rm=TRUE)))
  OngoingRels.lastyear.Male.lastAgeGap <- data.table(ddply(OngoingRels.lastyear, .(IDm), .fun = AgeGap.last))
  Partners.lastyear.Male <- daply(OngoingRels.lastyear, .(IDm), nrow)
  OngoingRels.lastyear.Male.lastAgeGap$Partners.lastyear.Male <- Partners.lastyear.Male
  
  OngoingRels.lastyear.Female.maxAgeGap <- data.table(ddply(OngoingRels.lastyear, .(IDw), summarise, max.AgeGap = max(AgeGap, na.rm=TRUE)))
  OngoingRels.lastyear.Female.lastAgeGap <- data.table(ddply(OngoingRels.lastyear, .(IDw), .fun = AgeGap.last))
  Partners.lastyear.Female <- daply(OngoingRels.lastyear, .(IDw), nrow)
  OngoingRels.lastyear.Female.lastAgeGap$Partners.lastyear.Female <- Partners.lastyear.Female
  
  OngoingRels.lastyear.Male.max.last <- merge(OngoingRels.lastyear.Male.maxAgeGap, OngoingRels.lastyear.Male.lastAgeGap, by="IDm")
  OngoingRels.lastyear.Female.max.last <- merge(OngoingRels.lastyear.Female.maxAgeGap, OngoingRels.lastyear.Female.lastAgeGap, by="IDw")
  setnames(OngoingRels.lastyear.Male.max.last, "V1", "last.AgeGap")
  setnames(OngoingRels.lastyear.Female.max.last, "V1", "last.AgeGap")
  
  setkey(OngoingRels.lastyear.Male.max.last, IDm) # These are the men that had ongoing relationships in the year before the time point
  setkey(OngoingRels.lastyear.Female.max.last, IDw) # These are the women that had ongoing relationships in the year before the time point
  
  # We add the relevant DTforincidencecalc.x data to the OngoingRels.lastyear.x.max.last dataset
  readyforInc.Male <- merge(OngoingRels.lastyear.Male.max.last, DTforincidencecalc.Male)
  readyforInc.Female <- merge(OngoingRels.lastyear.Female.max.last, DTforincidencecalc.Female)
  readyforInc.Male$MultiplePartners.lastyear <- readyforInc.Male$Partners.lastyear.Male > 1
  readyforInc.Female$MultiplePartners.lastyear <- readyforInc.Female$Partners.lastyear.Female > 1
  
  
  incidence.Male.Age.AgeGap <- readyforInc.Male[ ,Incidence.Age.maxAgeGap := sum(eventname.incident=="transmission", na.rm = TRUE) / sum(PY),by = "Age,max.AgeGap"] # And we calculate HIV Incidence, by age group and age gap 
  incidence.Female.Age.AgeGap <- readyforInc.Female[ ,Incidence.Age.maxAgeGap := sum(eventname.incident=="transmission", na.rm = TRUE) / sum(PY),by = "Age,max.AgeGap"] # And we calculate HIV Incidence, by age group and age gap 
  incidence.Male.Age.AgeGap <- readyforInc.Male[ ,Incidence.Age.lastAgeGap := sum(eventname.incident=="transmission", na.rm = TRUE) / sum(PY),by = "Age,last.AgeGap"] # And we calculate HIV Incidence, by age group and age gap 
  incidence.Female.Age.AgeGap <- readyforInc.Female[ ,Incidence.Age.lastAgeGap := sum(eventname.incident=="transmission", na.rm = TRUE) / sum(PY),by = "Age,last.AgeGap"] # And we calculate HIV Incidence, by age group and age gap 
  # incidence.Male.Age.AgeGap and incidence.Female.Age.AgeGap can be used to make contour plots
  # of the real HIV incidence, by age group and age gap.
  
  # Now let's prepare a dataset for Cox proportional hazard regression
  # Each individual at risk of HIV acquisition at the beginning of the window is a row.
  # survival time is either time of HIV infection, death, or end of the window (righ-censored).
  # In the simple analysis, we treat death as a censoring event. One should treat it as a 
  # competing risk, but the rate of non-AIDS death is so low that this is unlikely to
  # introduce much bias.
  readyforCox.Male <- subset(readyforInc.Male, select=c("IDm", "PY", "Age", "max.AgeGap", "last.AgeGap",
                                                        "Partners.lastyear.Male", "MultiplePartners.lastyear",
                                                        "eventname.incident", "eventname.dead", "InfectTime"))
  
  # Preparing data for "normal" coxph model with righ censored data
  readyforCox.Male$event <- NA
  readyforCox.Male$time <- NA
  readyforCox.Male$time2 <- NA
  readyforCox.Male$event[is.na(readyforCox.Male$eventname.incident)] <- 0 # Event is 2 for all censored observations: natural deaths and completed windows without infection.
  readyforCox.Male$time[is.na(readyforCox.Male$eventname.incident)] <- readyforCox.Male$PY[is.na(readyforCox.Male$eventname.incident)] # Event is 2 for all censored observations: natural deaths and completed windows without infection.
  readyforCox.Male$event[!is.na(readyforCox.Male$eventname.incident)] <- 1 # Event is 3 for all interval censored observations: HIV infection.
  readyforCox.Male$time[!is.na(readyforCox.Male$eventname.incident)] <- readyforCox.Male$InfectTime[!is.na(readyforCox.Male$eventname.incident)] - endtime + window
  readyforCox.Male$time2[!is.na(readyforCox.Male$eventname.incident)] <- readyforCox.Male$InfectTime[!is.na(readyforCox.Male$eventname.incident)] - endtime + window
  # Preparing data for "MIICD" coxph model with interval censored data
#   The data must contain at last two columns: left and right. For interval censored data, the left
#   and the right columns indicates lower and upper bounds of intervals respectively. Inf in the right
#   column stands for right censored observations.
  readyforCox.Male$left <- NA
  readyforCox.Male$right <- Inf # For right censored observations
  readyforCox.Male$left[is.na(readyforCox.Male$eventname.incident)] <- readyforCox.Male$time[is.na(readyforCox.Male$eventname.incident)] # For right censored observations
  readyforCox.Male$left[!is.na(readyforCox.Male$eventname.incident)] <- floor(readyforCox.Male$InfectTime[!is.na(readyforCox.Male$eventname.incident)]) - endtime + window
  readyforCox.Male$right[!is.na(readyforCox.Male$eventname.incident)] <- readyforCox.Male$left[!is.na(readyforCox.Male$eventname.incident)] + 1 #HIV testing is assumed to happen with 1-year intervals
  
  
  readyforCox.Female <- subset(readyforInc.Female, select=c("IDw", "PY", "Age", "max.AgeGap", "last.AgeGap",
                                                        "Partners.lastyear.Female", "MultiplePartners.lastyear",
                                                        "eventname.incident", "eventname.dead", "InfectTime"))
#   readyforCox.Female$left <- 0
#   readyforCox.Female$right <- NA
#   readyforCox.Female$right[is.na(readyforCox.Female$eventname.incident)] <- Inf
#   readyforCox.Female$right[!is.na(readyforCox.Female$eventname.incident)] <- window
  
  
  # Preparing data for "normal" coxph model with righ censored data
  readyforCox.Female$event <- NA
  readyforCox.Female$time <- NA
  readyforCox.Female$time2 <- NA
  readyforCox.Female$event[is.na(readyforCox.Female$eventname.incident)] <- 0 # Event is 2 for all censored observations: natural deaths and completed windows without infection.
  readyforCox.Female$time[is.na(readyforCox.Female$eventname.incident)] <- readyforCox.Female$PY[is.na(readyforCox.Female$eventname.incident)] # Event is 2 for all censored observations: natural deaths and completed windows without infection.
  readyforCox.Female$event[!is.na(readyforCox.Female$eventname.incident)] <- 1 # Event is 3 for all interval censored observations: HIV infection.
  readyforCox.Female$time[!is.na(readyforCox.Female$eventname.incident)] <- readyforCox.Female$InfectTime[!is.na(readyforCox.Female$eventname.incident)] - endtime + window
  readyforCox.Female$time2[!is.na(readyforCox.Female$eventname.incident)] <- readyforCox.Female$InfectTime[!is.na(readyforCox.Female$eventname.incident)] - endtime + window
  # Preparing data for "MIICD" coxph model with interval censored data
  #   The data must contain at last two columns: left and right. For interval censored data, the left
  #   and the right columns indicates lower and upper bounds of intervals respectively. Inf in the right
  #   column stands for right censored observations.
  readyforCox.Female$left <- NA
  readyforCox.Female$right <- Inf # For right censored observations
  readyforCox.Female$left[is.na(readyforCox.Female$eventname.incident)] <- readyforCox.Female$time[is.na(readyforCox.Female$eventname.incident)] # For right censored observations
  readyforCox.Female$left[!is.na(readyforCox.Female$eventname.incident)] <- floor(readyforCox.Female$InfectTime[!is.na(readyforCox.Female$eventname.incident)]) - endtime + window
  readyforCox.Female$right[!is.na(readyforCox.Female$eventname.incident)] <- readyforCox.Female$left[!is.na(readyforCox.Female$eventname.incident)] + 1 #HIV testing is assumed to happen with 1-year intervals
  
  
  Harling.output <- list(incidence.Male.Age.AgeGap=incidence.Male.Age.AgeGap,
                         incidence.Female.Age.AgeGap=incidence.Female.Age.AgeGap,
                         readyforCox.Male=readyforCox.Male,
                         readyforCox.Female=readyforCox.Female)
  return(Harling.output)
}



incidenceheatmapdata.list <- function(datalist){ #DTP = datalist$ptable, DTE = datalist$etable, cfg = configfile){
  DTP <- datalist$ptable
  DTE <- datalist$etable
  DTL <- datalist$ltable
  #cfg <- modeloutput[["input.params"]][["config"]]
  output <- data.table()
  timestep.width <- 1
  decplaces <- decimalplaces(timestep.width)
  timesteps <- DTL$Time #seq(as.numeric(cfg["hivseed.time"]), as.numeric(cfg["population.simtime"]), by = timestep.width)
  for (time_i in head(timesteps, -1)){
    DTcases <- incidentcase(DT = DTE, interval = c(time_i, (time_i + timestep.width)))  # New infections that happened at time x with: time_i < x <= time_(i+1)
    DTnonHIVdeaths <- nonHIVdeath(DT = DTE, interval = c(time_i, (time_i + timestep.width))) # New non-HIV deaths that happened at time x with: time_i < x <= time_(i+1)
    DTcases$ID <- DTcases$p2ID
    DTcases$Gender <- DTcases$p2gender
    DTnonHIVdeaths$ID <- DTnonHIVdeaths$p1ID
    DTnonHIVdeaths$Gender <- DTnonHIVdeaths$p1gender
    
    # Personyears at risk during the interval
    DTalive.infected <- alive.infected(DT = DTP, time = time_i, dec.places = decplaces) # First we only take the data of people who were alive at time_i
    DTalive.uninfected <- DTalive.infected[Infected == "FALSE"] # We only keep those that were not infected at the start of the interval
    # A merge with etable will find those that died and/or got infected during the interval
    # To do this, we need to add a key to DTalive.uninfected
    DTPkeycols <- c("ID", "Gender")
    DTcaseskeycols <- c("ID", "Gender") # The ID and gender of the newly infected person  
    DTnonHIVdeathskeycols <- c("ID", "Gender") # The ID and gender of the newly died person (non-HIV death)
    setkeyv(DTalive.uninfected, DTPkeycols)
    setkeyv(DTcases, DTcaseskeycols)
    setkeyv(DTnonHIVdeaths, DTnonHIVdeathskeycols)
    DTforincidencecalc_a <- merge(DTalive.uninfected, DTcases, all.x=TRUE)
    DTforincidencecalc <- merge(DTforincidencecalc_a, DTnonHIVdeaths, all.x=TRUE, suffixes = c(".incident", ".dead"))
    DTforincidencecalc$eventtime.incident[is.na(DTforincidencecalc$eventtime.incident)] <- time_i + timestep.width
    DTforincidencecalc$eventtime.dead[is.na(DTforincidencecalc$eventtime.dead)] <- time_i + timestep.width
    
    DTforincidencecalc$PYend <- pmin(DTforincidencecalc$eventtime.incident, DTforincidencecalc$eventtime.dead)
    DTforincidencecalc$PY <- DTforincidencecalc$PYend - time_i
    
    Incidencetimestepdata <- DTforincidencecalc[ ,Incidence_G_A := sum(eventname.incident=="transmission", na.rm = TRUE) / sum(PY),by = "Gender,Age"] # And we calculate HIV Incidence, by gender and age group
    Incidencetimestepdata <- DTforincidencecalc[ ,PY_G_A := sum(PY),by = "Gender,Age"] # And we calculate PY, by gender and age group
    Incidencetimestepdata <- DTforincidencecalc[ ,Incidence_G := sum(eventname.incident=="transmission", na.rm = TRUE) / sum(PY),by = "Gender"] # And we calculate HIV Incidence, by gender and age group
    Incidencetimestepdata <- DTforincidencecalc[ ,PY_G := sum(PY),by = "Gender"] # And we calculate PY, by gender and age group
    Incidencetimestepdata <- DTforincidencecalc[ ,Incidence_All := sum(eventname.incident=="transmission", na.rm = TRUE) / sum(PY)] # And we calculate HIV Incidence, by gender and age group
    Incidencetimestepdata <- DTforincidencecalc[ ,PY_All := sum(PY)] # And we calculate PY, by gender and age group
    # (so we can be selective when plotting, and not plot if PY per age-gender group is too low)
    # setnames(Incidencetimestepdata, "V1", "Incidence") # Cases per personyear at risk
    Incidencetimestepdata <- cbind(time_i, Incidencetimestepdata)
    output <- rbind(output, Incidencetimestepdata)
  }
  incid.table <- unique(subset(output, select = c(time_i, Incidence_All, PY_All)))
  incid.table$cases <- incid.table$Incidence_All * incid.table$PY_All
  mean.incidence <- sum(incid.table$cases) / sum(incid.table$PY_All)
  cumul.incidence <- 1 - exp(-mean.incidence*(length(timesteps)-1))
  outputlist <- list(output=output,
                     incid.table=incid.table,
                     mean.incidence=mean.incidence,
                     cumul.incidence=cumul.incidence)
  return(outputlist)
}

population <- function(modeloutput){
  datalist <- readthedata(modeloutput)
  DTP <- datalist$ptable
  output <- data.table()
  timestep.width <- 1
  decplaces <- decimalplaces(timestep.width)
  timesteps <- seq(as.numeric(cfg["hivseed.time"]), as.numeric(cfg["population.simtime"]), by = timestep.width)
  for (time_i in head(timesteps, -1)){
    DTalive.infected <- alive.infected(DT = DTP, time = time_i, dec.places = decplaces) # First we only take the data of people who were alive at time_i
    popsize <- nrow(DTalive.infected)
    pointprevalence <- sum(DTalive.infected$Infected) / popsize
    popdata <- data.table(cbind(time_i, popsize, pointprevalence))
    output <- rbind(output, popdata)
  }
  #outputlist <- list(output=output)
  #return(outputlist)
  return(output)
}

population.list <- function(datalist){ #DTP = datalist$ptable, DTE = datalist$etable, cfg = configfile){
  DTP <- datalist$ptable
  DTL <- datalist$ltable
  timestep.width <- 1
  decplaces <- decimalplaces(timestep.width)
  timesteps <- DTL$Time #seq(as.numeric(cfg["hivseed.time"]), as.numeric(cfg["population.simtime"]), by = timestep.width)
  nrows.output <- length(timesteps) - 1
  output <- data.table(time = rep(-1, nrows.output), popsize = rep(-1, nrows.output), pointprevalence = rep(-1, nrows.output))
  for (time_i in head(timesteps, -1)){
    DTalive.infected <- alive.infected(DT = DTP, time = time_i, dec.places = decplaces) # First we only take the data of people who were alive at time_i
    popsize <- nrow(DTalive.infected)
    pointprevalence <- sum(DTalive.infected$Infected) / popsize
    ifelse(length(pointprevalence) > 0, sum(DTalive.infected$Infected) / popsize, -1) # pointprevalence of -1 if the pointprevalence does not exist
    #popdata <- data.table(cbind(time_i, popsize, pointprevalence))
    output[(time_i+1), ] <- data.table(time_i, popsize, pointprevalence)
  }
  #outputlist <- list(output=output)
  #return(outputlist)
  return(output)
}

prevalence.18.49.after.15.30.list <- function(datalist = datalist,
                                              timepoints = c(as.numeric(results$input.params$config["hivseed.time"])+15,
                                                             as.numeric(results$input.params$config["hivseed.time"])+30)){
  DTP <- datalist$ptable
  output <- data.table()
  for (time_i in timepoints){
    DTalive.infected <- alive.infected(DT = DTP, time = time_i, dec.places = 0) # First we only take the data of people who were alive at time_i
    DTalive.infected.18.49 <- subset(DTalive.infected, TOB <= time_i - 18 & TOB > time_i - 50)
    popsize <- nrow(DTalive.infected.18.49)
    pointprevalence <- sum(DTalive.infected.18.49$Infected) / popsize
    popdata <- data.table(cbind(time_i, popsize, pointprevalence))
    output <- rbind(output, popdata)
  }
  #outputlist <- list(output=output)
  #return(outputlist)
  return(output)
}

concurrency.list <- function(datalist){ #DTP = datalist$ptable, DTE = datalist$etable, cfg = configfile){
  DTP <- datalist$ptable
  DTL <- datalist$ltable
  output <- data.table()
  timestep.width <- 1
  decplaces <- decimalplaces(timestep.width)
  timesteps <- DTL$Time #seq(as.numeric(cfg["hivseed.time"]), as.numeric(cfg["population.simtime"]), by = timestep.width)
  for (time_i in timesteps){
    DTalive.infected <- alive.infected(DT = DTP, time = time_i, dec.places = decplaces) # First we only take the data of people who were alive at time_i
    DT.degree <- degreecalc(DT = DTalive.infected, DTR = datalist$rtable, timepoint = time_i)# get instantaneous degree of these people
    popsize <- nrow(DTalive.infected)
    concurr.pointprevalence <- sum(DT.degree$concurrent) / popsize
    concurrdata <- data.table(cbind(time_i, concurr.pointprevalence))
    output <- rbind(output, concurrdata)
  }
  #outputlist <- list(output=output)
  #return(outputlist)
  return(output)
}

lifetime.partners.list <- function(datalist = datalist, timepoint = targetwindow.end){
  DTP <- datalist$ptable
  DTL <- datalist$ltable
  DTalive.infected <- alive.infected(DT = DTP, time = timepoint, dec.places = 0) # First we only take the data of people who were alive at timepoint
  DT.degree <- degreecalc(DT = DTalive.infected, DTR = datalist$rtable, timepoint = timepoint)# get lifetime degree of these people
  return(DT.degree)
}

# Calculate point prevalence of concurrency,
# Lifetime number of sex partners
# Cross-sectional degree distribution

degreecalc <- function(DT = datalist$ptable, DTR = datalist$rtable, timepoint = time_i){
  DT$ID.Gender <- paste(DT$ID, DT$Gender, sep=".")
  setkey(DT, ID.Gender)
  
  OngoingRels <- subset(DTR, FormTime <= timepoint & DisTime > timepoint)
  
  men_instant_degree <- data.table(count(OngoingRels, vars = "IDm"))

  men_degree <- data.table(data.frame(table(DTR$IDm)))
  if (nrow(men_degree) > 0) {
    men_degree$ID.Gender <- paste(men_degree$Var1, 0, sep=".")
    men_instant_degree$ID.Gender <- paste(men_instant_degree$IDm, 0, sep=".")
    #print(names(men_degree))
    setnames(men_degree, c("Var1", "Freq"), c("IDm", "degree"))
    setkey(men_degree, ID.Gender)
    setnames(men_instant_degree, "freq", "instantdegree")
    setkey(men_instant_degree, ID.Gender)
    
    women_instant_degree <- data.table(count(OngoingRels$IDw))
    women_degree <- data.table(data.frame(table(DTR$IDw)))
    women_degree$ID.Gender <- paste(women_degree$Var1, 1, sep=".")
    women_instant_degree$ID.Gender <- paste(women_instant_degree$IDw, 1, sep=".")
    
    setnames(women_degree, c("Var1", "Freq"), c("IDw", "degree"))
    setkey(women_degree, ID.Gender)
    setnames(women_instant_degree, "freq", "instantdegree")
    setkey(women_instant_degree, ID.Gender)
    
    DT_degree_MW <- merge(DT, men_degree, all.x=TRUE)
    DT_degree_W <- merge(DT, women_degree, all.x=TRUE)
    DT_degree_MW$degree[DT_degree_MW$Gender==1] <-  DT_degree_W$degree[DT_degree_W$Gender==1]
    DT_degree_MW$degree[is.na(DT_degree_MW$degree)] <- 0
    DT_degree_MW[ , IDm:=NULL] # before we add the instant_degree variable, we need to remove ID.m
    DT_degrees_MW <- merge(DT_degree_MW, men_instant_degree, all.x=TRUE)
    DT_degrees_W <- merge(DT_degree_MW, women_instant_degree, all.x=TRUE)
    DT_degrees_MW$instantdegree[DT_degrees_MW$Gender==1] <-  DT_degrees_W$instantdegree[DT_degrees_W$Gender==1]
    DT_degrees_MW$instantdegree[is.na(DT_degrees_MW$instantdegree)] <- 0
    DT_degrees_MW[ , IDm:=NULL] # cleaning things up: we remove ID.m
    DT_degrees_MW$concurrent <- DT_degrees_MW$instantdegree > 1
    result <- DT_degrees_MW
  } else {
    result <- data.table(concurrent = NA)
  }
 
  return(result)
}



partnerturnoverrate <- function(datalist){
  rels = datalist$rtable[FormTime < targetwindow.end & FormTime > targetwindow.start]
  nrel <- nrow(rels)
  halfnpeople <- nrow(subset(datalist$ptable, TOB <= targetwindow.midpoint-18 &
                               TOB > targetwindow.midpoint - 50 &
                               TOD > targetwindow.midpoint))/2
  relsperpersonperyear <- nrel / halfnpeople / targetwindow.duration
  return(relsperpersonperyear)
}

agegap.concur <- function(datalist = datalist){
  OngoingRels.targetwindow.end <- datalist$rtable[FormTime <= targetwindow.end & DisTime > targetwindow.end]
#   setkey(OngoingRels, IDm)
#   RelCount <- data.table(data.frame(table(OngoingRels$IDm)))
#   setnames(RelCount, c("IDm", "degree.man"))
#   setkey(RelCount, IDm)
#   RelCount$IDm <- as.numeric(RelCount$IDm)
  male.degree.end <- data.table(count(OngoingRels.targetwindow.end, "IDm"))
  setnames(male.degree.end, "freq", "male.degree.end")
  OngoingRels.targetwindow.end <- merge(OngoingRels.targetwindow.end, male.degree.end, by="IDm", all.x=TRUE)
  
  # Now we need to know which of the relationships that were ongoing at the targetwindow.end started out as a monogamous relationship
  degree.relat.start <- rep(NA, nrow(OngoingRels.targetwindow.end))
  for (relat in 1:nrow(OngoingRels.targetwindow.end)){
    formtime <- as.numeric(subset(OngoingRels.targetwindow.end[relat, ], select="FormTime"))
    ID.m <- as.numeric(subset(OngoingRels.targetwindow.end[relat, ], select="IDm"))
    degree.relat.start[relat] <- nrow(datalist$rtable[FormTime <= formtime & DisTime > formtime & IDm==ID.m])
  }
  
  OngoingRels.targetwindow.end <- cbind(OngoingRels.targetwindow.end, degree.relat.start)
  
#   OngoingRels.targetwindow.start <- datalist$rtable[FormTime <= targetwindow.start & DisTime > targetwindow.start]
#   male.degree.start <- data.table(count(OngoingRels.targetwindow.start, "IDm"))
#   setnames(male.degree.start, "freq", "male.degree.start")
#   OngoingRels.targetwindow.end <- merge(OngoingRels.targetwindow.end, male.degree.start, by="IDm", all.x=TRUE)
  
  monogam.mean.agegap <- mean(OngoingRels.targetwindow.end$AgeGap[OngoingRels.targetwindow.end$degree.relat.start==1]) #Relationships that started out as monogamous
  concur.mean.agegap <- mean(OngoingRels.targetwindow.end$AgeGap[OngoingRels.targetwindow.end$degree.relat.start > 1]) #Relationships that started out as concurrent
  agegap.diff <- concur.mean.agegap - monogam.mean.agegap
  return(agegap.diff)
}

Partner25ormore.girls <- function(datalist = datalist){
  OngoingRels.End <- datalist$rtable[FormTime <= targetwindow.end & DisTime > targetwindow.end] # - 1]
  Girls15to20yo <- subset(datalist$ptable, Gender==1 & TOB > (targetwindow.end-20) & TOB <= (targetwindow.end-15))
  Girls15to20yo <- Girls15to20yo[TOB <= targetwindow.end & TOD > targetwindow.end]
  
  Girls15to20yo$IDw <- Girls15to20yo$ID
  setkey(OngoingRels.End, IDw)
  setkey(Girls15to20yo, IDw)
  Girls15to20yo.AgeGaps <- merge(Girls15to20yo, OngoingRels.End, by="IDw", all.x=TRUE)
  if (nrow(Girls15to20yo) == 0)
  {
    Fraction.25ormore <- NA
  } else {
    # We need to know if no relationships exist for a girl
    Girls15to20yo.NumberOfAgeGaps <- as.data.table(table(Girls15to20yo.AgeGaps$IDw)) # data.table(ddply(Girls15to20yo.AgeGaps, .(IDw), summarise, number.AG  = length(AgeGap)))
    setnames(Girls15to20yo.NumberOfAgeGaps, c("IDw", "number.AG"))
    Girls15to20yo.NumberOfAgeGaps$IDw <- as.numeric(Girls15to20yo.NumberOfAgeGaps$IDw) 
    setkey(Girls15to20yo.NumberOfAgeGaps, IDw)
    Girls15to20yo.AgeGaps <- merge(Girls15to20yo.AgeGaps, Girls15to20yo.NumberOfAgeGaps, by="IDw", all.x=TRUE)
    Girls15to20yo.maxAgeGaps <- data.table(ddply(Girls15to20yo.AgeGaps, .(IDw), summarise, max.AgeGap = max(AgeGap, na.rm=FALSE)))
    Girls15to20yo.maxAgeGaps$max.AgeGap[is.infinite(Girls15to20yo.maxAgeGaps$max.AgeGap)]<-NA
    setkey(Girls15to20yo.maxAgeGaps, IDw)
    Girls15to20yo <- merge(Girls15to20yo, Girls15to20yo.maxAgeGaps, by="IDw", all.x=TRUE)
    
    Girls15to20yo$PartnerMaxAge <- (targetwindow.end - Girls15to20yo$TOB) + Girls15to20yo$max.AgeGap
    Girls15to20yo$Had25yoPartner <- Girls15to20yo$PartnerMaxAge>=25
    Fraction.25ormore <- ifelse(length(Girls15to20yo$Had25yoPartner) > 0, sum(Girls15to20yo$Had25yoPartner, na.rm=TRUE) /  length(Girls15to20yo$Had25yoPartner), 0)
    #   sumcheck <- sum(Girls15to20yo$Had25yoPartner, na.rm=TRUE)
    #   nrowcheck <- nrow(Girls15to20yo)
    #   divisioncheck <- sum(Girls15to20yo$Had25yoPartner, na.rm=TRUE) /  nrow(Girls15to20yo)
    #   print(c(sumcheck, nrowcheck, divisioncheck))
  }
  return(Fraction.25ormore)
}

## Process indicators for age-mixing analysis

# Reproductive Number over time (R)
allRtimestepsdata.list <- function(datalist){
  DT <- datalist$ptable
  DTE <- datalist$etable
  DTL <- datalist$ltable
  output <- data.table()
  timestep.width <- 1
  decplaces <- decimalplaces(timestep.width)
  timesteps <- DTL$Time
  for (time_i in head(timesteps, -1)){
    DTnewlyinfected <- DT[InfectTime > time_i & InfectTime <= (time_i + timestep.width) & TOD <= max(timesteps)] # People who got infected during the time slice and had a completed infectious period 
    transmissions <- DTE[eventname == "transmission" & eventtime > time_i] # transmissions that could possibly be caused by the people in DTnewlyinfected
    transmissions$ID <- transmissions$p1ID
    transmissions$Gender <- transmissions$p1gender
    
    DTnewlyinfectedkeycols <- c("ID", "Gender")
    transmissionskeycols <- c("ID", "Gender") # The ID and gender of the receptors, acquiring infection from people in DTnewlyinfected  
    setkeyv(DTnewlyinfected, DTnewlyinfectedkeycols)
    setkeyv(transmissions, transmissionskeycols)
    DTforRcalc <- merge(DTnewlyinfected, transmissions, all.x=TRUE)
    DTforRcalc$ID.Gender <- paste(DTforRcalc$ID, DTforRcalc$Gender, sep=".")
    DTwithR <- DTforRcalc[ ,infcaused :=  sum(!is.na(eventtime)), by = "ID.Gender"] # number of infections caused for each person in DTnewlyinfected 
    DTwithR$R <- sum(DTwithR$infcaused) / nrow(DTwithR)
    Rtimestepdata <- data.frame(time_i = time_i, suminfcaused = sum(DTwithR$infcaused), R = sum(DTwithR$infcaused) / nrow(DTwithR))
    output <- rbind(output, Rtimestepdata)
  }
  return(output)
}

allRtimestepsdata <- function(modeloutput){ #(DT = datalist$ptable, DTE = datalist$etable, cfg = configfile){
  datalist <- readthedata(modeloutput)
  DT <- datalist$ptable
  DTE <- datalist$etable
  cfg <- modeloutput[["input.params"]][["config"]]
  output <- data.table()
  timestep.width <- 1
  decplaces <- decimalplaces(timestep.width)
  timesteps <- seq(as.numeric(cfg["hivseed.time"]), as.numeric(cfg["population.simtime"]), by = timestep.width)
  for (time_i in head(timesteps, -1)){
    DTnewlyinfected <- DT[InfectTime > time_i & InfectTime <= (time_i + timestep.width) & TOD <= as.numeric(cfg["population.simtime"])] # People who got infected during the time slice and had a completed infectious period 
    # small hack: DT[ID == 151]$TOD <- as.numeric(cfg["population.simtime"]) - 1
    transmissions <- DTE[eventname == "transmission" & eventtime > time_i] # transmissions that could possibly be caused by the pople in DTnewlyinfected
    transmissions$ID <- transmissions$p1ID
    transmissions$Gender <- transmissions$p1gender

    DTnewlyinfectedkeycols <- c("ID", "Gender")
    transmissionskeycols <- c("ID", "Gender") # The ID and gender of the receptors, acquiring infection from people in DTnewlyinfected  
    setkeyv(DTnewlyinfected, DTnewlyinfectedkeycols)
    setkeyv(transmissions, transmissionskeycols)
    DTforRcalc <- merge(DTnewlyinfected, transmissions, all.x=TRUE)
    DTforRcalc$ID.Gender <- paste(DTforRcalc$ID, DTforRcalc$Gender, sep=".")
    DTwithR <- DTforRcalc[ ,infcaused :=  sum(!is.na(eventtime)), by = "ID.Gender"] # number of infections caused for each person in DTnewlyinfected 
    DTwithR$R <- sum(DTwithR$infcaused) / nrow(DTwithR)
    Rtimestepdata <- data.frame(time_i = time_i, suminfcaused = sum(DTwithR$infcaused), R = sum(DTwithR$infcaused) / nrow(DTwithR))
    output <- rbind(output, Rtimestepdata)
  }
  return(output)
}


# Create matrix of the average number of people that a man who was infected while in age group i infects in age group j.
# Do the same for female-to-male transmission.
# From there we can calculate the average number of men that get infected while in age groups <= i,
# with the "origin" infection coming from a man who got infected while in age group i.
# If this average is <1 for all age groups and genders, we have an unsustainable transmission pattern.
# We base this matrix on new infection data from the first year, so that
# (1) we are certain that they have completed their infectious period by the end of the simulation, and
# (2) we have an approximation for R0 instead of R. We will need to do multiple runs,
# introducing HIV in sigle 1-year age groups, for all age groups.
# Runs will need to last long (60 years or so) to ensure that the full infectious period is obtained for all.

recycler <- function(datalist = datalist, timepoint = as.numeric(cfgABC["hivseed.time"])){
  R.agegrouped <- R.agegrouped.list(datalist = datalist, timepoint = timepoint)
  DT <- datalist$ptable
  DThivseeds <- DT[InfectType == 0]
  # Find the number of male and female seeds for all 1-year age groups
  DThivseeds$age.floor <- floor(timepoint - DThivseeds$TOB)
  seeds.gender.agegrouped <- count(DThivseeds, vars = c("Gender", "age.floor"))
  seeds.gender.agegrouped.zero.supplement <- data.frame(Gender = rep(c(0,1), each=106^2), age.floor = rep(15:120, each=106, 2), freq = 0)
  seeds.gender.agegrouped <- merge.data.frame(seeds.gender.agegrouped, seeds.gender.agegrouped.zero.supplement, all.y = TRUE, by = c("Gender", "age.floor"))
  seeds.gender.agegrouped$freq.x[is.na(seeds.gender.agegrouped$freq.x)] <- 0
  seeds.gender.agegrouped$n <- seeds.gender.agegrouped$freq.x + seeds.gender.agegrouped$freq.y # Number of seeds of that age group and gender
  seeds.gender.agegrouped <- seeds.gender.agegrouped[ order(-seeds.gender.agegrouped$Gender), ]
  
  R.agegrouped$n.seeds <- seeds.gender.agegrouped$n
  R.agegrouped$R <- R.agegrouped$freq / R.agegrouped$n.seeds
  
  F2M.matrix <- matrix(head(R.agegrouped$R, length(R.agegrouped$R)/2), ncol=sqrt(length(R.agegrouped$R)/2), byrow=T)
  rownames(F2M.matrix) <- colnames(F2M.matrix) <- 15:120
  M2F.matrix <- matrix(tail(R.agegrouped$R, length(R.agegrouped$R)/2), ncol=sqrt(length(R.agegrouped$R)/2), byrow=T)
  rownames(M2F.matrix) <- colnames(M2F.matrix) <- 15:120
  
  F2M.matrix[is.na(F2M.matrix)] <- 0 # For local testing only
  M2F.matrix[is.na(M2F.matrix)] <- 0 # For local testing only
  
  F2F <- rep(NA, nrow(F2M.matrix))
  for(i in 1:nrow(F2M.matrix)){
    Female.le.i.infected <- 0
    for(j in 1:ncol(F2M.matrix)){
      Female.le.i.infected <- Female.le.i.infected + F2M.matrix[i,j] * sum(M2F.matrix[j, 1:i])
    }
    F2F[i] <- Female.le.i.infected
  }
  
  M2M <- rep(NA, nrow(F2M.matrix))
  for(i in 1:nrow(M2F.matrix)){
    Male.le.i.infected <- 0
    for(j in 1:ncol(M2F.matrix)){
      Male.le.i.infected <- Male.le.i.infected + M2F.matrix[i,j] * sum(F2M.matrix[j, 1:i])
    }
    M2M[i] <- Male.le.i.infected
  }
  # Which age groups can reclyce themselves?
  Male.recycle.success <- M2M >=1
  Male.agegroups.recycled <- colnames(F2M.matrix)[Male.recycle.success]
  Male.recycled.any <- sum(Male.recycle.success) > 0
  Female.recycle.success <- F2F >=1
  Female.agegroups.recycled <- colnames(F2M.matrix)[Female.recycle.success]
  Female.recycled.any <- sum(Female.recycle.success) > 0
  
  recycler.result <- list(data = R.agegrouped, M2M = M2M, F2F = F2F,
                          Male.agegroups.recycled = Male.agegroups.recycled,
                          Female.agegroups.recycled = Female.agegroups.recycled,
                          Male.recycled.any = Male.recycled.any,
                          Female.recycled.any = Female.recycled.any)
  return(recycler.result)
}
 
R.agegrouped.list <- function(datalist, timepoint = as.numeric(cfgABC["hivseed.time"])){
  DT <- datalist$ptable
  DTE <- datalist$etable
  DTnewlyinfected <- DT[InfectType == 0] # People who got infected during the time slice and had a completed infectious period
  if (sum(is.infinite(DTnewlyinfected$TOD)) > 0)
    warning("not all infectious periods are complete")
  transmissions <- DTE[eventname == "transmission" & eventtime > timepoint] # transmissions that could possibly be caused by the people in DTnewlyinfected
  transmissions$ID <- transmissions$p1ID
  transmissions$Gender <- transmissions$p1gender
  
  DTnewlyinfectedkeycols <- c("ID", "Gender")
  transmissionskeycols <- c("ID", "Gender") # The ID and gender of the receptors, acquiring infection from people in DTnewlyinfected  
  setkeyv(DTnewlyinfected, DTnewlyinfectedkeycols)
  setkeyv(transmissions, transmissionskeycols)
  DTforRcalc <- merge(DTnewlyinfected, transmissions, all.x=FALSE) # List of all transmission events caused by the seed infections
  DTforRcalc$p2age.floor <- floor(DTforRcalc$p2age) # Age bins of the people who got infected caused by seed infections
  DTforRcalc$p1age.floor <- floor(DTforRcalc$p1age - DTforRcalc$eventtime + timepoint) # Age bins of the seed infections (age at which they got infected)
  
  R.agegrouped <- count(DTforRcalc, vars = c("p2gender", "p1age.floor", "p2age.floor"))
  R.agegrouped.zero.supplement <- data.frame(p2gender = rep(c(0,1), each=106^2), p1age.floor = rep(15:120, each=106, 2), p2age.floor = rep(15:120, 2*106), freq = 0)
  R.agegrouped <- merge.data.frame(R.agegrouped, R.agegrouped.zero.supplement, all.y = TRUE, by = c("p2gender", "p1age.floor", "p2age.floor"))
  R.agegrouped$freq.x[is.na(R.agegrouped$freq.x)] <- 0
  R.agegrouped$freq <- R.agegrouped$freq.x + R.agegrouped$freq.y
  return(R.agegrouped)
}





transmissionmap <- function(modeloutput){ #(DT = datalist$ptable, DTE = datalist$etable, cfg = configfile){
  # requires the shape package
  # xlim <- c(0, 15)
  # ylim <- c(0, 108)
  # plot(0, type = "n", xlim = xlim, ylim = ylim)
  # First we insert diagonal "arrows" (age and time living with HIV)
  datalist <- readthedata(modeloutput)
  DT <- datalist$ptable
  DTE <- datalist$etable
  cfg <- modeloutput[["input.params"]][["config"]]
  Infecteds <- DT[InfectTime < Inf]
  timeatinf <- Infecteds$InfectTime
  ageatinf <- - Infecteds$TOB + Infecteds$InfectTime
  timeatdeath <- pmin(as.numeric(cfg["population.simtime"]), Infecteds$TOD)
  ageatdeath <- ageatinf + (timeatdeath - timeatinf)
  # Arrows(timeatinf, ageatinf, timeatdeath, ageatdeath, arr.length = 0.02, code = 2,
  #        arr.type = "T", arr.col = Infecteds$Gender + 2, col = Infecteds$Gender + 2)
  # Next we insert vertical dashed arrows (age and time at transmission)
  DTE <- datalist$etable
  DTincident <- DTE[eventname == "transmission"]
  timeattransm <- DTincident$eventtime
  ageoftransmitter <- DTincident$p1age
  genderoftransmitter <- DTincident$p1gender
  ageofreceptor <- DTincident$p2age + ((DTincident$p1age - DTincident$p2age) > 0) - ((DTincident$p1age - DTincident$p2age) < 0)
  # Arrows(timeattransm, ageoftransmitter, timeattransm, ageofreceptor, arr.length = 0.1, code = 2,
  #        arr.type = "triangle", col = DTincident$p1gender + 2)
  result <- list(livingwithHIVdata = data.frame(timeatinf = timeatinf,
                                                ageatinf = ageatinf,
                                                timeatdeath = timeatdeath,
                                                ageatdeath = ageatdeath),
                 transmissiondata = data.frame(timeattransm = timeattransm,
                                               ageoftransmitter = ageoftransmitter,
                                               genderoftransmitter = genderoftransmitter,
                                               ageofreceptor = ageofreceptor))
  return(result)
}


agemixing <- function(modeloutput){ #(DT = datalist$ptable, DTR = datalist$rtable){
  datalist <- readthedata(modeloutput)
  DT <- datalist$ptable
  DTR <- datalist$rtable
  #cfg <- modeloutput[["input.params"]][["config"]]
  DTRandDT_a <- merge(data.frame(DTR), data.frame(DT), by.x = "IDm", by.y = "ID")
  DTRandDT_b <- merge(DTRandDT_a, data.frame(DT), by.x = "IDw", by.y = "ID", suffixes = c(".m", ".w"))
  DTRandDT_b$relID <- factor(1:nrow(DTRandDT_b))
  setnames(DTRandDT_b, c("IDw", "IDm"), c("ID.w", "ID.m"))
  DTRlong <- reshape(data = DTRandDT_b,
                     idvar = "relID",
                     varying = names(DTRandDT_b)[c(1:2, 6:(ncol(DTRandDT_b)-1))], #, 6:ncol(DTRandDT_b))],
                     timevar = "g",
                     direction = "long")
  DTRlong$g.ID <- as.character(paste(DTRlong$g, DTRlong$ID, sep="."))
  
  DTRandDT_m <- merge(data.frame(DTR), data.frame(DT[Gender==0]), by.x = "IDm", by.y = "ID")
  DTRandDT_m$AgeMaleatForm <- - DTRandDT_m$TOB + DTRandDT_m$FormTime
  DTRandDT_m$AgeMaleatForm0 <- DTRandDT_m$AgeMaleatForm - min(DTRandDT_m$AgeMaleatForm)
  DTRandDT_m$AgeFemaleatForm <- DTRandDT_m$AgeMaleatForm - DTRandDT_m$AgeGap  
                
  # Analysis of age mixing pattern, using the DTR_andDT_m dataset, NOT the DTRlong dataset
  # Mean age gap, between-, and within-subject variance of age gaps
  
  meanagegap <- mean(DTRandDT_m$AgeGap)
  varagegap <- var(DTRandDT_m$AgeGap) #[DTRlong$g == "w"]) # We only want to count each relationship once
  # First we fit a linear mixed effects model
  MixingModel <- lme(AgeGap~1, random=~1|IDm, data = DTRandDT_m,
                     control=lmeControl(returnObject=TRUE, opt = "optim"))
  VarianceEstimates <- VarCorr(MixingModel)
  BetweenVar <- as.numeric(VarianceEstimates[1])
  WithinVar <- as.numeric(VarianceEstimates[2])
  # Then we fit a non-linear mixed effects model: within-subject variance does not increase linearly with age
  #Subset datasets for men as "respondents" reporting on all their relationships
  DTRmen <- DTRlong[DTRlong$g=="m",]
  DTRmen$AgeMaleatForm <- - DTRmen$TOB + DTRmen$FormTime
  DTRmen$AgeFemaleatForm <- DTRmen$AgeMaleatForm - DTRmen$AgeGap
  DTRmen$AgeMaleatForm0 <- DTRmen$AgeMaleatForm - min(DTRmen$AgeMaleatForm)
  
  men.lme <- lme(AgeFemaleatForm ~ AgeMaleatForm0,
                 data = DTRandDT_m,
                 control=lmeControl(returnObject=TRUE),
                 random = ~1 | IDm, 
                 method = "REML",
                 weight = varPower(value = 0.5, form = ~AgeMaleatForm0 + 1))
  #The formula to calculate the weights for variance is |v|^(2*t)
  #Age can't be at 0 in varPower formula because then it will evaluate to 0 variance for the first level (18 year olds)
  
  
  slope <- summary(men.lme)$tTable[2, 1]
  
  #add the predicted values to the dataset
  #Level 0 mean population level predictions
  DTRandDT_m$pred <- predict(men.lme, DTRandDT_m, level = 0)

  #Need to calculate a vector of within-subject variances because the variances should be different for each age
  #Create a vector of the variance weights. This is based upon the formula for varPower |v|^(2*t)
  #(2*t) is the power
  #|v| is the absolute value of each of the ages from the age vector
  powerm <- attributes(men.lme$apVar)$Pars["varStruct.power"]
  within.var.weights <- (1 + DTRandDT_m$AgeMaleatForm0)^powerm
  #Extracting the residual variance (within-subject variance) for youngest men (0— baseline)
  within.var.base <- men.lme$sigma^2

  #Now create a vector that contains the residual variances for all ages
  # DTRandDT_m$within.var <- within.var.base * within.var.weights

  #Extract the variance of the random intercept # This is the between-subject variance
  # DTRandDT_m$between.var <- getVarCov(men.lme)[1,1] # called "interceptvar" in LNS analysis

  #Caluclate the prediction stardard errors
  # DTRandDT_m$predict.sd <- sqrt(DTRandDT_m$between.var + DTRandDT_m$within.var) # called "se2" in LNS analysis
  
  result <- list(AAD = mean(DTRandDT_m$AgeGap),
                 VAD = var(DTRandDT_m$AgeGap),
                 men.lme = men.lme,
                 slope = slope,
                 powerm = attributes(men.lme$apVar)$Pars["varStruct.power"],
                 within.var.base = men.lme$sigma^2,
                 within.var.weights = within.var.weights,
                 WVAD = within.var.base * within.var.weights,
                 Pred = predict(men.lme, DTRandDT_m, level = 0),
                 BVAD = getVarCov(men.lme)[1,1],
                 agescatterdata = data.frame(DTRandDT_m))
  return(result)
}

agemixing.list <- function(datalist){ #(DT = datalist$ptable, DTR = datalist$rtable){
  DT <- datalist$ptable
  DTR <- datalist$rtable
  #cfg <- modeloutput[["input.params"]][["config"]]
  DTRandDT_a <- merge(data.frame(DTR), data.frame(DT), by.x = "IDm", by.y = "ID")
  DTRandDT_b <- merge(DTRandDT_a, data.frame(DT), by.x = "IDw", by.y = "ID", suffixes = c(".m", ".w"))
  DTRandDT_b$relID <- factor(1:nrow(DTRandDT_b))
  setnames(DTRandDT_b, c("IDw", "IDm"), c("ID.w", "ID.m"))
  DTRlong <- reshape(data = DTRandDT_b,
                     idvar = "relID",
                     varying = names(DTRandDT_b)[c(1:2, 6:(ncol(DTRandDT_b)-1))], #, 6:ncol(DTRandDT_b))],
                     timevar = "g",
                     direction = "long")
  DTRlong$g.ID <- as.character(paste(DTRlong$g, DTRlong$ID, sep="."))
  
  DTRandDT_m <- merge(data.frame(DTR), data.frame(DT[Gender==0]), by.x = "IDm", by.y = "ID")
  DTRandDT_m$AgeMaleatForm <- - DTRandDT_m$TOB + DTRandDT_m$FormTime
  DTRandDT_m$AgeMaleatForm0 <- DTRandDT_m$AgeMaleatForm - 18 #DTRandDT_m$AgeMaleatForm - min(DTRandDT_m$AgeMaleatForm)
  DTRandDT_m$AgeFemaleatForm <- DTRandDT_m$AgeMaleatForm - DTRandDT_m$AgeGap  
  
  DTRandDT_m <- subset(DTRandDT_m,
                 FormTime >=targetwindow.start &
                   FormTime < targetwindow.end &
                   AgeMaleatForm >= 18 & # summary(men$age) # in men aged 18 to 49
                   AgeMaleatForm < 50)
  
  # Analysis of age mixing pattern, using the DTR_andDT_m dataset, NOT the DTRlong dataset
  # Mean age gap, between-, and within-subject variance of age gaps
  
#   meanagegap <- mean(DTRandDT_m$AgeGap)
#   varagegap <- var(DTRandDT_m$AgeGap) #[DTRlong$g == "w"]) # We only want to count each relationship once
#   # First we fit a linear mixed effects model
#   MixingModel <- lme(AgeGap~1, random=~1|IDm, data = DTRandDT_m,
#                      control=lmeControl(returnObject=TRUE, opt = "optim"))
#   VarianceEstimates <- VarCorr(MixingModel)
#   BetweenVar <- as.numeric(VarianceEstimates[1])
#   WithinVar <- as.numeric(VarianceEstimates[2])
#   # Then we fit a non-linear mixed effects model: within-subject variance does not increase linearly with age
#   #Subset datasets for men as "respondents" reporting on all their relationships
#   DTRmen <- DTRlong[DTRlong$g=="m",]
#   DTRmen$AgeMaleatForm <- - DTRmen$TOB + DTRmen$FormTime
#   DTRmen$AgeFemaleatForm <- DTRmen$AgeMaleatForm - DTRmen$AgeGap
#   DTRmen$AgeMaleatForm0 <- DTRmen$AgeMaleatForm - min(DTRmen$AgeMaleatForm)
  
  men.lme <- lme(AgeFemaleatForm ~ AgeMaleatForm0,
                 data = DTRandDT_m,
                 control=lmeControl(maxIter=200, returnObject=TRUE, opt = "optim"),
                 random = ~1 | IDm, 
                 method = "REML",
                 weight = varPower(value = 0.5, form = ~AgeMaleatForm0 + 1))
  #The formula to calculate the weights for variance is |v|^(2*t)
  #Age can't be at 0 in varPower formula because then it will evaluate to 0 variance for the first level (18 year olds)
  
  #add the predicted values to the dataset
  #Level 0 mean population level predictions
  DTRandDT_m$pred <- predict(men.lme, DTRandDT_m, level = 0)
  
  #Need to calculate a vector of within-subject variances because the variances should be different for each age
  #Create a vector of the variance weights. This is based upon the formula for varPower |v|^(2*t)
  #(2*t) is the power
  #|v| is the absolute value of each of the ages from the age vector
  powerm <- attributes(men.lme$apVar)$Pars["varStruct.power"]
  within.var.weights <- (1 + DTRandDT_m$AgeMaleatForm0)^powerm
  #Extracting the residual variance (within-subject variance) for youngest men (0— baseline)
  within.var.base <- men.lme$sigma^2
  
  #Now create a vector that contains the residual variances for all ages
  # DTRandDT_m$within.var <- within.var.base * within.var.weights
  
  #Extract the variance of the random intercept # This is the between-subject variance
  # DTRandDT_m$between.var <- getVarCov(men.lme)[1,1] # called "interceptvar" in LNS analysis
  
  #Caluclate the prediction stardard errors
  # DTRandDT_m$predict.sd <- sqrt(DTRandDT_m$between.var + DTRandDT_m$within.var) # called "se2" in LNS analysis
  
  
  
  
  result <- list(AAD = mean(DTRandDT_m$AgeGap),
                 VAD = sd(DTRandDT_m$AgeGap),
                 men.lme = men.lme,
                 powerm = as.numeric(attributes(men.lme$apVar)$Pars["varStruct.power"]),
                 slope = summary(men.lme)$tTable[2, 1],
                 within.var.base = men.lme$sigma^2,
                 within.var.weights = within.var.weights,
                 WVAD = within.var.base * within.var.weights,
                 Pred = predict(men.lme, DTRandDT_m, level = 0),
                 BVAD = getVarCov(men.lme)[1,1],
                 agescatterdata = data.frame(DTRandDT_m))
  return(result)
}

agemixing.df.maker <- function(datalist){ #(DT = datalist$ptable, DTR = datalist$rtable){
  DT <- datalist$ptable
  DTR <- datalist$rtable
  #cfg <- modeloutput[["input.params"]][["config"]]
  DTRandDT_a <- merge(data.frame(DTR), data.frame(DT), by.x = "IDm", by.y = "ID")
  DTRandDT_b <- merge(DTRandDT_a, data.frame(DT), by.x = "IDw", by.y = "ID", suffixes = c(".m", ".w"))
  if (nrow(DTRandDT_b) > 0) {
    DTRandDT_b$relID <- factor(1:nrow(DTRandDT_b))
    setnames(DTRandDT_b, c("IDw", "IDm"), c("ID.w", "ID.m"))
    DTRlong <- reshape(data = DTRandDT_b,
                       idvar = "relID",
                       varying = names(DTRandDT_b)[c(1:2, 6:(ncol(DTRandDT_b)-1))], #, 6:ncol(DTRandDT_b))],
                       timevar = "g",
                       direction = "long")
    DTRlong$g.ID <- as.character(paste(DTRlong$g, DTRlong$ID, sep="."))
    
    DTRandDT_m <- merge(data.frame(DTR), data.frame(DT[Gender==0]), by.x = "IDm", by.y = "ID")
    DTRandDT_m$AgeMaleatForm <- - DTRandDT_m$TOB + DTRandDT_m$FormTime
    DTRandDT_m$AgeMaleatForm0 <- DTRandDT_m$AgeMaleatForm - 18 #DTRandDT_m$AgeMaleatForm - min(DTRandDT_m$AgeMaleatForm)
    DTRandDT_m$AgeFemaleatForm <- DTRandDT_m$AgeMaleatForm - DTRandDT_m$AgeGap  
    
    DTRandDT_m <- subset(DTRandDT_m,
                         FormTime >=targetwindow.start &
                           FormTime < targetwindow.end &
                           AgeMaleatForm >= 18 & # summary(men$age) # in men aged 18 to 49
                           AgeMaleatForm < 50)
    
    agemixing.men.df <- data.frame(DTRandDT_m)
  } else {
    agemixing.men.df <- data.frame(AgeGap = 0)
#     IDm = integer(0),
#                                    IDw = integer(0),
#                                    FormTime = numeric(0),
#                                    DisTime = numeric (0),
#                                    AgeGap = numeric(0),
#                                    Gender = integer(0),
#                                    TOB = numeric(0),
#                                    TOD = numeric(0),
#                                    IDF = integer(0),
#                                    IDM = integer(0),
#                                    TODebut = numeric(0),
#                                    FormEag = numeric(0),         "InfectTime"      "InfectOrigID"   
#                                    [15] "InfectType"      "log10SPVL"       "TreatTime"       "XCoord"          "YCoord"          "AIDSDeath"       "AgeMaleatForm"  
#                                    [22] "AgeMaleatForm0"  "AgeFemaleatForm")
  }
  return(agemixing.men.df)
}

agemixing.lme.fitter <- function(am.df){
  men.lme <- lme(AgeFemaleatForm ~ AgeMaleatForm0,
                 data = am.df,
                 control=lmeControl(maxIter=10000, returnObject=TRUE, opt = "optim"),
                 random = ~1 | IDm, 
                 method = "REML",
                 weight = varPower(value = 0.5, form = ~AgeMaleatForm0 + 1))
}




# allRtimestepsdata <- function(DT = datalist$ptable, DTE = datalist$etable, cfg = configfile){
#   DTE$ID <- DTE$
#   output <- data.table()
#   
#     transmissions <- DTE[eventname == "transmission" & eventtime > time_i] # transmissions that could possibly be caused by the pople in DTnewlyinfected
#     transmissions$ID <- transmissions$p1ID
#     transmissions$Gender <- transmissions$p1gender
#     
#     DTnewlyinfectedkeycols <- c("ID", "Gender")
#     transmissionskeycols <- c("ID", "Gender") # The ID and gender of the receptors, acquiring infection from people in DTnewlyinfected  
#     setkeyv(DTnewlyinfected, DTnewlyinfectedkeycols)
#     setkeyv(transmissions, transmissionskeycols)
#     DTforRcalc <- merge(DTnewlyinfected, transmissions, all.x=TRUE)
#     DTforRcalc$ID.Gender <- paste(DTforRcalc$ID, DTforRcalc$Gender, sep=".")
#     DTwithR <- DTforRcalc[ ,infcaused :=  sum(!is.na(eventtime)), by = "ID.Gender"] # number of infections caused for each person in DTnewlyinfected 
#     DTwithR$R <- sum(DTwithR$infcaused) / nrow(DTwithR)
#     Rtimestepdata <- data.frame(time_i = time_i, suminfcaused = sum(DTwithR$infcaused), R = sum(DTwithR$infcaused) / nrow(DTwithR))
#     output <- rbind(output, Rtimestepdata)
#   }
#   return(output)
# }

################
# Create cumulative network dataset
################
cumulnetwork <- function(DTR = datalist$rtable, cfg = configfile, windowwidth = 20){
  endwindow <- 30 # as.numeric(cfg["hivseed.time"]) + windowwidth
  startwindow <- 0 # endwindow - windowwidth
  Rels <- DTR[FormTime >=startwindow & FormTime < endwindow]
  rdata <- data.frame(
    tails = paste0(Rels$IDm, "m"),
    heads = paste0(Rels$IDw, "w"),
    form = Rels$FormTime,
    diss = Rels$DisTime,
    durat = pmin(endwindow, Rels$DisTime) - Rels$FormTime,
    stringsAsFactors=FALSE
  )
  cn <- network(rdata[ , 1:2], directed = FALSE, matrix.type = "edgelist") # cumulative network
  cnvertexnames <- network.vertex.names(cn)
  maleindices <- grep("m", cnvertexnames)
  Sex <- rep(0, network.size(cn))
  Sex[maleindices] <- 1
  cnseedindex <- grep(SeedID, cnvertexnames)
  Sex[cnseedindex] <- 2
  VertexSize <- rep(1, network.size(cn))
  VertexSize[cnseedindex] <- 2
  
  set.vertex.attribute(cn, attrname="Sex", value=Sex, v=seq_len(network.size(cn)))
  #cnplot<- plot(cn, label = "", vertex.col = Sex)
  # Saving the plot
  png(file="CN.png", height = 1800, width = 1800, res=300)
  cnplot <- plot(cn, label = "", vertex.col = Sex, vertex.cex = VertexSize)
  dev.off()
  
  
}



################
# Create cross-sectional network dataset
################
csnetwork <- function(DTR = datalist$rtable, cfg = configfile, time = as.numeric(cfg["hivseed.time"])){
  Rels <- DTR[FormTime <= time & DisTime > time]
  rdata <- data.frame(
    tails = paste0(Rels$IDm, "m"),
    heads = paste0(Rels$IDw, "w"),
    form = Rels$FormTime,
    diss = Rels$DisTime,
    durat = pmin(endwindow, Rels$DisTime) - Rels$FormTime,
    stringsAsFactors=FALSE
  )
  csn <- network(rdata[ , 1:2], directed = FALSE, matrix.type = "edgelist") # cumulative network
  csnvertexnames <- network.vertex.names(csn)
  maleindices <- grep("m", csnvertexnames)
  SexCSN <- rep(0, network.size(csn))
  SexCSN[maleindices] <- 1
  csnseedindex <- grep(SeedID, csnvertexnames)
  SexCSN[csnseedindex] <- 2
  VertexSize <- rep(1, network.size(csn))
  VertexSize[csnseedindex] <- 2
  set.vertex.attribute(csn, attrname="Sex", value=SexCSN, v=seq_len(network.size(csn)))
  #csnplot <- plot(csn, label = "", vertex.col = SexCSN)
  # Saving the plot
  png(file="CSN.png", height = 1800, width = 1800, res=300)
  plot(csn, label = "", vertex.col = SexCSN, vertex.cex = VertexSize)
  dev.off()
  
}

# Overlaying cn with csn
# matchingindices <- pmatch(csnvertexnames, cnvertexnames)
# csncoord <- cnplot[matchingindices, ]
# png(file="CSN_CS.png", height = 1800, width = 1800, res=300)
# plot(csn, label = "", vertex.col = SexCSN, coord = csncoord, vertex.cex = VertexSize)
# dev.off()


################
# Create potential transmission network (directed network) (=PTN)
################
# The potential transmission routes, starting from the seed infection
ptnetwork <- function(DT = datalist$ptable, DTR = datalist$rtable, cfg = configfile){
  
  DTR[, IDm := paste0(DTR[, IDm], "m")]
  DTR[, IDw := paste0(DTR[, IDw], "w")]
  Seed <- datalist$ptable[InfectTime == as.numeric(cfg["hivseed.time"])]
  gender_suffix <- c("m", "w")
  SeedID <- paste0(Seed$ID, gender_suffix[Seed$Gender + 1])
  
  PTN <- data.table()
  # Now we run through the partners of the seed and their partners and so on, while NEW relationships are being formed
  DTRfr <- data.frame(DTR)
  # Rels of seed
  SeedRels <- DTRfr[, (Seed$Gender + 1)] == SeedID
  AfterSeedTimeRels <- DTR$DisTime > as.numeric(cfg["hivseed.time"])
  RelsList <- DTR[SeedRels & AfterSeedTimeRels, ] # This is the list of ongoing Rels of the seed, after HIV was introduced
  
  iteration <- 0
  while (nrow(RelsList) > 0){
    iteration <- iteration + 1
    RelsListfr <- data.frame(RelsList)
    RelsListCopy <- data.table(RelsListfr)
    # Now we list the partners of the seed
    if (iteration == 1){
      PartnerColumnIndex <- abs(Seed$Gender - 1) + 1
      
      if (Seed$Gender == 1){
        setcolorder(RelsListCopy, c("IDw", "IDm", "FormTime", "DisTime", "AgeGap"))
      }
    } else {
      PartnerColumnIndex <- - PartnerColumnIndex + 3 # 1 becomes 2 and 2 becomes 1
      if (PartnerColumnIndex == 1){
        setcolorder(RelsListCopy, c("IDw", "IDm", "FormTime", "DisTime", "AgeGap"))
      }
    }
    RelsListCopyfr <- data.frame(RelsListCopy)
    if (iteration == 1){
    PTN <- as.matrix(RelsListCopyfr)#[, 1:2])
    } else {
      PTN <- rbind(as.matrix(PTN), as.matrix(RelsListCopyfr))#[, 1:2]))
    }
    
    PartnerIDs <- unique(RelsListfr[, PartnerColumnIndex])
    FirstRelIndex <- match(PartnerIDs, RelsListfr[, PartnerColumnIndex])
    StartPotentialTransm <- pmax(as.numeric(cfg["hivseed.time"]), RelsListfr[FirstRelIndex, "FormTime"])
    
    # Now list the relationships that these partners formed AFTER they formed the relationship with the seed that could have led to transmission
    RelsList <- data.table()
    for (i in (1:length(PartnerIDs))){
      PartnerRels <- DTRfr[, PartnerColumnIndex] == PartnerIDs[i]
      AfterRels <- DTRfr[, "FormTime"] > StartPotentialTransm[i]
      RelsListi <- DTR[PartnerRels & AfterRels, ]
      RelsList <- rbind(RelsList, RelsListi)
    }
  }
  rPTNdata <- data.frame(
    tails = PTN[, 1],
    heads = PTN[, 2],
    form = as.numeric(PTN[, 3]),
    diss = as.numeric(PTN[, 4]),
    durat = pmin(as.numeric(cfg["population.simtime"]), as.numeric(PTN[, 4])) - as.numeric(PTN[, 3]),
    stringsAsFactors=FALSE
  )
  ptn <- network(rPTNdata[ , 1:2], directed = TRUE, matrix.type = "edgelist") # potential tranmission network
  ptnvertexnames <- network.vertex.names(ptn)
  ptnmaleindices <- grep("m", ptnvertexnames)
  ptnseedindex <- grep(SeedID, ptnvertexnames)
  SexPTN <- rep(0, network.size(ptn))
  SexPTN[ptnmaleindices] <- 1
  SexPTN[ptnseedindex] <- 2
  set.vertex.attribute(ptn, attrname="Sex", value=list(SexPTN), v=seq_len(network.size(ptn)))
  VertexSize <- rep(1, network.size(ptn))
  VertexSize[ptnseedindex] <- 2
  # Plotting the network
  ptnplot <- plot(ptn,
                  label = "",
                  vertex.col = SexPTN,
                  vertex.cex = VertexSize,
                  jitter=T)
  # Saving the plot
  png(file="PTN.png", height = 1800, width = 1800, res=300)
  plot(ptn,
       label = "",
       vertex.col = SexPTN,
       vertex.cex = VertexSize,
       jitter=T)

  dev.off()
  
}


# Overlaying cn with ptn
# matchingindicesPTN <- pmatch(ptnvertexnames, cnvertexnames)
# ptncoord <- cnplot[matchingindicesPTN, ]
# png(file="PTN_CN.png", height = 1800, width = 1800, res=300)
# plot(ptn, label = "", vertex.col = SexPTN, coord = ptncoord, vertex.cex = VertexSize)
# dev.off()







################
# Create transmission tree dataset
################
transmissiontree <- function(datalist = datalist, cfg = configfile){
  TE = datalist$etable[eventname=="transmission"] # TE = transmission events
  
  Seed <- datalist$ptable[InfectTime == as.numeric(cfg["hivseed.time"])]
  gender_suffix <- c("m", "w")
  SeedID <- paste0(Seed$ID, gender_suffix[Seed$Gender + 1])
  
  edata <- data.frame(
    tails = paste0(TE$p1ID, gender_suffix[TE$p1gender + 1]), #TE$p1name,
    heads = paste0(TE$p2ID, gender_suffix[TE$p2gender + 1]),
    time = TE$eventtime,
    tailname = TE$p1name,
    tailgender = TE$p1gender,
    tailage = TE$p1age,
    headname = TE$p2name,
    headgender = TE$p2gender,
    headage = TE$p2age,
    stringsAsFactors=FALSE
  )
  
# Before we turn this edge list into a network object, we need to add the seeding transmission event
  edataplusroot <- rbind(c("root", edata$tails[1], as.numeric(cfg["hivseed.time"]),
                   "root", NA, NA, edata$tailname[1], edata$tailgender[1],
                   edata$tailage[1] - edata$time[1] + as.numeric(cfg["hivseed.time"])),
                 edata)
#   
#   ntips <- length(unique(c(TE$p1name, TE$p2name)))
#   tip.labels <- unique(c(TE$p1name, TE$p2name))
#                      
#   ninternalnodes <- nrow(edataplusroot)

tn <- network(edata[ , 1:2], directed = TRUE, matrix.type = "edgelist")
tnvertexnames <- network.vertex.names(tn)
tnmaleindices <- grep("m", tnvertexnames)
tnseedindex <- grep(SeedID, tnvertexnames) # Temporary solution for SACEMA Research days
SexTN <- rep(0, network.size(tn))
SexTN[tnmaleindices] <- 1
SexTN[tnseedindex] <- 2
set.vertex.attribute(tn, attrname="Sex", value=list(SexTN), v=seq_len(network.size(tn)))
VertexSizeTN <- rep(1, network.size(tn))
VertexSizeTN[tnseedindex] <- 2

tn.undir <- network(edata, directed = FALSE, matrix.type = "edgelist", ignore.eval = FALSE)
tndata <- list(tn = tn,
               tn.vertex.col = SexTN,
               tn.vertex.cex = VertexSizeTN,
               tn.undir = tn.undir)

return(tndata)
# Plotting the network
# tnplot <- plot(tn,
#                 label = "",
#                 vertex.col = SexTN,
#                 vertex.cex = VertexSizeTN,
#                 jitter=T)
# # Saving the plot
# png(file="TN.png", height = 1800, width = 1800, res=300)
# plot(tn,
#      label = "",
#      vertex.col = SexTN,
#      vertex.cex = VertexSizeTN,
#      jitter=T)
# 
# dev.off()

}

postABCplot <- function(ABCfit = ABC_LenormandResult,
                        targetSS = sum_stat_obs,
                        preprior = preprior,
                        param.prior = simpact_prior){
  
  #outputvector <- c(Exp.Growth.rate, Log.PTR, Exp.hivprevs, Exp.aAG, Exp.sdAG, Exp.slope, median.age.men, median.age.women, Exp.powerm, Exp.wsdad, Exp.bsdad, Exp.concurr.pointprev)
  
  
  stats <- ABCfit$stats[, c(1,2,18:ncol(ABCfit$stats))]
  param <- ABCfit$param
  prior.df <- ldply(simpact_prior)
  prior.df$V4 <- preprior
  
  statsdiff <- sweep(stats, 2, targetSS)
  statsreldiff <- sweep(statsdiff, 2, targetSS, FUN = "/")
  av.normalised.distance <- rowSums(statsreldiff)
  av.abs.normalised.distance <- rowSums(abs(statsreldiff))
  
  #bestruns <- head(order(av.abs.normalised.distance),10)
  bestruns <- order(av.abs.normalised.distance)
  
  stats.df <- as.data.table(stats[bestruns,])
  stats.df$run <- 1: nrow(stats.df)
  stats.df.long <- melt(data = stats.df, measure.vars = 1:((ncol(stats.df))-1))
  stats.df.long$target <- targetSS
  stats.df.long <- na.omit(stats.df.long)
  
  param.df <- as.data.table(param[bestruns, ])
  param.df$run <- 1:nrow(param.df)
  param.df.long <- melt(data = param.df, measure.vars = 1:((ncol(param.df))-1))
  param.df.long$prior.low <- prior.df$V2
  param.df.long$prior.high <- prior.df$V3
  param.df.long$prior.est <- prior.df$V4
  
  plot.SS <- ggplot(data = stats.df.long, aes(x = value)) +
    geom_histogram(stat = "bin",
                   binwidth = 0.1) +
    facet_wrap(~ variable, nrow = 2, scales = "free") +
    ggtitle("Summary statistics")
  vline.data <- data.frame(target = targetSS, variable = levels(stats.df.long$variable))
  plot.SS <- plot.SS +
    geom_vline(aes(xintercept = target), colour = "green4", size = 2, vline.data)
  plot.SS
  
  plot.param <- ggplot(data = param.df.long, aes(x = value)) +
    geom_histogram() +
    facet_wrap(~ variable, nrow = 2, scales = "free") +
    ggtitle("Model parameters")
  vline.data <- data.frame(prior.low = as.numeric(prior.df$V2),
                           prior.high = as.numeric(prior.df$V3),
                           prior.est = prior.df$V4,
                           variable = levels(param.df.long$variable))
  plot.param <- plot.param +
    geom_vline(aes(xintercept = prior.low), colour = "blue2", vline.data) +
    geom_vline(aes(xintercept = prior.high), colour = "blue2", vline.data) +
    geom_vline(aes(xintercept = prior.est), colour = "blue2", size = 2, vline.data)
  plot.param
  multiplot(plot.SS, plot.param)
  
}

# Overlaying cn with tn
# matchingindicesTN <- pmatch(tnvertexnames, cnvertexnames)
# tncoord <- cnplot[matchingindicesTN, ]
# png(file="TN_CN.png", height = 1800, width = 1800, res=300)
# plot(tn, label = "", vertex.col = SexTN, coord = tncoord, vertex.cex = VertexSizeTN)
# dev.off()

######
# Last step: plot phylogeny of transmission tree
######

# pg <- phylo4(nj(geodist(tn)$gdist))
# apepg <- as(pg, "phylo")
# apepgrooted <- multi2di(apepg, random = FALSE)
# apepgrooted$edge.length <- c(30-10.18568, # 7-6
#                              15.05207-10.18568, # 7-8
#                              30-15.05207, # 8-5
#                              16.13587-15.05207, # 8-9
#                              22.99193-16.13587, # 9-10
#                              24.29339-22.99193, # 10-11
#                              30-24.29339, # 11-2
#                              30-24.29339, # 11-1
#                              30-22.99193, # 10-3
#                              30-16.13587  #  9-4
# )
# 
# apepgrooted$root.edge <- 10.18568-10
# apepgrooted$tip.label <- c("Woman 46",
#                            "Woman 37",
#                            "Man 23",
#                            "Man 4",
#                            "Man 2",
#                            "Woman 28")
# png(file="PG.png", height = 1800, width = 1800, res=300)
# plot(apepgrooted, root.edge = TRUE)
# dev.off()

# 
# 
#   g <- network(edata[ , 1:2], directed = TRUE, matrix.type = "edgelist")
# vertexnames <- network.vertex.names(g)
# plot(network(edata, directed=T, matrix.type = "edgelist"), label = vertexnames)
# 
# plot(g, label = vertexnames)
# 
#   g2 <- network(edataplusroot[ , 1:2], directed = TRUE, matrix.type = "edgelist")
# 
# 
# test <- phylo4(nj(geodist(g)$gdist))
# apetest <- as(test, "phylo")
# apetestrooted <- multi2di(apetest, random = FALSE)
# plot(apetestrooted)
# 
# phylo4(nj(geodist(g2)$gdist))


  # Now we add node attributes
  # Now we add edge attributes

  # Now we extract essential information from edata to create phylo tree.
  # Constraints to input for phylo4 tree constructor:
  
  #   if the tree has edge lengths defined, the number of edge lengths must match
  #   the number of edges; the number of tip labels must match the number of tips;
  #   in a tree with ntips tips and nnodes (total) nodes, nodes 1 to ntips must be
  #   tips if the tree is rooted, the root must be node number ntips+1 and the
  #   root node must be the first row of the edge matrix; tip labels, node labels,
  #   edge labels, edge lengths must have proper internal names (i.e. internal
  #   names that match the node numbers they document); tip and node labels must be
  #   unique.

  # First essential element is x, the matrix of edges.
#   nodes <- 1:(ntips+ninternalnodes)
#   rootnodenumber <- ntips+1
#   internalnodenumbers <- (rootnodenumber+1):max(nodes)
#   tipnumbers <- 1:ntips
#   
#   firstrow <- c(rootnodenumber, 1)
  


################
# Create phylogeny dataset
################


#   nrels[i] <- nrow(r)
#   cumulincid[i] <- sum(p$InfectTime>=60 & p$InfectTime<Inf) #cfg["hivseed.time"] & p$InfectTime<Inf)
#   # How many infections to place in the last 10 years of the simulation? (simulation years 50-60)
#   
#   allrpdata <- rbind(allrpdata, rp)
# } 
# 
# summarystats <- (data.frame(simID, meanagegap, varagegap, BetweenVar, WithinVar,
#                             nrels, cumulincid))
# summarystats
# }
  
  



# hivincdata <- data.frame(timesteps,incidence)
# ggplot(hivincdata,aes(x=timesteps, y=1000*incidence)) +
#   geom_line(colour = "darkred", size=2) +
#   guides(colour = FALSE) +
#   xlab("Simulation time") +
#   ylab("HIV incidence (per 100 PY)")



incidentcase <- function(DT, interval){ # arguments are the eventlog data.table and a vector of start and end times of the interval
  DTincident <- DT[eventname == "transmission" & eventtime > interval[1] & eventtime <= interval[2]]
}

nonHIVdeath <- function(DT, interval){ # arguments are the eventlog data.table and a vector of start and end times of the interval
  DTnonHIVdeath <- DT[eventname == "normalmortality" & eventtime > interval[1] & eventtime <= interval[2]]
}


alive.infected <- function(DT, time, dec.places){ # arguments are the personlog data.table and a point in time
  DTalive <- DT[TOB <= time & TOD > time]
  DTalive$Age <- time -round(ceiling(DTalive$TOB), dec.places) # Next we allocate them in discrete age bins with bin size as wide as timestep
  DTalive$Infected <- time >= DTalive$InfectTime # Now we allocate infection status to all people in our table of living people
  return(DTalive)
}

gg_color_hue <- function(n, l) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=l, c=100)[1:n]
}


decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

# ####
# # To create heatplot, we need to loop through time steps, and calculate age- and gender-specific prevalence and incidence rates.
# alltimestepsdata <- data.table()
# timesteps <- seq(as.numeric(cfg["hivseed.time"]), as.numeric(cfg["population.simtime"]), by=1)
# for (time_i in timesteps) {
#   alivetable <- alivetab(time_i)   # First we only take the data of people who were alive at time_i
#   alivetable$Age.1 <- time_i -round(alivetable$TOB, 0) # Next we allocate them in discrete age bins with bin size equal as width of timesteps
#   alivetable$Infected <- time_i >= alivetable$InfectTime # Now we allocate infection status to all people in our table of living people
#   alive_data.table <- data.table(alivetable) # We turn the dataset into a proper data.table
#   timestepdata <- alive_data.table[,sum(Infected) / nrow(alive_data.table),by="Gender,Age.1"] # And we calculate HIV prevalence, by gender and age group
#   setnames(timestepdata, "V1", "Prevalence")
#   timestepdata <- cbind(time_i, timestepdata)
#   alltimestepsdata <- rbind(alltimestepsdata, timestepdata)
# }
# alltimestepsdata <- cbind(i, alltimestepsdata) # i is the ID number of the simulation
# 
# 
# 
# # Test sqldf query
# # x is evaltime (e.g. at end of simulation)
# relsandconcurM <- function(x){
#   query <- paste0("SELECT IDm, rels,
#                   CASE WHEN rels > 1 THEN 1 ELSE 0 END concur
#                   FROM (
#                   SELECT IDm, COUNT(*) rels
#                   FROM r
#                   WHERE ", x, "> FormTime AND ", x, "< DisTime
#                   GROUP BY IDm
#                   )"
#                   )
#   return (sqldf(query))
# }
# relsandconcurW <- function(x){
#   query <- paste0("SELECT IDw, rels,
#                   CASE WHEN rels > 1 THEN 1 ELSE 0 END concur
#                   FROM (
#                   SELECT IDw, COUNT(*) rels
#                   FROM r
#                   WHERE ", x, "> FormTime AND ", x, "< DisTime
#                   GROUP BY IDw
#                   )"
#                   )
#   return (sqldf(query))
# }
# # selecting people who are alive at time x
# alive <- function(x){
#   query <- paste0("SELECT ID, TOB
#                   FROM p
#                   WHERE ", x, ">= TOB AND ", x, "< TOD
#                   "#                  GROUP BY ID"
#   )
#   return (sqldf(query))
# }
# # number of infections caused
# infcaused <- function(){
#   query <- paste0("SELECT InfectOrigID, COUNT(*) infcaused
#                   FROM p
#                   GROUP BY InfectOrigID"
#   )
#   return (sqldf(query))
# }
# # infected at evaltime x
# infected <- function(x){
#   query <- paste0("SELECT ID, log10SPVL
#                   FROM p
#                   WHERE ", x, ">= InfectTime"
#   )
#   return (sqldf(query))
# }
# 
# # chronically infected at evaltime x (past the acute phase already)
# chroninfected <- function(x){
#   query <- paste0("SELECT ID, log10SPVL
#                   FROM p
#                   WHERE ", x, ">= InfectTime + 0.25"
#   )
#   return (sqldf(query))
# }
# 
# # dead and infected at evaltime x
# deadandinfected <- function(x){
#   query <- paste0("SELECT *
#                   FROM p
#                   WHERE ", x, ">= InfectTime AND ", x, "> TOD"
#   )
#   return (sqldf(query))
# }
# 
# 
# # dead and infected with number of infections at evaltime x
# secondaryinfections <- function(completedinfperiod, infectionscaused){
#   query <- paste0("SELECT *
#                   FROM completedinfperiod
#                   LEFT OUTER JOIN infectionscaused
#                   ON completedinfperiod.ID = infectionscaused.InfectOrigID"
#   )
#   return (sqldf(query))
# }
# 
# # alive and infected at evaltime x
# aliveandinfected <- function(infected, alive){
#   query <- paste0("SELECT ID, log10SPVL
#                   FROM infected
#                   JOIN alive
#                   USING (ID)"
#   )
#   return (sqldf(query))
# }
# 
# # alive and chronically infected at evaltime x
# aliveandchroninfected <- function(chroninfected, alive){
#   query <- paste0("SELECT ID, log10SPVL
#                   FROM chroninfected
#                   JOIN alive
#                   USING (ID)"
#   )
#   return (sqldf(query))
# }
# 
# # number of relationships per person in pop
# relsandpop <- function(relsMandW, pop){
#   query <-paste0("SELECT *
#                  FROM pop
#                  LEFT OUTER JOIN relsMandW
#                  USING (ID)"
#   )
#   return (sqldf(query))
# }
# 
# # Adding TOB to r, so that we can plot age mixing pattern
# rwithIDmdata <- function(r, p){
#   query <-paste0("SELECT *
#                  FROM r
#                  JOIN p
#                  ON r.IDm = p.ID"
#   )
#   return (sqldf(query))
# }
# 
# # Adding TOB to r, so that we can plot age mixing pattern
# pGender <- function(p, vIDs){
#   query <-paste0("SELECT ID, Gender
#                  FROM p
#                  JOIN e
#                  ON p.ID = vIDs.ID"
#   )
#   return (sqldf(query))
# }
# 
# alive <- function(x){
#   query <- paste0("SELECT ID, TOB
#                   FROM p
#                   WHERE ", x, ">= TOB AND ", x, "< TOD
#                   "#                  GROUP BY ID"
#   )
#   return (sqldf(query))
# }
# 
# alivetab <- function(x){
#   query <- paste0("SELECT *
#                   FROM DT
#                   WHERE ", x, ">= TOB AND ", x, "< TOD
# "#                  GROUP BY ID"
#   )
#   return (sqldf(query))
# }
