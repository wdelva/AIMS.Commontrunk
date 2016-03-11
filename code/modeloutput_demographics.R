# modeloutput.demographics produces a vector of summary statistics that will be compared to the target statistics.
# The summary statistics are:
# population growth rate (0 < target < 0.02)
# partner turnover rate (0.2 < target < 1 relationship per person per year)

modeloutput.all <- function(simpact.results = results, vsc = 0){
  source("/Users/wimdelva/Documents/AIMS Essays/2016/AIMS.Commontrunk/code/postsimscript.R")
  # Reading the data
  datalist <- readthedata(simpact.results)
  # Calculating population growth rate over the last 10 years of the simulation period
  start.popsize <- tail(datalist$ltable$PopSize, 11)[1]
  final.popsize <- tail(datalist$ltable$PopSize, 1)
  Exp.Growth.rate <- exp(log(final.popsize/start.popsize) / 10) #(length(datalist$ltable$PopSize) - 1))
  # Calculating partner turnover rate over the last year of the simulation period
  Number.relations <- nrow(subset(datalist$rtable, FormTime > tail(datalist$ltable$Time, 2)[1]))
  Population.size <- mean(tail(datalist$ltable$PopSize, 2))
  Log.PTR <- ifelse(Number.relations > 0, log((Number.relations) / (Population.size / 2)), log((Number.relations+0.000000001) / (Population.size / 2))) # To avoid a -Inf for Log.PTR
  #Calculating HIV prevalence 0, 5, 10, ... years after HIV is introduced
  hivtimes <- seq(from = simpact.results$input.params$config$hivseed.time,
                  to = simpact.results$input.params$config$population.simtime,
                  by = 1) # So now it goes from 10 to 45, i.e. a vector of length 36
  Exp.hivprevs <- rep(NA, length(hivtimes))
  for (i in 1:length(hivtimes)){
    hiv.time <- hivtimes[i]
    ptable.alive.18.50 <- subset(datalist$ptable, TOB <= hiv.time - 18 & TOB > hiv.time - 50 & TOD > hiv.time)
    ptable.alive.18.50$Infected <- hiv.time >= ptable.alive.18.50$InfectTime
    popsize <- nrow(ptable.alive.18.50)
    Exp.hivprevs[i] <- ifelse(popsize > 0,
                              exp(sum(ptable.alive.18.50$Infected) / popsize),
                              exp(0))
  }
  #Calculating mean and standard deviation of age gaps (aAG and sdAG)
  targetwindow.start <- simpact.results$input.params$config$population.simtime - 3
  targetwindow.end <- simpact.results$input.params$config$population.simtime
  ptable.targetwindow <- datalist$ptable[TOB < targetwindow.end & TOD > targetwindow.start]
  rtable.targetwindow <- datalist$rtable[FormTime < targetwindow.end & DisTime > targetwindow.start]
  datalist.targetwindow <- list(ptable=ptable.targetwindow, rtable=rtable.targetwindow)
  am.df <- agemixing.df.maker(datalist.targetwindow)
  Exp.aAG <- ifelse(dim(am.df)[1] > 1, exp(mean(am.df$AgeGap)), -99)
  Exp.sdAG <- ifelse(dim(am.df)[1] > 1, exp(sd(am.df$AgeGap)), exp(0))
  #Calculating other age-mixing summary statistics from lme fit
  men.lme <- tryCatch(agemixing.lme.fitter(am.df), error = agemixing.lme.errFunction) # Returns an empty list if the lme model can't be fitted
  Exp.slope <- ifelse(length(men.lme) > 0, exp(summary(men.lme)$tTable[2, 1]), exp(-9))
  median.age.men <- as.numeric(summary(am.df$AgeMaleatForm)["Median"])
  median.age.women <- as.numeric(summary(am.df$AgeFemaleatForm)["Median"])
  Exp.powerm <- ifelse(length(men.lme) > 0, exp(as.numeric(summary(men.lme)$modelStruct$varStruct)), -99)
  Exp.wsdad <- ifelse(length(men.lme) > 0, exp(men.lme$sigma), -99)
  Exp.bsdad <- ifelse(length(men.lme) > 0, exp(sqrt(getVarCov(men.lme)[1,1])), -99)
  #Calculating point prevalence of concurrency at the end of the simulation
  concurr.data <- concurrency.list(datalist)
  Exp.concurr.pointprev <- exp(tail(concurr.data$concurr.pointprevalence, 1))
  outputvector <- c(Exp.Growth.rate, Log.PTR, Exp.hivprevs, Exp.aAG, Exp.sdAG, Exp.slope, median.age.men, median.age.women, Exp.powerm, Exp.wsdad, Exp.bsdad, Exp.concurr.pointprev)
  return(outputvector)
}

VSC.modeloutput.demographics <- function(simpact.results = results){
  source("/user/data/gent/vsc400/vsc40070/simpact-test/VSC.postsimscript.R")
  # Reading the data
  datalist <- readthedata(simpact.results)
  # Calculating population growth rate over the last 10 years of the simulation period
  start.popsize <- tail(datalist$ltable$PopSize, 11)[1]
  final.popsize <- tail(datalist$ltable$PopSize, 1)
  Exp.Growth.rate <- exp(log(final.popsize/start.popsize) / 10) #(length(datalist$ltable$PopSize) - 1))
  # Calculating partner turnover rate over the last year of the simulation period
  Number.relations <- nrow(subset(datalist$rtable, FormTime > targetwindow.start))
  Population.size <- mean(tail(datalist$ltable$PopSize, 2))
  Log.PTR <- log((Number.relations+0.0000001) / (Population.size / 2)) # To avoid a -Inf for Log.PTR
  outputvector <- Exp.Growth.rate # c(Exp.Growth.rate, Log.PTR)
  return(outputvector)
}
  



modeloutput.agemixing <- function(simpact.results = results){
  datalist <- readthedata(simpact.results)
  if (nrow(datalist$rtable) < nrow(datalist$ptable)) {
    outputvector <- rep(NA, 13)
  } else {
    ptable.targetwindow <- datalist$ptable[TOB < targetwindow.end & TOD > targetwindow.start]
    rtable.targetwindow <- datalist$rtable[FormTime < targetwindow.end & DisTime > targetwindow.start]
    datalist.targetwindow <- list(ptable=ptable.targetwindow, rtable=rtable.targetwindow)
    
    # Let's calculate each summary statistic separately.
    # First, the age-mixing summary statistics
    am.df <- agemixing.df.maker(datalist.targetwindow)
    men.lme <- tryCatch(agemixing.lme.fitter(am.df), error = agemixing.lme.errFunction) # Returns an empty list if the lme model can't be fitted
    # We also need to calculate point prevalence of concurrency at the end of the simulation
    concurr.data <- concurrency.list(datalist)
    # And HIV prevalence
    prev.18.50 <- prevalence.18.49.after.15.30.list(datalist = datalist,
                                                    timepoints = c(as.numeric(simpact.results$input.params$config["hivseed.time"])+15,
                                                                   as.numeric(simpact.results$input.params$config["hivseed.time"])+30))
    
    
    # 1. aad
    aad <- mean(am.df$AgeGap)
    # 2. sdad
    sdad <- sd(am.df$AgeGap)
    # 3. slope
    slope <- ifelse(length(men.lme) > 0, summary(men.lme)$tTable[2, 1], 0.5)
    # 4. median.age.men
    median.age.men <- as.numeric(summary(am.df$AgeMaleatForm)["Median"])
    # 5. median.age.women
    median.age.women <- as.numeric(summary(am.df$AgeFemaleatForm)["Median"])
    # 6. powerm
    powerm <- ifelse(length(men.lme) > 0, as.numeric(attributes(men.lme$apVar)$Pars["varStruct.power"]), 0.5)
    # 7. wsdad
    wsdad <- ifelse(length(men.lme) > 0, men.lme$sigma, 0.5)
    # 8. bsdad
    bsdad <- ifelse(length(men.lme) > 0, getVarCov(men.lme)[1,1], 0.5)
    # 9. prev.18.50.after.15
    prev.18.50.after.15 <- prev.18.50$pointprevalence[1]
    # 10. prev.18.50.after.30
    prev.18.50.after.30 <- prev.18.50$pointprevalence[2]
    # 11. relsperpersonperyear
    relsperpersonperyear <- partnerturnoverrate(datalist) # based on relationships inside the target window (see function code)
    # 12. concurrency.pointprev
    concurrency.pointprev <- tail(concurr.data$concurr.pointprevalence, 1)
    # 13. Fraction.25ormore
    Fraction.25ormore <- Partner25ormore.girls(datalist)
    
    #   am <- tryCatch(agemixing.list(datalist.targetwindow), error = agemixing.errorFunction) # storing  age mixing data
    #   inc <- incidenceheatmapdata.list(datalist) # storing ALL HIV incidence data
    #   # rsteps <- allRtimestepsdata.list(datalist) Still a bug in this.
    #   pop <- population.list(datalist)
    #   prev.18.50 <- prevalence.18.49.after.15.30.list(datalist = datalist,
    #                                     timepoints = c(as.numeric(cfgABC["hivseed.time"])+15,
    #                                                    as.numeric(cfgABC["hivseed.time"])+30))
    #   #recycle <- recycler(datalist = datalist, timepoint = as.numeric(cfgABC["hivseed.time"]))
    #   
    # 
    #     
    #   if (is.null(am$powerm)){
    #     am$powerm <- NaN
    #     am$WVAD <- NaN
    #   }
    #   
    #   if (length(am$powerm) == 0)
    #   {
    #     am$powerm <- 3e10
    #   }
    #   
    # 
    #   median.age.men=as.numeric(summary(am$agescatterdata$AgeMaleatForm)["Median"])
    #   median.age.women=as.numeric(summary(am$agescatterdata$AgeFemaleatForm)["Median"])
    #   aad=am$AAD
    #   sdad=sqrt(am$VAD)
    #   bsdad=sqrt(am$BVAD)
    #   powerm=as.numeric(am$powerm)
    #   slope=am$slope
    #   wsdad=sqrt(am$within.var.base)
    #   wsdad.min=min(sqrt(am$WVAD))
    #   wsdad.max=max(sqrt(am$WVAD))
    #   wsdad.mean=mean(sqrt(am$WVAD))
    #   
    #   if (length(am)==0)
    #   {
    #     median.age.men <- 0 # something undesirable, but not outrageously off the scale bad
    #     median.age.women <- 0
    #     aad <- 0
    #     sdad <- 0
    #     bsdad <- 0
    #     powerm <- 0
    #     slope <- 0
    #     wsdad <- 0
    #     wsdad.min <- 0
    #     wsdad.max <- 0
    #     wsdad.mean <- 0
    #   }
    #   
    #   
    #   
    #   concurr.data <- concurrency.list(datalist)
    #   concurrency.pointprev <- tail(concurr.data$concurr.pointprevalence, 1)
    #   
    #   #agegap.diff <- agegap.concur(datalist)
    #   Fraction.25ormore <- Partner25ormore.girls(datalist)
    #   prev.18.50.after.15 <- prev.18.50$pointprevalence[1]
    #   prev.18.50.after.30 <- prev.18.50$pointprevalence[2]
    # #   starttime <- 20
    # #   window <- 10
    # #   # Harling.list <- incidence.Harling(datalist = datalist, endtime = as.numeric(cfgABC["population.simtime"]), window = window)
    # #   Harling.list <- incidence.Harling(datalist = datalist, endtime = starttime+window, window = window)
    # #   
    # # 
    # #   m.df <- Harling.list$readyforCox.Male
    # #   #m.df$left <- m.df$InfectTime - starttime #runif(nrow(m.df))*window#0
    # #   #m.df$left[m.df$event==0] <- m.df$PY[m.df$event==0] # For right-censored observations
    # #   #m.df$right <- Inf # For right-censored observations
    # #   #m.df$right[m.df$event==3] <- window/2#window # For interval-censored HIV infection events
    # #   ##m.df$event <- 3 # Now we overwrite the event variable: it must be 3 for all observations
    # #   m.df$AgeBin <- m.df$Age > 40
    # #   m.surv <- Surv(time = m.df$time, event = m.df$event, type = "right")
    # #   m.surv.age <- survfit(m.surv ~ AgeBin, data = m.df)
    # #   plot(m.surv.age)
    # #   
    # #   survreg(m.surv ~ max.AgeGap, data = m.df)
    # #   survreg(m.surv ~ last.AgeGap, data = m.df)
    # #   coxph(m.surv ~ max.AgeGap, data = m.df)
    # #   coxph(m.surv ~ last.AgeGap, data = m.df)
    # # 
    # #   coxph(m.surv ~ max.AgeGap + Age, data = m.df)
    # #   
    # #   # For women
    # #   f.df <- Harling.list$readyforCox.Female
    # #   #f.df$left <- f.df$InfectTime - starttime #runif(nrow(m.df))*window#0
    # #   #f.df$left[f.df$event==0] <- f.df$PY[f.df$event==0] # For right-censored observations
    # #   f.df$AgeBin <- f.df$Age > 40
    # #   f.surv <- Surv(time = f.df$time, event = f.df$event, type = "right")
    # #   f.surv.age <- survfit(f.surv ~ AgeBin, data = f.df)
    # #   plot(f.surv.age)
    # #   
    # #   survreg(f.surv ~ max.AgeGap, data = f.df)
    # #   survreg(f.surv ~ last.AgeGap, data = f.df)
    # #   coxph(f.surv ~ max.AgeGap, data = f.df)
    # #   coxph(f.surv ~ last.AgeGap, data = f.df)
    # #   
    # #   coxph(f.surv ~ max.AgeGap + Age + max.AgeGap*Age + MultiplePartners.lastyear, data = f.df)
    # #   
    # #   
    # #   m.misurv <- MI.surv(m=20, data=m.df, conf.int = TRUE, alpha = 0.05)
    # #   plot(m.misurv)
    # # 
    # #   # Now we use augmented data and multiple imputation to do coxph modelling with interval-censored data:
    # #   m.MIICD.cox.Partners.lastyear.Male <- MIICD.coxph(formula = ~ Partners.lastyear.Male, k = 10, m = 10, data = m.df, method = "ANDA", verbose = TRUE)
    # #   f.MIICD.cox.Partners.lastyear.Female <- MIICD.coxph(formula = ~ Partners.lastyear.Female, k = 10, m = 10, data = f.df, method = "ANDA", verbose = TRUE)
    # #   
    # # # And the Harling-like models we are most interested in:
    # #   f.MIICD.cox.Harling <- MIICD.coxph(formula = ~ Age + max.AgeGap + Partners.lastyear.Female, k = 10, m = 10, data = f.df, method = "ANDA", verbose = TRUE)
    # #   f.MIICD.cox.Harling <- MIICD.coxph(formula = ~ Age + last.AgeGap + Partners.lastyear.Female, k = 10, m = 10, data = f.df, method = "ANDA", verbose = TRUE)
    # #   
    # #   ########
    # #   # The data we prepared with tmerge (TestSurv.df)
    # #   ########
    # #   TestSurv.df <- outcome.datamaker(datalist = datalist,
    # #                                    endfollowup = as.numeric(cfgABC["population.simtime"]),
    # #                                    hivseed.time = as.numeric(cfgABC["hivseed.time"]))
    # #   # testing whether dataset was constructed correctly
    # #   head(subset(TestSurv.df, TOB > -5 & TOB <25 & !is.na(AgeGapMRP)), 10)[, c(1:3, 9, 19:25)]
    # #   subset(DTR, IDm==5891)
    # #   # Now let's assume the data are simply right-censored
    # #   TestSurv.df$Age <- TestSurv.df$survey.number - TestSurv.df$TOB
    # #   TestSurv.df <- subset(TestSurv.df, Age >15)
    # #   TestSurv.df$Age0 <- TestSurv.df$Age - min(TestSurv.df$Age)
    # #   TestSurv.df$Age0.sq <- TestSurv.df$Age0^2
    # #   #TestSurv.df$AgeGapMRP.sq <- TestSurv.df$AgeGapMRP^2
    # #   
    # #   f.TestSurv.df <- subset(TestSurv.df, Gender==1 & Age0 < 5) # only women between 15 and 20 yo
    # #   f.TestSurv.df$AgeGapMRP5 <- cut(f.TestSurv.df$AgeGapMRP, breaks = c(-15, 0, 5, 10, 30), right = FALSE) # [0,5) etc.
    # #   f.TestSurv.df$AgeGapMRP5 <- relevel(f.TestSurv.df$AgeGapMRP5, ref = "[0,5)")
    # #   f.TestSurv.df$AgeGapMRP2.5 <- cut(f.TestSurv.df$AgeGapMRP, breaks = c(-15, -2.5, 2.5, 10, 30), right = FALSE) # [-2.5, 2.5) etc.
    # #   f.TestSurv.df$AgeGapMRP2.5 <- relevel(f.TestSurv.df$AgeGapMRP2.5, ref = "[-2.5,2.5)")
    # #   f.testsurv.object <- Surv(f.TestSurv.df$tstart, f.TestSurv.df$tstop, f.TestSurv.df$hiv.postest==3)
    # #   # AgeGap as a continuous variable
    # #   fit1 <- coxph(f.testsurv.object ~ Age0+Age0.sq+Age0*AgeGapMRP+AgeGapMRP + strata(survey.number), data = f.TestSurv.df) 
    # #   baseline<-basehaz(fit1, centered=FALSE)
    # #   #draw baseline hazard (that's male)
    # #   plot(baseline$time, baseline$hazard, type='l',main="Hazard rates") 
    # #   #draw female hazard
    # #   lines(baseline$time, exp(0.8245)*baseline$hazard, col="blue") 
    # #   
    # #   
    # #   # AgeGap as categrotical variable [0,5) etc.
    # #   coxph(f.testsurv.object ~ Age0+Age0.sq+AgeGapMRP5 + strata(survey.number), data = f.TestSurv.df)
    # #   # AgeGap as categrotical variable [-2.5, 2.5) etc.
    # #   coxph(f.testsurv.object ~ Age0+Age0.sq++AgeGapMRP2.5 + strata(survey.number), data = f.TestSurv.df)
    # # 
    # #   coxph(f.testsurv.object ~ Age0+Age0.sq+Age0*AgeGapMRP+AgeGapMRP + strata(factor(survey.number)), data = f.TestSurv.df)
    # #   
    # #  # Preparing f.TestSurv.df for use with MIICD.coxph and MIICD.crreg:
    # #   # we treat each observation as a right censored one, except the ones that led to hiv.postest===3
    # #   TestSurv.df$left <- TestSurv.df$tstart # TestSurv.df$tstop - TestSurv.df$tstart
    # #   TestSurv.df$right <- TestSurv.df$tstop # Inf
    # #   TestSurv.df$left[TestSurv.df$hiv.postest!=3] <- TestSurv.df$tstop[TestSurv.df$hiv.postest!=3]
    # #   TestSurv.df$right[TestSurv.df$hiv.postest!=3] <- Inf
    # #   TestSurv.df$origin <- 0
    # #   TestSurv.df$origin[TestSurv.df$hiv.postest!=3] <- TestSurv.df$tstart[TestSurv.df$hiv.postest!=3]
    # # 
    # #   
    # #   
    # #   
    # #   # The data must contain at last four columns.
    # #   # One named left, one named right, the name of the 3^rd is indicated by the status parameter and one for the covariate to be tested.
    # #   # For interval censored data, the left and right columns indicates the lower and the upper bounds of the intervals respectively.
    # #   # Inf in the right column stands for right censored observations.
    # #   # When an observation is right censored, the status column must contain the censor indicator specified by cens.code.
    # #   # The transition of interest must be precised by the trans parameter.
    # #   
    # #   
    # #   f.MIICD.cox.Harling <- MIICD.coxph(formula = ~ Age0 + AgeGapMRP, k = 10, m = 10, data = f.TestSurv.df, method = "ANDA", verbose = TRUE)
    # #   
    # #   MIICD.crreg(formula, k, m, status, trans, cens.code, data, method = c("PMDA",
    # #                                                                         "ANDA"), verbose = FALSE)
    # #   
    # #   
    # # TestSurv.df$event <- NA
    # #   TestSurv.df$time <- NA
    # #   TestSurv.df$time2 <- NA
    # #   TestSurv.df$event[is.na(readyforCox.Male$eventname.incident)] <- 0 # Event is 2 for all censored observations: natural deaths and completed windows without infection.
    # #   TestSurv.df$time[is.na(readyforCox.Male$eventname.incident)] <- readyforCox.Male$PY[is.na(readyforCox.Male$eventname.incident)] # Event is 2 for all censored observations: natural deaths and completed windows without infection.
    # #   TestSurv.df$event[!is.na(readyforCox.Male$eventname.incident)] <- 1 # Event is 3 for all interval censored observations: HIV infection.
    # #   TestSurv.df$time[!is.na(readyforCox.Male$eventname.incident)] <- readyforCox.Male$InfectTime[!is.na(readyforCox.Male$eventname.incident)] - endtime + window
    # #   TestSurv.df$time2[!is.na(readyforCox.Male$eventname.incident)] <- readyforCox.Male$InfectTime[!is.na(readyforCox.Male$eventname.incident)] - endtime + window
    # #   
    # #   
    # #   
    # #   f.surv <- Surv(time = f.df$time, event = f.df$event, type = "right")
    # #   
    # #   
    # # 
    # #   
    # #   
    # #   m.misurv <- MI.surv(m=1, data=TestSurv.df, conf.int = TRUE, alpha = 0.05)
    # #   plot(m.misurv)
    # #   
    # #   # Now we use augmented data and multiple imputation to do coxph modelling with interval-censored data:
    # #   m.MIICD.cox.Partners.lastyear.Male <- MIICD.coxph(formula = ~ Partners.lastyear.Male, k = 10, m = 10, data = m.df, method = "ANDA", verbose = TRUE)
    # #   f.MIICD.cox.Partners.lastyear.Female <- MIICD.coxph(formula = ~ Partners.lastyear.Female, k = 10, m = 10, data = f.df, method = "ANDA", verbose = TRUE)
    # #   
    # #   # And the Harling-like models we are most interested in:
    # #   f.MIICD.cox.Harling <- MIICD.coxph(formula = ~ Age + max.AgeGap + Partners.lastyear.Female, k = 10, m = 10, data = f.df, method = "ANDA", verbose = TRUE)
    # #   f.MIICD.cox.Harling <- MIICD.coxph(formula = ~ Age + last.AgeGap + Partners.lastyear.Female, k = 10, m = 10, data = f.df, method = "ANDA", verbose = TRUE)
    # #   
    # #   
    # #   
    # #     
    # #   m.Surv <- Surv(time = m.df$left, time2 = m.df$right, type="interval2")
    # #   Fit <- intcox(m.Surv ~ Age, data=m.df)
    # #   
    # #   
    # #   survreg(m.Surv ~ AgeBin, data = m.df)
    # #   
    # #   
    # # 
    # #   m.MIICD.cox.Partners.lastyear.Male <- MIICD.coxph(formula = ~ Partners.lastyear.Male, k = 10, m = 10, data = m.df, method = "ANDA", verbose = TRUE)
    # #   
    # #   plot(m.misurv$est$time, m.misurv$est$surv)
    # #   m.Surv <- Surv(time = m.df$time, time2 = m.df$time2, event = m.df$event, type="interval")
    # #   survreg(m.Surv ~ AgeBin, data = m.df)
    # #   
    # #   plot(m.Surv)
    # #   m.cox <- survfit(m.Surv ~ 1, data = m.df)
    # #   plot(m.cox.Age)
    # #   
    # #   print(c("slope = ", slope))
    # #   print(c("powerm = ", powerm))
    # #   print(c("wsdad = ", wsdad))
    # #   print(c("bsdad = ", bsdad))
    # #   print(c("relsperpersonperyear =", relsperpersonperyear))
    
    outputvector <- c(aad,
                      sdad,
                      slope,
                      median.age.men,
                      median.age.women,
                      powerm,
                      wsdad,
                      bsdad,
                      prev.18.50.after.15,
                      prev.18.50.after.30,
                      relsperpersonperyear,
                      concurrency.pointprev,
                      Fraction.25ormore)
  }
  return(outputvector)
}
