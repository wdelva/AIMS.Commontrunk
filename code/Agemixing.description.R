# ==================
# Describe age-mixing pattern in LNS data
# ==================

# loading the essential libraries
library(nlme)
library(dplyr) # Data management and manipulation
library(magrittr) # Allows me to use pipe operator: %>%
library(GGally)

# loading the 50 datasets with imputed data, for each sex
load(mlong) # contains dataset "longm" which has 51 versions of the male dataset. The first 1688 rows are the original dataset, the next are the versions with imputations
load(mwide)
load(flong)
load(fwide)


# cycling through the orginal AND imputed datasets
imputations.min.max <- range(as.numeric(as.character(longm$.imp)))
imputations <- seq(from = imputations.min.max[1], to = imputations.min.max[2])
nrows <- length(imputations)

# preparing dataset
agemixing.features <- data.frame(aad.m = rep(NA, nrows),
                                 sdad.m = rep(NA, nrows),
                                 slope.m = rep(NA, nrows),
                                 power.m = rep(NA, nrows),
                                 wisd.m = rep(NA, nrows),
                                 bisd.m = rep(NA, nrows),
                                 median.maleage.m = rep(NA, nrows),
                                 median.femaleage.m = rep(NA, nrows))


for(i in 0:50) {
  
  #Smaller datasets
  men <- longm %>%
    filter(.imp == i)%>%
    select(id, age, age0, agep, missagep)
  
  if (i==0) { # because the lme models only fit for complete datasets
    men <- men[men$missagep=="Observed", ]
  }
  
  women <-longw %>%
    filter(.imp == i)%>%
    select(id, age, age0, agep, missagep)
  
  if (i==0) { # because the lme models only fit for complete datasets
    women <- women[women$missagep=="Observed", ]
  }
  
  #Now fit separate models for men and women
  #varPower variance structure for each model
  #The formula to calculate the weights for variance is |v|^(2*t)
  #Age can't be at 0 in varPower formula because then it will evaluate to 0 variance for the first level (18 year olds)
  m1m <- lme(agep ~ age0, 
             data = men,
             control=lmeControl(returnObject=TRUE),#, maxIter = 5000, msMaxIter = 5000),
             random = ~1 | id,  
             method = "REML",
             weights = varPower(value = 0.7, form = ~age0 + 1))
  
#   m1w <- lme(agep ~ age0, 
#              data = women,
#              control=lmeControl(returnObject=TRUE),
#              random = ~1 | id,
#              method = "REML",
#              weights = varPower(value = 0.5, form = ~age0 + 1))
  
  aad.m = mean(men$age - men$agep)
  sdad.m = sd(men$age - men$agep)
  slope.m <- summary(m1m)$tTable[2, 1]
  power.m <- as.numeric(summary(m1m)$modelStruct$varStruct)
  wisd.m <- m1m$sigma # within-individual standard deviation (for an 18 year old man)
  bisd.m <- sqrt(getVarCov(m1m)[1,1]) # between-individual standard deviation 
  median.maleage.m <- median(men$age)
  median.femaleage.m <- median(men$agep)
  
  # inserting features in dataframe
  agemixing.features[(i+1), ] <- c(aad.m, sdad.m, slope.m, power.m, wisd.m, bisd.m, median.maleage.m, median.femaleage.m)
  print(i)
}
agemixing.features$original <- factor(c(1, rep(0, 50)))


# Visualise distribution of features
ggscatmat(data = agemixing.features,
          columns = 1:6,
          color = "original")

# Get mean values of all features
colMeans(agemixing.features[, 1:8])
