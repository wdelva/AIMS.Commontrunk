# ==================
# Set up directories 
# ==================

# Define the project path to files
pp.0 <- "/Users/wimdelva/Documents/AIMS Essays/2016/AIMS.Commontrunk/" # Local environment. You can change this line to match your own directory
pp.1 <- "/user/data/gent/vsc400/vsc40070/simpact-test/" # VSC environment. Do not edit this line. 
ncluster.0 <- 8 # local environment has 8 cores
ncluster.1 <- 24 # nodes on VSC have 24 cores

project.paths <- c(pp.0, pp.1)
nclusters <- c(ncluster.0, ncluster.1)




