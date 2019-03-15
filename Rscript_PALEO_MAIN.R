##########################
# This master script writes all the MP files for the wooly mammoth model
#  Authors: Kevin Shoemaker and Damien Fordham
# 
#  22 May 2015 -- started scripting
#  11 June 2015 -- added dispersal correction
#  30 June 2016 -- cleaned up for initiation of model runs 
##########################


#######################
#   CLEAR THE WORKSPACE
#######################

rm(list=ls())

#########################
#   PRELIMINARY: SET PROJECT DIRECTORIES
#########################

KEVIN = FALSE # FALSE  # TRUE
DAMIEN = FALSE
HRA_LAB = FALSE # TRUE
KEVIN_LAB = TRUE

####################
#   PRELIMINARY: SET GLOBAL PARAMETERS  (USER SPECIFIED PARAMS)
####################

NREPS <- 2     # number of samples to draw from the multivariate uniform prior (LHS) per niche breadth value 
##NicheBreadths = c(30,40,50,60,70,80,90,100)   
NicheBreadths = c(70)                         # NOTE: changed to loop through single NB for testing

####################
#   PRELIMINARY: LOAD FUNCTIONS
#################### 

if(KEVIN) CODE_DIRECTORY <- "C:\\Users\\Kevin\\GIT\\paleo-models-sandbox"       # code directory should be your local copy of the GitHub repository   
if(HRA_LAB) CODE_DIRECTORY <- "C:\\Users\\Akcakaya\\Desktop\\Mammoth Model\\paleo-models-sandbox"
if(DAMIEN) CODE_DIRECTORY <-  "C:\\Users\\Damien Fordham\\Documents\\GitHub\\paleo-models-sandbox"
if(KEVIN_LAB) CODE_DIRECTORY <- "E:\\GIT\\paleo-models-sandbox"

setwd(CODE_DIRECTORY)
source("Rfunctions_PALEO_UTILITY.r")     # Load all functions for Paleo project
source("Rfunction_PALEO_MAKEMPs.r")
source("Rfunction_PALEO_RUNMPs.r")
source("Rfunction_PALEO_EXTRACTRESULTS.r")

####################
#  PRELIMINARY: SET UP WORKSPACE AND LOAD PACKAGES
####################

SetUpWorkspace()   # function loads packages and sets up the workspace...  
num_cores <- detectCores() - 1   # for setting up cluster... leave one core free for windows background processes?

####################
#  PRELIMINARY: DO DISPERSAL PRECALCULATIONS
#          for each cell, determine which neighboring cells are in each distance bin... 
####################

DistBins <- DispersalPreCalculations()

#####################
#  READ THE MP TEMPLATE FOR BUILDING MORE MP FILES

template <- ReadMPTemplate()

#############################
#       STEP 1. GENERATE THE MP FILES (parallel)
#############################

nb=70  # for testing...
for(nb in NicheBreadths){     # Loop through niche breadths

  NicheBreadth <- nb   # set the current niche breadth
  
  #####################
  #  SAMPLE FROM LATIN HYPERVOLUME
  #        result is a data frame with each row representing the key parameters for an mp file
  
  dir <- sprintf("%s\\Sample_%s",KCH_DIRECTORY,nb)
  if(nb==NicheBreadths[1]){
    masterDF <- MakeLHSSamples(nicheBreadthDir=dir,NicheBreadth)    
  } else{
    masterDF <- rbind(masterDF,MakeLHSSamples(dir,nb))
  }
  
  registerDoParallel(cores=num_cores)    # make the cluster
  
  #######################
    ## objects to export to each node in the cluster
  
  functionlist <- c('mp.read', 'mp.write')
  filelist <- c('MP_DIRECTORY','template','GENTIME','humanArrival.df','EXE_DIRECTORY','DLL_FILENAME',
                 'dispersalFunc.df','DistClasses','NPOPS','DistBins','masterDF','NicheBreadth')
  
  objectlist <- c(functionlist,filelist)   # full list of objects to export
  #packagelist <- c()
  
  all.mps <- foreach(i = 1:nrow(masterDF),
                       .export=objectlist,
                       #.packages = c("R2WinBUGS"),
                       .errorhandling=c("pass")
  ) %dopar% {   
    MakeMPfile(f=i,masterDF=masterDF,NicheBreadth=NicheBreadth,suspendtime=i)  # KTS: added and removed "newkch" argument to function, which uses existing KCH files instead of making new ones...
  }
  
}   # end loop through niche breadths

#  WRITE MODEL SPECS/METADATA
setwd(MP_DIRECTORY)
## name file for LHS parameters 
write.csv(masterDF,"masterDF.csv",row.names=F)


#############################
#       STEP 2. RUN THE MP FILES (in parallel)
#############################

#  READ MODEL SPECS/METADATA

###   read in the LHS draws/MP file specs
setwd(MP_DIRECTORY)
masterDF <- read.csv("masterDF.csv",header=T,stringsAsFactors = F)
nfiles <- nrow(masterDF)
 #MPsToRun <- as.character(masterDF$MPFilename)    # list of all MP files to run

nb=70  # for testing...
for(nb in NicheBreadths){     # Loop through niche breadths
  
  NicheBreadth <- nb   # set the current niche breadth
  
  registerDoParallel(cores=num_cores)    # make the cluster
  
  #######################
  ## objects to export to each node in the cluster
  
  functionlist <- c()   # 'mp.read', 'mp.write'
  filelist <- c('masterDF','NicheBreadth','MP_DIRECTORY') #,'template','GENTIME','humanArrival.df','EXE_DIRECTORY','DLL_FILENAME','dispersalFunc.df','DistClasses','NPOPS','DistBins',
  
  objectlist <- c(functionlist,filelist)   # full list of objects to export
  #packagelist <- c()
  
  all.mps <- foreach(i = 1:nrow(masterDF),
                     .export=objectlist,
                     #.packages = c("R2WinBUGS"),
                     .errorhandling=c("pass")
  ) %dopar% {   
    RunMPfile(f=i,masterDF=masterDF,NicheBreadth=NicheBreadth,ramasloc=METAPOP_LOCATION,suspendtime=i)
  }
  
}   # end loop through niche breadths



#############################
#       STEP 3. EXTRACT RESULTS FROM THE MP FILES (in parallel)
#############################

#  READ MODEL SPECS/METADATA

###   read in the LHS draws/MP file specs
setwd(MP_DIRECTORY)
masterDF <- read.csv("masterDF.csv",header=T,stringsAsFactors = F)
nfiles <- nrow(masterDF)

nb=70  # for testing...
for(nb in NicheBreadths){     # Loop through niche breadths
  
  NicheBreadth <- nb   # set the current niche breadth
  
  cl <- makeCluster(num_cores,outfile="LOG.TXT")
  registerDoParallel(cl=cl)    # make the cluster
  
  #######################
  ## objects to export to each node in the cluster
  
  functionlist <- c()   # , 'mp.write'
  filelist <- c('masterDF','NicheBreadth','NPOPS','TIMESTEPS','MP_DIRECTORY','GridCellAttributes')  #'MP_DIRECTORY','template','GENTIME','humanArrival.df','EXE_DIRECTORY','DLL_FILENAME','dispersalFunc.df','DistClasses','NPOPS','DistBins',
  
  objectlist <- c(functionlist,filelist)   # full list of objects to export
  #packagelist <- c()
  
  
  all.mps <- foreach(i = c(1), #1:nrow(masterDF),
                     .export=objectlist,
                     .packages =c("data.table"),     #c("sp","adehabitatHR","geosphere"),
                     .errorhandling=c("pass")
  ) %dopar% {   
    #ExtractMPresults(f=i,masterDF=masterDF,NicheBreadth=NicheBreadth,doMCP=TRUE,suspendtime=i)
    ExtractMPresults2(f=i,masterDF=masterDF,NicheBreadth=NicheBreadth,doMCP=TRUE,suspendtime=i)
  }
  
}   # end loop through niche breadths








################
# EXTRACT HARVEST RESULTS
################

template <- ReadMPTemplate()   # read the MP template

nb=70  # for testing...
for(nb in NicheBreadths){     # Loop through niche breadths
  
  NicheBreadth <- nb   # set the current niche breadth
  
  registerDoParallel(cores=num_cores)    # make the cluster
  
  #######################
  ## objects to export to each node in the cluster
  
  functionlist <- c('mp.read', 'mp.write')
  filelist <- c('MP_DIRECTORY','template','GENTIME','humanArrival.df','EXE_DIRECTORY','DLL_FILENAME',
                'dispersalFunc.df','DistClasses','NPOPS','DistBins','masterDF','NicheBreadth')
  
  objectlist <- c(functionlist,filelist)   # full list of objects to export
  #packagelist <- c()
  
  all.files <- foreach(i = 1:nrow(masterDF),
                     .export=objectlist,
                     #.packages = c("R2WinBUGS"),
                     .errorhandling=c("pass")
  ) %dopar% {   
    GetHarvest(f=i,masterDF=masterDF,NicheBreadth=NicheBreadth)  
  }
  
}   # end loop through niche breadths





