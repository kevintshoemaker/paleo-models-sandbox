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

<<<<<<< HEAD
DAMIEN = TRUE
=======
KEVIN = TRUE # FALSE  # TRUE
DAMIEN = FALSE
>>>>>>> origin/master
HRA_LAB = FALSE # TRUE

####################
#   PRELIMINARY: SET GLOBAL PARAMETERS  (USER SPECIFIED PARAMS)
####################

NREPS <- 5     # number of samples to draw from the multivariate uniform prior (LHS) per niche breadth value 
##NicheBreadths = c(30,40,50,60,70,80,90,100)   
NicheBreadths = c(40)                         # NOTE: changed to loop through single NB for testing

####################
#   PRELIMINARY: LOAD FUNCTIONS
#################### 

if(KEVIN) CODE_DIRECTORY <<- "C:\\Users\\Kevin\\GIT\\paleo-models-sandbox"       # code directory should be your local copy of the GitHub repository   
if(HRA_LAB) CODE_DIRECTORY <<- "C:\\Users\\Akcakaya\\Desktop\\Mammoth Model\\paleo-models-sandbox"
if(DAMIEN) CODE_DIRECTORY <<-  "C:\\Users\\Damien Fordham\\Documents\\GitHub\\paleo-models-sandbox"

setwd(CODE_DIRECTORY)
source("Rfunctions_PALEO_UTILITY.r")     # Load all functions for Paleo project
source("Rfunction_PALEO_MAKEMPs.r")
#source("Rfunction_PALEO_RUNMPs.r")

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

<<<<<<< HEAD
##NicheBreadths = c(30,40,50,60,70,80,90,100)
NicheBreadths = c(40)

FirstNicheBreadth = NicheBreadths[1]

## set first niche breadth to loop over
nb=40

for(nb in NicheBreadths){
=======
nb=40  # for testing...
for(nb in NicheBreadths){     # Loop through niche breadths
>>>>>>> origin/master
  
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
    MakeMPfile(f=i,masterDF=masterDF,NicheBreadth=NicheBreadth)
  }
  
}   # end loop through niche breadths

<<<<<<< HEAD
=======



#############################
#       STEP 2. RUN THE MP FILES (in parallel)
#############################







>>>>>>> origin/master
