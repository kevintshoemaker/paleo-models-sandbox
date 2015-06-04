##########################
# This script writes all the MP files for the wooly mammoth model
#  Authors: Kevin Shoemaker and Damien Fordham
# 
#  22 May 2015 -- started scripting

#TODO
 

######################
#   CLEAR THE WORKSPACE
######################

rm(list=ls())


######################
# SET PROJECT DIRECTORIES
######################
	  
KEVIN = TRUE
DAMIEN = FALSE

if(KEVIN) BASE_DIRECTORY <- "C:\\Users\\Kevin\\Dropbox\\Damien Fordham\\Mammoth Model"      
	                                                                            # code directory should be your local copy of the GitHub repository
if(KEVIN) CODE_DIRECTORY <- "C:\\Users\\Kevin\\GIT\\paleo-models-sandbox"           # directory for storing the r functions and scripts

DATA_DIRECTORY <- paste(BASE_DIRECTORY,"\\data",sep="")                         # directory for storing relevant data (CSV files)
if(is.na(file.info(DATA_DIRECTORY)[1,"isdir"])) dir.create(DATA_DIRECTORY)

MP_DIRECTORY <- paste(BASE_DIRECTORY,"\\mpfiles",sep="")                        # directory for storing and running MP files 
if(is.na(file.info(MP_DIRECTORY)[1,"isdir"])) dir.create(MP_DIRECTORY)

FIGURES_DIRECTORY <- paste(BASE_DIRECTORY,"\\figures\\raw",sep="")              # directory for storing raw figures 
if(is.na(file.info(FIGURES_DIRECTORY)[1,"isdir"])) dir.create(FIGURES_DIRECTORY)

HS_DIRECTORY <- paste(BASE_DIRECTORY,"\\hs",sep="")                             # directory for storing habitat suitability files (used for creating KCH files )
if(is.na(file.info(HS_DIRECTORY)[1,"isdir"])) dir.create(HS_DIRECTORY)

RESULTS_DIRECTORY <- paste(BASE_DIRECTORY,"\\results",sep="")                             # directory for storing habitat suitability files (used for creating KCH files )
if(is.na(file.info(RESULTS_DIRECTORY)[1,"isdir"])) dir.create(RESULTS_DIRECTORY)


####################
#  GLOBAL PARAMETERS  (USER SPECIFIED PARAMS)
####################

NREPS <- 50     # number of samples to draw from the multivariate uniform prior (LHS) 
TIMESTEPS <- 2468 #1000
NPOPS <- 7000
GENTIME <- 25
MP_TEMPLATE <- "Mammoth_7000_demo.mp"      # template MP file
VERBOSE=FALSE

      # save global params
setwd(MP_DIRECTORY)
save(NREPS,TIMESTEPS,NPOPS,GENTIME,MP_TEMPLATE,file="GlobalParams.RData")


####################
#  LOAD FUNCTIONS
#################### 

setwd(CODE_DIRECTORY)
source("Rfunctions_PaleoModels.r")

#############################
#  LOAD PACKAGES
#############################

loadPackage("fields")     # note: loadPackage should install the package from CRAN automatically if it is not already installed
loadPackage("lhs")
loadPackage("RCurl")

                 # Read in the functions for reading and writing MP files, from Matt Lammens

baseurl = "https://raw.githubusercontent.com/mlammens/SACode_Sandbox/master/"				 
source_github(baseurl,"mp.write.r")
source_github(baseurl,"mp.read.r")
source_github(baseurl,"metapopversion.r")
source_github(baseurl,"fill.matrix.df.r")
				 

###############################
#        LOAD DATA AND OTHER RELEVANT FILES
###############################

setwd(DATA_DIRECTORY)
humanArrival.df <- read.csv("Human_arrival_7000cells.csv",header=T)  # for each grid cell (population), earliest and latest human arrival dates (years before present))
dispersalFunc.df <- read.csv("dispersalfunctions.csv",header=T)      # links LHS samples to b and Dmax params

allHSfiles <- list.files(path=HS_DIRECTORY,pattern=".hs")            # char vector of all HS files

###############################
#        SPECIFY PARAMETER RANGES
################################

LHSParms <- list()    # initialize the container for parameter bounds

	 #### HS files (habitat suitability- ranging from zero to 1)
	       ### HS files have the pattern: pop1_rep1.hs
setwd(HS_DIRECTORY)
nHS <- length(list.files(path=HS_DIRECTORY,pattern="pop1_"))   #number of replicate KCH files  	
LHSParms <- specifyLHSParam(paramslist=LHSParms,name="HS",type="CAT",lb=1,ub=nHS)

     ####  RMAX
LHSParms <- specifyLHSParam(LHSParms,"RMAX",type="CONT",lb=1.28,ub=6.84)

     #### TEMPORAL VARIABILITY
LHSParms <- specifyLHSParam(LHSParms,"SD",type="CONT",lb=0,ub=0.175)

     ### ALLEE EFFECT    0 to 100 individuals...  
LHSParms <- specifyLHSParam(LHSParms,"ALLEE",type="CONT",lb=0,ub=100)   
   
     ### DENSITY PER 100X100 GRID CELL
LHSParms <- specifyLHSParam(LHSParms,"DENSITY",type="CONT",lb=1000,ub=40000)   

     ### DISPERSAL #1: PROPORTION DISPERSING (used to compute "a" parameter.  a = p/2)
LHSParms <- specifyLHSParam(LHSParms,"DISP1",type="CONT",lb=0.05,ub=0.25)

     ### DISPERSAL #2: 95% QUANTILE ON DISTANCE MOVED  [will translate into b and DMax according to Resit's pre-calculations]
LHSParms <- specifyLHSParam(LHSParms,"DISP2",type="CONT",lb=100,ub=500) 
  
  ### HUNTING  (proportion harvested per generation)
LHSParms <- specifyLHSParam(LHSParms,"HUNT",type="CONT",lb=0,ub=0.05)  
  
  #### HUMAN ARRIVAL (0 represents lower bound on a per-population basis, 1 represents upper bound)
nHS <- length(list.files(path=HS_DIRECTORY,pattern="pop1"))   #number of replicate KCH files  	
LHSParms <- specifyLHSParam(LHSParms,"HUMAN",type="CAT",lb=0,ub=1)  

##################
##### GENERATE LATIN HYPERCUBE SAMPLE

nVars <- length(names(LHSParms))

LHS <- randomLHS(NREPS, nVars )   # generate multiple samples from paramter space according to a LHS sampling scheme

masterDF <- as.data.frame(LHS)    #  storage container (data frame) to record relevant details for each MP file. Rows:MP file/LHS samples. Cols: relevant variables
	  

	  ### translate raw lhs samples into desired parameter space
colnames(masterDF) <- names(LHSParms)
parm=1
for(parm in 1:nVars){
  if(LHSParms[[parm]]$type=="CONT"){
     masterDF[,parm] <- LHSParms[[parm]]$lb + LHS[,parm]*(LHSParms[[parm]]$ub-LHSParms[[parm]]$lb)
  }
  if(LHSParms[[parm]]$type=="CAT"){
     masterDF[,parm] <- ceiling(LHSParms[[parm]]$lb + LHS[,parm]*(LHSParms[[parm]]$ub-LHSParms[[parm]]$lb))
  }
}

masterDF$dispndx <- as.numeric(cut(masterDF$DISP2,breaks=c(dispersalFunc.df$Distl,Inf)))   

#############################
#       GENERATE THE MP FILES
#############################

#### read in template MP file
setwd(MP_DIRECTORY)
mpfilein <- "Mammoth_7000_demo.mp"
template <- mp.read(mpfilein)   # read in template


###########################################
         ## BUILD THE MP FILES

masterDF$filenames <- character(1)                              

for(f in 1:NREPS){
  filename <- sprintf("LHS_Sample%s.mp",f)
  masterDF$filenames[f] <- filename 

     ### write out the KCH files for this model run (and get init abundances)   NOTE: KCH files are stored in the MP_DIRECTORY
  initabunds <- writeKCH(HSnum = masterDF$HS[f], density = masterDF$DENSITY[f])

     #### extract the main population arrays (to be edited)
  poplist1 <- template$mp.file$PopList
  poplist2 <- template$mp.file$PopData_df 

  poplist2$KchangeSt <- as.character(poplist2$KchangeSt)   

  for(i in 1:NPOPS){

              # SET GROWTH RATE
    poplist1[[i]]$MaxR <- masterDF$RMAX[f]   
    poplist2$MaxR[i] <-  masterDF$RMAX[f]  

              # SET KCH FILE
    KCHfilename <- sprintf("pop%s_rep%s.kch",i,masterDF$HS[f])
    poplist1[[i]]$KchangeSt <- KCHfilename
    poplist2$KchangeSt[i] <-  KCHfilename
 
              # SET INIT ABUND
    poplist1[[i]]$InitAbund <- round(initabunds[i])
    poplist2$InitAbund[i] <- round(initabunds[i])
  
              # SET ALLEE EFFECT
    poplist1[[i]]$localthr <- round(masterDF$ALLEE[f])
    poplist2$localthr[i] <- round(masterDF$ALLEE[f])

 }

        # start the new mp file
  version <- template$version
  newmp <- template$mp.file
  newmp$PopList <- poplist1
  newmp$PopData_df <- poplist2

  newmp$UseDispDistFunc <- TRUE
  newmp$UseCorrDistFunc <- TRUE

          # SET DISPERSAL
  params <- template$mp.file$DispDistFunc
  params[1] <- masterDF$DISP1[f]/2                                # "a" parameter describes the proportion staying/leaving a given population
  params[2] <- dispersalFunc.df$b[masterDF$dispndx[f]]            # "b" parameter describes the distance moved
  params[4] <- dispersalFunc.df$Dmax[masterDF$dispndx[f]]         # "Dmax" parameter limits the maximum distance moved.
  newmp$DispDistFunc <- params

         # SET MANAGEMENT ACTIONS (HARVEST
  
  newmp$NPopManage <- NPOPS                                       # num of management (harvest) actions equals the number of populations
  harvestdf <- newmp$PopManageProp                   
  mpstarttime <- newmp$MaxDur*GENTIME                             # for translating between human arrival (ybp) and generations in the model
  
  for(i in 1:NPOPS){
    if(i!=1) harvestdf <- rbind(harvestdf,newmp$PopManageProp) 
    harvstarttime <- humanArrival.df$year_min[i] + masterDF$HUMAN[f]*(humanArrival.df$year_max[i]-humanArrival.df$year_min[i])
    harvestdf$begin.time[i] <- round((mpstarttime - harvstarttime)/GENTIME)   # set start time for harvest
    harvestdf$end.time[i] <- newmp$MaxDur
    harvestdf$from.pop[i] <- i                                                # set population ID for management actions
    harvestdf$to.pop[i] <- i
    harvestdf$proportion <- round(masterDF$HUNT[f],5)                         # set proportion harvested
  }
  
  newmp$PopManageProp <- harvestdf                                # append all harvest-related parameters
  
         # SET INIT ABUND (again)
  newmp$StInit[,1] <- round(initabunds)                           # not sure which one is used, so set both at pop and global levels.

        # SET TEMP VAR
 newmp$SDMatr[[1]]$Matr[1,1] <- masterDF$SD[f]                    # set temporal variability



###############
      # write the new mp file
  setwd(MP_DIRECTORY)                                               
  mp.write(newmp,version,masterDF$filenames[f])                   # write MP file
 
}     # close loops through scenarios


##########################
  ####           SAVE METADATA FOR EACH MP FILE TO FILE

setwd(MP_DIRECTORY)
write.csv(masterDF,"masterDF.csv",row.names=F)


#######################
#######################





