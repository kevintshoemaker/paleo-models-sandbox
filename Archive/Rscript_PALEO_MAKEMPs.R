##########################
# This script writes all the MP files for the wooly mammoth model
#  Authors: Kevin Shoemaker and Damien Fordham
# 
#  22 May 2015 -- started scripting
#  11 June 2015 -- added dispersal correction

#TODO
 

#######################
#   CLEAR THE WORKSPACE
#######################

rm(list=ls())

#########################
# SET PROJECT DIRECTORIES
#########################
	  
KEVIN = FALSE  #TRUE
DAMIEN = FALSE
HRA_LAB = TRUE

####################
#  GLOBAL PARAMETERS  (USER SPECIFIED PARAMS)
####################

NREPS <- 5     # number of samples to draw from the multivariate uniform prior (LHS) 

####################
#  LOAD FUNCTIONS
#################### 

setwd(CODE_DIRECTORY)
source("Rfunctions_PaleoModels.r")     # Load all functions for Paleo project

####################
#  SET UP WORKSPACE AND LOAD PACKAGES

SetUpWorkspace()   # function loads packages and sets up the workspace...  
	
###############################
#        SPECIFY PARAMETER RANGES
################################

LHSParms <- list()    # initialize the container for parameter bounds

	 #### HS files (habitat suitability- ranging from zero to 1)
	 ###  HS files have the pattern: pop1_rep1.hs
setwd(KCH_DIRECTORY)
nHS <- length(list.files(path=KCH_DIRECTORY,pattern="Pop_1"))   #number of replicate KCH files  	
LHSParms <- specifyLHSParam(paramslist=LHSParms,name="HS",type="CAT",lb=1,ub=nHS)

     ####  RMAX
LHSParms <- specifyLHSParam(LHSParms,"RMAX",type="CONT",lb=1.28,ub=6.84)

     #### TEMPORAL VARIABILITY
LHSParms <- specifyLHSParam(LHSParms,"SD",type="CONT",lb=0,ub=0.175)

     ### ALLEE EFFECT    0 to 100 individuals...  
LHSParms <- specifyLHSParam(LHSParms,"ALLEE",type="CONT",lb=0,ub=100)   
   
     ### DENSITY PER 100kmX100km GRID CELL
LHSParms <- specifyLHSParam(LHSParms,"DENSITY",type="CONT",lb=1000,ub=40000)   

     ### DISPERSAL #1: PROPORTION DISPERSING (used to compute "a" parameter.  a = p/2)
LHSParms <- specifyLHSParam(LHSParms,"DISP1",type="CONT",lb=0.05,ub=0.25)

     ### DISPERSAL #2: MAXIMUM DISTANCE MOVED  [will translate into b and DMax according to Resit's pre-calculations]
LHSParms <- specifyLHSParam(LHSParms,"DISP2",type="CONT",lb=100,ub=500) 
  
  ### HUNTING  (proportion harvested per year)   # NOTE: this is specified on an annual time scale (not generational). 
LHSParms <- specifyLHSParam(LHSParms,"HUNT",type="CONT",lb=0,ub=0.05)

  ### DENSITY DEPENDENCE ON HARVEST (y intercept of the harvest rate/abundance relationship)

HUNTDD_posvals <- c(seq(-0.05,0,length=10),seq(0.1,1,length=10),seq(1.1,2,length=10))
LHSParms <- specifyLHSParam(LHSParms,"HUNTDD",type="CAT",lb=1,ub=length(HUNTDD_posvals))
  
  #### HUMAN ARRIVAL (0 represents lower bound on a per-population basis, 1 represents upper bound)
LHSParms <- specifyLHSParam(LHSParms,"HUMAN",type="CONT",lb=0,ub=1)    # NOTE: was "CAT" not sure why

##################
##### GENERATE LATIN HYPERCUBE SAMPLE

nVars <- length(names(LHSParms))  

LHS <- randomLHS(NREPS, nVars )   # generate multiple samples from parameter space according to a LHS sampling scheme

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

masterDF$HUNTDD <- HUNTDD_posvals[masterDF$HUNTDD]    # value input into RAMAS DLL for density dependent hunting 

#############################
# DISPERSAL PRE-CALCULATIONS
#############################
 
#?distMeeus
#?distm

#distMeeus(GridCellAttributes[1,c(1,2)],GridCellAttributes[2,c(1,2)])/1000

                                           # NPOPS*NPOPS matrix of pairwise distances between grid centroids (km) [uses geosphere package]
distmatrix <- distm(GridCellAttributes[,c(1,2)],GridCellAttributes[,c(1,2)],fun=distMeeus)/1000

    # Distance classes for computing dispersal matrix (in km)
MinDist <- 100
MaxDist <- 500
BinSize <- 20
DistClasses <- c(1,seq(MinDist,MaxDist,BinSize))

DistBins <- list()
for(p in 1:NPOPS){
                     # Determine the distance between this population and all others, in km
  DistBins[[p]] <- list()
  for(b in 2:length(DistClasses)){
    inBin <- which( (distmatrix[,p]>DistClasses[b-1]) & (distmatrix[,p]<=DistClasses[b]) )
    if(length(inBin)>0){
      DistBins[[p]][[b-1]] <- inBin
    }else{
      DistBins[[p]][[b-1]] <- NA
    }
  }
}

rm(distmatrix)   # save memory and remove the distance matrix...

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

   ## check to make sure the number of populations is correct- if not, correct it!
NPOPSold <- nrow(template$mp.file$PopData_df)
if(NPOPSold<NPOPS){
  #### extract the main population arrays (to be edited)
  poplist1 <- template$mp.file$PopList
  poplist2 <- template$mp.file$PopData_df 
  # start the new mp file
  version <- "60" #template$version
  newmp <- template$mp.file
  
  newmp$PopData_df$name <- as.character(newmp$PopData_df$name)
  
  templatePopLine <- poplist2[1,]
  templatePopList <- poplist1[[1]]
  
  newIDs <- paste("Pop",(NPOPSold+1):NPOPS)  
  counter=1 
  i=(NPOPSold+1)
  for(i in (NPOPSold+1):NPOPS){
    newmp$PopList[[i]] <- templatePopList 
    newmp$PopList[[i]]$name <- newIDs[counter]
    newmp$PopList[[i]]$X_coord <- GridCellAttributes$x.cord_modified[i]
    newmp$PopList[[i]]$Y_coord <- GridCellAttributes$y.cord[i]
    newmp$StInit <- rbind(newmp$StInit,newmp$StInit[1,])
    newmp$PopData_df <- rbind(newmp$PopData_df,templatePopLine)
    newmp$PopData_df[i,]$name <- newIDs[counter]
    newmp$PopData_df[i,]$X_coord <- GridCellAttributes$x.cord_modified[i]
    newmp$PopData_df[i,]$Y_coord <- GridCellAttributes$y.cord[i]  
    counter=counter+1
  }
  # tail(newmp$PopData_df)
    # write the MP file
  
  newmp$PopData_df$X_coord <- GridCellAttributes$x.cord_modified
  newmp$PopData_df$Y_coord <- GridCellAttributes$y.cord 
  
  setwd(MP_DIRECTORY)                                               
  mp.write(newmp,version,"NewMammothTemplate.mp")                   # write MP file
  MP_TEMPLATE <- "NewMammothTemplate.mp" 
  
  #### read in template MP file
  setwd(MP_DIRECTORY)
  mpfilein <- MP_TEMPLATE  #"Mammoth_7000_demo.mp"
  template <- mp.read(mpfilein)   # read in template
}



f=2
for(f in 1:NREPS){
  filename <- sprintf("LHS_Sample%s.mp",f)
  masterDF$filenames[f] <- filename 
  
     ### write out the KCH files for this model run (and get init abundances)   NOTE: KCH files are stored in the MP_DIRECTORY
     # initabunds <- writeKCH(HSnum = masterDF$HS[f], density = masterDF$DENSITY[f],areas=GridCellAttributes$Area2)     ### NOTE: takes a long time to run, but doesn't use up much memory
  initabunds <- writeKCH(Samplenum=40,Filenum=1,density=masterDF$DENSITY[f])
  
     #### extract the main population arrays (to be edited)
  poplist1 <- template$mp.file$PopList
  poplist2 <- template$mp.file$PopData_df 

  poplist2$KchangeSt <- as.character(poplist2$KchangeSt)   

  poplist2$UD1 <- round(masterDF$HUNT[f],3)
  poplist2$UD2 <- round(masterDF$HUNTDD[f],3)
  poplist2$UD3 <- 0
  
  i=1
  for(i in 1:NPOPS){

              # SET GROWTH RATE
    poplist1[[i]]$MaxR <- round(masterDF$RMAX[f],3)   
    poplist2$MaxR[i] <-  round(masterDF$RMAX[f],3)  

              # SET KCH FILE
    KCHfilename <- sprintf("pop_%s.kch",i)
    poplist1[[i]]$KchangeSt <- KCHfilename
    poplist2$KchangeSt[i] <-  KCHfilename
 
              # SET INIT ABUND
    poplist1[[i]]$InitAbund <- round(initabunds[i])
    poplist2$InitAbund[i] <- round(initabunds[i])
  
              # SET ALLEE EFFECT
    poplist1[[i]]$localthr <- round(masterDF$ALLEE[f])
    poplist2$localthr[i] <- round(masterDF$ALLEE[f])
    
    poplist1[[i]]$IncludeInSum <- TRUE
    poplist2$IncludeInSum[i] <- TRUE
    
              # SET INIT HARVEST
    mpstarttime <- newmp$MaxDur*GENTIME
    harvstarttime <- humanArrival.df$year_min[i] + masterDF$HUMAN[f]*(humanArrival.df$year_max[i]-humanArrival.df$year_min[i])
    poplist2$UD3[i] <- round((mpstarttime - harvstarttime)/GENTIME)

  }
  

        # start the new mp file
  version <- template$version
  newmp <- template$mp.file
  newmp$PopList <- poplist1
  newmp$PopData_df <- poplist2

  newmp$UseDispDistFunc <- FALSE #TRUE    # now we need to use a dispersal matrix, unfortunately...  
  newmp$UseCorrDistFunc <- TRUE
  
          # make sure env stoch is lognormal
  newmp$Distrib <- "Lognormal"
  
          # set the Dll
  newmp$DDforAllPop <- "UD"   # set to user-defined
  newmp$UserDllFileName <- sprintf("%s\\%s",EXE_DIRECTORY,DLL_FILENAME)
  
          # SET DISPERSAL
  params <- template$mp.file$DispDistFunc
  params[1] <- masterDF$DISP1[f]/2                                # "a" parameter describes the proportion staying/leaving a given population
  params[2] <- dispersalFunc.df$b[masterDF$dispndx[f]]            # "b" parameter describes the distance moved
  params[4] <- dispersalFunc.df$Dmax[masterDF$dispndx[f]]         # "Dmax" parameter limits the maximum distance moved.
  newmp$DispDistFunc <- round(params,3)

  # ensure that max dist is 95% quantile of dispersal kernel
  # compute the dispersal at each distance bin. 
  # re-weight so that only X% leaves the focal population
  # fill in the dispersal kernel.
  
  # b parameter and a parameter must be set so that 95% dispersal quantile at Dmax AND exactly Disp1 % of individs staying in the focal population
  #  instead, just use the params from above, and spread out as necessary

  dispForDist <- DispersalKernelFunc(DistClasses[-1],params)
  dispMatrix <- matrix(0,nrow=NPOPS,ncol=NPOPS) 
  dispForDist <- dispForDist[-which(dispForDist==0)]  
  i=1
  for(i in 1:NPOPS){
    allIndices <- numeric(0) 
    dispRates <- numeric(0)
    j=1
    for(j in 1:length(dispForDist)){
      if(!is.na(DistBins[[i]][[j]][1])){ 
        allIndices <- c(allIndices,DistBins[[i]][[j]])
        dispRates <- c( dispRates,rep(dispForDist[j],times=length(DistBins[[i]][[j]])) )  
      }
    }
    multiplier = params[1]/sum(dispRates)
    dispRates = dispRates * multiplier   # ensure that total # leaving the population is held constant.
       #  sum(dispRates)
    dispMatrix[,i][allIndices] <- dispRates     
  }
  newmp$DispMatr <- round(dispMatrix,3)    # fill in the dispersal matrix.
  

         # SET MANAGEMENT ACTIONS (HARVEST)  NOTE: this has changed!!! Now in DLL
  
  newmp$NPopManage <- NPOPS                                       # num of management (harvest) actions equals the number of populations

  newmp$NPopManage <- 0
  
  # harvestdf <- newmp$PopManageProp                   
  # mpstarttime <- newmp$MaxDur*GENTIME                             # for translating between human arrival (ybp) and generations in the model
  # 
  # i=1
  # for(i in 1:NPOPS){
  #   if(i!=1) harvestdf <- rbind(harvestdf,newmp$PopManageProp) 
  #   harvstarttime <- humanArrival.df$year_min[i] + masterDF$HUMAN[f]*(humanArrival.df$year_max[i]-humanArrival.df$year_min[i])
  #   harvestdf$begin.time[i] <- round((mpstarttime - harvstarttime)/GENTIME)   # set start time for harvest
  #   harvestdf$end.time[i] <- newmp$MaxDur
  #   harvestdf$from.pop[i] <- i                                                # set population ID for management actions
  #   harvestdf$to.pop[i] <- i
  #   harvestdf$proportion <- round(masterDF$HUNT[f],5)                         # set proportion harvested
  # }
  # 
  # newmp$PopManageProp <- harvestdf                                # append all harvest-related parameters
  # 
         # SET INIT ABUND (again)
  newmp$StInit[,1] <- round(initabunds)                           # not sure which one is used, so set both at pop and global levels.

        # SET TEMP VAR
  newmp$SDMatr[[1]]$Matr[1,1] <- round(masterDF$SD[f],3)                    # set temporal variability

  # time steps
  newmp$MaxDur <- TIMESTEPS   
  
   # replicates
  newmp$MaxRep <- 1
  
   # correlation (change from default, which is in km)
  newmp$CorrDistFunc[2] <- 8
  
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





