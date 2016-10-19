
###############
#### UTILITY FUNCTIONS FOR PALEO MODELS
###############




###########
## FUNCTION "MakeLHSSamples"
##
## Samples from the uniform LHS and translates into desired parameter space
##  Returns a master data frame that will specify all MP files to be run
###########


###############################
#        SPECIFY PARAMETER RANGES
################################

MakeLHSSamples <- function(nicheBreadthDir,NicheBreadth){
  
  LHSParms <- list()    # initialize the container for parameter bounds
  
  setwd(nicheBreadthDir)   
  nCUTS <- length(grep("file_",list.dirs(path=nicheBreadthDir,recursive=F)))   #number of "cuts"  pattern="File_" 	
  LHSParms <- specifyLHSParam(paramslist=LHSParms,name="CUTS",type="CAT",lb=1,ub=nCUTS)
  
  ####  RMAX
  LHSParms <- specifyLHSParam(LHSParms,"RMAX",type="CONT",lb=1.28,ub=6.84)
  
  #### TEMPORAL VARIABILITY
  LHSParms <- specifyLHSParam(LHSParms,"SD",type="CONT",lb=0,ub=0.175)
  
  ### ALLEE EFFECT    0 to 100 individuals...  
  LHSParms <- specifyLHSParam(LHSParms,"ALLEE",type="CONT",lb=0,ub=500)       # KTS: changed to 500
  
  ### DENSITY PER 100kmX100km GRID CELL
  LHSParms <- specifyLHSParam(LHSParms,"DENSITY",type="CONT",lb=625,ub=10000)     # was 1000 to 40000  
  
  ### DISPERSAL #1: PROPORTION DISPERSING (used to compute "a" parameter.  a = p/2)
  LHSParms <- specifyLHSParam(LHSParms,"DISP1",type="CONT",lb=0.05,ub=0.25)
  
  ### DISPERSAL #2: MAXIMUM DISTANCE MOVED  [will translate into b and DMax according to Resit's pre-calculations]
  LHSParms <- specifyLHSParam(LHSParms,"DISP2",type="CONT",lb=100,ub=500) 
  
  ### HUNTING  (proportion harvested per year)   # NOTE: this is specified on an annual time scale (not generational). 
  LHSParms <- specifyLHSParam(LHSParms,"HARV",type="CONT",lb=0,ub=0.05)
  
  ### DENSITY DEPENDENCE ON HARVEST (y intercept of the harvest rate/abundance relationship)
  
  #  # HUNTDD_posvals <- c(seq(-0.05,0,length=10),seq(0.1,1,length=10),seq(1.1,2,length=10))
  # LHSParms <- specifyLHSParam(LHSParms,"HUNTDD",type="CONT",lb=-0.05,ub=2)
  
  LHSParms <- specifyLHSParam(LHSParms,"HARVZ",type="CONT",lb=1,ub=2)
  
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
  
  #masterDF$HUNTDD <- HUNTDD_posvals[masterDF$HUNTDD]    # value input into RAMAS DLL for density dependent hunting 

  masterDF$NicheBreadth <- NicheBreadth   # set the niche breadth
  
  #######  Add the MP filename to the masterDF data frame
  
  masterDF$MPFilename <- ""
  for(i in 1:nrow(masterDF)){
    masterDF$MPFilename[i] <- sprintf("NicheBreadth%s_LHS_Sample%s.mp",NicheBreadth,i)
  }
  
  setwd(MP_DIRECTORY)
  ## name file for LHS parameters 
  write.csv(masterDF,sprintf("masterDF_prelim%s.csv",NicheBreadth),row.names=F)
  
  return(masterDF)
}

###########
## FUNCTION "ReadMPTemplate"
##
## Reads in a template MP file for later use...
###########

ReadMPTemplate <- function(){
  #### read in template MP file
  setwd(MP_DIRECTORY)
  mpfilein <- MP_TEMPLATE
  template <- mp.read(mpfilein)   # read in template
  return(template)
}

###########
## FUNCTION "DispersalPreCalculations"
##
## Loads all required packages for paleo project
###########


#############################
# DISPERSAL PRE-CALCULATIONS (only needs to be done once)
#############################


DispersalPreCalculations <- function(){
  #?distMeeus
  #?distm
  
  #distMeeus(GridCellAttributes[1,c(1,2)],GridCellAttributes[2,c(1,2)])/1000
  
  # NPOPS*NPOPS matrix of pairwise distances between grid centroids (km) [uses geosphere package]
  distmatrix <- distm(GridCellAttributes[,c(1,2)],GridCellAttributes[,c(1,2)],fun=distMeeus)/1000
  
  # Distance classes for computing dispersal matrix (in km)
  MinDist <- 100
  MaxDist <- 500
  BinSize <- 20
  DistClasses <<- c(1,seq(MinDist,MaxDist,BinSize))
  
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

  return(DistBins)
}

###########
## FUNCTION "LoadPackages"
##
## Loads all required packages for paleo project
###########

LoadPackages <- function(){
  #############################
  #  LOAD PACKAGES
  #############################
  
  loadPackage("fields")     # note: loadPackage should install the package from CRAN automatically if it is not already installed
  loadPackage("lhs")
  loadPackage("RCurl")
  loadPackage("geosphere")
  loadPackage("parallel")
  loadPackage("doParallel")
  loadPackage("adehabitatHR")
  
  # Read in the functions for reading and writing MP files, from Matt Lammens
  
  baseurl = "https://raw.githubusercontent.com/mlammens/SACode_Sandbox/master/"				 
  source_github(baseurl,"mp.write.r")
  source_github(baseurl,"mp.read.r")
  source_github(baseurl,"metapopversion.r")
  source_github(baseurl,"fill.matrix.df.r")
  source_github(baseurl,"mp.read.results.r")

}


###########
## FUNCTION "LoadData"
##
## Loads all required data for paleo project
###########

LoadData <- function(){
  ###############################
  #        LOAD DATA AND OTHER RELEVANT FILES
  ###############################
  
  setwd(DATA_DIRECTORY)
  
  ##########
  # Human arrival dates for each population
  #humanArrival.df <- read.csv("Human_arrival_7000cells.csv",header=T)  # for each grid cell (population), earliest and latest human arrival dates (years before present))
  humanArrival.df <<- read.csv("Mammoth&Human_estimates_v3.csv",header=T)
  humanArrival.df <<- humanArrival.df[order(humanArrival.df$Cell.ID),]    # make sure it's ordered the same as the RAMAS model
  
  humanArrival.df$year_min <<- humanArrival.df$Mean.timing.of.human.arrival-2*humanArrival.df$Sd.timing.of.Human.arrival
  humanArrival.df$year_max <<- humanArrival.df$Mean.timing.of.human.arrival+2*humanArrival.df$Sd.timing.of.Human.arrival
  
  ##########
  # For computing real dispersal params 
  dispersalFunc.df <<- read.csv("dispersalfunctions.csv",header=T)      # links LHS samples to b and Dmax params
  
  #########
  # Filenames of all HS files (precursors of KCH files)
  #allHSfiles <- list.files(path=HS_DIRECTORY,pattern=".hs")            # char vector of all HS files (not needed??)
  
  #########
  # Area and coordinates of all grid cells
  #GridCellAttributes <- read.csv("Lattice100km_7000Cells.csv",header=T)  studysite_xy_modified2.csv
  GridCellAttributes <<- read.csv("studysite_xy_modified2.csv",header=T)  
  
  #GridCellAttributes$RamasID <- as.numeric(gsub("[^0-9]", "", unlist(GridCellAttributes$ID)))
  GridCellAttributes$RamasID <<- GridCellAttributes$ID
  GridCellAttributes <<- GridCellAttributes[order(GridCellAttributes$RamasID),]    # make sure it's ordered the same as the RAMAS model
  
  NPOPS <<- nrow(GridCellAttributes)  # ensure that this global parameter is correct
  
  ##########
  # Compute area of each grid cell in km2
  GridCellAttributes$Area2 <<- 0
  df <- data.frame(x=numeric(4),y=0)
  i=1
  for(i in 1:nrow(GridCellAttributes)){	
    xval <- GridCellAttributes$x.cord[i]
    yval <- GridCellAttributes$y.cord[i]
    xvlb <- xval - 0.5
    xvub <- xval + 0.5
    yvlb <- yval - 0.5
    yvub <- yval + 0.5
    df$x <- c(xvlb,xvub,xvub,xvlb)
    df$y <- c(yvub,yvub,yvlb,yvlb)
    #df <- SpatialPoints(df,proj4string=proj)
    GridCellAttributes$Area2[i] <<- areaPolygon(df)/1e6    # area of the MCP, in km2
  }
  
  # write to file, with updated areas...
  setwd(DATA_DIRECTORY)
  write.csv(GridCellAttributes,"GridCellAttributes.csv")
  
}

###########
## FUNCTION "SetUpWorkspace"
##
## Sets up the R workspace with necessary global variables
###########

SetUpWorkspace <- function(){
  if(KEVIN) BASE_DIRECTORY <<- "C:\\Users\\Kevin\\Dropbox\\Damien Fordham\\Mammoth Model"  
  if(DAMIEN) BASE_DIRECTORY <<- "E:\\Mammoth model"
  if(HRA_LAB) BASE_DIRECTORY <<- "C:\\Users\\Akcakaya\\Desktop\\Mammoth Model"
  if(KEVIN_LAB) BASE_DIRECTORY <<- "E:\\Dropbox\\Damien Fordham\\Mammoth Model"
  
  
  ######## BELOW DOESNT NEED TO BE CHANGED
  DATA_DIRECTORY <<- paste(BASE_DIRECTORY,"\\data",sep="")                         # directory for storing relevant data (CSV files)
  if(is.na(file.info(DATA_DIRECTORY)[1,"isdir"])) dir.create(DATA_DIRECTORY)
  
  MP_DIRECTORY <<- paste(BASE_DIRECTORY,"\\mpfiles",sep="")                        # directory for storing and running MP files 
  if(is.na(file.info(MP_DIRECTORY)[1,"isdir"])) dir.create(MP_DIRECTORY)
  
  EXE_DIRECTORY <<- paste(BASE_DIRECTORY,"\\EXEs",sep="")                        # directory for storing DLLs and EXE files 
  if(is.na(file.info(EXE_DIRECTORY)[1,"isdir"])) dir.create(EXE_DIRECTORY)
  
  FIGURES_DIRECTORY <<- paste(BASE_DIRECTORY,"\\figures\\raw",sep="")              # directory for storing raw figures 
  if(is.na(file.info(FIGURES_DIRECTORY)[1,"isdir"])) dir.create(FIGURES_DIRECTORY)
  
  HS_DIRECTORY <<- paste(BASE_DIRECTORY,"\\hs",sep="")                             # directory for storing habitat suitability files (used for creating KCH files )
  if(is.na(file.info(HS_DIRECTORY)[1,"isdir"])) dir.create(HS_DIRECTORY)
  
  KCH_DIRECTORY <<- paste(BASE_DIRECTORY,"\\KCHFiles",sep="")                            # directory for storing habitat suitability files
  if(is.na(file.info(KCH_DIRECTORY)[1,"isdir"])) dir.create(KCH_DIRECTORY)
  
  RESULTS_DIRECTORY <<- paste(BASE_DIRECTORY,"\\results",sep="")                             # directory for storing habitat suitability files (used for creating KCH files )
  if(is.na(file.info(RESULTS_DIRECTORY)[1,"isdir"])) dir.create(RESULTS_DIRECTORY)
  
  ############
  ## GLOBAL PARAMETERS

  TIMESTEPS <<- 3201 # set according to length of the simulation 
  GENTIME <<- 25
  MP_TEMPLATE <<- "NewMammothTemplate.mp"      # template MP file
  DLL_FILENAME <<- "Mammoth2.dll"    # KTS: changed to "Mammoth2.dll" from "Mammoth.dll"
  VERBOSE=FALSE
  
  # save global params
  setwd(DATA_DIRECTORY)
  save(NREPS,TIMESTEPS,GENTIME,MP_TEMPLATE,file="GlobalParams.RData")
  
  LoadPackages()  # load all packages
  LoadData()      # load all data
  
}  


###########
## FUNCTION "LoadPackages"
##
## Basic dispersal function in Ramas...
###########

negexp <- function(distance,alpha, beta) {
   alpha*exp((-1*(distance))/beta)   
}


###########
## FUNCTION "LoadPackages"
##
## Visualize basic dispersal function in Ramas...
###########

visualizeDispersalFunction <- function(alpha, beta) {
 graphics.off()
 a = alpha      # dispersal parameters
 b = beta      
 curve(negexp(x,a,b),col="red",lwd=2,0,7)
}

###########
## FUNCTION "specifyLHSParam"
##
## Information necessary to translate standard uniform LHS sample into parameters of interest for paleo project 
###########

specifyLHSParam <- function(paramslist,name,type,lb,ub){
    newlist <- paramslist
	eval(parse(text=sprintf("newlist$%s <- list()",name)))
	eval(parse(text=sprintf("newlist$%s$type <- \"%s\"",name,type)))
	eval(parse(text=sprintf("newlist$%s$lb <- %s",name,lb)))
	eval(parse(text=sprintf("newlist$%s$ub <- %s",name,ub))) 	
	return(newlist)
}




###########
## FUNCTION "writeKCH"
##
## Selects raw KCH files that Damien built (0 to 1), translates to real abundance, c
## computes initial abundance, and moves KCH files into the necessary folder for processing (same folder where mp file is located)
###########

writeKCH <- function(toFolder=thisFolder,NicheBreadth=40,Filenum = 1,density = 1000){
	initabund <- numeric(NPOPS)
	directoryname = sprintf("%s\\Sample_%s\\file_%s",KCH_DIRECTORY,NicheBreadth,Filenum)
	setwd(directoryname)
	i=1
	for(i in 1:NPOPS){
	  KCHfilename = sprintf("pop_%s.kch",i)
	  setwd(directoryname)
	  K = read.table(KCHfilename,header=F)[,1]
	  K=K*density   # convert to abundance
	  initabund[i] = K[1]
	  setwd(toFolder)
	  write.table(K,file=KCHfilename,row.names=F,col.names=F)
	}
	return(initabund)
}

###########
## FUNCTION "loadPackage"
##
## GENERIC FUNCTION FOR INSTALLING/LOADING PACKAGES FROM CRAN
##
###########

loadPackage <- function(pkg){

  if(pkg %in% rownames(installed.packages()) == FALSE) {suppressMessages(suppressWarnings(install.packages(pkg)))}
  eval(parse(text=sprintf("suppressMessages(suppressWarnings(require(%s)))",pkg)), envir= .GlobalEnv)

}

###########
## FUNCTION "loadPackage"
##
## GENERIC FUNCTION FOR INSTALLING/LOADING SOURCE CODE FROM GITHUB
##
###########

source_github <- function(baseurl,scriptname) {
  # load package
   suppressMessages(suppressWarnings(require(RCurl)))
 
  # read script lines from website
  url <- sprintf("%s%s",baseurl,scriptname)
  script <- getURL(url, ssl.verifypeer = FALSE)
  
  script <- gsub("\r\n", "\n", script)     # get rid of carriage returns (not sure why this is necessary...)
 
  # parse lines and evaluate in the global environement
  eval(parse(text = script), envir= .GlobalEnv)
}


###########
## FUNCTION "DispersalKernelFunc"
##
## FUNCTION FOR CONVERTING DISTANCES TO DISPERSAL RATES USING RAMAS DISPERSAL PARAMS
##
###########

 DispersalKernelFunc <- function(dist,params){   
    ifelse(dist<=params[4],params[1]*exp((-dist^params[3])/params[2]),0)
  }
 
#########################
 ### DEPRECATED FUNCTIONS
#########################
 
 #### function to set up and write KCH files
 #      and... determine initial abundance for each population
 #      remember- we will need to scale by area... 
 
 # writeKCH <- function(HSnum = 1,density = 1000,areas=rep(1000,times=NPOPS)){
 # 	initabund <- numeric(NPOPS)
 # 	for(i in 1:NPOPS){
 # 	  HSfilename = sprintf("pop%s_rep%s.hs",i,HSnum)
 # 	  KCHfilename = sprintf("pop%s_rep%s.kch",i,HSnum)
 # 	  setwd(HS_DIRECTORY)
 # 	  hs = read.table(HSfilename,header=F)[,1]
 # 	  initabund[i] = hs[1]*density*(areas[i]/10000)
 # 	  K = round(hs*density*(areas[i]/10000))
 # 	  setwd(MP_DIRECTORY)
 # 	  write.table(K,file=KCHfilename,row.names=F,col.names=F)
 # 	}
 # 	return(initabund)
 # }




