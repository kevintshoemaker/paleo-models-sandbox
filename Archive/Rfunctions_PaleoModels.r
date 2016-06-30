
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
  
  # Read in the functions for reading and writing MP files, from Matt Lammens
  
  baseurl = "https://raw.githubusercontent.com/mlammens/SACode_Sandbox/master/"				 
  source_github(baseurl,"mp.write.r")
  source_github(baseurl,"mp.read.r")
  source_github(baseurl,"metapopversion.r")
  source_github(baseurl,"fill.matrix.df.r")
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
  GridCellAttributes$Area2 <- 0
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
    GridCellAttributes$Area2[i] <- areaPolygon(df)/1e6    # area of the MCP, in km2
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
  if(DAMIEN) BASE_DIRECTORY <<- "E:\\Mammoth"
  
  if(HRA_LAB) BASE_DIRECTORY <<- "C:\\Users\\Akcakaya\\Desktop\\Mammoth Model"
  
  if(KEVIN) CODE_DIRECTORY <<- "C:\\Users\\Kevin\\GIT\\paleo-models-sandbox"       # code directory should be your local copy of the GitHub repository   
  if(HRA_LAB) CODE_DIRECTORY <<- "C:\\Users\\Akcakaya\\Desktop\\Mammoth Model\\paleo-models-sandbox"
  
  
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
  
  # TODO: automate the KCH selection...
  KCH_DIRECTORY <<- paste(BASE_DIRECTORY,"\\KCHFiles",sep="") #\\Sample_40\\file_1"                             # directory for storing habitat suitability files
  if(is.na(file.info(KCH_DIRECTORY)[1,"isdir"])) dir.create(KCH_DIRECTORY)
  
  RESULTS_DIRECTORY <<- paste(BASE_DIRECTORY,"\\results",sep="")                             # directory for storing habitat suitability files (used for creating KCH files )
  if(is.na(file.info(RESULTS_DIRECTORY)[1,"isdir"])) dir.create(RESULTS_DIRECTORY)
  
  ############
  ## GLOBAL PARAMETERS

  TIMESTEPS <<- 3201 #2468 #1000
  #NPOPS <- 7000
  GENTIME <<- 25
  MP_TEMPLATE_OLD <<- "Mammoth_7000_demo.mp"      # template MP file
  DLL_FILENAME <<- "Mammoth.dll"
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

writeKCH <- function(Samplenum = 40,Filenum = 1,density = masterDF$DENSITY[f]){
	initabund <- numeric(NPOPS)
	directoryname = sprintf("%s\\Sample_%s\\file_%s",KCH_DIRECTORY,Samplenum,Filenum)
	setwd(directoryname)
	i=1
	for(i in 1:NPOPS){
	  KCHfilename = sprintf("pop_%s.kch",i)
	  setwd(directoryname)
	  K = read.table(KCHfilename,header=F)[,1]
	  K=K*density
	  initabund[i] = K[1]
	  setwd(MP_DIRECTORY)
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




