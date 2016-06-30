##########################
# This script writes all the MP files for the wooly mammoth model
#  Authors: Kevin Shoemaker and Damien Fordham
#
# Summary: ##### NOW: Run the mp files, record the predictor variables
#     and response variables. 
 
#  3 June 2015 -- started scripting
#  15 June 2015 -- tested script, loaded to GitHub

#TODO  
   #* deal with burnin period?

######################
#   CLEAR THE WORKSPACE
######################

rm(list=ls())


#########################
# SET PROJECT DIRECTORIES
#########################

KEVIN = FALSE  #TRUE
DAMIEN = FALSE
HRA_LAB = TRUE

if(KEVIN) BASE_DIRECTORY <- "C:\\Users\\Kevin\\Dropbox\\Damien Fordham\\Mammoth Model"  
if(DAMIEN) BASE_DIRECTORY <- "E:\\Mammoth"

if(HRA_LAB) BASE_DIRECTORY <- "C:\\Users\\Akcakaya\\Desktop\\Mammoth Model"

if(KEVIN) CODE_DIRECTORY <- "C:\\Users\\Kevin\\GIT\\paleo-models-sandbox"       # code directory should be your local copy of the GitHub repository   
if(HRA_LAB) CODE_DIRECTORY <- "C:\\Users\\Akcakaya\\Desktop\\Mammoth Model\\paleo-models-sandbox"


######## BELOW DOESNT NEED TO BE CHANGED
DATA_DIRECTORY <- paste(BASE_DIRECTORY,"\\data",sep="")                         # directory for storing relevant data (CSV files)
if(is.na(file.info(DATA_DIRECTORY)[1,"isdir"])) dir.create(DATA_DIRECTORY)

MP_DIRECTORY <- paste(BASE_DIRECTORY,"\\mpfiles",sep="")                        # directory for storing and running MP files 
if(is.na(file.info(MP_DIRECTORY)[1,"isdir"])) dir.create(MP_DIRECTORY)

EXE_DIRECTORY <- paste(BASE_DIRECTORY,"\\EXEs",sep="")                        # directory for storing DLLs and EXE files 
if(is.na(file.info(EXE_DIRECTORY)[1,"isdir"])) dir.create(EXE_DIRECTORY)

FIGURES_DIRECTORY <- paste(BASE_DIRECTORY,"\\figures\\raw",sep="")              # directory for storing raw figures 
if(is.na(file.info(FIGURES_DIRECTORY)[1,"isdir"])) dir.create(FIGURES_DIRECTORY)

HS_DIRECTORY <- paste(BASE_DIRECTORY,"\\hs",sep="")                             # directory for storing habitat suitability files (used for creating KCH files )
if(is.na(file.info(HS_DIRECTORY)[1,"isdir"])) dir.create(HS_DIRECTORY)

# TODO: automate the KCH selection...
KCH_DIRECTORY <- paste(BASE_DIRECTORY,"\\KCHFiles",sep="") #\\Sample_40\\file_1"                             # directory for storing habitat suitability files
if(is.na(file.info(KCH_DIRECTORY)[1,"isdir"])) dir.create(KCH_DIRECTORY)

RESULTS_DIRECTORY <- paste(BASE_DIRECTORY,"\\results",sep="")                             # directory for storing habitat suitability files (used for creating KCH files )
if(is.na(file.info(RESULTS_DIRECTORY)[1,"isdir"])) dir.create(RESULTS_DIRECTORY)

######################
# LOAD GLOBAL PARAMS
######################

setwd(MP_DIRECTORY)
load("GlobalParams.RData")

####################
#  LOAD PROJECT-SPECIFIC FUNCTIONS
#################### 

setwd(CODE_DIRECTORY)
source("Rfunctions_PaleoModels.r")

#############################
#  LOAD DATA
#############################

	#########
	# Area and coordinates of all grid cells
	
setwd(DATA_DIRECTORY)
GridCellAttributes <- read.csv("GridCellAttributes.csv",header=T)

 #GridCellAttributes$RamasID <- as.numeric(gsub("[^0-9]", "", unlist(GridCellAttributes$ID)))
GridCellAttributes <- GridCellAttributes[order(GridCellAttributes$RamasID),]    # make sure it's ordered the same as the RAMAS model


#############################
#  LOAD PACKAGES
#############################

loadPackage("fields")     # note: loadPackage should install the package from CRAN automatically if it is not already installed
loadPackage("lhs")
loadPackage("RCurl")
loadPackage("geosphere")
loadPackage("adehabitatHR")

                 # Read in the functions for reading and writing MP files, from Matt Lammens

baseurl = "https://raw.githubusercontent.com/mlammens/SACode_Sandbox/master/"				 
source_github(baseurl,"mp.write.r")
source_github(baseurl,"mp.read.r")
source_github(baseurl,"mp.read.results.r")
source_github(baseurl,"metapopversion.r")
source_github(baseurl,"fill.matrix.df.r")


###########################################
#  READ MODEL SPECS/METADATA
###########################################

     ###   read in the LHS draws/MP file specs
setwd(MP_DIRECTORY)
masterDF <- read.csv("masterDF.csv",header=T)
nfiles <- nrow(masterDF)

head(masterDF)

MPsToRun <- as.character(masterDF$filenames)    # list of all MP files to run

#########################################################################
#  RUN ALL MP FILES AND STORE/SUMMARIZE THE RESULTS
#########################################################################
# Response variables to record
	# Grid cell abundance over time    
	# Timing of grid cell extinction    
	# Population size over time              
	# Time of population level extinction
	# Sum cell occupancy over time
	# MCP over time

setwd(MP_DIRECTORY)
f=1
for(f in 1:nfiles){     # loop through MP files...

### FOR TESTING
#MPcon <- file("Mammoth_100_testing.mp",'r')
#NPOPS <- 100

	############
	# SET UP STORAGE STRUCTURES FOR KEY RESULTS METRICS
	############
	SimInfo <- masterDF[f,]
	name <- gsub(x=MPsToRun[f],pattern=".mp",replacement="",fixed=T)
    
	SimInfo <- as.list(SimInfo)   # coerce to list         

	  ## add field (matrix) for storing abundance over time for each grid cell
	SimInfo$PopAbund <- matrix(0,nrow=NPOPS,ncol=TIMESTEPS)

	  ## add field (vector) for storing the final occupied year for each grid cell
	SimInfo$FinalYear <- numeric(NPOPS)
	
	  ## add field (vector) for storing global abundance over time
	SimInfo$TotAbund <- numeric(TIMESTEPS)
	
	  ## add field (scalar) for storing the year of global extinction
	SimInfo$ExtinctionYear <- 0

	  ## add field (list) for storing MCP for each year
	SimInfo$MCPs <- list() 

	  ## add field (vector) for storing total range area (MCP) over time
	SimInfo$RangeArea <- numeric(TIMESTEPS)
	
	   ## add field (scalar) for storing occupancy (number of occupied cells) over time
	SimInfo$CellsOccupied <- numeric(TIMESTEPS)
	
	   ## add field (scalar) for storing occupancy (number of occupied cells) over time
	SimInfo$AreaOccupied <- numeric(TIMESTEPS)     

	###########################
	# RUN MODEL (takes some time!)
	###########################

	command <- sprintf("START /WAIT \"MammothModelTest\" \"C:\\Program Files (x86)\\RAMASGIS\\Metapop.exe\" \"%s\" /RUN=YES",MPsToRun[f])
	# cat(command)

	setwd(MP_DIRECTORY)
	write.table(command,file="Temp.bat",sep="",quote=FALSE,row.names=FALSE,col.names=FALSE)
	system("Temp.bat")    # invoke the batch file to run the MP file...

	if(VERBOSE) cat(paste("just finished running file number",f,"\n",sep=" "))

	###########################
	# GET RESULTS
	###########################

	### step 1: set up the MP file connection 

	MPcon <- file(MPsToRun[f], 'r')

		# while loop: find the string "Pop. 1" [indicates the beginning of the population results]
	stringToFind <- "Pop. 1"
	basendx <- 0
	CHUNKSIZE <- 1000				      
	while (length(input <- readLines(MPcon, n=CHUNKSIZE)) > 0){    # read in chunks until population results are found

		temp <- grep(stringToFind,input)
		if(length(temp)>0){ 
			ndx <- basendx + temp
			pushBack(input[(temp):CHUNKSIZE],MPcon)       # reset the file to where Pop. 1 began
			break
		} 
		basendx <- basendx + CHUNKSIZE

	}   # end while loop

	     # read in the population abundances over time
	for(pop in 1:NPOPS){
		stringToFind <- sprintf("Pop. %s",pop)
		temp <- readLines(MPcon,1)
		if(temp!=stringToFind){
			print(paste("ERROR!","Population #",pop))
			break
		}
		input <- readLines(MPcon, n=TIMESTEPS)
		Nvec <- sapply( strsplit(input, " "), function(t) as.numeric(t[1]))
		eval(parse(text=sprintf("SimInfo$PopAbund[%s,] <- Nvec",pop)))	   # RESULT: POP ABUNDANCE             
	}   # end loop through pops 

	closeAllConnections()

	# RESULT: TOTAL ABUNDANCE
	SimInfo$TotAbund <- apply(SimInfo$PopAbund,2,sum)
	
	# RESULT: GLOBAL EXTINCTION YEAR
      SimInfo$ExtinctionYear <- TIMESTEPS
	if(SimInfo$TotAbund[TIMESTEPS]<1) SimInfo$ExtinctionYear <- min(which(SimInfo$TotAbund==0))

	# RESULT: FINAL OCCUPIED YEAR FOR EACH POPULATION
	SimInfo$FinalYear <- apply(SimInfo$PopAbund,1,function(t) ifelse(sum(t)>0,max(which(t>0)),NA))
		
		# XY coords of all occupied sites
            occndx <- sapply(as.data.frame(SimInfo$PopAbund),function(t) which(t>1))  # indices of occupied populations for each year
            proj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
            for(t in 1:(SimInfo$ExtinctionYear-1)){
               ndx <- occndx[[t]]

               # RESULT: MCP OVER TIME

               if(length(ndx)>=5){       # At least 5 relocations are required to fit an home range   
                 df <- data.frame(x=numeric(length(ndx)),y=0)
                 df$x <- GridCellAttributes$x.cord[ndx]        # x and y coordinates for all occupied grid cells
                 df$y <- GridCellAttributes$y.cord[ndx]
                 df <- SpatialPoints(df,proj4string=proj)
                 MCP <- mcp(df)
                 MCPpoly <- MCP@polygons[[1]]@Polygons[[1]]@coords
                 MCParea <- areaPolygon(MCPpoly)/1e6    # area of the MCP, in km2
                 SimInfo$MCPs[[t]] <- MCP            # store the MCP as SpatialPolygonsDataFrame object
                 SimInfo$RangeArea[t] <- MCParea     # store MCP area for each year of the simulation
              }
              # RESULT: CELLS OCCUPIED OVER TIME
              SimInfo$CellsOccupied[t] <- length(ndx)

              # RESULT: OCCUPIED AREA OVER TIME AND CELLS OCCUPIED OVER TIME
              areaVec <-  GridCellAttributes$Area[ndx]
              SimInfo$AreaOccupied[t] <- sum(GridCellAttributes$Area2[ndx])

            }		
	
	

	####################
	# SAVE RESULTS TO HARD DISK AND REMOVE FROM RAM
	####################
	setwd(RESULTS_DIRECTORY)
	filename <- sprintf("%s.RData",name)
      eval(parse(text=sprintf("%s <- SimInfo",name)))
	eval(parse(text=sprintf("save(%s,file=filename)",name)))   # save to disk
	eval(parse(text=sprintf("rm(%s)",name)))   # remove from memory   
							
}  # end loop through files...



###############  END SCRIPT
############
########
####
#






