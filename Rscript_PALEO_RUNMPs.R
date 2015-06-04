##########################
# This script writes all the MP files for the wooly mammoth model
#  Authors: Kevin Shoemaker and Damien Fordham
#
# Summary: ##### NOW: Run the mp files, record the predictor variables
#     and response variables. 
 
#  3 June 2015 -- started scripting

#TODO  
   #* deal with burnin period?
   # test for run that goes extinct
   # for spatial analysis, need template raster of population names...
 

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
#  LOAD PACKAGES
#############################

loadPackage("fields")     # note: loadPackage should install the package from CRAN automatically if it is not already installed
loadPackage("lhs")
loadPackage("RCurl")

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

	
MPsToRun <- as.character(masterDF$filenames)


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

for(f in 1:nfiles){     # loop through MP files...

	############
	# SET UP STORAGE STRUCTURES FOR KEY RESULTS METRICS
	############
	temp <- masterDF[f,]
	name1 <- gsub(x=MPsToRun[f],pattern=".mp",replacement="",fixed=T)
	eval(parse(text=sprintf("%s <- as.list(temp)",name1)))         # coerce to list

	  ## add field (data frame) for storing abundance over time for each grid cell
	eval(parse(text= sprintf("%s$PopAbund <- matrix(0,nrow=NPOPS,ncol=TIMESTEPS)",name1) ))

	  ## add field (list) for storing the year(s) of extinction for each grid cell [not needed?]
	#eval(parse(text= sprintf("%s$PopExtinctionYears <- list()",name1) ))

	  ## add field (list) for storing the final occupied year for each grid cell
	eval(parse(text= sprintf("%s$FinalYear <- numeric(NPOPS)",name1) ))

	  ## add field (vector) for storing global abundance over time
	eval(parse(text= sprintf("%s$TotAbund <- numeric(TIMESTEPS)",name1) ))

	  ## add field (scalar) for storing the year of global extinction
	eval(parse(text= sprintf("%s$ExtinctionYear <- 0",name1) )) 

	   ## add field (scalar) for storing occupancy (number of occupied cells) over time
	eval(parse(text= sprintf("%s$CellsOccupied <- numeric(TIMESTEPS)",name1) ))     

	  ## add field (vector) for storing total range area (MCP) over time
	eval(parse(text= sprintf("%s$TotAbund <- numeric(TIMESTEPS)",name1) ))


	###########################
	# RUN MODEL (takes a while)
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

	## MP <- mp.read(MPsToRun[f])   # read in template     # NOTE: Matt's mp.read function can't handle these big files!


      # new strategy: find "Pop. 1", "Pop. 2" etc and read in chunks from there... 

	  ### step 1: find the 

      MPcon <- file(MPsToRun[f], 'r')

### FOR TESTING
MPcon <- file("Mammoth_100_testing.mp",'r')
NPOPS <- 100


					# while loop: find the string "Pop. 1"
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

      for(pop in 1:NPOPS){
            stringToFind <- sprintf("Pop. %s",pop)
            temp <- readLines(MPcon,1)
            if(temp!=stringToFind){
              print(paste("ERROR!","Population #",pop))
              break
            }
      	input <- readLines(MPcon, n=TIMESTEPS)
            Nvec <- sapply( strsplit(input, " "), function(t) as.numeric(t[1]))
            eval(parse(text=sprintf("%s$PopAbund[%s,] <- Nvec",name1,pop)))	   # RESULT: POP ABUNDANCE             

	}   # end loop through pops 

	closeAllConnections()

	# RESULT: TOTAL ABUNDANCE
	eval(parse(text= sprintf("%s$TotAbund <- apply(%s$PopAbund,2,sum)",name1,name1) ))

	# RESULT: GLOBAL EXTINCTION YEAR
	min(which(ta==0))
	eval(parse(text= sprintf("%s$ExtinctionYear <- min(which(%s$TotAbund==0))",name1,name1) ))

	# RESULT: LOCAL EXTINCTION YEAR(S)  [not really needed??]
	#min(which(ta==0))
      #for(pop in 1:NPOPS){
	#  eval(parse(text= sprintf("%s$PopExtinctionYears$Pop%s <- which(%s$PopAbund[pop,]==0)",name1,pop,name1) ))
	#}

	# RESULT: FINAL OCCUPIED YEAR FOR EACH POPULATION
	eval(parse(text= sprintf("%s$FinalYear <- apply(%s$PopAbund,1,function(t) max(which(t>0)))",name1,name1) ))




	####################
    	# SAVE RESULTS TO HARD DISK AND REMOVE FROM RAM
	####################
	setwd(RESULTS_DIRECTORY)
      filename <- sprintf("%s.RData",name1)
	eval(parse(text=sprintf("save(%s,file=filename)",name1)))   # save to disk
      eval(parse(text=sprintf("rm(%s)",name1)))   # remove from memory   

}  # end loop through files...



write.csv(allVars,"AllVars.csv",row.names=F)


###############
############
########
####
#




###############
######  Calculate mean inter-colony distance...

df <- read.csv("AllVars.csv",header=T)

names(df)

dispToDist <- function(disp){
  ifelse(disp>0.0001, (log(func[1])-log(disp))*func[2], 
      15)
} 


meannndist <- numeric(nrow(df))
meandisp <- numeric(nrow(df)) 
meanconn <- numeric(nrow(df))
for(i in 1:length(filenames)){
   temp <- mp.read(as.character(df$filenames[i]))   # read in the file 

   
   #  names(temp$mp.file)
   func <- temp$mp.file$DispDistFunc
   # curve(ifelse(x<10,func[1]*exp(-1*(x^func[3])/func[2]),0),0,15)
   # curve((log(func[1])-log(x))*func[2],0,0.1)
   DispMatr <- temp$mp.file$DispMatr
   
   temp_nnd <- apply(DispMatr,2,function(t) min(dispToDist(t))) 
   meannndist[i] <- mean(temp_nnd)

   meandisp[i] <- mean(apply(DispMatr,2,max))

   meanconn[i] <- mean(apply(DispMatr,2,function(t) length(which(t>=0.005))))
   closeAllConnections()
}

df$nn_dist <- meannndist
df$mean_disp <- meandisp
df$mean_conn <- meanconn

write.csv(df,"AllVars2.csv",row.names=F)



