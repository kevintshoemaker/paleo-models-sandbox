##########################
# This script extracts results from a single MP file generated for the paleo project 
#  Authors: Kevin Shoemaker and Damien Fordham
#
# Summary: 
#  3 June 2015 -- started scripting
#  15 June 2015 -- tested script, loaded to GitHub
#  1 July 2016  -- modified for parallelization

#TODO  
##########################


#########################################################################
#  EXTRACT RESULTS FROM MP FILE
#########################################################################
# Response variables to record
# Grid cell abundance over time    
# Timing of grid cell extinction    
# Population size over time              
# Time of population level extinction
# Sum cell occupancy over time
# MCP over time (maybe save for later)


ExtractMPresults <- function(f=i,masterDF=masterDF,NicheBreadth=NicheBreadth){
  
  ## set up the new folder to store the MP file and associated KCH files... (specifies the niche breadth)
  thisFolder <- sprintf("%s\\Sample_%s\\LHS_Sample%s",MP_DIRECTORY,NicheBreadth,f)
  
  ############
  # SET UP STORAGE STRUCTURES FOR KEY RESULTS METRICS
  ############
  SimInfo <- masterDF[f,]
  SimInfo <- as.list(SimInfo)   # coerce to list 
  MPFilename <- as.character(SimInfo$MPFilename)
  name <- gsub(x=MPFilename,pattern=".mp",replacement="",fixed=T)
  
  
  ## add field (matrix) for storing abundance over time for each grid cell
  SimInfo$PopAbund <- matrix(0,nrow=NPOPS,ncol=TIMESTEPS)
  
  ## add field (vector) for storing the final occupied year for each grid cell
  SimInfo$FinalYear <- numeric(NPOPS)
  
  ## add field (vector) for storing global abundance over time
  SimInfo$TotAbund <- numeric(TIMESTEPS)
  
  ## add field (scalar) for storing the year of global extinction
  SimInfo$ExtinctionYear <- 0
  
  ## add field (list) for storing MCP for each year
  #SimInfo$MCPs <- list() 
  
  ## add field (vector) for storing total range area (MCP) over time
  #SimInfo$RangeArea <- numeric(TIMESTEPS)
  
  ## add field (scalar) for storing occupancy (number of occupied cells) over time
  #SimInfo$CellsOccupied <- numeric(TIMESTEPS)
  
  ## add field (scalar) for storing occupancy (number of occupied cells) over time
  #SimInfo$AreaOccupied <- numeric(TIMESTEPS)     
  
  ###########################
  # GET RESULTS
  ###########################
  
  ### step 1: set up the MP file connection 
  setwd(thisFolder)
  MPcon <- file(MPFilename, 'r')
  
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
  SimInfo$ExtinctionYear <- NA
  if(SimInfo$TotAbund[TIMESTEPS]<1) SimInfo$ExtinctionYear <- min(which(SimInfo$TotAbund==0))
  
  # RESULT: FINAL OCCUPIED CELL(s)
  SimInfo$FinalOccCell <- NA
  if(SimInfo$TotAbund[TIMESTEPS]<1){ 
    SimInfo$FinalOccCell <- which(SimInfo$PopAbund[,(SimInfo$ExtinctionYear-1)]>0)
  }else{
    SimInfo$FinalOccCell <- NA #which(SimInfo$PopAbund[,(TIMESTEPS-1)]>0)
  }
  
  
  # RESULT: FINAL OCCUPIED YEAR FOR EACH POPULATION
  SimInfo$FinalYear <- apply(SimInfo$PopAbund,1,function(t) ifelse(sum(t)>0,max(which(t>0)),NA))
  
  SimInfo$PopAbund <- NULL
  
  ####################
  # SAVE RESULTS TO HARD DISK AND REMOVE FROM RAM 
  ####################
  setwd(thisFolder)
  filename <- sprintf("%s.RData",name)
  eval(parse(text=sprintf("%s <- SimInfo",name)))
  eval(parse(text=sprintf("save(%s,file=filename)",name)))   # save to disk
  eval(parse(text=sprintf("rm(%s)",name)))   # remove from memory   
  
  #}  # end loop through files...
  
  
}  # end function 'ExtractMPresults'   




