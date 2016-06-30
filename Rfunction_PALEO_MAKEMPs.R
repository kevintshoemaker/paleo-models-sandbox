##########################
#  FUNCTION MakeMPs
#
#  This function writes a single MP file for the wooly mammoth model
#  Authors: Kevin Shoemaker and Damien Fordham
# 
#  22 May 2015 -- started scripting
#  11 June 2015 -- added dispersal correction
#  30 June 2016 -- cleaned up for initiation of model runs 
##########################




#############################
#       GENERATE THE MP FILES: now parallelizable...
#############################

MakeMPfile <- function(f=1,masterDF=masterDF,NicheBreadth=40){
    
    ###########################################
             ## BUILD THE MP FILES
    # thisLine <- masterDF[] 
    # masterDF$filenames <- character(1)      
    
    # f=2
    # for(f in 1:NREPS){
    filename <- sprintf("NicheBreadth%s_LHS_Sample%s.mp",NicheBreadth,f)
    
       ## set up the new folder to store the MP file and associated KCH files... (specifies the niche breadth)
    thisFolder <- sprintf("%s\\Sample_%s\\LHS_Sample%s",MP_DIRECTORY,NicheBreadth,f)
    if(is.na(file.info(thisFolder)[1,"isdir"])) dir.create(thisFolder,recursive = T)
    
    # masterDF$filenames[f] <- filename 
    
       ### write out the KCH files for this model run (and get init abundances)   NOTE: KCH files are stored in the MP_DIRECTORY
       # initabunds <- writeKCH(HSnum = masterDF$HS[f], density = masterDF$DENSITY[f],areas=GridCellAttributes$Area2)     ### NOTE: takes a long time to run, but doesn't use up much memory
    initabunds <- writeKCH(toFolder=thisFolder,NicheBreadth=40,Filenum=masterDF$CUTS[f],density=masterDF$DENSITY[f])
    
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
      mpstarttime <- template$mp.file$MaxDur*GENTIME
      harvstarttime <- humanArrival.df$year_min[i] + masterDF$HUMAN[f]*(humanArrival.df$year_max[i]-humanArrival.df$year_min[i])
      poplist2$UD3[i] <- round((mpstarttime - harvstarttime)/GENTIME) + 3  # KTS: added "+3" to make human effect have a time lag of 100 years 

    }
    
  
          # start the new mp file
    version <- 60 #template$version
    newmp <- template$mp.file
    newmp$PopList <- poplist1
    newmp$PopData_df <- poplist2
  
    newmp$UseDispDistFunc <- FALSE #TRUE    # now we need to use a dispersal matrix, unfortunately...  
    newmp$UseCorrDistFunc <- TRUE
    
            # make sure env stoch is lognormal
    newmp$Distrib <- "Lognormal"
    
            # set the Dll
    newmp$DDforAllPop <- "UD"   # set to user-defined
    newmp$UserDllFileName <- sprintf("%s\\%s",EXE_DIRECTORY,DLL_FILENAME)   # TODO: change this? 
    
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
    newmp$CorrDistFunc[1] <- 0.99 
    
  ###############
        # write the new mp file
    setwd(thisFolder)                                               
    mp.write(newmp,version,filename)                   # write MP file
     
    #}     # close loops through scenarios
    
    
      ### return something...
    return(filename)
    
    #######################
    #######################

}


# ##########################
# ####           SAVE METADATA FOR EACH MP FILE TO FILE
# 
# setwd(MP_DIRECTORY)
# write.csv(masterDF,"masterDF.csv",row.names=F)
