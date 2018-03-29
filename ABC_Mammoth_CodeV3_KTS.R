##########################
# This master script performs Approximate Bayesian Analysis
#     algorithms and goodness-of-fit testing for the Paleo simulations
#         USES THE 'ABC' PACKAGE (for post-hoc elimination-style ABC)
#
#  Authors: Kevin Shoemaker and Damien Fordham
# 
#  22 May 2017 -- started scripting
#
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
#   PRELIMINARY: LOAD FUNCTIONS
#################### 

if(KEVIN) CODE_DIRECTORY <- "C:\\Users\\Kevin\\GIT\\paleo-models-sandbox"       # code directory should be your local copy of the GitHub repository   
if(HRA_LAB) CODE_DIRECTORY <- "C:\\Users\\Akcakaya\\Desktop\\Mammoth Model\\paleo-models-sandbox"
if(DAMIEN) CODE_DIRECTORY <-  "C:\\Users\\Damien Fordham\\Documents\\GitHub\\paleo-models-sandbox"
if(KEVIN_LAB) CODE_DIRECTORY <- "E:\\GIT\\paleo-models-sandbox"

setwd(CODE_DIRECTORY)
source("Rfunctions_PALEO_UTILITY.r")     # Load all functions for Paleo project

####################
#  PRELIMINARY: SET UP WORKSPACE AND LOAD PACKAGES
####################

SetUpWorkspace()   # function loads packages and sets up the workspace...  

####################
# LOAD ADDITIONAL PACKAGES
###################

##ABC code
library(abc)
library(raster)
#library(ggmap)
library(mapproj)
library(rworldmap)
library(abctools)    # for fancy calibration of the posterior - for later? 


####################
# LOAD FUNCTIONS
####################

# prepare data for ABC

prepareForABC <- function(data=all,subset=c(1,3),standardized=T){
  par.sim <<- data[,estimable.params]     # samples run from parameter space 
  if(standardized){ 
    maxstats <- apply(data[,test.statistics[subset],drop=FALSE],2,function(t) max(t,na.rm=T))
    for(i in test.statistics[subset]){
      data[,i] <- data[,i]/maxstats[i]
    }
    stat.obs <<- test.statistics.obs[test.statistics[subset]]/maxstats
  } else{
    stat.obs <<- test.statistics.obs[test.statistics[subset]]
  }
  stat.sim <<- data[,test.statistics[subset],drop=FALSE]
}


## plot priors for all parameters (to check for uniformity) 

visualize.priors <- function(){
  par(mfrow=c(2, 3),ask=TRUE)
  
  for(i in 1:length(estimable.params)){
    hist(all[,estimable.params[i]], breaks = 50,main=estimable.params[i])
  }
}


## note: this function only works with rejection method...
visualize.posteriors <- function(){
  par(mfrow=c(2, 3),ask=TRUE)
  posterior_ndx <- as.numeric(rownames(as.matrix(rej$ss))) 
  all_post <- all[posterior_ndx,]
  for(i in 1:length(estimable.params)){
    hist(all_post[,estimable.params[i]], breaks = 50,main=estimable.params[i],freq=FALSE)
    lines(density(all_post[,estimable.params[i]]),col="blue",lwd=2)
  }
}



visualize <- function(target=target){
  
  par(ask=TRUE)
  
  #if(toPlot=="ABUNDANCE"){
    matches <- as.numeric(unlist(regmatches(target, gregexpr("[[:digit:]]+", target))))
    nichebreadth <- matches[1]
    sample <- matches[2]
    
    
    filename <- sprintf("NicheBreadth%i_TotAbund_output_alldata.csv",nichebreadth)
    
    setwd(SIM_RESULTS_DIR)
    
    #getwd()
    #list.files()
    
    
    abund <- read.csv(filename, header = T,stringsAsFactors = FALSE)    # if(!"abund"%in%ls(name=.GlobalEnv))  
  
    rownames(abund) <- abund[,1]
  
    #head(abund[,2000:2010])
    
    time<- rev(seq(from =0, to =80000, by = 25))
    
    timeall <- time[2361:2928]
    timebreak1 <- time[2708:2928]
    timebreak2 <- time[2361:2707]
    
    ## Provide row id of simulation of interest
    abundall <- as.numeric(abund[target,2362:2929])     
    abundbreak1 <- as.numeric(abund[target,2709:2929])
    abundbreak2 <- as.numeric(abund[target,2362:2708])
    
    #graphics.off()
    temp = sprintf("abundfig_%s.pdf",filename)
    pdf(file = temp,width = 5,height=5)
    plot(timeall, abundall, pch = 16, cex = 0.2, col = "blue", main = "Change in Ne", 
         ylim=c(0,min(2000000,max(abundall+10000))),xlab = "Time (year)", ylab = "Abundance")   #  ,
    
    break1coef <- round(lm(abundbreak1~timebreak1)$coefficients[2],2) 
    break2coef <- round(lm(abundbreak2~timebreak2)$coefficients[2],2)
    
    abline(lm(abundbreak1~timebreak1))
    abline(lm(abundbreak2~timebreak2))
    
    
    ## Extinction presumed to have occurred by 3000 years ago; and only interested in records after 21,000 years ago
    time2 <- (eff_popsize[61:420,1])
    abund2 <-(eff_popsize[61:420,2])
    data2 <- cbind(time2, abund2)
    
    points(time2, abund2, cex = 0.1, col = "red",pch=20)
  #}
    dev.off()
  
  
  #if(toPlot=="EXTINCTION_DATE"){
    matches <- as.numeric(unlist(regmatches(target, gregexpr("[[:digit:]]+", target))))
    nichebreadth <- matches[1]
    sample <- matches[2]
    
    colors <- c("red","orange","yellow","green","blue","purple","violet","black")
    cutoffs <- c(100000,70000,40000,20000,15000,10000,5000,4000,0)
    
    #colname <- sprintf("NicheBreadth%i_LHS_Sample%i.mp",nichebreadth,sample)
    
    filename <- sprintf("NicheBreadth%i_FinalYear_output_alldata_formatted.csv",nichebreadth)
    
    setwd(SIM_RESULTS_DIR)
    
    extinct <- read.csv(filename, header = T,stringsAsFactors = FALSE)    # if(!"abund"%in%ls(name=.GlobalEnv))  
    
    #rownames(extinct) <-  
    
    xcords <- extinct[,1]
    ycords <- extinct[,2]
    head(extinct[,1:10])
    
    if(any(xcords<0)) xcords[which(xcords<0)] <- 180+(180 - abs(xcords[which(xcords<0)]))
    
    thissample <- extinct[,target]
    
    time<- rev(seq(from =0, to =80000, by = 25))
    
    ext.time <- time[thissample]
    
    ext.time.cols <- 9-as.numeric(cut(ext.time,cutoffs))
    
    #exttime_scaled <- ext.time/80000
    
    cols <- colors[ext.time.cols] #colorRamp(c("blue","red"), bias = 1, space = c("rgb"),
             #         interpolate = c("linear"), alpha = FALSE)(exttime_scaled)
    
    #graphics.off()
    #map <- get_map(location = 'Asia', zoom = 3)
    newmap <- getMap(resolution = "low")
    #ggmap(map)
    
    
    
    par(mfrow=c(2,1))
    plot(newmap,xlim = c(20, 200),
         ylim = c(30, 80),
         asp = 1,
         main="simulated"
    )
    points(xcords,ycords,pch=20,cex=0.01,col=cols)
    
    ### location of final extinction
    ndx <- which(ext.time==min(ext.time,na.rm=T))
    xcords_lastext <- extinct[ndx,1]
    
    if(any(xcords_lastext<0)) xcords_lastext[which(xcords_lastext<0)] <- 180+(180 - abs(xcords_lastext[which(xcords_lastext<0)]))
    
    ycoord_lastext <- extinct[ndx,2]
    
    # plot(newmap,xlim = c(50, 219),
    #      ylim = c(30, 71),
    #      asp = 1,
    #      main="real"
    # )
    points(xcords_lastext,ycoord_lastext,pch="X",cex=2,col="black")
    
    
    plot(newmap,xlim = c(20, 200),
         ylim = c(30, 80),
         asp = 1,
         main="real"
    )
    
    ext.time.cols <- 9-as.numeric(cut(extinct_date$Ext,cutoffs))
    
    #exttime_scaled <- ext.time/80000
    
    cols <- colors[ext.time.cols]
    #cols <- colorRamp(c("blue","red"), bias = 1, space = c("rgb"),
     #                 interpolate = c("linear"), alpha = FALSE)(extinct_date$Ext/80000)

    points(extinct_date$Mammoth_estimates.CI.Barnoski._Long,extinct_date$Mammoth_estimates.CI.Barnoski._Lat,pch=20,cex=0.01,col=cols)    
    

    points(last.locations[5,],pch="X",cex=2,col="yellow")
    
    
    
    
  #}
  #graphics.off()
  
}


###############
# MAKE MOVIE   (NOTE: need ffmpeg and imagemagick installed)
###############

# for testing
target = "NicheBreadth70_LHS_Sample1.mp"

MakeMovie <- function(target){
  
  setwd(SIM_RESULTS_DIR)
  
  ############
  # SET UP STORAGE STRUCTURES FOR KEY RESULTS METRICS
  ############
  
  ###########################
  # GET RESULTS
  ###########################
  
  result = tryCatch({   ## try to catch errors! 
    ### step 1: set up the MP file connection 
    
    MPcon <- file(target, 'r')
    
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
    
    PopMat <- array(0,dim=c(NPOPS,TIMESTEPS))
    
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
      PopMat[pop,] <- Nvec
      #eval(parse(text=sprintf("SimInfo$PopAbund[%s,] <- Nvec",pop)))	   # RESULT: POP ABUNDANCE             
    }   # end loop through pops 
    
    if(isOpen(MPcon)){ 
      close.connection(MPcon)
      rm("MPcon")
    }
    
    if(exists("MPcon")){
      rm("MPcon")
    }
  }, warning = function(w){
    as.character(w)
  }, error = function(e){
    setwd(thisFolder)
    filename <- sprintf("%s.RData",name)
    eval(parse(text=sprintf("%s <- e",name)))
    eval(parse(text=sprintf("save(%s,file=filename)",name)))   # save to disk
    eval(parse(text=sprintf("rm(%s)",name)))   # remove from memory
    as.character(e)
  }, finally = {
    SimInfo$PopAbund <- NULL   ## remove from memory
    if(exists("MPcon")){
      close.connection(MPcon)
      rm("MPcon")
    }
  })   # end tryCatch

  
  thisMOVIE_DIR <- sprintf("%s\\%s",getwd(),gsub("\\.mp","",target))
  if(is.na(file.info(thisMOVIE_DIR)[1,"isdir"])) dir.create(thisMOVIE_DIR)
  
  setwd(thisMOVIE_DIR)
  
  #### remove figures that are not part of this simulation?
  
  fileyrs <- as.numeric(unlist(regmatches(list.files(), gregexpr("[[:digit:]]+", list.files()))))
  notthissim <- fileyrs > TIMESTEPS
  toremove <- list.files()[notthissim]
  
  if(any(notthissim)) file.remove(toremove)
  
  ################
  # MAKE PLOTS
  
  width = 800
  height= 300

  t = 1
  counter=1
  for(t in seq(1,TIMESTEPS,4)){
    file = sprintf("AbundanceMap_year%04d.tif",counter)
    tiff(file, width=width,height=height)
    
    #colors <- c("red","orange","yellow","green","blue","purple","violet","black","black")
    colors <- viridis::plasma(9)
    
    
    
    cutoffs <- c(100000,10000,5000,2500,1000,500,100,10,0,-1)
    
    #rownames(extinct) <-  
    
    xcords <- GridCellAttributes$x.cord
    ycords <- GridCellAttributes$y.cord
    head(extinct[,1:10])
    
    if(any(xcords<0)) xcords[which(xcords<0)] <- 180+(180 - abs(xcords[which(xcords<0)]))
    
    time<- rev(seq(from =0, to =80000, by = 25))
    
    abundnow <- PopMat[,t]
    
    abund.cols <- 10-as.numeric(cut(abundnow,cutoffs))
    
    cols <- colors[abund.cols] #colorRamp(c("blue","red"), bias = 1, space = c("rgb"),
    #         interpolate = c("linear"), alpha = FALSE)(exttime_scaled)
    
    newmap <- getMap(resolution = "low")
    
    plot(newmap,xlim = c(20, 200),
         ylim = c(30, 80),
         asp = 1,
         main=sprintf("%s years bp",time[t])
    )
    nonzero <- abundnow>0
    points(xcords[nonzero],ycords[nonzero],pch=20,cex=0.01,col=cols[nonzero])
  
   
    dev.off() 
    counter=counter+1
  }
  
  ## NOTE: need command line like this: ffmpeg -f image2 -framerate 2 -i AbundanceMap_year%03d.tif -s 500x500 test.avi -y
  
  # MAKING THE REAL MOVIE HERE! USE IMAGE MAGICK AND FFMPEG SOFTWARE  (https://blogazonia.wordpress.com/2016/01/19/making-a-movie-with-r/)
  
  
  # create the movie
  cmd_abundmov <- paste0("ffmpeg -f image2 -framerate 2 -i AbundanceMap_year%04d.tif -s 800x300 ", 
                         sprintf("%s\\AbundanceMovie.avi",thisMOVIE_DIR)," -y")
  
  
  #sink(tempfile())
  system(cmd_abundmov,ignore.stdout = T,ignore.stderr = T)
  
  #sink()
}  ### end function "MakeMovie"




####################
# LOAD SIMULATIONs
#
# specify directory that contains all simulation results and parameter values, prepared for ABC
###################

KEVIN_LAB <- TRUE # FALSE
DAMIEN <- FALSE # TRUE 

if(KEVIN_LAB) BASE_RESULTS_DIR <- "E:\\MammothResults"
if(DAMIEN) BASE_RESULTS_DIR <- "F:\\MammothResultsV2"

SIM_RESULTS_DIR <- sprintf("%s\\ABCAnalysisInputs\\",BASE_RESULTS_DIR)
TARGETS_DIR <- sprintf("%s\\ABCAnalysisTargets",BASE_RESULTS_DIR)
GENETIC_TARGETS_DIR <- sprintf("%s\\Genetic",TARGETS_DIR)
FOSSIL_TARGETS_DIR <- sprintf("%s\\Fossil",TARGETS_DIR)   
INPUTS_DIR <- sprintf("%s\\ABCAnalysis",BASE_RESULTS_DIR)

setwd(INPUTS_DIR)


dat1 <- read.csv("NicheBreadth40_ABC_data.csv", header = T,stringsAsFactors = FALSE)
dat2 <- read.csv("NicheBreadth50_ABC_data.csv", header = T,stringsAsFactors = FALSE)
dat3 <- read.csv("NicheBreadth60_ABC_data.csv", header = T,stringsAsFactors = FALSE)
dat4 <- read.csv("NicheBreadth70_ABC_data.csv", header = T,stringsAsFactors = FALSE)
dat5 <- read.csv("NicheBreadth80_ABC_data.csv", header = T,stringsAsFactors = FALSE)
dat6 <- read.csv("NicheBreadth90_ABC_data.csv", header = T,stringsAsFactors = FALSE)

all <- rbind(dat1,dat2,dat3,dat4,dat5,dat6)

nrow(all)   # 89999 simulations

head(all)   # check the combined dataframe 

## Generate a version with no NA values (some analyses do not allow NA values)

all_nona <- na.omit(all)   ### remove all na observations


### Generate a version with no persistent simulations (remove all sims that never went extinct or extinct in first time step)

keep <- !is.na(all$extinction_yr)
length(which(keep))   # 20126 simulations remaining
nrow(all)-length(which(keep))   # 20087 sims did not undergo global extinction within the simulation period

all_extinct <- all[keep,]     # remove the ca. 20k simulations that were not plausible.


# explore any remaining NAs in the data frame, after removal of the non-extinctions...
remaining.nas <- all_extinct[which(is.na(all_extinct),arr.ind=TRUE)[,1],]

nrow(remaining.nas)     # 306 simulations with NAs for dist_centroid and dist_fossil
remaining.nas$dist_fossil
remaining.nas$dist_centroid
tail(remaining.nas)



####################
# LOAD REAL-WORLD DATA 

# targets for matching with simulated data- for ABC
####################

setwd(GENETIC_TARGETS_DIR)
  
    # load data on effective population size from genetic data
eff_popsize <- read.csv("Mammoth_Ne.csv", header = T)


setwd(FOSSIL_TARGETS_DIR)

  # load data on extinction date from the fossil record -- INTERPOLATED  STATISTICAL DATA
extinct_date <- read.csv("Barnoski_matched_estimates_fossils.csv", header = T)

extinct_date$Mammoth_estimates.CI.Barnoski._Long[which(extinct_date$Mammoth_estimates.CI.Barnoski._Long<0)]<- 180+(180 - abs(extinct_date$Mammoth_estimates.CI.Barnoski._Long[which(extinct_date$Mammoth_estimates.CI.Barnoski._Long<0)]))


  # Generate a version that has extinction date only for sites with fossil records (not interpolated)
extinct_date_fossilsites <- extinct_date[extinct_date$Fossil.Data==1,]

    # final locations of extinction, observed
last.locations <- extinct_date[which(extinct_date$Ext==min(extinct_date$Ext)),c(1,2)]

names(last.locations) <- c("X","Y")

hist(extinct_date$Ext)

####################
# OBSERVED STATISTICS FOR MATCHING WITH SIMULATION RESULTS
###################

test.statistics.obs <- c(64.0765, 11.4448, 0, 0, 3993, 0, 0)

names(test.statistics.obs) <- c("abundcoef1","abundcoef2","extinctpattern_all","extinctpattern_fossilsites","extinction_yr","dist_centroid","dist_fossil")


test.statistics.obs

####################
# SET UP VARIABLEs FOR ABC

# All variable (estimable) parameters and all summary statistics for matching with data
###################

      #### What is the multivariate parameter space we are attempting to estimate?

estimable.params <- c(
  "RMAX",      
  "SD",        # env. stochasticity
  "ALLEE",     
  "DENSITY",   # max density of mammoths
  "DISP1",     # determines rate of staying/leaving a given cell
  "DISP2",     # max dispersal distance
  "HARV",      # maximal harvest rate
  "HARVZ",     # regulates shape of functional response
  "HUMAN",
  "summ_precip_median",
  "summ_precip_range",
  "jan_temp_median",
  "jan_temp_range",
  "jul_temp_median",
  "jul_temp_range"
)

estimable.params

   #### Which test statistics should we use? 

test.statistics <- c(
  "abundcoef1",    # slope of abundance decline, most recent
  "abundcoef2",     # slope of abuncance decline, older
  "extinctpattern_all",
  "extinctpattern_fossilsites",
  "extinction_yr",
  "dist_centroid",
  "dist_fossil"
)

test.statistics


### NOTES: which summary statistics to use??

## use just coefficient 1 (abundance by itself) 
## timing of extinction is critical
## play around with both distance measures but don't use both at the same time because meaningless...
## pattern by itself
## pattern with extinction?   bring NAs in?


#####################
# PERFORM BASIC VISUALIZATIONS AND DATA CHECKS
#####################
 
### VISUALIZE TEST STATISTICS  (with targets overlaid)

par(mfrow=c(2,3))
for(i in test.statistics){
  hist(all[,i],main=i)
  abline(v=test.statistics.obs[i],col="green",lwd=2)
}


### VISUALIZE PRIORS (check for uniformity- don't need to run this every time)

visualize.priors()

####################
# SELECT SUMMARY STATISTICS
####################

?selectsumm

test <- prepareForABC(all_nona,subset=c(1:7),standardized = TRUE)

selectstats <- abctools::selectsumm(test.statistics.obs, par.sim, 
                                     stat.sim, ssmethod =mincrit,  # AS.select,    #  
                                     final.dens = FALSE)

# which summary statistics were selected by the "selectsumm" algorithm?   # # given by selectstats$best
colnames(selectstats$best) <- test.statistics

selectstats$best

#### make final selections...

test.statistics


   ## here we subset the targets to the best set identified above. This time extinction year is best
prepareForABC(all,subset=c(4,5,6),standardized=T)


head(par.sim)
head(stat.sim)   # make sure the right summary statistics are being used
stat.obs


#####################
# PERFORM PRELIMINARY ABC, VISUALIZATIONS AND DATA CHECKS
#####################

### Using the rejection Method to find the best fits


#####
# ABC analysis
#####
rej <- abc(target=stat.obs, param=par.sim, sumstat=stat.sim, tol=.005, 
           method ="neuralnet", numnet=20, sizenet=10,trace=T,maxit=1000,lambda = c(0.0001,0.001,0.01))   # rejection  # loclinear



str(rej)

# NOTE: neural net method seems to perform the best...

   # explore the abc object...
head(rej$adj.values)
head(rej$unadj.values)
rej$weights

rej$dist    # this is the euclidean distance between simulation results and observations 

######     
 # summary and diagnostic plots
######

plot(rej,param=par.sim)     # visualize all priors and posteriors

rej$ss  # all simulation runs that were kept for the posterior

    ### credible intervals
summary(rej)



    ### cross validation (can take a while with neural net)

### the cross validation does NOT use the observed targets. Instead, it takes simulated data, one at a time, as
#    the "observed" targets, and uses the remaining simulations to perform an ABC analysis with the new
#    "observed" summary statistics. 
#    this way, we can look at the overall quality of the ABC (ability to estimate selected parameters) 
#    if the algorithm is good (summary stats are informative) then the parameters estimated by ABC should
#    match the parameters of the simulation that in fact generated the psuedo-observed data... 

cv.res.reg <- cv4abc(data.frame(Harv=par.sim$HARV,Harvz=par.sim$HARVZ,Jantemp=par.sim$jan_temp_median), stat.sim,
                     nval=20, tols=c(0.01,0.005),   # method="rejection")
                     method ="neuralnet", numnet=20, sizenet=10,trace=F,maxit=1000,lambda = c(0.0001,0.001,0.01)) 
plot(cv.res.reg)

summary(cv.res.reg)   # good summary of model performance



#############################
# VISUALIZE ABUNDANCE TRAJECTORIES AND EXTINCTION TIMES
#############################

### find the best simulations to explore further

if(length(test.statistics)>1){
  ndx <- as.numeric(rownames(rej$ss))
  if(is.null(ndx)){
    intersect(rej$adj.values,par.sim)
  }
  #ranking <- order(apply(as.matrix(rej$ss,min),1,sum))
  ranking <- order(as.vector(rej$weights))
  ndx <- ndx[ranking]
  ndx_min <- as.numeric(names(which.min(apply(as.matrix(rej$ss,min),1,sum))))
}else{
  ndx <- as.numeric(names(rej$ss))
  ranking <- order(as.vector(rej$ss))
  ndx <- ndx[ranking]
  #ndx_min <- as.numeric(names(which(rej$ss == min(rej$ss), arr.ind = TRUE)))[1]
  #mins <- apply(as.matrix(rej$ss,min),2,min)
}



index <- ndx[1]   # one of a set of all models used in the posterior
index <- ndx_min   # the best fit model

target <- all[index,]$model

#toPlot <-  "EXTINCTION_DATE"  # "ABUNDANCE"     # 




### NOTE on color scheme: red is >70k time to extinction. 
                      # orange is from 40 to 70k
                      # yellow is from 20 to 40k
                      # green is from 15 to 20k
                      # blue is from 10 to 15k
                      # purple is from 5 to 10k
                      # indigo is from 4-5k
                      # black is from 3-4k

graphics.off()
visualize(target)



#####################
# VISUALIZE POSTERIORS
#####################


    # only works for elimination method
 # visualize.posteriors()



#####################
# from ABC package

## Note: cv4postpr is only for model selection...
          # also "postpr" is only for model selection
          # use "cv4abc" function for our purposes


## Note: the neural network method seems to be the best- local linear does not work at all!


##### TESTS

# does removing nas do anything to the results? it shouldn't!
# does standardizing the targets have an effect?





######


function (target, param, sumstat, tol, method, hcorr = TRUE, 
          transf = "none", logit.bounds = c(0, 0), subset = NULL, kernel = "epanechnikov", 
          numnet = 10, sizenet = 5, lambda = c(1e-04, 0.001, 0.01), 
          trace = FALSE, maxit = 500, ...) 
{
  call <- match.call()
  if (missing(target)) 
    stop("'target' is missing")
  if (missing(param)) 
    stop("'param' is missing")
  if (missing(sumstat)) 
    stop("'sumstat' is missing")
  if (!is.matrix(param) && !is.data.frame(param) && !is.vector(param)) 
    stop("'param' has to be a matrix, data.frame or vector.")
  if (!is.matrix(sumstat) && !is.data.frame(sumstat) && !is.vector(sumstat)) 
    stop("'sumstat' has to be a matrix, data.frame or vector.")
  if (missing(tol)) 
    stop("'tol' is missing")
  if (missing(method)) 
    stop("'method' is missing with no default")
  if (!any(method == c("rejection", "loclinear", "neuralnet", 
                       "ridge"))) {
    stop("Method must be rejection, loclinear, or neuralnet or ridge")
  }
  if (method == "rejection") 
    rejmethod <- TRUE
  else rejmethod <- FALSE
  if (!any(kernel == c("gaussian", "epanechnikov", "rectangular", 
                       "triangular", "biweight", "cosine"))) {
    kernel <- "epanechnikov"
    warning("Kernel is incorrectly defined. Setting to default (Epanechnikov)")
  }
  if (is.data.frame(param)) 
    param <- as.matrix(param)
  if (is.data.frame(sumstat)) 
    sumstat <- as.matrix(sumstat)
  if (is.list(target)) 
    target <- unlist(target)
  if (is.vector(sumstat)) 
    sumstat <- matrix(sumstat, ncol = 1)
  if (length(target) != dim(sumstat)[2]) 
    stop("Number of summary statistics in 'target' has to be the same as in 'sumstat'.")
  nss <- length(sumstat[1, ])
  cond1 <- !any(as.logical(apply(sumstat, 2, function(x) length(unique(x)) - 
                                   1)))
  if (cond1) 
    stop("Zero variance in the summary statistics.")
  ltransf <- length(transf)
  if (is.vector(param)) {
    numparam <- 1
    param <- matrix(param, ncol = 1)
  }
  else numparam <- dim(param)[2]
  for (i in 1:ltransf) {
    if (sum(transf[i] == c("none", "log", "logit")) == 0) {
      stop("Transformations must be none, log, or logit.")
    }
    if (transf[i] == "logit") {
      if (logit.bounds[i, 1] >= logit.bounds[i, 2]) {
        stop("Logit bounds are incorrect.")
      }
    }
  }
  if (rejmethod) {
    if (!all(transf == "none")) {
      warning("No transformation is applied when the simple rejection is used.", 
              call. = F)
    }
    transf[1:numparam] <- "none"
  }
  else {
    if (numparam != ltransf) {
      if (length(transf) == 1) {
        transf <- rep(transf[1], numparam)
        warning("All parameters are \"", transf[1], "\" transformed.", 
                sep = "", call. = F)
      }
      else stop("Number of parameters is not the same as number of transformations.", 
                sep = "", call. = F)
    }
  }
  gwt <- rep(TRUE, length(sumstat[, 1]))
  gwt[attributes(na.omit(sumstat))$na.action] <- FALSE
  if (missing(subset)) 
    subset <- rep(TRUE, length(sumstat[, 1]))
  gwt <- as.logical(gwt * subset)
  if (!length(colnames(param))) {
    warning("No parameter names are given, using P1, P2, ...")
    paramnames <- paste("P", 1:numparam, sep = "")
  }
  else paramnames <- colnames(param)
  if (!length(colnames(sumstat))) {
    warning("No summary statistics names are given, using S1, S2, ...")
    statnames <- paste("S", 1:nss, sep = "")
  }
  else statnames <- colnames(sumstat)
  scaled.sumstat <- sumstat
  for (j in 1:nss) {
    scaled.sumstat[, j] <- normalise(sumstat[, j], sumstat[, 
                                                           j][gwt])
  }
  for (j in 1:nss) {
    target[j] <- normalise(target[j], sumstat[, j][gwt])
  }
  sum1 <- 0
  for (j in 1:nss) {
    sum1 <- sum1 + (scaled.sumstat[, j] - target[j])^2
  }
  dist <- sqrt(sum1)
  dist[!gwt] <- floor(max(dist[gwt]) + 10)
  nacc <- ceiling(length(dist) * tol)
  ds <- sort(dist)[nacc]
  wt1 <- (dist <= ds)
  aux <- cumsum(wt1)
  wt1 <- wt1 & (aux <= nacc)
  if (kernel == "gaussian") {
    wt1 <- rep(TRUE, length(dist))
  }
  for (i in 1:numparam) {
    if (transf[i] == "log") {
      if (min(param[, i]) <= 0) {
        cat("log transform: values out of bounds - correcting...")
        x.tmp <- ifelse(param[, i] <= 0, max(param[, 
                                                   i]), param[, i])
        x.tmp.min <- min(x.tmp)
        param[, i] <- ifelse(param[, i] <= 0, x.tmp.min, 
                             param[, i])
      }
      param[, i] <- log(param[, i])
    }
    else if (transf[i] == "logit") {
      if (min(param[, i]) <= logit.bounds[i, 1]) {
        x.tmp <- ifelse(param[, i] <= logit.bounds[i, 
                                                   1], max(param[, i]), param[, i])
        x.tmp.min <- min(x.tmp)
        param[, i] <- ifelse(param[, i] <= logit.bounds[i, 
                                                        1], x.tmp.min, param[, i])
      }
      if (max(param[, i]) >= logit.bounds[i, 2]) {
        x.tmp <- ifelse(param[, i] >= logit.bounds[i, 
                                                   2], min(param[, i]), param[, i])
        x.tmp.max <- max(x.tmp)
        param[, i] <- ifelse(param[, i] >= logit.bounds[i, 
                                                        2], x.tmp.max, param[, i])
      }
      param[, i] <- (param[, i] - logit.bounds[i, 1])/(logit.bounds[i, 
                                                                    2] - logit.bounds[i, 1])
      param[, i] <- log(param[, i]/(1 - param[, i]))
    }
  }
  ss <- sumstat[wt1, ]
  unadj.values <- param[wt1, ]
  statvar <- as.logical(apply(cbind(sumstat[wt1, ]), 2, function(x) length(unique(x)) - 
                                1))
  cond2 <- !any(statvar)
  if (cond2 && !rejmethod) 
    stop("Zero variance in the summary statistics in the selected region. Try: checking summary statistics, choosing larger tolerance, or rejection method.")
  if (rejmethod) {
    if (cond2) 
      warning("Zero variance in the summary statistics in the selected region. Check summary statistics, consider larger tolerance.")
    weights <- rep(1, length = sum(wt1))
    adj.values <- NULL
    residuals <- NULL
    lambda <- NULL
  }
  else {
    if (cond2) 
      cat("Warning messages:\nStatistic(s)", statnames[!statvar], 
          "has/have zero variance in the selected region.\nConsider using larger tolerance or the rejection method or discard this/these statistics, which might solve the collinearity problem in 'lsfit'.\n", 
          sep = ", ")
    if (kernel == "epanechnikov") 
      weights <- 1 - (dist[wt1]/ds)^2
    if (kernel == "rectangular") 
      weights <- dist[wt1]/ds
    if (kernel == "gaussian") 
      weights <- 1/sqrt(2 * pi) * exp(-0.5 * (dist/(ds/2))^2)
    if (kernel == "triangular") 
      weights <- 1 - abs(dist[wt1]/ds)
    if (kernel == "biweight") 
      weights <- (1 - (dist[wt1]/ds)^2)^2
    if (kernel == "cosine") 
      weights <- cos(pi/2 * dist[wt1]/ds)
    if (method == "loclinear") {
      fit1 <- lsfit(scaled.sumstat[wt1, ], param[wt1, ], 
                    wt = weights)
      pred <- t(structure(cbind(fit1$coefficients)[fit1$qr$pivot, 
                                                   ], names = names(fit1$coefficients))) %*% c(1, 
                                                                                               target)
      pred <- matrix(pred, ncol = numparam, nrow = sum(wt1), 
                     byrow = TRUE)
      residuals <- param[wt1, ] - t(t(structure(cbind(fit1$coefficients)[fit1$qr$pivot, 
                                                                         ], names = names(fit1$coefficients))) %*% t(cbind(1, 
                                                                                                                           scaled.sumstat[wt1, ])))
      residuals <- cbind(residuals)
      the_m <- apply(residuals, FUN = mean, 2)
      residuals <- sapply(1:numparam, FUN = function(x) {
        residuals[, x] - the_m[x]
      })
      pred <- sapply(1:numparam, FUN = function(x) {
        pred[, x] + the_m[x]
      })
      sigma2 <- apply(as.matrix(residuals), FUN = function(x) {
        sum((x)^2 * weights)/sum(weights)
      }, MARGIN = 2)
      aic <- sum(wt1) * sum(log(sigma2)) + 2 * (nss + 1) * 
        numparam
      bic <- sum(wt1) * sum(log(sigma2)) + log(sum(wt1)) * 
        (nss + 1) * numparam
      if (hcorr == TRUE) {
        fit2 <- lsfit(scaled.sumstat[wt1, ], log(residuals^2), 
                      wt = weights)
        auxaux <- t(structure(cbind(fit2$coefficients)[fit2$qr$pivot, 
                                                       ], names = names(fit2$coefficients))) %*% c(1, 
                                                                                                   target)
        pred.sd <- sqrt(exp(auxaux))
        pred.sd <- matrix(pred.sd, nrow = sum(wt1), ncol = numparam, 
                          byrow = T)
        pred.si <- t(t(structure(cbind(fit2$coefficients)[fit2$qr$pivot, 
                                                          ], names = names(fit2$coefficients))) %*% t(cbind(1, 
                                                                                                            scaled.sumstat[wt1, ])))
        pred.si <- sqrt(exp(pred.si))
        adj.values <- pred + (pred.sd * residuals)/pred.si
        residuals <- (pred.sd * residuals)/pred.si
      }
      else {
        adj.values <- pred + residuals
      }
      colnames(adj.values) <- colnames(unadj.values)
      lambda <- NULL
    }
    if (method == "neuralnet") {
      linout <- TRUE
      param.mad <- c()
      for (i in 1:numparam) {
        param.mad[i] <- mad(param[, i][gwt])
        param[, i] <- normalise(param[, i], param[, i][gwt])
      }
      lambda <- sample(lambda, numnet, replace = T)
      fv <- array(dim = c(sum(wt1), numparam, numnet))
      pred <- matrix(nrow = numparam, ncol = numnet)
      for (i in 1:numnet) {
        fit1 <- nnet(scaled.sumstat[wt1, ], param[wt1, 
                                                  ], weights = weights, decay = lambda[i], size = sizenet, 
                     trace = trace, linout = linout, maxit = maxit, 
                     ...)
        cat(i)
        fv[, , i] <- fit1$fitted.values
        pred[, i] <- predict(fit1, data.frame(rbind(target)))
      }
      cat("\n")
      pred.med <- apply(pred, 1, median)
      pred.med <- matrix(pred.med, nrow = sum(wt1), ncol = numparam, 
                         byrow = T)
      fitted.values <- apply(fv, c(1, 2), median)
      residuals <- param[wt1, ] - fitted.values
      if (hcorr == TRUE) {
        pred2 <- matrix(nrow = numparam, ncol = numnet)
        fv2 <- array(dim = c(sum(wt1), numparam, numnet))
        for (i in 1:numnet) {
          fit2 <- nnet(scaled.sumstat[wt1, ], log(residuals^2), 
                       weights = weights, decay = lambda[i], size = sizenet, 
                       trace = trace, linout = linout, ...)
          cat(i)
          fv2[, , i] <- fit2$fitted.values
          pred2[, i] <- predict(fit2, data.frame(rbind(target)))
        }
        cat("\n")
        pred.sd <- sqrt(exp(apply(pred2, 1, median)))
        pred.sd <- matrix(pred.sd, nrow = sum(wt1), ncol = numparam, 
                          byrow = T)
        fv.sd <- sqrt(exp(apply(fv2, c(1, 2), median)))
        adj.values <- pred.med + (pred.sd * residuals)/fv.sd
        residuals <- (pred.sd * residuals)/fv.sd
      }
      else {
        adj.values <- pred.med + residuals
      }
      colnames(adj.values) <- colnames(unadj.values)
      for (i in 1:numparam) {
        adj.values[, i] <- adj.values[, i] * param.mad[i]
      }
    }
    if (method == "ridge") {
      param.mad <- c()
      for (i in 1:numparam) {
        param.mad[i] <- mad(param[, i][gwt])
        param[, i] <- normalise(param[, i], param[, i][gwt])
      }
      numnet <- length(lambda)
      fv <- array(dim = c(sum(wt1), numparam, numnet))
      pred <- matrix(nrow = numparam, ncol = numnet)
      mataux <- sqrt(diag(weights))
      paramaux <- as.matrix(mataux %*% param[wt1, ])
      scaledaux <- mataux %*% scaled.sumstat[wt1, ]
      for (parcur in (1:numparam)) {
        fit1 <- lm.ridge(paramaux[, parcur] ~ scaledaux, 
                         lambda = lambda)
        for (i in 1:numnet) {
          fv[, parcur, i] <- drop(cbind(1, scaled.sumstat[wt1, 
                                                          ]) %*% (rbind(coef(fit1))[i, ]))
          pred[parcur, i] <- drop(c(1, target) %*% (rbind(coef(fit1))[i, 
                                                                      ]))
        }
      }
      pred.med <- apply(pred, 1, median)
      pred.med <- matrix(pred.med, nrow = sum(wt1), ncol = numparam, 
                         byrow = T)
      fitted.values <- apply(fv, c(1, 2), median)
      residuals <- param[wt1, ] - fitted.values
      if (hcorr == TRUE) {
        pred2 <- matrix(nrow = numparam, ncol = numnet)
        fv2 <- array(dim = c(sum(wt1), numparam, numnet))
        for (parcur in (1:numparam)) {
          lresidaux <- (mataux %*% (log(residuals[, parcur]^2)))
          fit2 <- lm.ridge(lresidaux ~ scaledaux, lambda = lambda)
          for (i in 1:numnet) {
            fv2[, parcur, i] <- drop(cbind(1, scaled.sumstat[wt1, 
                                                             ]) %*% (rbind(coef(fit2))[i, ]))
            pred2[parcur, i] <- drop(c(1, target) %*% 
                                       (rbind(coef(fit2))[i, ]))
          }
        }
        cat("\n")
        pred.sd <- sqrt(exp(apply(pred2, 1, median)))
        pred.sd <- matrix(pred.sd, nrow = sum(wt1), ncol = numparam, 
                          byrow = T)
        fv.sd <- sqrt(exp(apply(fv2, c(1, 2), median)))
        adj.values <- pred.med + (pred.sd * residuals)/fv.sd
        residuals <- (pred.sd * residuals)/fv.sd
      }
      else {
        adj.values <- pred.med + residuals
      }
      colnames(adj.values) <- colnames(unadj.values)
      for (i in 1:numparam) {
        adj.values[, i] <- adj.values[, i] * param.mad[i]
      }
    }
  }
  if (numparam == 1) {
    unadj.values <- matrix(unadj.values, ncol = 1)
    if (method != "rejection") {
      adj.values <- matrix(adj.values, ncol = 1)
      residuals <- matrix(residuals, ncol = 1)
    }
  }
  for (i in 1:numparam) {
    if (transf[i] == "log") {
      unadj.values[, i] <- exp(unadj.values[, i])
      adj.values[, i] <- exp(adj.values[, i])
    }
    else if (transf[i] == "logit") {
      unadj.values[, i] <- exp(unadj.values[, i])/(1 + 
                                                     exp(unadj.values[, i]))
      unadj.values[, i] <- unadj.values[, i] * (logit.bounds[i, 
                                                             2] - logit.bounds[i, 1]) + logit.bounds[i, 1]
      adj.values[, i] <- exp(adj.values[, i])/(1 + exp(adj.values[, 
                                                                  i]))
      adj.values[, i] <- adj.values[, i] * (logit.bounds[i, 
                                                         2] - logit.bounds[i, 1]) + logit.bounds[i, 1]
    }
  }
  abc.return(transf, logit.bounds, method, call, numparam, 
             nss, paramnames, statnames, unadj.values, adj.values, 
             ss, weights, residuals, dist, wt1, gwt, lambda, hcorr, 
             aic, bic)
}

