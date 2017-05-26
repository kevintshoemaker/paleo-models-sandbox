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

rm(list=ls())

####################
# LOAD PACKAGES
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


####################
# LOAD SIMULATIONs
#
# specify directory that contains all simulation results and parameter values, prepared for ABC
###################

KEVIN_LAB <- TRUE
DAMIEN <- FALSE 

if(KEVIN_LAB) BASE_RESULTS_DIR <- "E:\\MammothResults"

SIM_RESULTS_DIR <- sprintf("%s\\ABCAnalysisInputs\\",BASE_RESULTS_DIR)
TARGETS_DIR <- sprintf("%s\\ABCAnalysisTargets",BASE_RESULTS_DIR)
GENETIC_TARGETS_DIR <- sprintf("%s\\Genetic",TARGETS_DIR)
FOSSIL_TARGETS_DIR <- sprintf("%s\\Fossil",TARGETS_DIR)   
INPUTS_DIR <- sprintf("%s\\ABCAnalysis\\MergedInputs",BASE_RESULTS_DIR)

setwd(INPUTS_DIR)


dat1 <- read.csv("NicheBreadth40_ABC_data_revised.csv", header = T,stringsAsFactors = FALSE)
dat2 <- read.csv("NicheBreadth50_ABC_data_revised.csv", header = T,stringsAsFactors = FALSE)
dat3 <- read.csv("NicheBreadth60_ABC_data_revised.csv", header = T,stringsAsFactors = FALSE)
dat4 <- read.csv("NicheBreadth70_ABC_data_revised.csv", header = T,stringsAsFactors = FALSE)
dat5 <- read.csv("NicheBreadth80_ABC_data_revised.csv", header = T,stringsAsFactors = FALSE)
dat6 <- read.csv("NicheBreadth90_ABC_data_revised.csv", header = T,stringsAsFactors = FALSE)

all <- rbind(dat1,dat2,dat3,dat4,dat5,dat6)

nrow(all)   # 74931 simulations

head(all)   # check the combined dataframe 

## Generate a version with no NA values (some analyses do not allow NA values)

all_nona <- na.omit(all)   ### remove all na observations


### Generate a version with no persistent simulations (remove all sims that never went extinct or extinct in first time step

keep <- !is.na(all$extinction_yrs)
length(which(keep))   # 20126 simulations remaining
nrow(all)-length(which(keep))   # 54805 sims did not undergo global extinction within the simulation period

all_extinct <- all[keep,]     # remove the ca. 55k simulations that were not plausible.


  # explore any remaining NAs in the data frame, after removal of the non-extinctions...
remaining.nas <- all_extinct[which(is.na(all_extinct),arr.ind=TRUE)[,1],]

nrow(remaining.nas)     # 240 simulations with NAs for dist_centroid and dist_fossil
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

  # load data on extinction date from the fossil record, only fossil locations... INTERPOLATED
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

test.statistics.obs <- c(64.0765, 11.4448, 0, 0, 0, 0)

names(test.statistics.obs) <- c("abundcoef1","abundcoef2","extinctpattern","extinction_yrs","dist_centroid","dist_fossil")

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
  "extinctpattern",
  "extinction_yrs",
  "dist_centroid",
  "dist_fossil"
)

test.statistics



### NOTES: which summary statistics to use??

## use just coefficient 1 (abundance by itself) 
## timing of extinction is critical
## play around with both distances but don't use both...
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

prepareForABC(all_nona,subset=c(1:6),standardized = TRUE)

selectstats <- abctools::selectsumm(test.statistics.obs, par.sim, 
                                     stat.sim, ssmethod =mincrit,  # AS.select,    #  
                                     final.dens = FALSE)

# which summary statistics were selected by the "selectsumm" algorithm?   # # best subset is variables 1 and 3?

colnames(selectstats$best) <- test.statistics

selectstats$best

#### make final selections...

test.statistics


#prepareForABC(all,subset=c(1,2,4,6))


   ## here we subset the targets to the best set identified above
prepareForABC(all,subset=c(5,6),standardized=T)


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
rej <- abc(target=stat.obs, param=par.sim, sumstat=stat.sim, tol=.01, 
           method ="neuralnet", numnet=20, sizenet=10,trace=T,maxit=1000,lambda = c(0.0001,0.001,0.01))   # rejection  # loclinear

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

cv.res.reg <- cv4abc(data.frame(Density=par.sim$DENSITY,Allee=par.sim$ALLEE,Jultemp=par.sim$jul_temp_median), stat.sim,
                     nval=20, tols=c(0.01,0.05),   # method="rejection")
                     method ="neuralnet", numnet=10, sizenet=10,trace=F,maxit=500,lambda = c(0.0001,0.001,0.01)) 
plot(cv.res.reg)

summary(cv.res.reg)   # good summary of model performance



#############################
# VISUALIZE ABUNDANCE TRAJECTORIES AND EXTINCTION TIMES
#############################

### find the best simulations to explore further

if(length(test.statistics)>1){
  ndx <- as.numeric(rownames(rej$ss))
  ranking <- order(apply(as.matrix(rej$ss,min),1,sum))
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


