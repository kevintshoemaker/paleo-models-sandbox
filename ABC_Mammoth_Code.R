##########################
# This master script performs Approximate Bayesian Analysis
#     algorithms and goodness-of-fit testing for the Paleo simulations
#         USES THE ABC PACKAGE (for post-hoc elimination-style ABC)
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
library(ggmap)
library(mapproj)
library(rworldmap)
library(abctools)    # for fancy calibration of the posterior - for later? 

####################
# LOAD SIMULATION RESULTS
###################


wd <- "E:\\Dropbox\\Damien Fordham\\Mammoth Model\\Kevin\\MergedInputs\\"

wd_results <- "E:\\MammothResults\\ABCAnalysisInputs\\"

setwd(wd)
dat1 <- read.csv(sprintf("%sNicheBreadth40_ABC_data_revised.csv",wd), header = T,stringsAsFactors = FALSE)
dat2 <- read.csv(sprintf("%sNicheBreadth50_ABC_data_revised.csv",wd), header = T,stringsAsFactors = FALSE)
dat3 <- read.csv(sprintf("%sNicheBreadth60_ABC_data_revised.csv",wd), header = T,stringsAsFactors = FALSE)
dat4 <- read.csv(sprintf("%sNicheBreadth70_ABC_data_revised.csv",wd), header = T,stringsAsFactors = FALSE)
dat5 <- read.csv(sprintf("%sNicheBreadth80_ABC_data_revised.csv",wd), header = T,stringsAsFactors = FALSE)
dat6 <- read.csv(sprintf("%sNicheBreadth90_ABC_data_revised.csv",wd), header = T,stringsAsFactors = FALSE)

all <- rbind(dat1,dat2,dat3,dat4,dat5,dat6)
nrow(all)

head(all)


### remove all sims that never went extinct or extinct in first time step

#keep <- !is.na(all$extinction_yrs)

#length(which(keep))

#all <- all[keep,]     # remove the ca. 55k simulations that were not plausible. 

  # explore any remaining NAs in the data frame

#remaining.nas <- all[which(is.na(all),arr.ind=TRUE),]

#nrow(remaining.nas)     # 480 simulations with NAs for dist_centroid and dist_fossil

####################
# LOAD REAL WORLD DATA
####################

    # load data on effective population size from genetic data
gendat <- read.csv("E:\\MammothResults\\ABCAnalysisTargets\\Genetic\\Mammoth_Ne.csv", header = T)

   # load data on extinction date from the fossil record, interpolated to entire range
#extdat <- read.csv("E:\\MammothResults\\ABCAnalysisTargets\\Fossil\\interprolated_estimate_xyz.csv", header = T)

  # load data on extinction date from the fossil record, only fossil locations... 
extdat <- read.csv("E:\\MammothResults\\ABCAnalysisTargets\\Fossil\\Barnoski_matched_estimates_fossils.csv", header = T)
extdat_fossilsites <- extdat[extdat$Fossil.Data==1,]

last.locations <- extdat[which(extdat$Ext==min(extdat$Ext)),c(1,2)]

names(last.locations) <- c("X","Y")




## plot out the real-world data



####################
# OBSERVED STATISTICS FOR MATCHING WITH SIMULATION RESULTS
###################

test.statistics.obs <- c(64.0765, 11.4448, 0, 0, 0, 0)

names(test.statistics.obs) <- c("abundcoef1","abundcoef2","extinctpattern","extinction_yrs","dist_centroid","dist_fossil")

test.statistics.obs

####################
# SET UP VARIABLEs FOR ABC
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


## use just coefficient 1 (abundance by itself) 

## timing of extinction is critical

## play around with both distances but don't use both...

## pattern by itself

## pattern with extinction?   bring NAs in?





#####################
# VISUALIZE PRIORS (check for uniformity- don't need to run this every time)
#####################

## plot parameters 
visualize.priors <- function(){
  par(mfrow=c(2, 3),ask=TRUE)
  
  for(i in 1:length(estimable.params)){
    hist(all[,estimable.params[i]], breaks = 50,main=estimable.params[i])
  }
}
  
# visualize.priors()


#####################
# PERFORM BASIC VISUALIZATIONS AND DATA CHECKS
#####################
 


   # prepare data for ABC

prepareForABC <- function(data=all,subset=c(1,3)){
  par.sim <<- data[,estimable.params]     # samples run from parameter space 
  stat.sim <<- scale(data[,test.statistics[subset],drop=FALSE])    # statistics to match with the real data
  stat.obs <<- test.statistics.obs[test.statistics[subset]]
}

prepareForABC(all,subset=c(1:6))


head(stat.sim)


####################
# SELECT SUMMARY STATISTICS
####################


### remove all na observations

all_nona <- na.omit(all)

prepareForABC(all_nona,subset=c(1:6))

selectstats <- abctools::selectsumm(test.statistics.obs, par.sim, 
                                     stat.sim, ssmethod =AS.select,
                                     final.dens = FALSE)

selectstats

?selectsumm




prepareForABC(all,subset=c(1,3))   # go back to using all observations.


# best subset is variables 1 and 3...


#####################
# PERFORM PRELIMINARY ABC, VISUALIZATIONS AND DATA CHECKS
#####################

### Using the rejection Method to find the best fits


#####
# ABC analysis
#####
rej <- abc(target=stat.obs, param=par.sim, sumstat=stat.sim, tol=.01, method ="neuralnet")   # rejection  # loclinear


head(rej$adj.values)
head(rej$unadj.values)
rej$weights

rej$dist

######     
 # summary and diagnostic plots
######

plot(rej,param=par.sim)

rej$ss

    ### credible intervals
summary(rej)



    ### cross validation (can take a while with neural net)

cv.res.reg <- cv4abc(data.frame(Density=par.sim$DENSITY,Allee=par.sim$ALLEE), stat.sim,
                     nval=20, tols=c(0.005,0.01), method="rejection")
plot(cv.res.reg)

summary(cv.res.reg)   # good summary of model performance




### find the best simulations to explore further

if(length(test.statistics)>1){
  ndx <- as.numeric(rownames(rej$ss))
  mins <- apply(as.matrix(rej$ss,min),2,min)
}else{
  ndx <- as.numeric(names(rej$ss))
  ndx_min <- as.numeric(names(which(rej$ss == min(rej$ss), arr.ind = TRUE)))[1]
  mins <- apply(as.matrix(rej$ss,min),2,min)
}




#############################
# VISUALIZE ABUNDANCE TRAJECTORIES
#############################

index <- ndx[1]
index <- ndx_min

target <- all[index,]$model
#toPlot <- "ABUNDANCE"
toPlot <- "EXTINCTION_DATE"

visualize <- function(target=target,toPlot=toPlot){
  
  if(toPlot="EXTINCTION_DATE"){
    matches <- as.numeric(unlist(regmatches(target, gregexpr("[[:digit:]]+", target))))
    nichebreadth <- matches[1]
    sample <- matches[2]
    
    
    if(toPlot=="ABUNDANCE") filename <- sprintf("NicheBreadth%i_TotAbund_output_alldata.csv",nichebreadth)
    
    setwd(wd_results)
    getwd()
    
    list.files()
    
    
    abund <- read.csv(filename, header = T)    # if(!"abund"%in%ls(name=.GlobalEnv))  
    
    head(abund[,2000:2010])
    
    time<- rev(seq(from =0, to =80000, by = 25))
    
    timeall <- time[2361:2928]
    timebreak1 <- time[2708:2928]
    timebreak2 <- time[2361:2707]
    
    ## Provide row id of simulation of interest
    abundall <- as.numeric(abund[sample,2362:2929])     
    abundbreak1 <- as.numeric(abund[sample,2709:2929])
    abundbreak2 <- as.numeric(abund[sample,2362:2708])
    
    graphics.off()
    plot(timeall, abundall, pch = 16, cex = 0.2, col = "blue", main = "Change in Ne", 
         ylim=c(0,min(1000000,max(abundall+10000))),xlab = "Time (year)", ylab = "Abundance")   #  ,
    
    break1coef <- round(lm(abundbreak1~timebreak1)$coefficients[2],2) 
    break2coef <- round(lm(abundbreak2~timebreak2)$coefficients[2],2)
    
    abline(lm(abundbreak1~timebreak1))
    abline(lm(abundbreak2~timebreak2))
  
    
    ## Extinction presumed to have occurred by 3000 years ago; and only interested in records after 21,000 years ago
    time2 <- (gendat[61:420,1])
    abund2 <-(gendat[61:420,2])
    data2 <- cbind(time2, abund2)
    
    points(time2, abund2, cex = 0.1, col = "red",pch=20)
  }
  

  if(toPlot="EXTINCTION_DATE"){
    matches <- as.numeric(unlist(regmatches(target, gregexpr("[[:digit:]]+", target))))
    nichebreadth <- matches[1]
    sample <- matches[2]
    
    colname <- sprintf("NicheBreadth%i_LHS_Sample%i.mp",nichebreadth,sample)
    
    filename <- sprintf("NicheBreadth%i_FinalYear_output_alldata_formatted.csv",nichebreadth)
    
    setwd(wd_results)
    
    extinct <- read.csv(filename, header = T)    # if(!"abund"%in%ls(name=.GlobalEnv))  
    
    xcords <- extinct[,1]
    ycords <- extinct[,2]
    head(extinct[,1:10])
    
    thissample <- extinct[,colname]
    
    time<- rev(seq(from =0, to =80000, by = 25))
    
    ext.time <- time[thissample]
    
    exttime_scaled <- ext.time/80000
    
    cols <- colorRamp(c("blue","red"), bias = 1, space = c("rgb"),
                      interpolate = c("linear"), alpha = FALSE)(exttime_scaled)
    
    
    #map <- get_map(location = 'Asia', zoom = 3)
    newmap <- getMap(resolution = "low")
    #ggmap(map)
    
    
    
    par(mfrow=c(2,1))
    plot(newmap,xlim = c(20, 180),
         ylim = c(30, 71),
         asp = 1,
         main="simulated"
         )
    points(xcords,ycords,pch=20,cex=0.01,col=cols)
    
    
    plot(newmap,xlim = c(20, 180),
         ylim = c(30, 71),
         asp = 1,
         main="real"
    )
    cols <- colorRamp(c("blue","red"), bias = 1, space = c("rgb"),
                      interpolate = c("linear"), alpha = FALSE)(extdat$Ext/80000)
    
    points(extdat$Mammoth_estimates.CI.Barnoski._Long,extdat$Mammoth_estimates.CI.Barnoski._Lat,pch=20,cex=0.01,col=cols)    
    
    
    points(last.locations)
    
    
     ### location of final extinction
    ndx <- which(ext.time==min(ext.time,na.rm=T))
    xcords_lastext <- extinct[ndx,1]
    ycoord_lastext <- extinct[ndx,2]
    
    plot(newmap,xlim = c(50, 219),
         ylim = c(30, 71),
         asp = 1,
         main="real"
    )
    points(xcords_lastext,ycoord_lastext,pch=20,cex=3,col="black")
    points(last.locations)
    
  }  
  
}




visualize(target,toPlot)





#####################
# VISUALIZE POSTERIORS
#####################

visualize.posteriors <- function(){
  par(mfrow=c(2, 3),ask=TRUE)
  
  posterior_ndx <- as.numeric(rownames(as.matrix(rej$ss))) 
  all_post <- all[posterior_ndx,]
  for(i in 1:length(estimable.params)){
    hist(all_post[,estimable.params[i]], breaks = 50,main=estimable.params[i])
  }
}

visualize.posteriors()




?abc










all$





## ABC with local linear regression correction without/with correction
## for heteroscedasticity
##

lin <- abc(target=stat.obs, param=par.sim, sumstat=stat.sim, tol=.1, hcorr = FALSE, method = "loclinear", transf=c("none","log"))



linhc <- abc(target=stat.obs, param=par.sim, sumstat=stat.sim, tol=.1, method =
               "loclinear", transf=c("none","log"))




#####################
# from ABC package

## Note: cv4postpr is only for model selection...
          # also "postpr" is only for model selection
          # use "cv4abc" function for our purposes


## Note: the neural network method seems to be the best- local linear does not work at all!



library(abc)
data(human)
cv.modsel <- cv4postpr(models, stat.3pops.sim, nval=50, tol=.01, method="mnlogistic")   # Not applicable
plot(cv.modsel)    # not applicable

stat.italy.sim <- subset(stat.3pops.sim, subset=models=="bott")    # 
cv.res.reg <- cv4abc(data.frame(Na=par.italy.sim[,"Ne"]), stat.italy.sim,
                       nval=200, tols=c(.005,.001), method="loclinear")
plot(cv.res.reg, caption="Ne")

res <- abc(target=stat.voight["italian",], param=data.frame(Na=par.italy.sim [, "Ne"]),
             + sumstat=stat.italy.sim, tol=0.005, transf=c("log"), method="neuralnet")
plot(res, param=par.italy.sim [, "Ne"])













##### TESTS

# does removing nas do anything to the results? it shouldn't!


