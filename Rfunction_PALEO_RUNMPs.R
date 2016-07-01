##########################
# This script runs a single MP file generated for the paleo project 
#  Authors: Kevin Shoemaker and Damien Fordham
#
# Summary: 
#  3 June 2015 -- started scripting
#  15 June 2015 -- tested script, loaded to GitHub
#  1 July 2016  -- modified for parallelization

#TODO  
##########################


RunMPfile <- function(f,masterDF=masterDF,NicheBreadth=NicheBreadth){

  #setwd(MP_DIRECTORY)
  #  for(f in 1:nfiles){     # loop through MP files...  [now paralellizable]
    filename <- as.character(masterDF$MPFilename[f])
  
         ## set up the new folder to store the MP file and associated KCH files... (specifies the niche breadth)
    thisFolder <- sprintf("%s\\Sample_%s\\LHS_Sample%s",MP_DIRECTORY,NicheBreadth,f)
  

  	###########################
  	# RUN MODEL (takes some time!)
  	###########################
  
  	command <- sprintf("START /WAIT \"MammothModelTest\" \"C:\\Program Files\\RAMAS GIS 6\\Metapop.exe\" \"%s\" /RUN=YES",filename)
  	# cat(command)
  
  	setwd(thisFolder)
  	write.table(command,file="Temp.bat",sep="",quote=FALSE,row.names=FALSE,col.names=FALSE)
  	system("Temp.bat")    # invoke the batch file to run the MP file...
  
  	#if(VERBOSE) cat(paste("just finished running file number",f,"\n",sep=" "))

}

###############  END FUNCTION 'RunMPfile'
############
########
####
#






