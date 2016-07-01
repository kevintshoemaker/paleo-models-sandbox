### test for running multiple MP files...

foldernames <- c("C:\\MPTEST1","C:\\MPTEST2")

filename <- "TEST.mp"

RunMPFile <- function(i){
  ## set up the new folder to store the MP file and associated KCH files... (specifies the niche breadth)
  thisFolder <- foldernames[i]
  
  
  ###########################
  # RUN MODEL (takes some time!)
  ###########################
  
  command <- sprintf("START /WAIT \"MammothModelTest\" \"C:\\Program Files\\RAMAS GIS 6\\Metapop.exe\" \"%s\" /RUN=YES",MPsToRun[f])
  # cat(command)
  
  setwd(thisFolder)
  write.table(command,file="Temp.bat",sep="",quote=FALSE,row.names=FALSE,col.names=FALSE)
  system("Temp.bat")    # invoke the batch file to run the MP file...

}
#if(VERBOSE) cat(paste("just finished running file number",f,"\n",sep=" "))


all.mps <- foreach(i = 1:nrow(masterDF),
                   .export=objectlist,
                   #.packages = c("R2WinBUGS"),
                   .errorhandling=c("pass")
) %dopar% {   
  RunMPfile(f=i,masterDF=masterDF,NicheBreadth=NicheBreadth)
}

