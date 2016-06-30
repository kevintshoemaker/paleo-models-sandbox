
#################
### This script takes the old mammoth template and creates a new template with the correct
### number of populations and parameters that don't need to change from model to model
#################

#### read in template MP file
setwd(MP_DIRECTORY)

MP_TEMPLATE_OLD <- "Mammoth_7000_demo.mp"

mpfilein <- MP_TEMPLATE_OLD
template <- mp.read(mpfilein)   # read in template

## check to make sure the number of populations is correct- if not, correct it!
NPOPSold <- nrow(template$mp.file$PopData_df)
if(NPOPSold<NPOPS){
  #### extract the main population arrays (to be edited)
  poplist1 <- template$mp.file$PopList
  poplist2 <- template$mp.file$PopData_df 
  # start the new mp file
  version <- "60" #template$version
  newmp <- template$mp.file
  
  newmp$PopData_df$name <- as.character(newmp$PopData_df$name)
  
  templatePopLine <- poplist2[1,]
  templatePopList <- poplist1[[1]]
  
  newIDs <- paste("Pop",(NPOPSold+1):NPOPS)  
  counter=1 
  i=(NPOPSold+1)
  for(i in (NPOPSold+1):NPOPS){
    newmp$PopList[[i]] <- templatePopList 
    newmp$PopList[[i]]$name <- newIDs[counter]
    newmp$PopList[[i]]$X_coord <- GridCellAttributes$x.cord_modified[i]
    newmp$PopList[[i]]$Y_coord <- GridCellAttributes$y.cord[i]
    newmp$StInit <- rbind(newmp$StInit,newmp$StInit[1,])
    newmp$PopData_df <- rbind(newmp$PopData_df,templatePopLine)
    newmp$PopData_df[i,]$name <- newIDs[counter]
    newmp$PopData_df[i,]$X_coord <- GridCellAttributes$x.cord_modified[i]
    newmp$PopData_df[i,]$Y_coord <- GridCellAttributes$y.cord[i]  
    counter=counter+1
  }
  # tail(newmp$PopData_df)
  # write the MP file
  
  newmp$PopData_df$X_coord <- GridCellAttributes$x.cord_modified
  newmp$PopData_df$Y_coord <- GridCellAttributes$y.cord 
  
  setwd(MP_DIRECTORY)                                               
  mp.write(newmp,version,"NewMammothTemplate.mp")                   # write MP file
  MP_TEMPLATE <- "NewMammothTemplate.mp" 
  
  #### read in template MP file
  setwd(MP_DIRECTORY)
  mpfilein <- MP_TEMPLATE  #"Mammoth_7000_demo.mp"
  template <- mp.read(mpfilein)   # read in template
}
