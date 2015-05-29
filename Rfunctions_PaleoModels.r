
######  Dispersal function...
negexp <- function(distance,alpha, beta) {
   alpha*exp((-1*(distance))/beta)   
}

visualizeDispersalFunction <- function(alpha, beta) {
 graphics.off()
 a = alpha      # dispersal parameters
 b = beta      
 curve(negexp(x,a,b),col="red",lwd=2,0,7)
}

specifyLHSParam <- function(paramslist,name,type,lb,ub){
    newlist <- paramslist
	eval(parse(text=sprintf("newlist$%s <- list()",name)))
	eval(parse(text=sprintf("newlist$%s$type <- \"%s\"",name,type)))
	eval(parse(text=sprintf("newlist$%s$lb <- %s",name,lb)))
	eval(parse(text=sprintf("newlist$%s$ub <- %s",name,ub))) 	
	return(newlist)
}


#### function to set up and write KCH files
#      and... determine initial abundance for each population
#      remember- we will need to scale by area... 

writeKCH <- function(HSnum = 1,density = 1000){
	initabund <- numeric(NPOPS)
	for(i in 1:NPOPS){
	  HSfilename = sprintf("pop%s_rep%s.hs",i,HSnum)
	  KCHfilename = sprintf("pop%s_rep%s.kch",i,HSnum)
	  setwd(HS_DIRECTORY)
	  hs = read.table(HSfilename,header=F)[,1]
	  initabund[i] = hs[1]*density
	  K = round(hs*density)
	  setwd(MP_DIRECTORY)
	  write.table(K,file=KCHfilename,row.names=F,col.names=F)
	}
	return(initabund)
}



########################################
## GENERIC FUNCTION FOR INSTALLING/LOADING PACKAGES FROM CRAN
########################################

loadPackage <- function(pkg){

  if(pkg %in% rownames(installed.packages()) == FALSE) {suppressMessages(suppressWarnings(install.packages(pkg)))}
  eval(parse(text=sprintf("suppressMessages(suppressWarnings(require(%s)))",pkg)), envir= .GlobalEnv)

}


##########################################
## GENERIC FUNCTION FOR PULLING SOURCE CODE FROM GITHUB OR OTHER URL
##########################################

source_github <- function(baseurl,script) {
  # load package
   suppressMessages(suppressWarnings(require(RCurl)))
 
  # read script lines from website
  url <- sprintf("%s%s",baseurl,script)
  script <- getURL(url, ssl.verifypeer = FALSE)
 
  # parse lines and evaluate in the global environement
  eval(parse(text = script), envir= .GlobalEnv)
}



 




