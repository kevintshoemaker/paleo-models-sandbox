######################
### KTS: Code for extracting more results from the mp files: tabled for now...


# XY coords of all occupied sites
occndx <- sapply(as.data.frame(SimInfo$PopAbund),function(t) which(t>1))  # indices of occupied populations for each year
proj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")


# RESULT: MCP OVER TIME
for(t in 1:(SimInfo$ExtinctionYear-1)){
  ndx <- occndx[[t]]
  if(length(ndx)>=5){       # At least 5 relocations are required to fit an home range   
    df <- data.frame(x=numeric(length(ndx)),y=0)
    df$x <- GridCellAttributes$x.cord[ndx]        # x and y coordinates for all occupied grid cells
    df$y <- GridCellAttributes$y.cord[ndx]
    df <- SpatialPoints(df,proj4string=proj)
    MCP <- mcp(df)
    MCPpoly <- MCP@polygons[[1]]@Polygons[[1]]@coords
    MCParea <- areaPolygon(MCPpoly)/1e6    # area of the MCP, in km2
    SimInfo$MCPs[[t]] <- MCP            # store the MCP as SpatialPolygonsDataFrame object
    SimInfo$RangeArea[t] <- MCParea     # store MCP area for each year of the simulation
  }
  # RESULT: CELLS OCCUPIED OVER TIME
  SimInfo$CellsOccupied[t] <- length(ndx)
  
  # RESULT: OCCUPIED AREA OVER TIME AND CELLS OCCUPIED OVER TIME
  areaVec <-  GridCellAttributes$Area[ndx]
  SimInfo$AreaOccupied[t] <- sum(GridCellAttributes$Area2[ndx])

}		
