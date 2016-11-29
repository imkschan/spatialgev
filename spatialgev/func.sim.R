getGridIndex <- function(size,length){
	tmp <- seq(0,length,length/size)
	gridIndex <- tmp[-(length(tmp))] + length/size/2
	if (length <= 1) gridIndex <- (seq(0,length*10-1,length*10/size) + length*10/size/2)/10
	# site <- cbind(gridIndex,gridIndex)
	# colnames(site)<-c("lon","lat")
	return(gridIndex)
}

# getIRGridIndex <- function(sizeY,sizeX,length){
# 	tmpX <- seq(0,length,length/sizeX)
# 	tmpY <- seq(0,length,length/sizeY)
# 	gridIndexX <- tmpX[-(length(tmpX))] + length/sizeX/2
# 	gridIndexY <- tmpY[-(length(tmpY))] + length/sizeY/2
# 	if (length <= 1) gridIndexX <- (seq(0,length*10-1,length*10/sizeX) + length*10/sizeX/2)/10
# 	if (length <= 1) gridIndexY <- (seq(0,length*10-1,length*10/sizeY) + length*10/sizeY/2)/10
# 	# site <- cbind(gridIndex,gridIndex)
# 	# colnames(site)<-c("lon","lat")
# 	return(cbind(gridIndexY, gridIndexX))
# }

genRandomSite <- function(n.site,length){
    	# set.seed(seed)
    	site <- matrix(runif(2*n.site,0,length),ncol=2)
    	colnames(site)<-c("lon","lat")
    	sim.site <- site
    }

getexp <- function(site,var,scale){
	nsite <- nrow(site)
	tmpdist <- as.matrix(dist(site,diag=TRUE,upper=TRUE))
	getexp <- var * exp(-tmpdist/scale)
	return(getexp)
}

getexpFromDist <- function(sitedist,var,scale){
	tmp <- var * exp(-sitedist/scale)
	return(tmp)
}

getdist <- function(site){
	tmp <- as.matrix(dist(site,diag=TRUE,upper=TRUE))
	return(tmp)
}