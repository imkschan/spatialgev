loadenvf <- function(filename,chosen=FALSE){
	cat('Path:',filename,'is loading... ')
	raw <- read.csv(filename,header=TRUE,skip=5,as.is=TRUE)
	loc <- data.matrix(raw[1:2,2:ncol(raw)])	# loc id names
	val <- data.matrix(raw[4:nrow(raw),2:ncol(raw)])	# pollutant levels
	dt <- as.POSIXct(raw[4:nrow(raw),1])	# date time
	# chose only stations having data since the first day  AND  missing value less than 15%
	if (chosen==TRUE) {
		fR <- firstRecord(val,dt)
		chosen.name <- names(fR[which(yday(fR)==1 & year(fR)==2013)])	
		if (is.null(chosen.name)==FALSE) val <- val[,colnames(val) %in% chosen.name]
		if (is.null(chosen.name)==FALSE) loc <- loc[,colnames(loc) %in% chosen.name]
		delmiss <- which(missall(val)*100>15)
		val <- val[,-delmiss]
		loc <- loc[,-delmiss]
	}
	# replace invalid data by NA
	val[val==-99999] <- NA
	val[val<0] <- NA
	cat('done!\n')
	return(list(val=val,loc=loc))
}

# loadrawenvf <- function(filename){
# 	cat('Path:',filename,'is loading... ')
# 	raw <- read.csv(filename,header=TRUE,skip=5,as.is=TRUE)
# 	loc <- data.matrix(raw[1:2,2:ncol(raw)])	# loc id names
# 	val <- data.matrix(raw[4:nrow(raw),2:ncol(raw)])	# pollutant levels
# 	dt <- as.POSIXct(raw[4:nrow(raw),1])	# date time
# 	# chose only stations having data since the first day  AND  missing value less than 15%
# 	if (chosen==TRUE) {
# 		fR <- firstRecord(val,dt)
# 		chosen.name <- names(fR[which(yday(fR)==1 & year(fR)==2013)])	
# 		if (is.null(chosen.name)==FALSE) val <- val[,colnames(val) %in% chosen.name]
# 		if (is.null(chosen.name)==FALSE) loc <- loc[,colnames(loc) %in% chosen.name]
# 		delmiss <- which(missall(val)*100>15)
# 		val <- val[,-delmiss]
# 		loc <- loc[,-delmiss]
# 	}
# 	# replace invalid data by NA
# 	val[val==-99999] <- NA
# 	val[val<0] <- NA
# 	cat('done!\n')
# 	return(list(val=val,loc=loc))
# }


firstRecord <- function(data,dt,print=FALSE){
	n.col <- ncol(data)
	tmp <- .POSIXct(character(n.col))
	for (i in 1:n.col){
		tmp[i] <- dt[which(data[,i]>-1)[1]]
	}
	names(tmp) <- colnames(data)
	return(tmp)
}

commonenvf <- function(dl){
	nn <- Reduce(intersect,list(colnames(dl[[1]][[1]]),colnames(dl[[2]][[1]]),colnames(dl[[3]][[1]]),colnames(dl[[4]][[1]]),colnames(dl[[5]][[1]])))
	nn <-nn[-which(nn=="CN_1360A.1")]   # delete repeated data
	for (i in 1:5){
		for (j in 1:2){
			dl[[i]][[j]] <- dl[[i]][[j]][,colnames(dl[[i]][[j]]) %in% nn]
		}
	}
	return(dl)
}

loadnames <- function(dl,index){
	tmp.names<-lapply(dl,function(x) colnames(x[[2]]))
	rawn <- read.xlsx('data/stations2014.xlsx',sheet=index)
	tmptn <- tmp.names[[index]]
	# locnames <- array(dim=length(tmptn))
	# locnames <- rawn[which(rawn[,8] %in% tmp.names[[index]]),4]
	locnames <- sapply(1:length(tmptn), function(x) rawn[which(rawn[,8] == tmptn[x]),4])
	return(locnames)
}

which.dl.na <- function(dl){
	tt <- dim(dl)[1]
	nsite <- dim(dl)[3]
	wm <- list()
	for (i in 1:tt){
		tmpls <- list()
		for (j in 1:nsite){
			tmp <- which(is.na(dl[i,,j]))
			tmpls <- cbind(tmpls,list(tmp))
		}
		wm <- cbind(wm,list(tmpls))
	}
	return(wm)
}

which.lstdl.na <- function(mat){
	# tt <- dim(dl)[1]
	nsite <- dim(mat)[2]
	wm <- list()
	# for (i in 1:tt){
		tmpls <- list()
		for (j in 1:nsite){
			tmp <- which(is.na(mat[,j]))
			tmpls <- cbind(tmpls,list(tmp))
		}
		wm <- cbind(wm,list(tmpls))
	# }
	return(wm)
}
