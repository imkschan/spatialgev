which.na <- function(data){
	tmp <- which(is.na(data))
	return(tmp)
}

ListfirstRecord <- function(data,dt,print=FALSE) return(apply(data,2,function(x) dt[which(x>-1)[1]]))

getdm <- function(val,dt,na){
	if (!(na %in% c('all.in','del.all','na.fail'))) stop("Wrong na.action!")
	fac<-as.Date(dt[hour(dt)==0])
	tmp.val <- fac
	tmp.index <- numeric(length(fac))
	
	if (na=='all.in') val[is.na(val)] <- -1
	
	for (i in 1:ncol(val)){
	c.max <- as.numeric(tapply(val[,i],as.Date(dt),max))
		c.index <- as.numeric(tapply(val[,i],as.Date(dt), which.max))-1  # hours (minus 1)

		if (na=='all.in') c.max[c.max==-1] <- NA
		if (na=='all.in') c.index[which(c.index==0 & is.na(c.max)==TRUE)] <- NArm
		if (na=='del.all') c.index[which(is.na(c.max)==TRUE)] <- NA

		tmp.val <- cbind(tmp.val,c.max)
		tmp.index <- cbind(tmp.index,c.index)
	}
	tmp.val <- tmp.val[,-1]
	tmp.index <- tmp.index[,-1]
	colnames(tmp.val) <- colnames(val)
	colnames(tmp.index) <- colnames(val)
	return(list(dm=tmp.val,hours=tmp.index))
}


getwm_oneyear <- function(data,dt,year){
	ind <- which(dt==dt[year(dt)==year])
	val <- data[ind,]
	dw <- dt[ind]
	val <- val[1:8736,]
	dw <- dt[1:8736]
	nsite <- ncol(data)
	
	tmp<-ceiling(yday(dw) / 7)
	tmp.val <- 1:52
	
	for (i in 1:nsite){
		c.max <- as.numeric(tapply(val[,i],tmp,max))
		tmp.val <- cbind(tmp.val,c.max)
	}	
	tmp.val <- tmp.val[,-1]
	colnames(tmp.val) <- colnames(val)
	return(tmp.val)
}


print.n <- function(obj) {
	names(obj) <- NULL
	print(obj)
}

miss.all <- function(data) return(as.matrix(apply(data,2,function(x) length(x[x==-99999])/length(x))))

missall <- function(data) {
	tmp <- numeric(ncol(data))
	for (i in 1:ncol(data)){
		x <- data[1:nrow(data),i]
		tmp[i] <- length(x[x==-99999]) / length(x)
	}
	return(as.data.frame(tmp))
}

negall <- function(data) {
	tmp <- numeric(ncol(data))
	for (i in 1:ncol(data)){
		x <- data[1:nrow(data),i]
		tmp[i] <- length(x[which(x<0)])
	}
	return(as.data.frame(tmp))
}

missall.h <- function(data,dt) {
	tmp <- matrix(nrow=24,ncol=ncol(data))
	for (i in 1:ncol(data)){
		for (h in 0:23){
			ind<-which(hour(dt)==h)
			x <- data[ind,i]
			tmp[(h+1),i] <- length(x[x==-99999]) / length(x)
		}
	}
	return(as.data.frame(tmp))
}


firstIndex <- function(data) return(as.matrix(apply(data,2,function(x) which(x>-1)[1])))
missingStart <- function(data) {
	tI <- firstIndex(data)
	tmp <- numeric(ncol(data))
	for (i in 1:ncol(data)){
		x <- data[tI[i]:nrow(data),i]
		tmp[i] <- length(x[x==-99999]) / length(x)
	}
	return(as.data.frame(tmp))
}

print.missing <- function(data){
	return(cbind(missingStart(data)*100,missing1990(data)*100))
}

plotdensity <- function(x) plot(density(x))

plotline <- function(x,...) plot(x,type='l',...)

plot3line <- function(y){
	df <- data.frame(x=1:length(y[,1]),y1=y[,1],y2=y[,2],y3=y[,3])
	p1<-ggplot(df, aes(x)) + geom_line(aes(y=y1), colour="red")   # first layer
	p2<-ggplot(df, aes(x)) + geom_line(aes(y=y2), colour="blue")   # first layer
	p3<-ggplot(df, aes(x)) + geom_line(aes(y=y3), colour="green")   # first layer
	grid.arrange(p1, p2, p3, nrow = 3)
}

plot5line <- function(y){
	df <- data.frame(x=1:length(y[,1]),y1=y[,1],y2=y[,2],y3=y[,3],y4=y[,4],y5=y[,5])
	p1<-ggplot(df, aes(x)) + geom_line(aes(y=y1), colour="red")   # first layer
	p2<-ggplot(df, aes(x)) + geom_line(aes(y=y2), colour="blue")   # first layer
	p3<-ggplot(df, aes(x)) + geom_line(aes(y=y3), colour="green")   # first layer
	p4<-ggplot(df, aes(x)) + geom_line(aes(y=y4), colour="brown")   # first layer
	p5<-ggplot(df, aes(x)) + geom_line(aes(y=y5), colour="purple")   # first layer
	grid.arrange(p1, p2, p3, p4, p5, nrow = 5)
}

plot5line2 <- function(y){
	df <- data.frame(x=1:length(y[,1]),NO2=y[,1],SO2=y[,2],PM10=y[,3],PM25=y[,4],O3=y[,5])
	p1<-ggplot(df, aes(x)) + geom_line(aes(y=NO2))   # first layer
	p2<-ggplot(df, aes(x)) + geom_line(aes(y=SO2))   # first layer
	p3<-ggplot(df, aes(x)) + geom_line(aes(y=PM10))   # first layer
	p4<-ggplot(df, aes(x)) + geom_line(aes(y=PM25))   # first layer
	p5<-ggplot(df, aes(x)) + geom_line(aes(y=O3))   # first layer
	grid.arrange(p1, p2, p3, p4, p5, nrow = 5)
}


lh <- function(x) lapply(x,head)
llh <- function(x) lapply(x,function(y) lapply(y,head))

multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  require(grid)

  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots == 1) {
    print(plots[[1]])

  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

plotline.imput <- function(dl,wna,i,j,takelog=FALSE){
	# data <- dl[[i]][[1]][,j]
	data <- dl[i,,j]
	if (takelog==TRUE) data <- log(data)
	plotline(data);points(x=wna[[i]][[j]],y=data[wna[[i]][[j]]],col='red')
}

p2P <- function(param, d){
    P <- diag(1, nrow = d)
    P[lower.tri(P)] <- param
    P <- P + t(P)
    diag(P) <- rep.int(1, d)
    return(P)
}

cor2cov = function(corMat, varVec) {
  # test the input
  if (!is.matrix(corMat)) stop("'corMat must be a matrix")
  n = nrow(corMat)
  if (ncol(corMat) != n) stop("'corMat' must be square")
  if (mode(corMat) != "numeric") stop("'corMat must be numeric")
  if (mode(varVec) != "numeric") stop("'varVec must be numeric")
  if (!is.null(dim(varVec))) {
    if (length(dim(varVec)) != 2) stop("'varVec' should be a vector")
    if (any(dim(varVec)==1)) stop("'varVec' cannot be a matrix")
    varVec = as.numeric(varVec) # convert row or col matrix to a vector
  }
  if (!all(diag(corMat) == 1)) stop("correlation matrices have 1 on the diagonal")
  if (any(corMat < -1 | corMat > +1)) 
    stop("correlations must be between -1 and 1")
  if (any(varVec <= 0)) stop("variances must be non-negative")
  if (length(varVec) != n) stop("length of 'varMat' does not match 'corMat' size")

  # Compute the covariance
  sdMat = diag(sqrt(varVec))
  rtn = sdMat %*% corMat %*% t(sdMat)
  if (det(rtn)<=0) warning("covariance matrix is not positive definite")
  return(rtn)
}


ratiodiff <- function(up,down,takeexp=TRUE){
	tmp <- up - down
	if (takeexp==TRUE) tmp <- exp(tmp)
	if (is.nan(tmp)==TRUE) tmp <- 0
	return(tmp) 
}

ratiodiff <- function(up,down,takeexp=TRUE){
	tmp <- up - down
	if (takeexp==TRUE) tmp <- exp(tmp)
	if (is.nan(tmp)==TRUE) tmp <- 0
	return(tmp) 
}

