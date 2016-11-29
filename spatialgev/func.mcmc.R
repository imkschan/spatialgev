loglikgau <- function(loc,site,beta.loc,var,scale,neg=FALSE){
	Sigma <- getexp(site,var,scale)
	# tmp <- dmvnorm(loc,mean=rep(beta.loc,nrow(site)),sigma=Sigma,log=TRUE)
	tmp <- dmvnrm_arma(loc,mean=rep(beta.loc,nrow(site)),sigma=Sigma,log=TRUE)
	# tmp[which(is.in cfinite(tmp))] <-  -743
	# tmp[which(is.infinite(tmp))] <-  -100000
	if (neg==TRUE) tmp <- -tmp
	return(sum(tmp))
}

loglikgau1 <- function(loc,site,beta.loc,Sigma,neg=FALSE){
	# tmp <- dmvnorm(loc,mean=rep(beta.loc,nrow(site)),sigma=Sigma,log=TRUE)
	tmp <- dmvnrm_arma(loc,mean=rep(beta.loc,nrow(site)),sigma=Sigma,log=TRUE)
	# tmp[which(is.in cfinite(tmp))] <-  -743
	# tmp[which(is.infinite(tmp))] <-  -100000
	if (neg==TRUE) tmp <- -tmp
	return(sum(tmp))
}

loglikgau2 <- function(loc,site,sploc,Sigma,neg=FALSE){
	# tmp <- dmvnorm(loc,mean=rep(beta.loc,nrow(site)),sigma=Sigma,log=TRUE)
	tmp <- dmvnrm_arma(loc,mean=sploc,sigma=Sigma,log=TRUE)
	# tmp[which(is.in cfinite(tmp))] <-  -743
	# tmp[which(is.infinite(tmp))] <-  -100000
	if (neg==TRUE) tmp <- -tmp
	return(sum(tmp))
}

loglikgev <- function(data,loc,sc,sh,neg=FALSE){
	tmp<-dgev(x=data,loc=loc,scale=sc,shape=sh,log=TRUE)
	if (neg==TRUE) tmp <- -tmp
	return(tmp)
}


# loglikgau(raw.loc[,,k,r],raw.site[,,r],raw.beta.loc[k],raw.gp.var[k],raw.gp.scale[k])
loglikgev1 <- function(data,loc,sc,sh,neg=FALSE){
	tmp<-dgev(x=data,loc=loc,scale=sc,shape=sh,log=TRUE)
	# tmp[which(is.infinite(tmp))] <-  -743
	# tmp[which(is.infinite(tmp))] <-  -100000
	# if (is.infinite(tmp)==TRUE) tmp <- -743
	if (neg==TRUE) tmp <- -tmp
	return(sum(tmp))
}
# loglikgev1(raw[,,k,r],raw.loc[,,k,r],raw.sc[,,k,r],0)

getmatu <- function(mdata,mloc,msc,msh){
	n <- dim(mdata)[1]
	nsite <- dim(mdata)[2]
	tt <- dim(mdata)[3]
	u <- array(dim=c(n,nsite,tt))
	for (k in 1:tt){
		if (length(dim(msc))==3) u[,,k]<-pgev(mdata[,,k],loc=mloc[,,k],scale=msc[,,k],shape=msh[k])
		if (length(dim(msc))==0) u[,,k]<-pgev(mdata[,,k],loc=mloc[,,k],scale=msc[k],shape=msh[k])
	}
	matu <- matrix(u,ncol=tt)
	return(matu)
}

getmatucpp <- function(mdata,mloc,msc,msh){	
	matu <- matrix(mdata,ncol=3)
	return(matu)
}

logliknormcop <- function(matu,rho12,rho13,rho23,neg=FALSE){
	norm.cop <- normalCopula(dim=3,dispstr = "un")
	tmp <- loglikCopula(param=c(rho12, rho13, rho23),matu,norm.cop)	
	# tmp[which(is.infinite(tmp))] <-  -1000000
	if (neg==TRUE) tmp <- -tmp
	return(tmp)
}

logliknormcopG <- function(matu,x,neg=FALSE){
	tmp <- logliknormcop(matu=matu,rho12=x[1],rho13=x[2],rho23=x[3],neg=FALSE)
	return(tmp)
}

logliknormcopfull <- function(mdata,mloc,msc,msh,rho12,rho13,rho23,neg=FALSE){
	n <- dim(mdata)[1]
	nsite <- dim(mdata)[2]
	tt <- dim(mdata)[3]
	u <- array(dim=c(n,nsite,tt))
	for (k in 1:tt){
			u[,,k]<-pgev(mdata[,,k],loc=mloc[,,k],scale=msc[,,k],shape=msh[k])
	}
	matu <- matrix(u,ncol=3)
	norm.cop <- normalCopula(dim=3,dispstr = "un")
	tmp <- loglikCopula(param=c(rho12, rho13, rho23),matu,norm.cop)	
	if (is.infinite(tmp)==TRUE) tmp <- -743
	if (neg==TRUE) tmp <- -tmp
	return(tmp)
}

loglikgumbelcop <- function(matu,theta,neg=FALSE){
	gumbel.cop <- gumbelCopula(dim=3)
	tmp <- loglikCopula(param=theta,matu,gumbel.cop)	
	# tmp[which(is.infinite(tmp))] <-  -1000000
	if (neg==TRUE) tmp <- -tmp
	return(tmp)
}



getloccov <- function(loc.t,bin){
	nsite <- dim(loc.t)[3]
	tmpcov <- array(dim=c(length(bin),nsite,nsite))
	mucov <- array(dim=c(nsite,nsite))
	for (i in (min(bin)):(max(bin))){
		tmpcov[i-min(bin)+1,,] <- cov(loc.t[i,,])
	}
	for (i in 1:nsite){
		for (j in 1:nsite){
			mucov[i,j] <- mean(tmpcov[,i,j])
		}
	}
	return(mucov)
}


updatelcov <- function(covlt,lp,i){
	covlp <- cov(lp)
	if (i > 1){
		tmp <- (i-1)/i * covlt + 1/i * covlp
	} else if (i == 1){
		tmp <- covlp
	}
	return(tmp)
}

# for (i in (min(bin2)):(max(bin2))){
# 	knloc[,,1] <- updatelcov(knloc[,,1],loc.t[i,,,1],i-(min(bin2))+1)
# }


# p2P(c(0.4,0.7,0.2),3)

# benchmark(fgau=loglikgau(raw.loc[,,k,r],raw.site[,,r],raw.beta.loc[k],raw.gp.var[k],raw.gp.scale[k]),fmardata=loglikgev1(raw[,,k,r],raw.loc[,,k,r],raw.sc[,,k,r],0),fcop=logliknormcopfull(raw[,,,r],raw.loc[,,,r],raw.sc[,,,r],raw.beta.sh,0.4,0.7,0.2))
# benchmark(fgau=loglikgau(raw.loc[,,k,r],raw.site[,,r],raw.beta.loc[k],raw.gp.var[k],raw.gp.scale[k]),fmardata=loglikgev1(raw[,,k,r],raw.loc[,,k,r],raw.sc[,,k,r],0),fmatu=getmatu(raw[,,,r],raw.loc[,,,r],raw.sc[,,,r],raw.beta.sh),fcop=logliknormcop(matu,0.4,0.7,0.2))

# logliknormcop(raw[,,,r],raw.loc[,,,r],raw.sc[,,,r],raw.beta.sh,0.4,0.7,0.2)

# loglikgau(raw.loc[,,k,r],raw.site[,,r],raw.beta.loc[k],raw.gp.var[k],raw.gp.scale[k]) + 
# 	loglikgau(raw.loc[,,k,r],raw.site[,,r],raw.beta.loc[k],raw.gp.var[k],raw.gp.scale[k]) + 
# 	loglikgau(raw.loc[,,k,r],raw.site[,,r],raw.beta.loc[k],raw.gp.var[k],raw.gp.scale[k]) +
# 	loglikgev1(raw[,,k,r],raw.loc[,,k,r],raw.sc[,,k,r],0) + 
# 	loglikgev1(raw[,,k,r],raw.loc[,,k,r],raw.sc[,,k,r],0) + 
# 	loglikgev1(raw[,,k,r],raw.loc[,,k,r],raw.sc[,,k,r],0) + 
# 	logliknormcop(raw[,,,r],raw.loc[,,,r],raw.sc[,,,r],raw.beta.sh,0.4,0.7,0.2)