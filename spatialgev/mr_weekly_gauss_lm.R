options(digits=5)
# library(CDVine)
# library(evd)
library(openxlsx)
library(lubridate)
# library(SpatialExtremes)
# library(ggmap)
# library(colorRamps)
# library(mtsdi)

startTimeHere <- Sys.time()
# if (getwd()=="/Users/kschan" || getwd()=="/home/kschan" ) setwd("~/research/thesis-cpp/proj3/")
source('loadlibrary.R')
load('dlirwexp.RData')
library(SpatialExtremes)

savefit <- TRUE
fname <- 'real_gauss_2013_lm.RData'
yr <- 2013

all.loc.t <- FALSE

totalType <- tt <- dim(dli1)[1]
nsite <- dim(dli1)[3]

n <- 52
dimcop <- 5
lengthcop <- 10
r <- 1

mdata <- array(dim=c(n,nsite,tt))
for (k in 1:tt){
	mdata[,,k] <- getwm_oneyear(dli1[k,,],dt,yr)
}
site <- t(siteloc)
site[,1]<-(site[,1]-22.2)*10
site[,2]<-(site[,2]-113.0)*10


step <- 30000
cutoff1 <- 3000; cutoff2 <- 12000; 
if (cutoff1==cutoff2) {
		bin1 <- bin2 <- 2:cutoff2; bin3 <- (cutoff2+1):step
	} else {
		if (cutoff2!=step){
			bin1 <- 1:cutoff1; bin2 <- (cutoff1+1):cutoff2; bin3 <- (cutoff2+1):step
		} else {
			bin1 <- 1:cutoff1; bin2 <- bin3 <- (cutoff1+1):cutoff2
		}
	}
if (step<=cutoff1) bin1<-bin2<-bin3<-1:step

bin <- list(bin1,bin2,bin3)

nparm <- 5 * totalType + dimcop
nblock <- 2 * totalType + 1

nb <- n
sinb <- n %/% nb

bloc <- sc <- sh <- array(data=NA,dim=c(step,totalType))
beta1 <- beta2 <- array(data=NA,dim=c(step,totalType))
copparm <- array(data=NA,dim=c(step,lengthcop))
if (all.loc.t==TRUE) loc.t <- array(data=NA,dim=c(step,n,nsite,totalType)) else loc.t <- array(data=NA,dim=c(1,n,nsite,totalType))
likloc.t <- array(data=NA,dim=c(step,totalType))

l.t <- array(data=NA,dim=c(step,nblock))
tlt <- array(data=NA,dim=c(step,1 * totalType + 1))
lt <- array(data=NA,dim=c(1 * totalType + 1))

loc.p <- array(data=NA,dim=c(n,nsite,totalType))
l.p <- array(data=NA,dim=nblock)
likloc.p <- array(data=NA,dim=totalType)

us <- array(data=NA,dim=c(step,nblock)); us[1,] <- rep(0, nblock)
al <- array(data=NA,dim=c(step,nblock)); al[1,] <- rep(0, nblock)
ar <- array(data=NA,dim=c(3,nblock))

uw <- array(data=NA,dim=c(step,totalType)); uw[1,] <- rep(0,totalType)	# for random effect
alw <- array(data=NA,dim=c(step,totalType)); alw[1,] <- rep(0,totalType)	# for random effect
arw <- array(data=NA,dim=totalType)
loc.sample <- array(dim=c(step,totalType))

# kn.bloc <- numeric(totalType)
kn.bbloc <- array(dim=c(3,3,totalType))
# kn.gsigma <- array(dim=c(2,2,totalType))
# kn.loc <- array(dim=c(nsite,nsite,totalType))
# kn.loc.up <- array(dim=c(nsite,nsite,totalType))
kn.gev <- array(dim=c(2,2,totalType))

kn.cop <- array(dim=c(lengthcop,lengthcop))

# loc.ac <- array(0,dim=tt)

form.loc <- loc ~ lat + lon
form.scale <- scale ~ 1
form.shape <- shape ~ 1

for (k in 1:tt){
	# fg <- fgev(mdata[,,k])
	fg <- fitspatgev(mdata[,,k],site,form.loc,form.scale,form.shape)
	
	bloc[1,k] <- fg$param[1]; beta1[1,k] <- fg$param[2]; beta2[1,k] <- fg$param[3]
	sc[1,k] <- fg$param[4]; sh[1,k] <- fg$param[5]
	# gvar[1,k] <- fg$std[1]
}
	# gvar[1,k] <- fg$std[1]
	# gvar[1,] <- c(0.075, 0.02, 0.025, 0.012, 0.011) * 1
	# gsc[1,] <- c(0.1, 0.08, 0.35, 0.15, 0.025) * 10


for (k in 1:totalType){
	# tmpSigma <- getexp(site,gvar[1,k],gsc[1,k])
	tmpmu <- bloc[1,k] + beta1[1,k] * site[,1] + beta2[1,k] * site[,2]
	loc.t[1,,,k] <- matrix(tmpmu,nrow=n,ncol=nsite,byrow=T)
	lt[k] <- llgev1(mdata[,,k],loc.t[1,,,k],sc[1,k],sh[1,k])
	# lt[k+totalType] <- loglikgau2(loc.t[1,,,k],site,tmpmu,tmpSigma)
	}
# tmp <- sapply(1:55,function(x) mean(loc.t[1,,x,1]))
# quilt.plot(site[,2],site[,1],tmp)

mtu <- pgev_mat(matrix(mdata,ncol=tt),matrix(loc.t[1,,,],ncol=tt),sc[1,],sh[1,])
if (length(mtu[is.nan(mtu)])>0){
	cat('mtu at t=1 has NaN!')
	break
}
tmpfit <- fit.gausscopula(mtu)$P
copparm[1,] <- tmpfit[lower.tri(tmpfit)]
# mtu2 <- getmatu(matrix(mdata,ncol=tt),matrix(loc.t[1,,,],ncol=tt),sc[1,],sh[1,])
tmp.sig <- p2P(copparm[1,],dimcop)
lt[6] <- sum(dcopula.gauss(mtu, tmp.sig, log=T))
sumlt <- sum(lt)
tlt[1,] <- lt

if (step > 1) for (i in 2:step){
	if (all.loc.t) loci <- i else loci <- 1
	if (all.loc.t) loci.1 <- i-1 else loci.1 <- 1
	if (i %% 100 == 0) {
		if(i <= cutoff1) bb <- 1 else if (i > cutoff2) bb <- 3 else bb <- 2
		ii <- i-1
		cat('\n',i,' ',sapply(1:nblock, function(x) mean(us[which(bin[[bb]]<=(ii)),x]<=al[which(bin[[bb]]<=(ii)),x])),'\n')
		# cat('',loc.ac,' ',loc.ac/(i*nb),'\n')
		cat('',bloc[ii,],' ',beta1[ii,],' ',beta2[ii,],'\n',sc[ii,],' ',sh[ii,],'\n',copparm[ii,], ' ',sum(tlt[ii,]), ' ',format(Sys.time(),usetz = TRUE),'\n')
	}
	for (k in 1:totalType){

		b <- k
		
		# if (i==2) kn.bloc <- (c(0.075, 0.05, 0.08, 0.07, 0.04)*0.2)^2
		if (i==2) kn.bbloc[,,k] <- (diag(c(bloc[1,k] * 0.1, beta1[1,k] * 0.3, beta2[1,k]) * 0.3)*1.0)^2
		# if (i==(cutoff1+1)) kn.bloc[k] <- var(bloc[bin1,b]) * 2
		if (i==(cutoff2+1)) kn.bbloc[,,k] <- var(cbind(bloc[bin2,k], beta1[bin2,k], beta2[bin2,k]) ) * 1.3
		# if (i==(cutoff3+1)) kn.bbloc[,,k] <- var(cbind(bloc[bin3,k], beta1[bin,k], beta2[bin2,k]) ) * 1.3
		# p.bloc <- mvrnorm(n=1, mu=bloc[i-1,k], Sigma=kn.bloc[k])
		p.bbloc <- mvrnormArma(n=1, mu=c(bloc[i-1,k],beta1[i-1,k],beta2[i-1,k]), Sigma=kn.bbloc[,,k])
		pmu <- p.bbloc[1] + p.bbloc[2]*site[,1] + p.bbloc[3] * site[,2]
		pmatmu <- matrix(pmu,nrow=n,ncol=nsite,byrow=T)
		# tmpmu <- bloc[i-1,k] + beta1[i-1,k] * site[,1] + beta2[i-1,k] * site[,2]

		# llgev1(mdata[,,k],loc.t[1,,,k],sc[1,k],sh[1,k])
		l.p[b] <- llgev1(mdata[,,k],pmatmu,sc[1,k],sh[1,k])
		tmplik <- llgev1(mdata[,,k],loc.t[1,,,k],sc[1,k],sh[1,k])
		al[i,b] <- rdiff(l.p[b],tmplik)
		us[i,b] <- runif(1)
		if (us[i,b] <= al[i,b]) {
			bloc[i,k] <- p.bbloc[1] 
			beta1[i,k] <- p.bbloc[2]
			beta2[i,k] <- p.bbloc[3]
			l.t[i,b] <- l.p[b]
			tmplik <- l.p[b]
			tmpmu <- pmu
			loc.t[1,,,k] <- pmatmu
		} else {
			bloc[i,k] <- bloc[i-1,k]; 
			beta1[i,k] <- beta1[i-1,k]
			beta2[i,k] <- beta2[i-1,k]
			l.t[i,b] <- tmplik
		}

		b <- k + totalType
		if (i==2) kn.gev[,,k] <- diag((c(sc[1,k] * 0.23,sh[1,k] * 0.5) * 0.1 )^2,nrow=2,ncol=2)
		if (i==(cutoff1+1)) {
			# kn.gev[,,k] <- var(cbind(sc[bin1,k], sh[bin1,k])) 
		}
		if (i==(cutoff2+1)) {
			kn.gev[,,k] <- var(cbind(sc[bin2,k], sh[bin2,k])) * 1.1
		}
		p.gev <- mvrnormArma(n=1,mu=c(sc[i-1,k],sh[i-1,k]),Sigma=kn.gev[,,k])
		if (p.gev[1] > 0){
			l.p[b] <- llgev1(mdata[,,k],loc.t[loci,,,k],p.gev[1],p.gev[2])
			tmplik <- llgev1(mdata[,,k],loc.t[loci,,,k],sc[i-1,k],sh[i-1,k])
		} else l.p[b] <- -Inf
		al[i,b] <- rdiff(l.p[b],tmplik)
		us[i,b] <- runif(1)
		if (us[i,b] <= al[i,b]) {
			sc[i,k] <- p.gev[1]; sh[i,k] <- p.gev[2]
			l.t[i,b] <- l.p[b]
		} else {
			sc[i,k] <- sc[i-1,k]; sh[i,k] <- sh[i-1,k]
			l.t[i,b] <- tmplik
		} 

	}

	b <- nblock
	matu <- pgev_mat(matrix(mdata,ncol=tt),matrix(loc.t[loci,,,],ncol=tt),sc[i,],sh[i,])
	if (lengthcop == 10){
		# if (i==2) kn.cop <- diag(x=(copparm[1,]*0.01)^2,nrow=lengthcop,ncol=lengthcop)
		if (i==2) kn.cop <- diag(x=(0.005)^2,nrow=lengthcop,ncol=lengthcop)
		# if (i==(cutoff1+1)) kn.cop[,,k] <- var(copparm[bin1,]) * 0.25
		if (i==(cutoff2+1)) kn.cop <- var(copparm[bin2,]) * 0.25
		p.cop <- mvrnormArma(n=1, mu=copparm[i-1,], Sigma=kn.cop)
		p.cop.sigma <- p2P(p.cop,dimcop)
		l.p[b] <- sum(dcopula.gauss(matu, p.cop.sigma, log=T))
		tmp.sigma <- p2P(copparm[i-1,],dimcop)
		tmplik <- sum(dcopula.gauss(matu, tmp.sigma, log=T))			
	} else if (lengthcop == 1){
		if (i==2) kn.cop <- 0.12 ^ 2
		# if (i==(cutoff1+1)) kn.cop[,,k] <- var(copparm[bin1,]) * 0.25
		if (i==(cutoff2+1)) kn.cop <- var(copparm[bin2,]) * 0.25
		p.cop <- mvrnorm(n=1, copparm[i-1,], Sigma=kn.cop)
		l.p[b] <- loglikgumbelcop(matu,p.cop)
		tmplik <- loglikgumbelcop(matu,copparm[i-1,])
	}

	al[i,b] <- rdiff(l.p[b],tmplik)
	us[i,b] <- runif(1)
	if (us[i,b] <= al[i,b]) {
		copparm[i,] <- p.cop; 
		l.t[i,b] <- l.p[b]
	} else {
		copparm[i,] <- copparm[i-1,]; 
		l.t[i,b] <- tmplik
	}

	# tmp <- sum(sapply(1:3, function(k) llgev1(mdata[,,k],loc.t[i,,,k],sc[i,k],sh[i,k]) + 
	# loglikgau(loc.t[i,,,k],site,bloc[i,k],gvar[i,k],gsc[i,k])))
	# mtu <- pgev_mat(matrix(mdata,ncol=3),matrix(loc.t[i,,,],ncol=3),sc[i,],sh[i,])
	# tlt[i] <- tmp + logliknormcopG(mtu,copparm[i,])

	for (k in 1:totalType){
		# tmpSigma <- getexp(site,gvar[i,k],gsc[i,k])
		# tmpmu <- bloc[i,k] + beta1[i,k] * site[,1] + beta2[i,k] * site[,2]
		tlt[i,k] <- llgev1(mdata[,,k],loc.t[loci,,,k],sc[i,k],sh[i,k])
		# tlt[i,k+totalType] <- loglikgau2(loc.t[loci,,,k],site,tmpmu,tmpSigma)
	}
	mtu <- pgev_mat(matrix(mdata,ncol=tt),matrix(loc.t[loci,,,],ncol=tt),sc[i,],sh[i,])
	tmp.sig <- p2P(copparm[i,],dimcop)
	tlt[i,6] <- sum(dcopula.gauss(mtu, tmp.sig, log=T))

	
}
# sapply(1:nblock,function(x) sapply(1:3,function(y) mean(us[bin[[y]],x]<=al[bin[[y]],x])))
# apply(loc.t,3,mean)
# apply(sc[bin[[3]],],2,mean)
# apply(sh[bin[[3]],],2,mean)
# apply(gvar[bin[[3]],],2,mean)
# apply(gsc[bin[[3]],],2,mean)
# apply(copparm[bin[[3]],],2,mean)
if (savefit==TRUE){
	save(n,yr,dimcop,lengthcop,tt,totalType,nsite,site,
		bloc,loc.t,sc,sh,copparm,
		l.t,likloc.t,us,al,uw,alw,bin,nblock,step,cutoff1,cutoff2,nparm,lt,tlt,file=fname)
	cat('\nSaved in ',fname,'\n')
}


cat('\nTime used:',difftime(Sys.time(),startTimeHere,units='mins'),'minutes \n')