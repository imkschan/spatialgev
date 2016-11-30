options(digits=5)

startTimeHere <- Sys.time()
# if (getwd()=="/Users/kschan" || getwd()=="/home/kschan" ) setwd("~/research/thesis-cpp/proj3-jasa/")

# load libraries
source('loadlibrary.R')
load('dlirwexp.RData')	# load dataset
library(SpatialExtremes)
library(openxlsx)
library(lubridate)
library(pscl)

yr <- 2015   # Year to be fit. 
savefit <- TRUE # save results or not
fname <- 'real_gauss_2015_bhm_kn.RData' # filename of saved output

#  define control parameters
all.loc.t <- TRUE
totalType <- tt <- dim(dli1)[1]
nsite <- dim(dli1)[3]	# number of sites
n <- 52	# number of weeks
dimcop <- 5	# dim of copula matrix
lengthcop <- 10	# length of copula parameters
r <- 1

# obtain weekly block maxima data
mdata <- array(dim=c(n,nsite,tt))
for (k in 1:tt){
	mdata[,,k] <- getwm_oneyear(dli1[k,,],dt,yr)
}
# rescale the site lon and lat
site <- t(siteloc)
site[,1]<-(site[,1]-22.2)*10
site[,2]<-(site[,2]-113.0)*10

# mcmc setup
step <- 30000	# no. of iteration
cutoff1 <- 3000; cutoff2 <- 12000 # no of burn-in. only cutoff2 is used here.
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

nparm <- 7 * totalType + dimcop
nblock <- 3 * totalType + 1

nb <- n
sinb <- n %/% nb


# define param vectors to save mcmc iterations
# bloc is beta_0 in the paper, sc is scale of gev, sh is shape of gev, 
# gvar is v^2, gsc is lambda, copparm are the vectorized copula parameters
bloc <- sc <- sh <- gvar <- gsc <- array(data=NA,dim=c(step,totalType))
beta1 <- beta2 <- array(data=NA,dim=c(step,totalType))
copparm <- array(data=NA,dim=c(step,lengthcop))
if (all.loc.t==TRUE) loc.t <- array(data=NA,dim=c(step,n,nsite,totalType)) else loc.t <- array(data=NA,dim=c(1,n,nsite,totalType))
likloc.t <- array(data=NA,dim=c(step,totalType))

# likelihood, U[0,1], acceptance ratios in MH algorithm
l.t <- array(data=NA,dim=c(step,nblock))
tlt <- array(data=NA,dim=c(step,2 * totalType + 1))
lt <- array(data=NA,dim=c(2 * totalType + 1))

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

#  kernal of the MH algorithm
kn.bbloc <- array(dim=c(3,3,totalType))
kn.gsigma <- array(dim=c(2,2,totalType))
kn.loc <- array(dim=c(nsite,nsite,totalType))
kn.loc.up <- array(dim=c(nsite,nsite,totalType))
kn.gev <- array(dim=c(2,2,totalType))
kn.cop <- array(dim=c(lengthcop,lengthcop))
load('kn2015.RData')

loc.ac <- array(0,dim=tt)

# set up initial guess of the mcmc 
form.loc <- loc ~ lat + lon
form.scale <- scale ~ 1
form.shape <- shape ~ 1
for (k in 1:tt){
	fg <- fitspatgev(mdata[,,k],site,form.loc,form.scale,form.shape)
	bloc[1,k] <- fg$param[1]; beta1[1,k] <- fg$param[2]; beta2[1,k] <- fg$param[3]
	sc[1,k] <- fg$param[4]; sh[1,k] <- fg$param[5]
}
gvar[1,] <- c(200,15,1500,700,2000)
gsc[1,] <- c(7, 5.5, 15, 18, 12)

for (k in 1:totalType){
	tmpSigma <- getexp(site,gvar[1,k],gsc[1,k])
	tmpmu <- bloc[1,k] + beta1[1,k] * site[,1] + beta2[1,k] * site[,2]
	loc.t[1,,,k] <- mvrnorm(n=n,mu=tmpmu,Sigma=tmpSigma/1000000)
	loc.sample[1,k] <- loc.t[1,1,1,k]
	lt[k] <- llgev1(mdata[,,k],loc.t[1,,,k],sc[1,k],sh[1,k])
	lt[k+totalType] <- loglikgau2(loc.t[1,,,k],site,tmpmu,tmpSigma) 
	}
mtu <- pgev_mat(matrix(mdata,ncol=tt),matrix(loc.t[1,,,],ncol=tt),sc[1,],sh[1,])
if (length(mtu[is.nan(mtu)])>0){
	cat('mtu at t=1 has NaN!')
	break
}
copparm[1,] <- c(0.089, 0.118, 0.119, 0.019, 0.020, 0.044, 0.032, 0.458, 0.108, 0.080)
tmp.sig <- p2P(copparm[1,],dimcop)
lt[11] <- sum(dcopula.gauss(mtu, tmp.sig, log=T))
sumlt <- sum(lt)
tlt[1,] <- lt


# start mcmc
cat('\nIt may take 60 mins or more ...')
if (step > 1) for (i in 2:step){
	if (all.loc.t) loci <- i else loci <- 1
	if (all.loc.t) loci.1 <- i-1 else loci.1 <- 1
	if (i %% 100 == 0) {
		if(i <= cutoff1) bb <- 1 else if (i > cutoff2) bb <- 3 else bb <- 2
		ii <- i-1
		bbin <- which(bin[[bb]]<=(ii))+min(bin[[bb]])-1
		# cat('\n',i,' ',sapply(1:nblock, function(x) mean(us[bbin,x]<=al[bbin,x])),' ',loc.ac/(i*nb),'\n')
		# cat('',bloc[ii,],' ',beta1[ii,],' ',beta2[ii,],'\n',gvar[ii,],' ',gsc[ii,],'\n',loc.sample[ii,],' ',sc[ii,],' ',sh[ii,],'\n',copparm[ii,], ' ',sum(tlt[ii,]), ' ',format(Sys.time(),usetz = TRUE),'\n')
		cat('\n',i,'/', step,' iterations have been done.')
	}
	for (k in 1:totalType){

		# MH with Gibbs step 1: sample beta_0, beta_1, beta_2
		b <- k
		tmpSigma.1 <- getexp(site,gvar[i-1,k],gsc[i-1,k])
		if (i==(cutoff2+1)) kn.bbloc[,,k] <- var(cbind(bloc[bin2,k], beta1[bin2,k], beta2[bin2,k]) ) * 1.3
		p.bbloc <- mvrnormArma(n=1, mu=c(bloc[i-1,k],beta1[i-1,k],beta2[i-1,k]), Sigma=kn.bbloc[,,k])
		pmu <- p.bbloc[1] + p.bbloc[2]*site[,1] + p.bbloc[3] * site[,2]
		tmpmu <- bloc[i-1,k] + beta1[i-1,k] * site[,1] + beta2[i-1,k] * site[,2]
		l.p[b] <- loglikgau2(loc.t[loci.1,,,k],site,pmu,tmpSigma.1)		
		tmplik <- loglikgau2(loc.t[loci.1,,,k],site,tmpmu,tmpSigma.1)
		al[i,b] <- rdiff(l.p[b],tmplik)
		us[i,b] <- runif(1)
		if (us[i,b] <= al[i,b]) {
			bloc[i,k] <- p.bbloc[1] 
			beta1[i,k] <- p.bbloc[2]
			beta2[i,k] <- p.bbloc[3]
			l.t[i,b] <- l.p[b]
			tmplik <- l.p[b]
			tmpmu <- pmu
		} else {
			bloc[i,k] <- bloc[i-1,k]; 
			beta1[i,k] <- beta1[i-1,k]
			beta2[i,k] <- beta2[i-1,k]
			l.t[i,b] <- tmplik
		}

		# MH with Gibbs step 2: sample v^2, lambda
		b <- k + totalType
		if (i==(cutoff2+1)) {
			 kn.gsigma[,,k] <- var(cbind(gvar[bin2,k], gsc[bin2,k]) ) * 0.15
		}
		p.gsigma <- mvrnormArma(n=1,mu=c(gvar[i-1,k],gsc[i-1,k]),Sigma=kn.gsigma[,,k])
		if (p.gsigma[1] > 0 && p.gsigma[2] > 0){
			pSigma <- getexp(site,p.gsigma[1],p.gsigma[2])
			l.p[b] <- loglikgau2(loc.t[loci.1,,,k],site,tmpmu,pSigma) + densigamma(p.gsigma[1],1,1)
			tmplik <- loglikgau2(loc.t[i-1,,,k],site,tmpmu,tmpSigma.1) + densigamma(gvar[i-1,k],1,1)
		} else l.p[b] <- -Inf
		al[i,b] <- rdiff(l.p[b],tmplik)
		us[i,b] <- runif(1)
		if (us[i,b] <= al[i,b]) {
			gvar[i,k] <- p.gsigma[1]; gsc[i,k] <- p.gsigma[2]
			l.t[i,b] <- l.p[b]
			tmpSigma <- pSigma
		} else {
			gvar[i,k] <- gvar[i-1,k]; gsc[i,k] <- gsc[i-1,k]
			l.t[i,b] <- tmplik
			tmpSigma <- tmpSigma.1
		}
	
		# MH with Gibbs step 3: sample latent location variables mu
		if (all.loc.t==TRUE){
			 if (i==(cutoff2+1)) kn.loc[,,k] <- getloccov(loc.t[,,,k],bin2) * 0.008
		} else {
			if ((i-1) %in% bin[[2]] && length(bin2)!=length(bin3)){
				kn.loc.up[,,k] <- updatelcov(kn.loc.up[,,k],loc.t[1,,,k],i-bin[[2]][1])
			}
			if (i==(cutoff2+1)) kn.loc[,,k] <- kn.loc.up[,,k] * 0.1
		}
		for (p in 1:nb){
			x <- p			
			p.loc <- mvrnormArma(n=1,mu=loc.t[loci.1,x,,k],Sigma=kn.loc[,,k])
			likloc.p <- llgev1(matrix(mdata[x,,k]),matrix(p.loc),sc[i-1,k],sh[i-1,k]) + 
							loglikgau2(p.loc,site,tmpmu,tmpSigma)
			tmplik <- llgev1(matrix(mdata[x,,k]),matrix(loc.t[loci.1,x,,k]),sc[i-1,k],sh[i-1,k]) +
							loglikgau2(t(loc.t[loci.1,x,,k]),site,tmpmu,tmpSigma)
			aa <- rdiff(likloc.p, tmplik)
			uu <- runif(1)

			if (uu <= aa){
				loc.t[loci,x,,k] <- p.loc
				likloc.t[i,k] <- likloc.p
				loc.ac[k] <- loc.ac[k] + 1
			} else {
				loc.t[loci,x,,k] <- loc.t[loci.1,x,,k]
				likloc.t[loci,k] <- likloc.t[loci.1,k]
			}
		}
		loc.sample[i,k] <- loc.t[loci,1,1,k]

		# MH with Gibbs step 4: sample sigma, xi
		b <- k + 2 * totalType
		if (i==(cutoff2+1)) {
			kn.gev[,,k] <- var(cbind(sc[bin2,k], sh[bin2,k])) * 0.3
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

	# MH with Gibbs step 5: sample Omega
	b <- nblock
	matu <- pgev_mat(matrix(mdata,ncol=tt),matrix(loc.t[loci,,,],ncol=tt),sc[i,],sh[i,])
	if (lengthcop == 10){
		if (i==(cutoff2+1)) kn.cop <- var(copparm[bin2,]) * 0.25
		p.cop <- mvrnormArma(n=1, mu=copparm[i-1,], Sigma=kn.cop)
		p.cop.sigma <- p2P(p.cop,dimcop)
		l.p[b] <- sum(dcopula.gauss(matu, p.cop.sigma, log=T))
		tmp.sigma <- p2P(copparm[i-1,],dimcop)
		tmplik <- sum(dcopula.gauss(matu, tmp.sigma, log=T))			
	} else if (lengthcop == 1){
		if (i==2) kn.cop <- 0.12 ^ 2
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

	for (k in 1:totalType){
		tmpSigma <- getexp(site,gvar[i,k],gsc[i,k])
		tmpmu <- bloc[i,k] + beta1[i,k] * site[,1] + beta2[i,k] * site[,2]
		tlt[i,k] <- llgev1(mdata[,,k],loc.t[loci,,,k],sc[i,k],sh[i,k])
		tlt[i,k+totalType] <- loglikgau2(loc.t[loci,,,k],site,tmpmu,tmpSigma)
	}
	mtu <- pgev_mat(matrix(mdata,ncol=tt),matrix(loc.t[loci,,,],ncol=tt),sc[i,],sh[i,])
	tmp.sig <- p2P(copparm[i,],dimcop)
	tlt[i,11] <- sum(dcopula.gauss(mtu, tmp.sig, log=T))

	
}
# sapply(1:nblock,function(x) sapply(1:3,function(y) mean(us[bin[[y]],x]<=al[bin[[y]],x])))
# apply(loc.t,3,mean)

if (savefit==TRUE){
	save(n,yr,dimcop,lengthcop,tt,totalType,nsite,site,mdata,
		bloc,loc.sample,loc.t,sc,sh,gvar,gsc,copparm,beta1,beta2,
		l.t,likloc.t,us,al,uw,alw,loc.ac,bin,nblock,step,cutoff1,cutoff2,nparm,lt,tlt,file=fname)
	cat('\nSaved in ',fname,'\n')
}

# Print table 5.
table5_mean <- matrix(NA, 12, 5)
table5_mean[1,] <- apply(bloc[bin[[3]],],2,mean)
table5_mean[2,] <- apply(beta1[bin[[3]],],2,mean)
table5_mean[3,] <- apply(beta2[bin[[3]],],2,mean)
table5_mean[4,] <- apply(sc[bin[[3]],],2,mean)
table5_mean[5,] <- apply(sh[bin[[3]],],2,mean)
table5_mean[6,] <- apply(gvar[bin[[3]],],2,function(x) mean(sqrt(x)))
table5_mean[7,] <- apply(gsc[bin[[3]],],2,mean)
table5_mean[8:12,] <- p2P(apply(copparm[bin[[3]],],2,mean),5)
colnames(table5_mean) <- c("NO2", "SO2", "PM10", "PM2.5", "O3")
rownames(table5_mean) <- c("beta_0", "beta_lat*", "beta_lon*", "sigma", "xi", "v", "lambda" ,"Omega_NO2:", "Omega_SO2:", "Omega_PM10:", "Omega_PM2.5:", "Omega_O3:")
print("posterior mean for table 5:")
print(table5_mean)

table5_lowerCI <- matrix(NA, 12, 5)
table5_lowerCI[1,] <- apply(bloc[bin[[3]],],2,function(x) quantile(x, 0.025))
table5_lowerCI[2,] <- apply(beta1[bin[[3]],],2,function(x) quantile(x, 0.025))
table5_lowerCI[3,] <- apply(beta2[bin[[3]],],2,function(x) quantile(x, 0.025))
table5_lowerCI[4,] <- apply(sc[bin[[3]],],2,function(x) quantile(x, 0.025))
table5_lowerCI[5,] <- apply(sh[bin[[3]],],2,function(x) quantile(x, 0.025))
table5_lowerCI[6,] <- apply(gvar[bin[[3]],],2,function(x) quantile(sqrt(x), 0.025))
table5_lowerCI[7,] <- apply(gsc[bin[[3]],],2,function(x) quantile(x, 0.025))
table5_lowerCI[8:12,] <- p2P(apply(copparm[bin[[3]],],2,function(x) quantile(x, 0.025)),5)
colnames(table5_lowerCI) <- c("NO2", "SO2", "PM10", "PM2.5", "O3")
rownames(table5_lowerCI) <- c("beta_0", "beta_lat*", "beta_lon*", "sigma", "xi", "v", "lambda" ,"Omega_NO2:", "Omega_SO2:", "Omega_PM10:", "Omega_PM2.5:", "Omega_O3:")
print("Lower 2.5% CI for table 5:")
print(table5_lowerCI)

table5_upperCI <- matrix(NA, 12, 5)
table5_upperCI[1,] <- apply(bloc[bin[[3]],],2,function(x) quantile(x, 0.025))
table5_upperCI[2,] <- apply(beta1[bin[[3]],],2,function(x) quantile(x, 0.025))
table5_upperCI[3,] <- apply(beta2[bin[[3]],],2,function(x) quantile(x, 0.025))
table5_upperCI[4,] <- apply(sc[bin[[3]],],2,function(x) quantile(x, 0.025))
table5_upperCI[5,] <- apply(sh[bin[[3]],],2,function(x) quantile(x, 0.025))
table5_upperCI[6,] <- apply(gvar[bin[[3]],],2,function(x) quantile(sqrt(x), 0.025))
table5_upperCI[7,] <- apply(gsc[bin[[3]],],2,function(x) quantile(x, 0.025))
table5_upperCI[8:12,] <- p2P(apply(copparm[bin[[3]],],2,function(x) quantile(x, 0.025)),5)
colnames(table5_upperCI) <- c("NO2", "SO2", "PM10", "PM2.5", "O3")
rownames(table5_upperCI) <- c("beta_0", "beta_lat*", "beta_lon*", "sigma", "xi", "v", "lambda" ,"Omega_NO2:", "Omega_SO2:", "Omega_PM10:", "Omega_PM2.5:", "Omega_O3:")
print("Upper 2.5% CI for table 5:")
print(table5_upperCI)





cat('\nTime used:',difftime(Sys.time(),startTimeHere,units='mins'),'minutes \n')
