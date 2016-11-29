startTime <- Sys.time()

if (getwd()=="/Users/kschan" || getwd()=="/home/kschan" ) setwd("~/research/thesis-cpp/proj3-jasa/")
source('loadlibrary.R')
library(maps)
library(fields)
library(rworldmap)
library(condMVNorm)
# load('dlirwexp.RData')
load('real_gauss_2015_bhm_kn.RData')
load('retlv_2015_bhm_kn.RData')



newmap <- getMap(resolution='high')

savefile <- FALSE
savefilename = 'pd_2015_bhm_kn.RData'

# the pollutants considered to draw the predictive distribution
# k = 1,2,3,4,5 to specify the pollutant type to be considered. 1: nitrogen dioxide, 2: sulfur dioxide, 3: PM10, 4: PM2.5, 5: ozone
k_choice <- 1


# init random number generator
tmp <- runif(n=100)

realsite <- cbind(site[,1]/10+22.2,site[,2]/10+113)
pd1 <- array(dim=c(step,n,npredsite,tt))

# calculate the quantities
for (i in bstep:estep){
	tmpU <- rCopula(1, normalCopula(param=copparm[i,],dim=5,dispstr='un'))
	if (i %% 250 == 1) cat(i,' ')
	for (j in 1:1){
		for (k in 1:5)	pd1[i,j,,k] <- qgev(tmpU[k], ploc[i,j,,k],sc[i,k],sh[i,k])
		}
	}

if (savefile==TRUE) save(newmap,spaceLengthX,spaceLengthY,gridResolutionX,gridResolutionY,
	npredsite,site,n,tt,giX,giY,predsite,ts,p.bin,s.bin,allsite,
	p1,p2,ploc,rl1,bstep,estep, file=savefilename)

k <- k_choice   # the pollutants to be drew

a <- numeric(npredsite)
for (s in 1:npredsite){
	a[s] <- median(pd1[bstep:estep,1,s,k])
}

# draw the plot
aa<-matrix(a,nrow=gridResolutionY,ncol=gridResolutionX)
image.plot(giX/10+113,giY/10+22.2,t(aa),col=terrain.colors(1000),xlab='lon',ylab='lat')
points(cbind(realsite[,2],realsite[,1]),pch=16,col='black',cex=1)
plot(newmap,xlim=c(113,114.6),ylim=c(22.2,23.4),add=TRUE)


cat('Time used:',difftime(Sys.time(),startTime,units='mins'),'minutes \n')