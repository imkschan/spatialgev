startTime <- Sys.time()

# if (getwd()=="/Users/kschan" || getwd()=="/home/kschan" ) setwd("~/research/thesis-cpp/proj3-jasa/")
source('loadlibrary.R')
library(maps)
library(fields)
library(rworldmap)
library(condMVNorm)

load('dlirwexp.RData')
load('real_gauss_2015_bhm_kn.RData')

# The pollutant considered to print the contour plot
# k = 1,2,3,4,5 to specify the pollutant type to be considered. 1: nitrogen dioxide, 2: sulfur dioxide, 3: PM10, 4: PM2.5, 5: ozone
k <- 5

# draw the map
newmap <- getMap(resolution='high')
spaceLengthY <- 12	# length of the space
spaceLengthX <- 16	# length of the space
gridResolutionX <- 16 * 2
gridResolutionY <- 12 * 2
npredsite <- gridResolutionX * gridResolutionY
site <- t(siteloc)
site[,1]<-(site[,1]-22.2)*10
site[,2]<-(site[,2]-113.0)*10
n <- 1
tt <- 5
giY <- getGridIndex(size=gridResolutionY, length=spaceLengthY) # + 22.2 * 10
giX <- getGridIndex(size=gridResolutionX, length=spaceLengthX) # + 113 * 10
predsite <- as.matrix(expand.grid(giY,giX))
ts <- nrow(predsite) + nrow(site)
s.bin <- 1:nrow(site)
p.bin <- (nrow(site)+1):ts
allsite <- rbind(site,predsite)
bstep <- min(bin[[3]])
estep <- max(bin[[3]])

# control parameters
icor20 <- array(dim=c(30000,1,npredsite,5))
sx <- 21
snx <- 20


# calculate the correlation
sitecor <- function(x,vec,gvar,gsc){
	tmp <- numeric(nrow(vec))
	for (y in 1:nrow(vec)){
		tmpdist <- sqrt( (x[1]-vec[y,1])^2 + (x[2]-vec[y,2])^2 )
		# tmp[y] <- (gvar * exp(-tmpdist/gsc)) / gvar
		tmp[y] <- exp(-tmpdist/gsc)
	}
	return(tmp)
}

cat('\nStart loading ...\n\n')
for (k in 1:5){
	cat('\n\nk =',k,'\ni = ')
	for (i in bstep:estep){
		if (i %% 1000 == 1) cat(i,' ')
		tmpSigma <- get_sitecor(site[sx,],predsite,gvar[i,k],gsc[i,k])
		icor20[i,1,,k] <- tmpSigma
		# ig.tmp1 <- matrix(tmp1,nrow=gridResolutionX,ncol=gridResolutionY)
	}
}

a <- numeric(npredsite)
for (s in 1:npredsite){
	a[s] <- mean(icor20[bstep:estep,1,s,k])
	# a[s] <- quantile(rl1[bstep:estep,1,s,k],0.975)-quantile(rl1[bstep:estep,1,s,k],0.025)
}
# # quilt.plot(predsite[,2],predsite[,1],a)


# draw the plot
realsite <- cbind(site[,1]/10+22.2,site[,2]/10+113)
aa<-matrix(a,nrow=gridResolutionY,ncol=gridResolutionX)
# image.plot(giX/10+113,giY/10+22.2,t(aa),col=terrain.colors(100), xlab='lon',ylab='lat',zlim=c(0,1))
contour(giX/10+113,giY/10+22.2,t(aa), xlab='lon',ylab='lat',labcex=1,lwd=1.2) # nlevels=5
points(cbind(realsite[,2],realsite[,1]),pch=16,col='black',cex=1)
plot(newmap,xlim=c(113,114.6),ylim=c(22.2,23.4),add=TRUE,border=grey(0.7))
