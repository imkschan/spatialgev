startTime <- Sys.time()

if (getwd()=="/Users/kschan" || getwd()=="/home/kschan" ) setwd("~/research/thesis-cpp/proj3/")
source('loadlibrary.R')
library(maps)
library(fields)
library(rworldmap)
load('dlirwexp.RData')
load('real_gauss_2015_bhm_kn.RData')

library(condMVNorm)

# parameters
newmap <- getMap(resolution='high')

savefile <- TRUE
savefilename = 'retlv_2015_bhm_kn.RData'

spaceLengthY <- 12	# length of the space
spaceLengthX <- 16	# length of the space

gridResolutionX <- 16 * 2
gridResolutionY <- 12 * 2
npredsite <- gridResolutionX * gridResolutionY

# init random number generator
tmp <- runif(n=100)

cat('\nThis program can take many hours ... \n')

# generate site 
cat('Generating sites ... ')

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


cat('OK! \n')
cat('Generating conditional latent locations ... ')
p1<-1/52; xtmp1 <- -log(1-p1)
p2<-1/520; xtmp2 <- -log(1-p2)
ploc <- array(dim=c(step,n,npredsite,tt))
rl1 <- array(dim=c(step,n,npredsite,tt))
bstep <- 12001
estep <- 30000
for (k in 1:5){
	cat('\n\nk =',k,'\ni = ')
	for (i in bstep:estep){
		if (i %% 25 == 1) cat(i,' ')
		tmpSigma <- getexpcpp(allsite,gvar[i,k],gsc[i,k])
		tmpmu <- bloc[i,k] + beta1[i,k] * allsite[,1] + beta2[i,k] * allsite[,2]
		for (j in 1:1){
			ploc[i,j,,k] <- rcmvnorm(n=1,mean=tmpmu,sigma=tmpSigma,dep=p.bin,given=s.bin,X=loc.t[i,j,,k],method='chol')
			rl1[i,j,,k] <- ploc[i,j,,k] - sc[i,k]/sh[i,k] * (1-xtmp1^(-sh[i,k]))
		}
	}
}

if (savefile==TRUE) save(newmap,spaceLengthX,spaceLengthY,gridResolutionX,gridResolutionY,
	npredsite,site,n,tt,giX,giY,predsite,ts,p.bin,s.bin,allsite,
	p1,p2,ploc,rl1,bstep,estep, file=savefilename)


cat('Time used:',difftime(Sys.time(),startTime,units='mins'),'minutes \n')