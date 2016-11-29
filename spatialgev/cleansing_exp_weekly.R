library(openxlsx)
library(lubridate)

startTime <- Sys.time()

if (getwd()=="/Users/kschan" || getwd()=="/home/kschan" ) setwd("~/research/thesis-cpp/proj3-jasa/")
source('func.R')
cat('Sourced func.R\n')
source('func.cleansing.R')
cat('Sourced func.cleansing.R\n')
source('func.sim.R')
cat('Sourced func.sim.R\n')
source('func.mcmc.R')
cat('Sourced func.mcmc.R\n')

tmp.dt <- read.csv('data/raw/DateTime-20130101-20151231.csv',header=FALSE,as.is=TRUE)  # read raw datetime
dt <- as.POSIXlt(as.character(tmp.dt[,1]),format='%d/%m/%Y %H:%M')	# extract date time into R format
dd <- as.Date(dt[hour(dt)==0])	# extract date from date time
mtemp <- read.csv('data/raw/MAVG_TEMP_13_15.csv',header=TRUE)	# read mean average temperature of HK

# create mean temp vector with index date
m2dtemp <- matrix(nrow=length(dd),ncol=3)	
for (i in 1:length(dd)){
	m2dtemp[i,] <- as.numeric(mtemp[which(year(dd[i])==mtemp[,1] & month(dd[i])==mtemp[,2]),])
}

# create mean temp vector with index date and hour
m2htemp <- matrix(nrow=length(dt),ncol=4)
for (i in 1:length(dd)){
	for (j in 1:24){
		m2htemp[(i-1)*24+j,] <- c(j-1,as.numeric(mtemp[which(year(dd[i])==mtemp[,1] & month(dd[i])==mtemp[,2]),]))
	}
}

# read raw pollutants data (chosen determine the conditions, see func.cleansing.R)
lstno2      <- loadenvf('data/raw/AQ_NO2-20130101-20151231.csv',chosen=TRUE)
lstso2     <- loadenvf('data/raw/AQ_SO2-20130101-20151231.csv',chosen=TRUE)
lstrsp     <- loadenvf('data/raw/AQ_RSPMC-20130101-20151231.csv',chosen=TRUE)
lstfsp      <- loadenvf('data/raw/AQ_FSPMC-20130101-20151231.csv',chosen=TRUE)
lstoz3      <- loadenvf('data/raw/AQ_O3-20130101-20151231.csv',chosen=TRUE)

lstdl <- list(no2=lstno2,so2=lstso2,rsp=lstrsp,fsp=lstfsp,oz3=lstoz3)
# reduce the data set into non-misaligned locations
lstdl <- commonenvf(lstdl)

dl <- array(dim=c(5,26280,55))
for (i in 1:5){
	dl[i,,] <- lstdl[[i]][[1]]
}
siteloc <- lstdl[[1]][[2]]
rownames(siteloc) <- c('lat','lon')

tt <- dim(dl)[1]
nsite <- dim(dl)[3]

pnames <- c('NO2','SO2','RSP','FSP','O3')
sn <- loadnames(lstdl,1)	# load station names

rm(lstdl,lstno2,lstso2,lstrsp,lstfsp,lstoz3)

wdlna <- which.dl.na(dl)

dli1 <- dl # prepare for imputed dl in level 1 (seasonality)
dlr1 <- dl # prepare for residual of imputed dl in level 1 (seasonality)

fit1 <- array(dim=c(25,nsite,tt))

for (i in 1:tt){
	for (j in 1:nsite){
		cat('i =',i,'  j =',j,'\n')
		data.e <- dl[i,,j]	# init data value (in original unit)
		data.e[data.e==0] <- 0.1    # AD HOC!!!!!!!!!!!!!!!!!!!!!
		data <- data.e		# work on log data.

		formula <- data ~
			I(m2htemp[,1]==1)+I(m2htemp[,1]==2)+I(m2htemp[,1]==3)+I(m2htemp[,1]==4)+
			I(m2htemp[,1]==5)+I(m2htemp[,1]==6)+I(m2htemp[,1]==7)+I(m2htemp[,1]==8)+
			I(m2htemp[,1]==9)+I(m2htemp[,1]==10)+I(m2htemp[,1]==11)+I(m2htemp[,1]==12)+
			I(m2htemp[,1]==13)+I(m2htemp[,1]==14)+I(m2htemp[,1]==15)+I(m2htemp[,1]==16)+
			I(m2htemp[,1]==17)+I(m2htemp[,1]==18)+I(m2htemp[,1]==19)+I(m2htemp[,1]==20)+
			I(m2htemp[,1]==21)+I(m2htemp[,1]==22)+I(m2htemp[,1]==23)+m2htemp[,4]		

		proceed <- 0
		for (r in 1:1000){
			if (proceed < 2){
				fit <- lm(formula=formula,na.action=na.exclude)
				nd <- model.frame(formula,na.action=na.pass)	# prepare for union of response and design matrix
				colnames(nd)<-c('data',paste('h',seq(1,23),sep=''),"temp")
				tmp <- predict(fit,nd)[wdlna[[i]][[j]]]
				fitc <- fit$coefficients
				if (r==1) fitc.o <- fitc
				if (r>1) if (all(abs(fitc-fitc.o)<1e-6)) proceed <- proceed + 1 else proceed <- 0
				# print.n(head(exp(tmp)))
				print.n(head(fit$coefficients))
				data[wdlna[[i]][[j]]] <- tmp
				fitc.o <- fitc
			}
		}
		dli1[i,,j] <- data
		dlr1[i,,j] <- resid(fit)
		fit1[,j,i] <- fitc
		cat('\n')
	}
}

dli2 <- dli1
dlr2 <- dlr1

save(dli1,dlr1,fit1,dt,dd,m2dtemp,m2htemp,dl,pnames,sn,wdlna,siteloc,file='dlirwexp.RData')
cat('Saved results in dlirwexp.RData.\n')
cat('Time used:',difftime(Sys.time(),startTime),'\n')
