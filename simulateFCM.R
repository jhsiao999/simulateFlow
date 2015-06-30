############################################################################
############################################################################
###                                                                      ###
### This file provides R code to reproduce a simulated                   ###
### flow cytometry data set that mimics a real flow                      ###
### data set.                                                            ###
###                                                                      ###
###                                                                      ###
### Law et al. "Mapping cell populations in flow cytometry               ###
### data for cross-sample comparison using the Friedman-Rafsky test      ###
### statistic as a distance measure", Cytometry A, 2015.                 ###
###                                                                      ###
############################################################################
############################################################################

## Use R package sn and MBESS to simulate multivariate skew t distributions
library(sn)
library(MBESS)


## Import a real flow cytometry data set
realData <- read.table("./sample7.txt", sep="\t", header=TRUE)


## Include 4 marker channels: CD14, CD23, CD3, and CD19, as well
## as the cell population identifying label
realData <- realData[,c(1:4,6)]
colnames(realData)[5] <- "id"


## Recode cell population ID numbers to a sequence of 
## 1 to 9
iiChange <- realData$id > 5
realData$id[iiChange] <- realData$id[iiChange] - 1


## Specificy mean and covariance parameters of each 
## cell populations. These parameters are also described
## in the paper. 

# CP1 covariance
covarr_pop1 <- .8*cov( realData[realData$id==1,1:4] )

# CP2 covariance
covarr_pop2 <- matrix(c(3500,0,0,0,0,18000,0,0,0,0,
                        6100,0,0,0,0,1200), 
                      nrow = 4, ncol = 4, byrow = FALSE)

# CP3 covariance
rho12 = 0.9
rho13 = 0.9
rho14 = 0.2
rho23 = 0.7
rho24 = 0.2
rho34 = 0.2
corrmat <- matrix(c(1,rho12,rho13,rho14,rho12,1,rho23,
	rho24,rho13,rho23,1,rho34,rho14,rho24,rho34,1), 4, 4)
covarr_pop3 <- cor2cov(corrmat, sd=.4*c(sqrt(3862),
		sqrt(2627),sqrt(4672),sqrt(9948)))

# CP4 covariance
covarr_pop4 <- cov( realData[realData$id==4,1:4] )

# CP5 covariance
rho12 = -0.1
rho13 = -.2
rho14 = .1
rho23 = 0.2
rho24 = -0.2
rho34 = -0.2
corrmat <- matrix(c(1,rho12,rho13,rho14,rho12,1,rho23,
                    rho24,rho13,rho23,1,rho34,rho14,rho24,rho34,1),4,4)
covarr_pop5 <- cor2cov(corrmat, sd=.4*c(sqrt(2300),sqrt(2700),
                                        sqrt(2300),sqrt(2600)))
# CP6 covariance
covarr_pop6 <- cov( realData[realData$id==6,1:4] )

# CP7 covariance
covarr_pop7 <- cov( realData[realData$id==7,1:4] )

# CP8 covariance
rho12 = .6
rho13 = .6
rho14 = .3
rho23 = .6
rho24 = .3
rho34 = .3
corrmat <- matrix(c(1,rho12,rho13,rho14,rho12,1,rho23,
	rho24,rho13,rho23,1,rho34,rho14,rho24,rho34,1),4,4)
covarr_pop8 <- cor2cov(corrmat, sd=.4*c(sqrt(2440),sqrt(2080),
	sqrt(2552),sqrt(1606)))

# CP9 covariance
rho12 = .6
rho13 = .5
rho14 = .5
rho23 = .7
rho24 = .4
rho34 = .4
corrmat <- matrix(c(1,rho12,rho13,rho14,rho12,1,rho23,
                    rho24,rho13,rho23,1,rho34,rho14,rho24,rho34,1),4,4)
covarr_pop9 <- cor2cov(corrmat, sd=.4*c(sqrt(9444),
                                        sqrt(12856),sqrt(9483),sqrt(5423)))


## Coerce the cell population paramters of mean, covariance, and skewness 
## into a list object 

CP1 <- list(xi = c(113,129,190,300),
		Omega = covarr_pop1,
		alpha = c(0,0,0,-4),df = 1 )
CP2 <- list(xi = c(68,585,95,497), 
		Omega = covarr_pop2, 
		alpha = c(0,-2.5,0,0), df = 3 )
CP3 <- list(xi = c(392,464,416,563), 
            Omega = covarr_pop3, 
            alpha = c(0,0,0,0), df = 3 )
CP4 <- list(xi = c(214,277,270,240), 
		Omega = .7*covarr_pop4, 
		alpha = c(0,0,-7,0), df = 3 )
CP5 <- list(xi = c(186,248,252,276), 
            Omega = covarr_pop5, 
            alpha = c(0,0,0,3), df = 2.5 )
CP6 <- list(xi = c(670,191,218,262), 
		Omega = .8*covarr_pop6, 
		alpha = c(-6,0,0,1.5), df = 3 )
CP7 <- list(xi = c(98,157,650,260), 
		Omega = .5*covarr_pop7, 
		alpha = c(0,-1,-8,-2), df = 3 )
CP8 <- list(xi = c(238,303,290,293), 
		Omega = covarr_pop8, 
		alpha = c(3,3,2,2), df = 3 )
CP9 <- list(xi = c(530,592,496,412), 
            Omega = 2*covarr_pop9, 
            alpha = c(0,0,0,0), df = 4 )

paramlist <- list(CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9)




