############################################################################
############################################################################
###                                                                      ###
### This file provides R code to reproduce the simulated                 ###
### flow cytometry data sets used in the evaluation of                   ###
### flowMap-FR under scenaiors of biological and technical               ###
### variations. The evaluation results are preseneted in                 ###             ###                                                                      ###
###                                                                      ###
### Hsiao et al. "Mapping cell populations in flow cytometry               ###
### data for cross-sample comparison using the Friedman-Rafsky test      ###
### statistic as a distance measure", Cytometry Part A, 2015.            ###
###                                                                      ###
############################################################################
############################################################################

library(sn)
library(MBESS)

#' Simulate multivariate skew t distribution
#' 
#' This is a modified version of the rmst function in the sn package.
#' The original version ran with an error.
#' 
#' @param n number of events or data points to be drawn from the 
#'        multivariate skew t distribution.
#' @param xi a vector of mean parameters.
#' @param Omega a matrix of covariance paramters.
#' @param df degrees of freedom.
rmst_ed <- function(n = 1, 
                    xi=rep(0,length(alpha)), 
                    Omega, alpha, df=Inf, dp=NULL)
{ 
  if(!(missing(alpha) & missing(Omega)) && !is.null(dp)) 
    stop("You cannot set both component parameters and dp")
  if(!is.null(dp)){
    if(!is.null(dp$xi)) xi <- dp$xi
    else
      if(!is.null(dp$beta)) xi <- as.vector(dp$beta)
      Omega <- dp$Omega
      alpha <- dp$alpha
      df <- dp$df
  }  
  d <- length(alpha)
  x <- if(df==Inf) 1 else rchisq(n,df)/df
  z <- rmsn(n, rep(0,d), Omega, alpha)
  y <- sweep(z/sqrt(x),2,STAT=xi,FUN="+")
  # y <- t(xi+ t(z/sqrt(x)))
  attr(y,"parameters") <- list(xi=xi, Omega=Omega, alpha=alpha, df=df)
  return(y)
}

#' Simulate a single cell population
#' 
#' @param nn number of events in the cell population.
#' @param xi mean a vector of mean parameters.
#' @param Omega a matrix of covariance paramters.
#' @param alpha a scalor vector of skewness paramter.
#' @param df degrees of freedom.
simOnePop <- function(nn,xi,Omega,alpha,df) {
  
  rsam <- rmst_ed(n=nn,xi = xi,Omega=Omega,alpha=alpha,df=df)
  IIout <- rowSums(rsam<0 | rsam>1000) >0
  rsam <- rsam[!IIout,]
  nn <- sum(IIout)
  
  while(nn > 0) {
    temp <- rmst_ed(n=nn,xi = xi,Omega=Omega,alpha=alpha,df=df)
    IIout <- rowSums(temp<0 | temp>1000) >0
    rsam2 <- temp[!IIout,]
    nn <- sum(IIout)
    rsam <- rbind(rsam,rsam2)
  }
  rsam <- data.frame(rsam,check.names=F,row.names=NULL)
  return(rsam)	
}





############################################################################
## Proportion changes
## Helper function to evaluate the scenario under which cell populations
## are changed in proportion and with fixed mean, covariance and 
## skewness parameters.

#' Simulate cell populations with changed proportions and compute
#' FR statistic between the original cell populations and the changed
#' cell populations
#' 
#' @param id identifying number of the cell population that is to be
#'        simulated in varying proportions.
#' @param simulatedData a simulated data set that contains cell population
#'        mimicing the a real flow cytometry data set.
#' @param sampleMethod method used to downsample the number of events in each 
#'        cell population comparison when computing FR statistic.
#' @param ndraws number of random samples that is to be drawn in 
#'        each cell population comparison in order to compute the FR statistic.
#' @param samToMatch a data.frame that contains the cell populations that are
#'        going to be compared with the changed cell populations.
#' @export propcase data.frames of the simulated cases.
#' @export propcase_res data.frames of the results of comparing simulated
#'         cell populations with the target cell population.
make.propcase <- function(id, 
  simulatedata, sampleMethod, ndraws, samToMatch) {
  # sameToMatch: the sample that serves as the reference
  prop_list <- c(1,10,seq(25,150,25))/100
  simprops <- round( 100*table(simulatedata$id)[id]*prop_list/length(simulatedata$id[id]),2)
  propcase <- simSam2(pop=id,param=simprops/100,sam1=simulatedata,
                      xi = paramlist[[id]]$xi, Omega = paramlist[[id]]$Omega, 
                      alpha = paramlist[[id]]$alpha, df = paramlist[[id]]$df) 
  names(propcase) <- simprops
  propcase_res <- lapply(1:length(simprops), function(i) 
    makeDistmat(list(propcase[[samToMatch]],propcase[[i]]),
                sampleMethod = sampleMethod, ndraws = ndraws) )
  list(propcase = propcase, propcase_res = propcase_res)
}


#' Simulate a single cell population under the scenaior of changed propotion.
#' 
#' This is a wrapper function of simOnePop. Simulate a series of cell populations
#' with changes in the number of events (proportional to the original cell population)
#' and fixed mean, covariance, and skewness parameters.
#' 
#' @param pop identifying number of the cell population to be simulated. 
#' @param param a vector that contains a series of number of events that is to be 
#'        drawn from the cell population.
#' @param sam1 the complete simulated data set.
#' @param xi mean a vector of mean parameters.
#' @param Omega a matrix of covariance paramters.
#' @param alpha a scalor vector of skewness paramter.
#' @param df degrees of freedom.
simSam2 <- function(pop,param,sam1,xi,Omega,alpha,df) {
  
  tempmat <- sam1[sam1$id!=pop,]
  nn <- nrow(sam1)
  
  sam2 <- lapply(1:length(param),function(i) {
    # i=1
    nsam <- round(param[i])	
    if (nsam==0) {
      tempmat
    } else if (nsam > 0) {
      mat <- simOnePop(n=nsam, xi=xi, Omega=Omega, 
                       alpha=alpha, df=df)
      mat$id <- rep(pop, nsam)
      colnames(mat)[1:4] <- colnames(tempmat)[1:4]
      rbind(mat, tempmat)
    }
  })
  return(sam2)
}




############################################################################
## Single marker shifts
## Helper function to evaluate the scenario under which cell populations
## are shifted in a single marker channel

#' Match cell populations under location shift in CD23
#' 
#' @param id identifying number of the cell population that is to be
#'        simulated in varying location shift paramters.
#' @param idToMatch identifying number of the cell population that is 
#'        to be compared with the shifted cell population. 
#' @param simulatedData a simulated data set that contains cell population
#'        mimicing the a real flow cytometry data set.
#' @param paramlist a list object that contains multivariate skew t 
#'        paramters of each cell population (mean, covariance, and skewness).
#' @param shifts amount of shift (in the unit of each marker channel).
#' @param ndraws number of random samples that is to be drawn in 
#'        each cell population comparison in order to compute the FR statistic.
#' @export dshiftcase data.frames of the simulated cases.
#' @export dshiftcase_res data.frames of the results of comparing simulated
#'         cell populations with the target cell population.
make.singleshift.CD23_constant <- function(id, 
  idToMatch, simulatedata, paramlist, shifts, ndraws ) {

  popsToKeep <- simulatedata[simulatedata$id!=id,]
  popsToSim <- simulatedata[simulatedata$id==id,]
  
  CD23params <- (shifts*IQR(simulatedata[simulatedata$id==id,"CD23"]))
  
  sam20 <- lapply(1:length(CD23params),function(j) {
    CD23param <- CD23params[j]
    mat <- sweep(popsToSim,2,STATS=c(0,CD23param,0,0,0),FUN="+")
    return(mat)
  })
  names(sam20) <- paste("CD23:",round(CD23params),sep="")
  
  dshiftcase_res <- lapply( 1:length(shifts),
     function(i) {
     res <- getFRest(XX1=simulatedata[simulatedata$id==idToMatch,],
              XX2=sam20[[i]],sampleMethod="proportional",
              sampleSize = 200, estStat="median",ndraws=ndraws )
               res <- res@ww   })
  names(dshiftcase_res) <- names(sam20)
  
  list(dshiftcase = sam20, dshiftcase_res = dshiftcase_res )
}


#' Match cell populations under location shift in CD3
#' 
#' @param id identifying number of the cell population that is to be
#'        simulated in varying location shift paramters.
#' @param idToMatch identifying number of the cell population that is 
#'        to be compared with the shifted cell population. 
#' @param simulatedData a simulated data set that contains cell population
#'        mimicing the a real flow cytometry data set.
#' @param paramlist a list object that contains multivariate skew t 
#'        paramters of each cell population (mean, covariance, and skewness).
#' @param shifts amount of shift (in the unit of each marker channel).
#' @param ndraws number of random samples that is to be drawn in 
#'        each cell population comparison in order to compute the FR statistic.
#' @export dshiftcase data.frames of the simulated cases.
#' @export dshiftcase_res data.frames of the results of comparing simulated
#'         cell populations with the target cell population.
make.singleshift.CD3_constant <- function(id,
  idToMatch, simulatedata, paramlist, shifts, ndraws) {
  
  popsToKeep <- simulatedata[simulatedata$id!=id,]
  popsToSim <- simulatedata[simulatedata$id==id,]
  
  CD3params <- (shifts*IQR(simulatedata[simulatedata$id==id,"CD3"]))
  
  sam20 <- lapply(1:length(CD3params),function(j) {
    CD3param <- CD3params[j]
    mat <- sweep(popsToSim,2,STATS=c(0,0,CD3param,0,0),FUN="+")
    return(mat)
  })
  names(sam20) <- paste("CD3:",round(CD3params),sep="")
  
  dshiftcase_res <- lapply(1:length(shifts),
                           function(i) {
           res <- getFRest(XX1=simulatedata[simulatedata$id==idToMatch,],                               XX2=sam20[[i]],sampleMethod="proportional",sampleSize=200,                                 estStat="median",ndraws=ndraws)
                             res <- res@ww
                           })
  names(dshiftcase_res) <- names(sam20)
  
  list(dshiftcase=sam20,dshiftcase_res=dshiftcase_res)
}


#' Match cell populations under location shift in CD14
#' 
#' @param id identifying number of the cell population that is to be
#'        simulated in varying location shift paramters.
#' @param idToMatch identifying number of the cell population that is 
#'        to be compared with the shifted cell population. 
#' @param simulatedData a simulated data set that contains cell population
#'        mimicing the a real flow cytometry data set.
#' @param paramlist a list object that contains multivariate skew t 
#'        paramters of each cell population (mean, covariance, and skewness).
#' @param shifts amount of shift (in the unit of each marker channel).
#' @param ndraws number of random samples that is to be drawn in 
#'        each cell population comparison in order to compute the FR statistic.
#' @export dshiftcase data.frames of the simulated cases.
#' @export dshiftcase_res data.frames of the results of comparing simulated
#'         cell populations with the target cell population.
make.singleshift.CD14_constant <- function(id,
  idToMatch, simulatedata, paramlist, shifts, ndraws ) {
  
  popsToKeep <- simulatedata[simulatedata$id!=id,]
  popsToSim <- simulatedata[simulatedata$id==id,]
  
  CD14params <- (shifts*IQR(simulatedata[simulatedata$id==id,"CD14"]))
  
  sam20 <- lapply(1:length(CD14params),function(j) {
    CD14param <- CD14params[j]
    mat <- sweep(popsToSim,2,STATS=c(CD14param,0,0,0,0),FUN="+")
    return(mat)
  })
  names(sam20) <- paste("CD14:",round(CD14params),sep="")
  
  dshiftcase_res <- lapply(1:length(shifts),
                           function(i) {
                   res <- getFRest(XX1=simulatedata[simulatedata$id==idToMatch,],                               XX2=sam20[[i]],sampleMethod="proportional",sampleSize=200,
                      estStat="median",ndraws=ndraws)
                   res <- res@ww
                   })
  names(dshiftcase_res) <- names(sam20)
  
  list(dshiftcase=sam20,dshiftcase_res=dshiftcase_res)
}



#' Match cell populations under location shift in CD19
#' 
#' @param id identifying number of the cell population that is to be
#'        simulated in varying location shift paramters.
#' @param idToMatch identifying number of the cell population that is 
#'        to be compared with the shifted cell population. 
#' @param simulatedData a simulated data set that contains cell population
#'        mimicing the a real flow cytometry data set.
#' @param paramlist a list object that contains multivariate skew t 
#'        paramters of each cell population (mean, covariance, and skewness).
#' @param shifts amount of shift (in the unit of each marker channel).
#' @param ndraws number of random samples that is to be drawn in 
#'        each cell population comparison in order to compute the FR statistic.
#' @export dshiftcase data.frames of the simulated cases.
#' @export dshiftcase_res data.frames of the results of comparing simulated
#'         cell populations with the target cell population.
make.singleshift.CD19_constant <- function(id,
  idToMatch, simulatedata, paramlist, shifts, ndraws) {

  popsToKeep <- simulatedata[simulatedata$id!=id,]
  popsToSim <- simulatedata[simulatedata$id==id,]
  
  CD19params <- (shifts*IQR(simulatedata[simulatedata$id==id,"CD19"]))
  
  sam20 <- lapply(1:length(CD19params),function(j) {
    CD19param <- CD19params[j]
    mat <- sweep(popsToSim,2,STATS=c(0,0,0,CD19param,0),FUN="+")
    return(mat)
  })
  names(sam20) <- paste("CD19:",round(CD19params),sep="")
  
  dshiftcase_res <- lapply(1:length(shifts),
                           function(i) {
          res <- getFRest(XX1=simulatedata[simulatedata$id==idToMatch,],
          XX2=sam20[[i]],sampleMethod="proportional",sampleSize=200,
          estStat="median",ndraws=ndraws)
          res <- res@ww
          })
  names(dshiftcase_res) <- names(sam20)
  
  list(dshiftcase=sam20,dshiftcase_res=dshiftcase_res)
}




############################################################################
## Inappropriate partitioning of the cell population
## Helper function to evaluate the scenario under which a cell population
## is inappropriately partitioned into two cell populations


#' Match and simulate the scenarior under which a cell population
#' is inappropriately divided into two cell populations.
#' 
#' @param id identifying number of the cell population that is to be
#'        simulated in varying location shift paramters.
#' @param idToMatch identifying number of the cell population that is 
#'        to be compared with the shifted cell population. 
#' @param simulatedData a simulated data set that contains cell population
#'        mimicing the a real flow cytometry data set.
#' @param sampleMethod the method used to downsample the number of events 
#'        in each cell population comparison when computing FR statistic.
#' @param ndraws number of random samples that is to be drawn in 
#'        each cell population comparison in order to compute the FR statistic.
#' @param which.park to include the events above or below the cutoffs.
#' 
#' @export partcase data.frames of the simulated cases.
#' @export partcase_res data.frames of the results of comparing simulated
#'         cell populations with the target cell population.
make.partcase <- function(id,
  idToMatch, simulatedata, sampleMethod, ndraws, which.part) {
  popsToKeep <- simulatedata[simulatedata$id!=id,]
  popsToSim <- simulatedata[simulatedata$id==id,]
  breaklist <- quantile(popsToSim$CD23,prob=seq(0.1,0.9,0.1))
  
  partcase <- lapply(1:length(breaklist),function(j) {
    nbreak <- breaklist[j]
    if (which.part=="upper") { mat <- subset(popsToSim,popsToSim$CD23>nbreak) }
    if (which.part=="lower") { mat <- subset(popsToSim,popsToSim$CD23<=nbreak) }
    mat <- rbind(mat,popsToKeep)
    mat <- mat[order(mat$id),]
    return(mat)
  })
  names(partcase) <- breaklist
  
  # compare every test sample (partitioned sample) with the reference sample
  res <- lapply(1:length(breaklist), function(i) 
    makeDistmat(list(simulatedata,partcase[[i]]),sampleMethod=sampleMethod,ndraws=ndraws))
  list(partcase = partcase,partcase_res=res)
}

