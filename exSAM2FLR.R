##~--------------------------------------------------------------------------
# Code to take the SAM assessment results from stockassessment.org (new TMB fits), 
# and get an FLstock object
# D.C.M.Miller
# (slightly modified by Anders Nielsen)
##~--------------------------------------------------------------------------


rm(list=ls())


library(stockassessment)
require(FLCore)

exportFLstock <- function(name){

  getCatch<-function(fit){
    aux <- fit$data$aux
    logobs <- fit$data$logobs
    .goget <- function(y, a) {
      ret <- exp(logobs[aux[, "fleet"] == 1 & aux[, "year"] == y & aux[, "age"] == a])
      ifelse(length(ret) == 0, 0, ret)
    }
    CW<-fit$data$catchMeanWeight
    tmp <- outer(rownames(CW), colnames(CW), Vectorize(.goget))
    dimnames(tmp)[[1]] <- rownames(CW)
    dimnames(tmp)[[2]] <- colnames(CW)
    return(tmp)
  }


  fit<-fitfromweb(name, character.only=TRUE)

  ages <- fit$conf$minAge:fit$conf$maxAge
  years <- fit$data$years
  meanFages <- fit$conf$fbarRange[1]:fit$conf$fbarRange[2]


  flq <- FLQuant(NA, dimnames = list(age = ages, year = years), quant='age')
  stk <- FLStock(stock.n = flq, name = name, desc = "FLStock_from_SAM")

  units(stk)[1:17]    <- as.list(c(rep(c("tonnes","thousands","kg"),4), rep("NA",2),"f",rep("NA",2)))
  # Mean F range
  range(stk)[c("minfbar","maxfbar")]    <- c(min(meanFages), max(meanFages))
  # Last age a plusgroup
  if(fit$conf$maxAgePlusGroup==1){
    stk  <- setPlusGroup(stk,stk@range["max"])
  }
    
  # add catches
  #catch.n(stk)[,ac(years[1]:(max(years)-1))] <- landings.n(stk)[,ac(years[1]:(max(years)-1))] <- tmpCat; rm(tmpCat)
  tmpCat <- t(getCatch(fit))
  tmpLF <- t(fit$data$landFrac)
  dms <- list(intersect(ac(ages),dimnames(tmpCat)[[1]]),intersect(years,dimnames(tmpCat)[[2]]))
  catch.n(stk)[dms[[1]],dms[[2]]] <- tmpCat[dms[[1]],dms[[2]]]
  catch.n(stk)[is.na(catch.n(stk))] <- 0
  landings.n(stk)[dms[[1]],dms[[2]]] <- tmpCat[dms[[1]],dms[[2]]] * tmpLF[dms[[1]],dms[[2]]]
  landings.n(stk)[is.na(landings.n(stk))] <- 0
  discards.n(stk)[] <- catch.n(stk) - landings.n(stk)
  rm(tmpCat, dms)

   #catch.wt(stk)[,ac(years[1]:(max(years)-1))] <- t(read.ices("cw.dat"))[-1,]
  tmpCwt <- t(fit$data$catchMeanWeight)
  dms <- list(intersect(ac(ages),dimnames(tmpCwt)[[1]]),intersect(years,dimnames(tmpCwt)[[2]]))
  catch.wt(stk)[dms[[1]],dms[[2]]] <- tmpCwt[dms[[1]],dms[[2]]]; rm(tmpCwt,dms)

  tmpLwt <- t(fit$data$landMeanWeight)
  dms <- list(intersect(ac(ages),dimnames(tmpLwt)[[1]]),intersect(years,dimnames(tmpLwt)[[2]]))
  landings.wt(stk)[dms[[1]],dms[[2]]] <- tmpLwt[dms[[1]],dms[[2]]]; rm(tmpLwt,dms)

  tmpDwt <- t(fit$data$catchMeanWeight)
  dms <- list(intersect(ac(ages),dimnames(tmpDwt)[[1]]),intersect(years,dimnames(tmpDwt)[[2]]))
  discards.wt(stk)[dms[[1]],dms[[2]]] <- tmpDwt[dms[[1]],dms[[2]]]; rm(tmpDwt,dms)

  discards(stk) <- computeDiscards(stk)
  landings(stk) <- computeLandings(stk)
  catch(stk) <- computeCatch(stk)

  # add bio
  tmp <- t(fit$data$propMat); dms <- list(intersect(ac(ages),dimnames(tmp)[[1]]),intersect(years,dimnames(tmp)[[2]]))
  mat(stk)[dms[[1]],dms[[2]]] <- tmp[dms[[1]],dms[[2]]]; rm(tmp,dms)
  tmp <- t(fit$data$natMor); dms <- list(intersect(ac(ages),dimnames(tmp)[[1]]),intersect(years,dimnames(tmp)[[2]]))
  m(stk)[dms[[1]],dms[[2]]] <- tmp[dms[[1]],dms[[2]]]; rm(tmp,dms)
  tmp <- t(fit$data$propF); dms <- list(intersect(ac(ages),dimnames(tmp)[[1]]),intersect(years,dimnames(tmp)[[2]]))
  harvest.spwn(stk)[dms[[1]],dms[[2]]] <- tmp[dms[[1]],dms[[2]]]; rm(tmp,dms)
  tmp <- t(fit$data$propM); dms <- list(intersect(ac(ages),dimnames(tmp)[[1]]),intersect(years,dimnames(tmp)[[2]]))
  m.spwn(stk)[dms[[1]],dms[[2]]] <- tmp[dms[[1]],dms[[2]]]; rm(tmp,dms)

  # Update stock and fisheries from SAM fit
  stock.n(stk)[] <-  t(ntable(fit))
  tmp <- t(fit$data$stockMeanWeight); dms <- list(intersect(ac(ages),dimnames(tmp)[[1]]),intersect(years,dimnames(tmp)[[2]]))
  stock.wt(stk)[dms[[1]],dms[[2]]] <- tmp[dms[[1]],dms[[2]]]; rm(tmp,dms)
  stock(stk)[] <- computeStock(stk)
  harvest(stk)[]<-t(faytable(fit))

  return(stk)
}

stk<-exportFLstock("gobherring_2018")

