# Stat Arb HW3 11.01.17
# Attempt #2

require(data.table)
require(lubridate)
require(caTools)
require(fastAdaboost)
require(quadprog)
# require(xlsx) #I have a problem with xlsx something to do with Java...
rm(list = ls())
gc()

setwd('***')

databaseMat <- R.matlab::readMat('database.mat')
for (i in names(databaseMat)){
  if(!i %in% c('allstocks', 'myday')) {assign(x = i, databaseMat[[i]])}
  else {assign(x = i, databaseMat[i])}
}

allstocks <- allstocks[[1]]
codesDS <- unlist(allstocks[1, , ])
# dimnames(allstocks)[[3]] <- codeDS
# industryDS <- unlist(allstocks[7, , ])
industryDS <- unlist(lapply(allstocks[7, , ], function(x) x[2]))

Date <- as.Date(unlist(myday), format = '%d-%b-%Y')
myDate <- data.table('Date' = lubridate::ymd(Date))
myDate[, YrMo := year(Date) + ((month(Date) - 1) / 12)]
myDate[, I := as.integer(rownames(.SD))]
 

# for (i in c('cap', "isactivenow", "mtbv", "price", "rec", "tcost", "tri", "volume")){
#   assign(x = i, cbind(Date, YrMo, get(i)))
#   setnames(get(i), c('Date', 'YrMo', codesDS))
# }

for (i in c('cap', "isactivenow", "mtbv", "price", "rec", "tcost", "tri", "volume")){
  # assign(x = i, cbind(Date, YrMo, get(i)))
  # setnames(get(i), codesDS)
  # colnames(get(i)) = codesDS
  assign(i, structure(get(i), dimnames = list(NULL, codesDS)))
  # "dimnames<-"(eval(i), list(NULL, codesDS))
  }


# rm(list = ls()[which(!(ls() %in% c('allstocks', 'DB', 'databaseMat', )))]) # Remove unneeded objects
# gc() # clears memory of deleted objects

# ret <- copy(tri)
# ret[, names(ret) := lapply(.SD, function(x) (x / shift(x) - 1))]
ret <- apply(tri, 2, function(x) (x / shift(x) - 1))
# ret <- as.matrix(ret)
# ret[, 3:568 := lapply(.SD, function(x) ifelse(is.na(x), 0, x)), .SDcols = 3:568]

# Industry dummy matrix for alphas:
R <- sapply(unique(industryDS), 
  function(x){
    as.numeric(industryDS == x)
  }
)




firstDays <- unlist(myDate[Date >= min(myDate$Date) + 365,
                    .SD[1, ], 
                    .(YrMo)][, .(I)])

activeUniv <- sapply(setNames(firstDays, firstDays), 
                     function(x) {which(unlist(isactivenow[x, ]) == 1)}, 
                     USE.NAMES = T)

empCovMats <- lapply(firstDays, function(x){
  # This is taking much longer than anticipated to run... fixed
  # x <- firstDays[1]
  activeUniv <- which(unlist(isactivenow[x, ]) == 1)
  
  pastYrUniv <- ret[myDate[I %between% c(x - 1 - 249, x - 1), I], activeUniv]
  pastYrUniv <- pastYrUniv[apply(pastYrUniv, 1, function(x) length(unique(x))>1),]
  
  # Omega is estimation error:
  omegaList <- lapply(as.data.frame(t(pastYrUniv)), #Notice the transpose to use lapply
                      function(x){
                        # browser()
                        x[is.na(x)] <- 0
                        ((x) %*% t(x))
                      })
  
  S_Mat <- Reduce('+', omegaList) / length(omegaList)
  # S_Mat <- apply(simplify2array(omegaList), 1:2, function(x) mean(x, na.rm = T)) #More correct but takes longer to run
  
  sq_omega <- mean(sapply(omegaList, function(x) {sum((x - S_Mat)^2, na.rm = T)})) 
  
  # delta is dispersion
  # delta plus estimation is total var
  # shrinkage target is the scaled identity matrix:
  avg_Var <- mean(apply(pastYrUniv, 2, function(x) var(x, na.rm = T)))
  scaledIden <- avg_Var * diag(length(activeUniv))
  sq_delta <- sum((S_Mat - scaledIden)^2)
  totalVar <- sq_omega + sq_delta
  
  shrinkBeta <- 1 - sq_omega / totalVar
  covMat <- var(pastYrUniv, na.rm = T)
  shrunkCovMat <- (1 - shrinkBeta) * scaledIden + shrinkBeta * covMat
  
  # Alphas:
  # Easier to implement in a for loop!!!
  
  
  cat('\r', as.character(myDate[x, Date]), paste(dim(pastYrUniv)))
  
  return(list('Beta' = shrinkBeta, 
              'covMat' = covMat, 
              'shrCovMat' = shrunkCovMat))
  
})

# Alpha - Momentum:
alphaMom_all <- apply(log(1 + ret), 
                  2, 
                  function(x) {
                    zoo::rollapply(x, 
                                   width = 250, 
                                   function(TS) sum(TS[1:230], na.rm =  T), 
                                   align = 'right', 
                                   fill = NA)
                  })

alphaMom_all <- exp(alphaMom_all) - 1

alphaMom_all <- apply(alphaMom_all, 2, function(x){
  x[x == 0] <- NA
  return(x)
})

# Alpha - Analyst Revision:
alphaRec_all <- apply(ifelse(rec>0, 1, ifelse(rec<0, -1, 0)), 2, function(x){
  zoo::rollapply(x, 
                 width = 45, 
                 function(TS){sum(TS, na.rm = T)}, 
                 align = 'right', 
                 fill = NA)
  
})
# Windsorization:

windsor <- function(vec, unit_var = F){
  # Takes only vectors
  sd_vec <- sd(vec, na.rm = T)
  vec <- vec - median(vec, na.rm = T)
  sd_3 <- 3 * sd_vec + median(vec, na.rm = T)
  sd_less_3 <- -3 * sd_vec + median(vec, na.rm = T)
  vec[vec > sd_3] <- sd_3
  vec[vec < sd_less_3] <- sd_less_3
  # vec <- vec - mean(vec, na.rm = T)
  if (unit_var == T) {vec <- (vec - median(vec)) / sd(vec, na.rm = T)}
  return(vec)
  }


# Alphas:

wghts <- function(len){
  linB <- -1/(sum(0:len)) #linear decay beta
  linA <- len / sum(0:len) #linear devay alpha
  
  wghts <- rev(linA + linB * (0:len))
  return(wghts[-1])
}

TriangularDecayMean <- function(MAT, wid, weights = NA){
  if (!is.na(weights[1])) {meanFN <- function(x){mean(x * weights, na.rm = T)}
  } else {meanFN <- function(x){mean(x, na.rm = T)}}
  
  retVal <- apply(MAT, 
                  MARGIN = 2, 
                  function(x) {
                    zoo::rollapply(x, 
                                   width = wid, 
                                   FUN = function(x) {meanFN(x)}, 
                                   fill = NA, 
                                   align = 'right')
                  })
  
  return(retVal)
}

CleanAndLag <- function(AlphaMat, mean_length = 1, mean_weights = NA, Windsorize = T, windsor_uni_var = F){
  # Need description for the function (Clean and lags the alpha matrix)
  if(is.na(mean_weights)) {
    mean_weights <- wghts(mean_length)
    }
  
  if(mean_length == 1) {
    AlphaMat <- TriangularDecayMean(AlphaMat, mean_length, mean_weights)
  }
 
   AlphaMat <- apply(AlphaMat, 2, shift)[-(1:mean_length), ] # The alpha is now lagged and the excess is dropped off
  
   if(Windsorize == T) {AlphaMat <- t(apply(AlphaMat, 1, windsor, unit_var = windsor_uni_var))}
  
  return(AlphaMat)
}

for (firstDate in firstDays[14]){
  subUniv <- activeUniv[[as.character(firstDate)]]
  # In-Sample:
  BTL <- 250 # Back Test Length
  ST_Len <- 21 # Mean revrersion triangular decay length
  
  meanRev_Univ <- ret[(firstDate - BTL - ST_Len):(firstDate - 1), subUniv]
  
  subR <- R[subUniv, unique(industryDS[subUniv])]
  
  alphaMeanRev <- meanRev_Univ %*% 
    (diag(length(subUniv)) - subR %*% solve(t(subR) %*% subR) %*% t(subR))
  
  alphaMeanRev <- CleanAndLag(alphaMeanRev, ST_Len)
  
  alphaMom <- alphaMom_all[(firstDate - BTL - 1):(firstDate - 1), subUniv]
  alphaMom <- CleanAndLag(alphaMom, 1)
  # sub_alphaMom <- t(apply(sub_alphaMom, 1, windsor))
  # sub_alphaMom <- apply(sub_alphaMom, 2, shift)[-1, ]# The alpha is now lagged and the excess is dropped off
  
  
  Rec_Len <- 45
  alphaRec <- alphaRec_all[(firstDate - BTL - Rec_Len):(firstDate - 1), subUniv] # Do not windsorize
  alphaRec <- CleanAndLag(alphaRec, Rec_Len)
  
  alphaMTB <- mtbv[(firstDate - BTL - 1):(firstDate - 1), subUniv]
  alphaMTB <- CleanAndLag(alphaMTB, 1)
  
  univRet <- ret[(firstDate - BTL):(firstDate - 1), subUniv]
  mktCap <- cap[(firstDate - BTL-1):(firstDate - 1), subUniv]
  laggedMktCap <- apply(mktCap, 2, shift)[-1, ]
  mktUnivRet <- apply(univRet * laggedMktCap, 1, sum, na.rm = T) * 
    1 / apply(laggedMktCap, 1, sum, na.rm = T)
  
  excUnivRet <- univRet - mktUnivRet
  
  # lmOut <- lm(as.vector(excUnivRet) ~ 
  #               as.vector(alphaMeanRev) + 
  #               as.vector(alphaMom) +
  #               as.vector(alphaRec) + 
  #               as.vector(alphaMTB))
  # 
  # lmOut <- lm(as.vector(excUnivRet) ~ 
  #               # as.vector(alphaMeanRev) + 
  #               as.vector(alphaMom) +
  #               # as.vector(alphaRec) + 
  #               as.vector(alphaMTB))
  # 
  # summary(lmOut)
  
  # xTrain <- cbind(as.vector(alphaMeanRev), 
  #                   as.vector(alphaMom), 
  #                   as.vector(alphaRec),  
  #                   as.vector(alphaMTB))
  # 
  # trainData <- data.frame(Y = factor(as.vector(excUnivRet)), 
  #                         A1 = as.vector(alphaMeanRev), 
  #                         A2 = as.vector(alphaMom), 
  #                         A3 = as.vector(alphaRec),  
  #                         A4 = as.vector(alphaMTB))
  # 
  # model <- LogitBoost(xlearn = xTrain, 
  #            ylearn =  as.vector(excUnivRet), nIter = 4)
  
  # LogitBoost can work on 10 levels prediction

  # Out of Sample:
  OOS_EndDate <- myDate[YrMo == myDate[512, YrMo], .SD[.N, I]]
  OOS_alphaMeanRev <- ret[(firstDate - ST_Len):(OOS_EndDate), subUniv]
  
  OOS_alphaMeanRev <- OOS_alphaMeanRev %*% 
    (diag(length(subUniv)) - subR %*% solve(t(subR) %*% subR) %*% t(subR))
  
  OOS_alphaMeanRev <- CleanAndLag(OOS_alphaMeanRev, ST_Len + 1, windsor_uni_var = T)
  
  
  OOS_alphaMom <- alphaMom_all[(firstDate):(OOS_EndDate), subUniv]
  OOS_alphaMom <- CleanAndLag(OOS_alphaMom, windsor_uni_var = T)
  
  # Rec_Len <- 45
  OOS_alphaRec <- alphaRec_all[(firstDate - Rec_Len):OOS_EndDate, subUniv] # Do not windsorize
  OOS_alphaRec <- CleanAndLag(OOS_alphaRec, mean_length = Rec_Len + 1, Windsorize = F)
  
  OOS_alphaMTB <- mtbv[(firstDate):OOS_EndDate, subUniv]
  OOS_alphaMTB <- CleanAndLag(OOS_alphaMTB, windsor_uni_var = T)
  
  # Note: the alphas are not unit variance which makes it incorrect to statically blend them...
  OOS_AggAlpha <- .5 * OOS_alphaMeanRev + .25 * OOS_alphaRec + .15 * OOS_alphaMTB + .1 * OOS_alphaMom
  
  OOS_RetMat <- ret[(firstDate + 1):(OOS_EndDate), subUniv]
  
  w <- rep(0, length(subUniv)) # w is the current portfolio
  
  # summary(lm(as.vector(OOS_RetMat) ~ as.vector(OOS_AggAlpha)))
  # Need to define mu - risk parameter
  mu <- .2
  H <- 2 * mu * rbind(cbind(empCovMats[[14]][['shrCovMat']], -empCovMats[[14]][['shrCovMat']]), 
                  cbind(-empCovMats[[14]][['shrCovMat']], empCovMats[[14]][['shrCovMat']]))
  
  lambda <- .2 # lambda is the trading penalty parameter
  g <- rbind(2 * mu * empCovMats[[14]][['shrCovMat']] %*% w - alpha + lambda * tau, 
             -2 * mu * empCovMats[[14]][['shrCovMat']] %*% w + alpha + lambda * tau)
  
  A <- rbind(cbind(t(subR), -t(subR)), 
             cbind(-t(subR), t(subR)))
  
  indLimit <- 100
  b <- rbind(indLimit * rep(1, ncol(subR)) - t(subR) %*% w, 
             -indLimit * rep(1, ncol(subR)) - t(subR) %*% w)
  # Note, I'm not aware of a function in R that nicely reformats matrices into one matrix.
  #   I passed on writing a function because I didn't see a way to save time by using it.
  
  # Country constraint omitted
  # C <- rbind(t(eqBetas, -t(eqBetas)))
  # 
  # d <- -t(eqBetas) %*% w
  ADV <-  volume[(firstDate - 21):(firstDate - 1), subUniv]# Average daily volume
  ADV <- apply(ADV, 2, mean, na.rm = T)
  
  maxTrade <- .01 * ADV
  maxPosition <- 150
  
  LB <- rep(0, 2 * length(subUniv)) # Lower bound
  UB <- c(pmax(0, pmin(maxTrade, maxPosition - w)), 
              pmax(0, pmin(maxTrade, maxPosition + w)))
  # theta = max trade size
  # pi = max position size
  
  # Quadratic Programming, here we go:
  QP_Amat <- rbind(A, 
                   matrix(1, ncol = 2 * length(w), nrow = 1), 
                   matrix(1, ncol = 2 * length(w), nrow = 1))
  QP_b0 <- rbind(-b, LB, -UB)
  solve.QP(Dmat = H, dvec = -g, Amat = )
  
}
  







# *******************************************************************************
###!!!!!! Everything below this line is essentially pseudocode for reference only
# *******************************************************************************
# # 3
# # Short Term Mean Reversion over 21 days:
# # This should only be evaluated per monthly universe
# ST_Period <- 21
# revWeights <- rev((1/11) - ((1/231) * (0:(21-1))))
# alpharev <- apply(ret[, -(1:2)], 
#                   2, 
#                   function(x) {
#                     zoo::rollapply(x, 
#                       width = 21, 
#                       FUN = function(x) {-mean(x * revWeights, na.rm = T)}, 
#                       fill = NA, 
#                       align = 'right')
#                   })
# 
# 
# 
# 
# 
# 
# 
# 
# # Alpha Rec:
# recWeights <- rev((1/23) - ((1/1035) * (0:(45-1))))
# alpharec <- apply(rec[, -(1:2)], 
#                   2, 
#                   function(x) {
#                     zoo::rollapply(x, 
#                                    width = 45, 
#                                    FUN = function(x) {-mean(x * recWeights, na.rm = T)}, 
#                                    fill = NA, 
#                                    align = 'right')
#                   })
# 
# # Alpha Val - take weighted values by day
# alphaval <- apply(mtbv[, -(1:2)], 
#                        1,
#                        function(x) x/sum(x))
# 
# # Alpha Momentum

# # quodprog
# mktRet <- 
# apply(isactivenow[, -(1:2)], 
#       1, 
#       function(x) {
#         activeUnivNum <- which(x == 1)
#       })
# mktRet <- sapply(isactivenow$Date, 
#                  function(x) {
#                    activeUniv <- unlist(isactivenow[Date == x, -(1:2)])
#                    activeUniv <- names(activeUniv)[which(activeUniv == 1)]
#                    capWeights <- unlist(cap[Date == x, activeUniv, with = F])
#                    capWeights <- capWeights / sum(capWeights)
#                    # print(x)
#                    sum(unlist(ret[, activeUniv, with = F]) * capWeights)
#                  })
# 
# activeUniv <- apply(isactivenow[, -(1:2)],
#       1, 
#       function(x) {
#         which(x == 1) + 2
#       })
# 
# capWghts <- sapply(1:nrow(cap), 
#                    function(x){
#                      activeCap <- unlist(cap[x, activeUniv[[x]], with = F])
#                      activeCap/sum(activeCap, na.rm = T)
#                    })
# 
# mktRet <- sapply(1:nrow(ret), 
#                  function(x){
#                    sum(unlist(ret[x, activeUniv[[x]], with = F]) *
#                          capWghts[[x]])
#                  })
#            
# univBetas <- lapply(ret[, -(1:2)], 
#                     function(returns){
#                       cat('I ')
#                       shift(zoo::rollapply(seq_along(returns),
#                                      width = 245, # Number of trading days for a year
#                                      FUN = function(DTS){
#                                       # lm(mktRet[DTS] ~ returns[DTS])$coef[2]
#                                       X <- cbind(rep(1, 245), mktRet[DTS])
#                                       # View(X)
#                                       # browser()
#                                       if(class(try(solve(m),silent=T))=="matrix"){return(NA)}
#                                       (solve(t(X) %*% X) %*% t(X) %*% returns[DTS])[2]
#                                      }, 
#                                      fill = NA, 
#                                      align = 'right'))
#                     })
# 
# univBetas2 <- (do.call(cbind, univBetas))
# # lm(ret[(1504-245):1504, 3][[1]]~mktRet[(1504-245):1504])$coef
# 
# # 4
# 
# countryDS <- unlist(allstocks[2, , ])
# countryDS <- unlist(lapply(allstocks[2, , ], function(x) x[2]))
#  
# matF <- matrix(0, length(codesDS), length(unique(countryDS)))
# 
# countryDummy <- sapply(countryDS,
#                    function(x) which(x == unique(countryDS))[1])
# 
# for (i in 1:nrow(matF)) {matF[i, countryDummy[i]] <- 1}
# 
# alphaRev2 <- as.matrix(apply(alpharev, 
#                    2, 
#                    function(x){
#                      x[which(is.na(x))] <- 0
#                      x
#                    }))
# 
# tcost2 <- as.matrix(apply(tcost, 
#                    2, 
#                    function(x){
#                      x[which(is.na(x))] <- .1
#                      x
#                    }))
# 
# ADV <- as.matrix(apply(copy(volume[, -(1:2)]),
#              2,
#               function(x){
#                 newX <- zoo::rollapply(x, 
#                                width = 245, 
#                                FUN = function(subX){
#                                  mean(subX, na.rm = T)
#                                }, 
#                                align = 'right', 
#                                fill = NA)
#                 newX[which(is.na(newX))] <- 0
#                 newX
#               }))
# 
# w <- matrix(rep(0, (ncol(ret)-2) * 246), ncol = length(codesDS), nrow = 246)
# 
# 
# for (i in 246){
#   fStar <- .1 # $100,000
#   rStar <- .3 # $300,000
#   theta <- pmax(.01 * ADV[i, ], .15) # S/B capped at 150000 
#   mu <- .4
#   lambda <- .4
#   tau <- as.numeric(tcost2[i, (activeUniv[[i]]-2)])
#   a <- alphaRev2[i, (activeUniv[[i]]-2)]
#   
#   Beta <- univBetas2[i, (activeUniv[[i]]-2)]
#   
#   sigMat <- cov(ret[(i-245):i, activeUniv[[i]], with = F])
#   
#   wghts <- w[i, (activeUniv[[i]]-2)]
#   
#   subR <- R[(activeUniv[[i]]-2), ]
#   
#   subF <- matF[(activeUniv[[i]]-2), ]
#   
#   # mat.u <- matrix(c(y, z), ncol = 1)
#   # Need to define sigMat
#   # Need to define w
#   PI <- 1 # Max position size
#   mat.H <- 2 * mu * cbind(rbind(sigMat, -sigMat), rbind(-sigMat, sigMat))
#   mat.g <- rbind(2 * mu * sigMat %*% wghts - a + lambda * tau, 
#                  -2 * mu * sigMat %*% wghts + a + lambda * tau)
#   
#   mat.A <- cbind(rbind(t(subR), -t(subR), t(subF), -t(subF)), 
#                  rbind(-t(subR), t(subR), -t(subF), t(subF)))
#   
#   
#   One <- matrix(rep(1, 40), ncol = 1)
#   mat.b <- rbind(rStar * One - t(subR) %*% wghts, 
#                  rStar * One + t(subR) %*% wghts, 
#                  fStar * One - t(subF) %*% wghts, 
#                  fStar * One + t(subF) %*% wghts)
#   
#   mat.C <- cbind(t(Beta), -t(Beta))
#   
#   mat.d <- -t(Beta) * w
#   
#   n <- ncol(ret)-2
#   LB <- matrix(rep(0, 2*n), ncol = 1)
#   
#   UB <- rbind(pmax(0, pmin(theta, PI - w)), 
#               pmax(0, pmin(theta, PI + w))) 
# }
