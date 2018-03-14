# Stat Arb - European Indices Constituents
#   Written with scalability and universe changes in mind
# author: Danny Ferdman
# begin date: 02/15/2018
# completion date: 3/12/2018

# Spent time making the backtest as realistic as possible.
#   ie: when a stock goes ex univ need to figure out how to start getting rid  of it 
#   as soon as possible. One stock had a good alpha while the ADV kept going down
#   and the optimizer didn't cut down the position. This created a large illiquid position.
# 10 day ADV constraint is not binding - it is working correctly I got stuck in a large position when
#   vol completely collapsed.

# Assume all trading occurs at the last second of every trading day but the book is updated a second later
#   and is applied over the following date for reporting and analysis.

# Some aspects of the code use slightly different indexing. Need to fix it eventually.
# As stock universe expands quad prog takes much longer to solve - main bottle neck.

# Code that was deemed too valuable to be shared publicly online is replaced with bracketed text.


require(data.table)
require(lubridate)
require(quadprog)
rm(list = ls())
gc()

setwd('***')

# Unpacking all the data modules from an array object
databaseMat <- R.matlab::readMat('database.mat')
for (i in names(databaseMat)){
  if(!i %in% c('allstocks', 'myday')) {assign(x = i, databaseMat[[i]])}
  else {assign(x = i, databaseMat[i])}
}

allstocks <- allstocks[[1]]
codesDS <- unlist(allstocks[1, , ])
industryDS <- unlist(lapply(allstocks[7, , ], function(x) x[2]))

Date <- as.Date(unlist(myday), format = '%d-%b-%Y')
myDate <- data.table('Date' = lubridate::ymd(Date))
myDate[, YrMo := year(Date) + ((month(Date) - 1) / 12)]
myDate[, I := as.integer(rownames(.SD))]
 

for (i in c('cap', "isactivenow", "mtbv", "price", "rec", "tcost", "tri", "volume")){
  assign(i, structure(get(i), dimnames = list(NULL, codesDS)))
  }

ret <- apply(tri, 2, function(x) (x / shift(x) - 1))

cleanRet <- ret # creating a returns matrix that doesn't have NaNs for attribution purposes
cleanRet[which(is.nan(cleanRet))] <- 0

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

empCovMats <- lapply(firstDays, function(x){ # Empirical convariance matrices
  
  activeUniv <- which(unlist(isactivenow[x, ]) == 1)
  
  pastYrUniv <- ret[myDate[I %between% c(x - 1 - 249, x - 1), I], activeUniv]
  pastYrUniv <- pastYrUniv[apply(pastYrUniv, 1, function(x) length(unique(x))>1),]
  
  # [Stable covariance estimation code omitted]
  covMat <- var(pastYrUniv, na.rm = T)
  
  cat('\r', as.character(myDate[x, Date]), paste(dim(pastYrUniv)))
  
  return(list('covMat' = covMat))
  
})

names(empCovMats) <- as.character(firstDays)

# Alpha Signals:
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

# Average Daily Volume (used for max trade and size constraints in the optimizer)
ADV_all <- apply(volume, 2, function(V){
  retVal <- zoo::rollapply(V, 
                 width = 21, 
                 mean, 
                 na.rm = T, 
                 fill = NA, 
                 align = 'right')
  retVal <- zoo::na.locf(retVal, na.rm = F)
  if(length(retVal) != length(V)){browser()}
  return(retVal)
})


# Creating market index returns based on monthly sub universes:====
first_ofMonth <- unlist(myDate[, .SD[1, ], 
                              .(YrMo)][, .(I)])

names(first_ofMonth) <- myDate[first_ofMonth, as.character(Date)]

mktIndexRet <- lapply(first_ofMonth, function(dayInd){
  # browser()
  activeNow <- as.logical(isactivenow[dayInd, ])
  
  dayInd_end <- first_ofMonth[which(first_ofMonth == dayInd) + 1] - 1
  if(which(first_ofMonth == dayInd) == length(first_ofMonth)) {
    dayInd_end <- nrow(myDate)
  }
  retMat <- ret[dayInd:dayInd_end, activeNow]
  
  capMat <- cap[(dayInd-1):(dayInd_end-1), activeNow]
  if(dayInd == 1){ capMat <- rbind(cap[1, activeNow], capMat)}
  
  retVal <- rowSums(retMat * capMat, na.rm = T) / rowSums(capMat, na.rm = T)
  cat(dayInd, '\n')
  return(retVal)
})

mktIndexRet <- unlist(mktIndexRet)

# Functions:====
# Windsorization:
windsor <- function(vec, unit_var = F){
  # Takes only vectors
  sd_vec <- sd(vec, na.rm = T)
  vec <- vec - median(vec, na.rm = T)
  sd_3 <- 3 * sd_vec + median(vec, na.rm = T)
  sd_less_3 <- -3 * sd_vec + median(vec, na.rm = T)
  vec[vec > sd_3] <- sd_3
  vec[vec < sd_less_3] <- sd_less_3
  if (unit_var == T) {vec <- (vec - median(vec, na.rm = T)) / sd(vec, na.rm = T)}
  return(vec)
  }

wghts <- function(len){
  linB <- -1/(sum(0:len)) #linear decay beta
  linA <- len / sum(0:len) #linear decay alpha
  
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

CleanAndDecay <- function(AlphaMat, mean_length = 1, mean_weights = NA, Windsorize = T, windsor_uni_var = F){
  # """
  # Clean the data and calculates a weighted mean (decay).
  # Mainly used for the alpha signals to make the signal usable.
  # """
  # Removed the lagging feature
  if(is.na(mean_weights)) {
    mean_weights <- wghts(mean_length)
    }
  
  if(mean_length != 1) {
    AlphaMat <- TriangularDecayMean(AlphaMat, mean_length, mean_weights)
    AlphaMat <- AlphaMat[-(1:(mean_length - 1)), ] # The excess from taking the mean is dropped off
  }
  
   if(Windsorize == T) {AlphaMat <- t(apply(AlphaMat, 1, windsor, unit_var = windsor_uni_var))}
   
   AlphaMat <- apply(AlphaMat, 2, function(x){
     x[is.na(x)] <- 0
     return(x)
   })
  
  return(AlphaMat)
}

subjectTo_Matrices <- function(argList) {
  # """
  # Creates a subject to condition matrix from a list to be
  # used in qudaratic programming.
  # Must supply a list of constraints for argList.
  # Note that incorporating regex is the way to go but I cannot spend much more time on it.
  # """
  
  # [Code omitted]
  
  QP_argVar <- do.call('rbind', lapply(matList, function(Item) Item[[1]]))
  QP_Bounds <- do.call('rbind', lapply(matList, function(Item) Item[[2]]))
  
  return(list(QP_argVar, QP_Bounds))
}

# Backtesting:====

backtest_Port <- cleanRet * NA
rownames(backtest_Port) <- as.character(myDate$Date)
for (firstDate in firstDays[-(1:3)]){
  Timer <- Sys.time()
  isFirst <- ifelse(any(is.na(backtest_Port[firstDate - 1, ])), # checks if this is the first iteration
                    T, 
                    F)
  if(isFirst == T){
    printCounter <- 0
    backtest_Port[firstDate - 1, ] <- 0
  }
  
  index <- as.character(firstDate)
  subUniv <- activeUniv[[index]]

  # [Factor blend code omitted]

  # Out of Sample:
  subR <- R[subUniv, unique(industryDS[subUniv])]
  OOS_EndDate <- myDate[YrMo == myDate[firstDate, YrMo], .SD[.N, I]]
  
  backtest_Port[firstDate:OOS_EndDate, ] <- 0 # Initiating the values from NA to 0
  
  ST_Len <- 21
  OOS_alphaMeanRev <- ret[(firstDate - ST_Len):(OOS_EndDate - 1), subUniv] # NB: not subtracting 1 from firstDate
  #   because the ST_Len takes care of it
  
  OOS_alphaMeanRev <- OOS_alphaMeanRev %*% 
    (diag(length(subUniv)) - subR %*% solve(t(subR) %*% subR) %*% t(subR)) # demeans by industry
  
  OOS_alphaMeanRev <- CleanAndDecay(OOS_alphaMeanRev, ST_Len, windsor_uni_var = T)
  
  
  OOS_alphaMom <- alphaMom_all[(firstDate - 1):(OOS_EndDate - 1), subUniv]
  OOS_alphaMom <- CleanAndDecay(OOS_alphaMom, windsor_uni_var = T)
  
  Rec_Len <- 45
  OOS_alphaRec <- alphaRec_all[(firstDate - Rec_Len):(OOS_EndDate - 1), subUniv] # Do not windsorize
  OOS_alphaRec <- CleanAndDecay(OOS_alphaRec, mean_length = Rec_Len, Windsorize = F)
  
  OOS_alphaMTB <- mtbv[(firstDate - 1):(OOS_EndDate - 1), subUniv]
  OOS_alphaMTB <- CleanAndDecay(OOS_alphaMTB, windsor_uni_var = T)
  
  # Note: the alphas are not unit variance which makes it incorrect to statically blend them...
  # This area is stated for future improvement.
  OOS_AggAlpha <- -.5 * OOS_alphaMeanRev + .25 * OOS_alphaRec - .15 * OOS_alphaMTB + .1 * OOS_alphaMom
  
  OOS_RetMat <- ret[(firstDate):(OOS_EndDate), subUniv]
  
  # Calculate equity betas over 250 days on a rolling basis:
  # NB: the calculation is also used on an out of sample basis. I'm not
  #   using the returns for the day on which the calculated beta is used.
  rowsSelect <- (firstDate - 250):(OOS_EndDate - 1)
  OOS_Betas <- zoo::rollapply(rowsSelect, 
                              width = 250, 
                              fill = NULL,
                              align = 'right', 
                              FUN = function(rows){
                                apply(ret[rows, subUniv], 2, 
                                      function(Y){
                                        notNA <- which(!is.na(Y))
                                        X <- cbind(1, mktIndexRet[rows][notNA])
                                        betas <- solve(t(X) %*% X) %*% t(X) %*% Y[notNA]
                                        return(betas[2])
                                      }
                                )
                              })
  
  # mu and lambda parameters should constantly be adjusted. Future improvement.
  mu <- 2
  
  H <- 2 * mu * rbind(cbind(empCovMats[[index]][['covMat']], -empCovMats[[index]][['covMat']]),
                  cbind(-empCovMats[[index]][['covMat']], empCovMats[[index]][['covMat']]))
  # [Stable H omitted]
  # T-Cost data is monthly averages
  tau <- tcost[firstDate - 1, subUniv] / 2 
  lambda <- 5 # lambda is the trading penalty parameter
  
  for (day in (1:nrow(OOS_AggAlpha) - 1)){
    # NB: day object starts from 0
    alpha <- OOS_AggAlpha[day + 1, ]

    # Check that ex universe holdings are liquidated:
    if(isFirst == F & any(backtest_Port[(firstDate + day - 1), -subUniv] != 0)){ 
      
      # NB: I decided to ID and work with the explicit non-zero ex-univ only instead of the entire ex Universe for 
      #   easier debugging
      exUnivNames <- names(which(backtest_Port[(firstDate + day - 1), -subUniv] != 0))
      exUnivMat <- backtest_Port[(firstDate + day - 1), exUnivNames]
      
      # Max trade and size, s/b the aggressive version of the one a few lines down
      maxTrade <- .01 * ADV_all[(firstDate + day - 1), exUnivNames]
      exUnivTrade <- -sign(exUnivMat) * pmin(maxTrade, abs(exUnivMat - 0))
      backtest_Port[(firstDate + day), exUnivNames] <- 
        backtest_Port[(firstDate + day - 1), exUnivNames] + exUnivTrade
    } else {exUnivNames <- 0}
    
    w <- backtest_Port[(firstDate + day - 1), subUniv] # Begining of day portfolio
    
    # quadprog takes the form: min(.5*b'*D*b -d'b) s.t 
    #                          A'b >= b0
    
    # Minimization setup with H/D defined outside of the loop since it doesn't change
    g <- cbind(2 * mu * empCovMats[[index]][['covMat']] %*% w - alpha + lambda * tau, 
               -2 * mu * empCovMats[[index]][['covMat']] %*% w + alpha + lambda * tau)
    
    # Constraints:
    # Industry constraint
    indLimit <- 100
    indCons <- list(pmin(-indLimit * rep(1, ncol(subR)) - t(subR) %*% w, -1e-04),
         '<', 
         cbind(t(subR), -t(subR)), 
         '<', 
         pmax(1e-04, indLimit * rep(1, ncol(subR)) - t(subR) %*% w))
    
    # Max trade and position size constraint
    maxTrade <- pmin(50, .01 * ADV_all[(firstDate + day - 1), subUniv])
    # maxPosition <- 150
    maxPosition <- pmin(150, 10 * .01 * ADV_all[(firstDate + day - 1), subUniv])
    # NB: the lower bound is a little less than zero because a zero bound throws an error
    #   in the solver when both bounds are the same, I think it is looking for meq option
    #   to be specified for the bound when both bounds are the same.
    tradeAndSize_Constraint <- list(c(pmin(-pmin(1e-04, maxPosition + w), maxTrade - 1e-04), 
                                      pmin(-pmin(1e-04, maxPosition - w), maxTrade - 1e-04)), 
                                    '<', 
                                    diag(2 * length(w)), 
                                    '<', 
                                    c(pmax(0, pmin(maxTrade, maxPosition - w)),
                                          pmax(0, pmin(maxTrade, maxPosition + w)))
                                    )
    
    # theta = max trade size
    # pi = max position size
    
    # Long-short portfolio constraint
    LS_Constraint <- list(-(sum(w) + sum(backtest_Port[(firstDate + day), exUnivNames])), 
                          '=', 
                          matrix(c(rep(1, length(w)), 
                                   rep(-1, length(w))), 
                                 nrow = 1))
    
    # Beta constraint
    mktNeutralCons <- list(-.05 * sum(abs(w)) - (w %*% OOS_Betas[day + 1, ]), 
                           '=', 
                           c(OOS_Betas[day + 1, ], 
                             -OOS_Betas[day + 1, ]))
    
    QP_List <- subjectTo_Matrices(list(LS_Constraint
                                       ,mktNeutralCons
                                        ,indCons
                                       , tradeAndSize_Constraint
                                       ))

    QP_Amat <- QP_List[[1]]
    QP_b0 <- QP_List[[2]] 

    
    # **************************************
    # Quadratic Programming, here we go:
    # quadprog takes the form: min(.5*b'*D*b -d'b) s.t 
    #                          A'b >= b0
    # **************************************
    
    sol1 <- solve.QP(Dmat = H + diag(ncol(H)) * 5e-4, dvec = -g, Amat = t(QP_Amat), bvec = QP_b0, meq = 2)
    
    optimW <- sol1$solution
    y <- optimW[1:(length(optimW)/2)]
    z <- optimW[-(1:(length(optimW)/2))]
    x <- w + y - z
    
    backtest_Port[firstDate + day, subUniv] <- x # Close of day portfolio/ begining of the subsequent day

    if(day != 0 & abs(sum(backtest_Port[firstDate + day,])) > 1) {
      message('error: sum of weights is not 0')
              browser()}
  }

  cat('\r', 
      paste(myDate[as.numeric(index), Date]), 
      ' ', 
      paste(round(Sys.time() - Timer, 2), attributes(Sys.time()-Timer)[1]))
  if(printCounter < 10) {cat('\n')
    } else if(printCounter == 10){cat('\n...\n')}
  
  printCounter <- printCounter + 1

}

# Results analysis and visualization:====

createAxis <- function(){
  # Labels the x-axis on the plot
  axis(side = 1, 
       at = which(names(OOS_cumNetRet) %in% as.character(myDate$Date[firstDays])), 
       labels = names(OOS_cumNetRet)[names(OOS_cumNetRet) %in% as.character(myDate$Date[firstDays])], 
       col.axis = 'blue', 
       las = 2, 
       hadj = .92)
}
# Extract the backtested days index, drop the first day which just consists of transaction costs.
OOS_dayIndex <- which(apply(backtest_Port, 1, function(Row){any(!is.na(Row))}) == T)[-1]

OOS_PortRet <- ((backtest_Port * cleanRet) %*% matrix(1, ncol = 1, nrow = ncol(cleanRet)))[OOS_dayIndex]
names(OOS_PortRet) <- rownames(backtest_Port)[OOS_dayIndex]


OOS_Trades <- rbind(backtest_Port[-1, ], 
                    matrix(NA, nrow = 1, ncol = ncol(backtest_Port))) - 
              backtest_Port

OOS_Trades <- OOS_Trades[OOS_dayIndex, ]
rownames(OOS_Trades) <- rownames(backtest_Port)[OOS_dayIndex]

OOS_TCost <- rowSums(abs(OOS_Trades) * tcost[OOS_dayIndex, ] / 2, na.rm = T) 
# mean(OOS_PortRet, na.rm = T); sd(OOS_PortRet, na.rm = T)
# mean(OOS_PortRet, na.rm = T)/ sd(OOS_PortRet, na.rm = T) *sqrt(252)
# plot(cumsum(na.omit(OOS_PortRet)), type = 'l') 

OOS_netPort <- OOS_PortRet - OOS_TCost # Out of sample net portfolio returns

OOS_cumNetRet <- cumsum(OOS_netPort) # Out of sample cumulative net returns

OOS_Size <- rowSums(abs(backtest_Port[OOS_dayIndex, ]))

mean(OOS_netPort, na.rm = T); sd(OOS_netPort, na.rm = T)
mean(OOS_netPort, na.rm = T)/ sd(OOS_netPort, na.rm = T) * sqrt(252)
# Sharpe ratio:
mean(OOS_netPort / (.2 * OOS_Size) - .05/252) / 
  sd(OOS_netPort / (.2 * OOS_Size)) * sqrt(252)
paste0('T-Stat: ', mean(OOS_netPort, na.rm = T) / 
         (sd(OOS_netPort, na.rm = T)/sqrt(length(na.omit(OOS_netPort)))))

cor(cbind(OOS_netPort/(.2*OOS_Size), mktIndexRet[OOS_dayIndex]))
  
# Drawdown
drawDownMat <- sapply(seq_along(OOS_cumNetRet), function(dayInd){
  rolling_draw <- OOS_cumNetRet[dayInd:length(OOS_cumNetRet)] - OOS_cumNetRet[dayInd]
  return(c(rep(NA, dayInd - 1), 
           rolling_draw))
  }
)

drawDowns <- list()

while(any(!is.na(drawDownMat))){
  # Identify the top 3 draw downs and record them in a list:
  currentDraw_Index <- which.min(apply(drawDownMat, 2, function(vec){
    if(any(!is.na(vec))){return(min(vec, na.rm = T))
    } else{return(10^6)}
  }))
  currentDraw_End <- which.min(drawDownMat[, currentDraw_Index])
  drawDowns[[i]] <- drawDownMat[currentDraw_Index:currentDraw_End, currentDraw_Index]
  if(length(drawDowns[[i]]) <= 2) break
  
  drawDownMat[, seq(currentDraw_Index,currentDraw_End)] <- NA
  drawDownMat[-(1:currentDraw_Index), 1:currentDraw_Index] <- NA
}

# Market draw down: #NB: copying this over is not efficient but I only anticipate doing this twice.
pltMktIndexRet <- c(0, na.omit(.2 * OOS_Size * mktIndexRet))
drawDownMat <- sapply(seq_along(pltMktIndexRet), function(dayInd){
  rolling_draw <- pltMktIndexRet[dayInd:length(pltMktIndexRet)] - pltMktIndexRet[dayInd]
  return(c(rep(NA, dayInd - 1), 
           rolling_draw))
})

mkt_drawDowns <- list()

while(any(!is.na(drawDownMat))){
  # Identify the top 3 draw downs and record them in a list:
  currentDraw_Index <- which.min(apply(drawDownMat, 2, function(vec){
    if(any(!is.na(vec))){return(min(vec, na.rm = T))
    } else{return(10^6)}
  }))
  currentDraw_End <- which.min(drawDownMat[, currentDraw_Index])
  mkt_drawDowns[[i]] <- drawDownMat[currentDraw_Index:currentDraw_End, currentDraw_Index]
  if(length(mkt_drawDowns[[i]]) <= 2) break
  
  drawDownMat[, seq(currentDraw_Index,currentDraw_End)] <- NA
  drawDownMat[-(1:currentDraw_Index), 1:currentDraw_Index] <- NA
}

# Create cumulative returns chart and comparison
plot(OOS_cumNetRet, 
     type = 'l', 
     xlab = "", 
     ylab = 'Net Cumulative Return (in \u20AC1,000)', 
     xaxt = 'n', 
     main = 'Strategy Cumulative Returns', 
     ylim = c(-750, max(OOS_cumNetRet) + 750))
polygon(x = c(seq_along(OOS_cumNetRet), rev(seq_along(OOS_cumNetRet))), 
        y = c(OOS_cumNetRet, rep(0, length(OOS_cumNetRet))), 
        col = 'darkolivegreen1', 
        border = NA)
lines(OOS_cumNetRet, 
      col = 'green4', 
      lwd = 1.25)
createAxis()

for(i in 1:3){
  # Shade the top 3 draw downs on the plot:
  dayIndex <- which(names(OOS_cumNetRet) %in% names(drawDowns[[i]]))
  polygon(x = unlist(mapply(list, range(dayIndex), range(dayIndex), SIMPLIFY = F)), # draws out the date arange
          y = c(0, min(OOS_cumNetRet[dayIndex] - 400), rev(c(0, min(OOS_cumNetRet[dayIndex] - 400)))),
          border = 'purple', 
          lty = 3, 
          lwd = .25)
  polygon(x = c(dayIndex, rev(dayIndex)), # fills in the draw down area
          y = c(OOS_cumNetRet[dayIndex] - 400, 
                rep(min(OOS_cumNetRet[dayIndex] - 400), length(dayIndex))),
          col = 'orange', border = 'purple')

}

mktPerf <- cumsum(c(0, na.omit(.2 * OOS_Size * mktIndexRet))) 
polygon(x = c(seq_along(OOS_cumNetRet), rev(seq_along(OOS_cumNetRet))), 
        y = c(pmin(0, mktPerf), rep(0, length(mktPerf))), 
        col = 'red', 
        border = NA)
polygon(x = c(seq_along(OOS_cumNetRet), rev(seq_along(OOS_cumNetRet))), 
        y = c(pmax(0, mktPerf), rep(0, length(mktPerf))), 
        col = 'gray83', 
        border = NA)
lines(mktPerf, 
      col = 'gray20', 
      lwd = 1.25)

legend('topleft',
       legend = c('Start Cumulative Net Return', 
                  'Market Index Cumulative Return', 
                  'Draw Down (top 3, offset on graph)'), 
       lty = c(1, 1, NA),
       lwd = c(2.75, 2.75, NA), 
       col = c('green4', 'gray20', 'purple'), 
       pt.cex = c(1, 1, 2), 
       pch = c(NA, NA, 22), 
       pt.bg = c(NA, NA, 'orange'), 
       cex = .75,
       text.width = 3.5,  
       bty = 'n')

pltOOS_Size <- na.omit(OOS_Size)
plot(na.omit(pltOOS_Size), 
     type = 'l', 
     xlab = "", 
     ylab = 'Total Book Size (in \u20AC1,000)', 
     xaxt = 'n', 
     main = 'Strategy Total Book Size')
createAxis()

polygon(x = c(seq_along(pltOOS_Size), rev(seq_along(pltOOS_Size))), 
        y = c(pltOOS_Size, rep(min(pltOOS_Size), length(pltOOS_Size))), 
        col = 'darkolivegreen1', 
        border = NA)

lines(pltOOS_Size)

plot(na.omit(apply(OOS_Trades, 1, function(S){
  if(any(!is.na(S))) {return(sum(abs(S), na.rm = T))
  } else {return(NA)}
})), type = 'p', ylab = 'Trade Size (\u20AC1000)')

# Comparative stats:

statsMat <- matrix(NA, nrow = 7, ncol = 2)
rownames(statsMat) <- c('Annual Sharpe Ratio', 
                        'Market correlation', 
                        'Max drawdown', 
                        'Longest drawdown', 
                        'Average trade size', 
                        'Skewness', 
                        'Kurtosis')
colnames(statsMat) <- c('Strategy', 'Market')

statsMat[1, 1] <- mean((OOS_netPort / (.2 * OOS_Size))[testedIndex][-1] - .05 / 252) / 
  sd((OOS_netPort / (.2 * OOS_Size))[testedIndex][-1]) * sqrt(252) # Sharpe ratio

statsMat[1, 2] <- mean(mktIndexRet[testedIndex]) / sd(mktIndexRet[testedIndex]) *sqrt(252) # Mkt SR

statsMat[2, ] <- cor(na.omit(cbind(OOS_netPort/(.2*OOS_Size), mktIndexRet)))[2, ]

statsMat[3, ] <- c(min(unlist(lapply(drawDowns, min))), 
                   min(unlist(lapply(mkt_drawDowns, min))))
                   
statsMat[4, ] <- c(max(unlist(lapply(drawDowns, length))), 
                   max(unlist(lapply(mkt_drawDowns, length))))

statsMat[5, 1] <- mean(rowSums(abs(OOS_Trades)), na.rm = T)

statsMat[6, ] <- c(moments::skewness((OOS_netPort / (.2 * OOS_Size))[testedIndex][-1]), 
                   moments::skewness(mktIndexRet[testedIndex]))

statsMat[7, ] <- c(moments::kurtosis((OOS_netPort / (.2 * OOS_Size))[testedIndex][-1]), 
                   moments::kurtosis(mktIndexRet[testedIndex]))

statsMat <- round(statsMat, 4)

# save.image(file='myEnvironment.RData')
