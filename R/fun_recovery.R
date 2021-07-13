#' Calculate recovery metrics from a time series with known disturbance date. The calcFrazier function derives the RRI, R80P and YrYr recovery indicators,
#' defined by Frazier et al. (2018). The indicators are originally developped for annual long-term time series of optical vegetation indices.
#' Yet, in order to be able to derive the indicators as well for dense and/or short time series, a modified version is suggested.
#' Here, the user can define the time period before, during and after the disturbance that is used to derive the indicators.
#' To reduce the interference of the seasonal pattern of dense time series, the chosen time period should cover blocks of n years.
#' (Frazier, R. J., Coops, N. C., Wulder, M. A., Hermosilla, T., & White, J. C. (2018). Analyzing spatial and temporal variability in short-term rates
#' of post-fire vegetation return from Landsat time series. Remote Sensing of Environment, 205, 32-45.)
#'
#' @param tsio vector of observations (time series with a fixed observation frequency)
#' @param tdist observation number of disturbance, indicating the timing of the disturbance
#' @param obspyr number of observations per year
#' @param nPre number of years prior to the disturbance used to calculate the pre-disturbance value
#' @param nDist number of years used to quantify the time series value during the disturbance
#' @param nPostMin the post-disturbance condition is quantified starting from nPostMin years after the disturbance
#' @param nPostMax max number of years after the disturbance used to quantify the post-disturbance condition
#'
#' @return a list containing the RRI recovery indicator, R80p recovery indicator and YrYr recovery indicator
#' @export
#'
calcFrazier <- function(tsio, tdist, obspyr, nPre, nDist, nPost, nPostStart, nDelta, nDeltaStart){
  # check if there are enough observations before and after the disturbance to calculate the metrics
  if ( (tdist > (nPre*obspyr) ) & ( (tdist + ((nPost+nPostStart)*obspyr) - 2) < length(tsio)  ) & ( sum(!is.na(tsio)) > 2) ) {
    # translate parameters to those needed for the recovery functions
    ys <- tsio# response values of time series
    ts <- seq(1, length(tsio))# observation numbers of time series
    # the observations during the perturbation
    if (nDist == 0 ){# if annual observatons or if duration of perturbation equals one time step
      tpert <- tdist
    }else if (obspyr == 1){
      tpert <- seq(tdist, tdist + nDist-1 )
    }else{
      tpert <- seq(tdist, tdist + nDist*obspyr - 1)
    }
    # the observations that represent the pre-disturbed state
    ts_pre <- seq(tdist - nPre*obspyr, tdist - 1)

    # the observations that represent the post-disturbed state
    ts_post <-  seq(tdist + (nPostStart*obspyr), tdist + ((nPostStart + nPost)*obspyr) - 1)

    # the observations that define the post-disturbance state for YrYr metric (relative to start disturbance)
    # the post-disturbance period for the YrYr metric is defined as the last year of the regular post-disturbance period
    deltat <- seq(nDeltaStart*obspyr,(nDelta+nDeltaStart)*obspyr-1)#ts_post[(length(ts_post) - obspyr+1):length(ts_post)]-tdist

    # Derive recovery indicators
    RRI <- rri(ts, ys, tpert, ts_pre, ts_post)
    R80P <- r80p(ts, ys, r = 0.8, ts_pre, ts_post)
    YrYr <- yryr(ts, ys, tpert, deltat)
    # make list of recovery indicators as output of the function
    lst <- list(RRI, R80P, YrYr)
    names(lst) <- c('RRI', 'R80P', 'YrYr')
    # give NA as output if not able to calculate the recovery indicators
  } else {
    lst <- list(NA, NA, NA)
    names(lst) <- c('RRI', 'R80P', 'YrYr')
  }
  lst
}
# calcFrazier <- function(tsio, tdist, obspyr, nPre, nDist, nPostMin, nPostMax){
#     # check if there are enough observations before and after the disturbance to calculate the metrics
#     if ( (tdist > (nPre*obspyr) ) & ( tdist < (length(tsio) - ((nPost+nPostStart)*obspyr) + 1) )  & ( sum(!is.na(tsio)) > 2) ) {
#       # translate parameters to those needed for the recovery functions
#       ys <- tsio# response values of time series
#       ts <- seq(1, length(tsio))# observation numbers of time series
#         # the observations during the perturbation
#       if (obspyr == 1 | nDist == 0 ){# if annual observatons or if duration of perturbation equals one time step
#         tpert <- seq(tdist, tdist + nDist*obspyr )
#       }else{
#         tpert <- seq(tdist, tdist + nDist*obspyr - 1)
#       }
#         # the observations that represent the pre-disturbed state
#       ts_pre <- seq(tdist - nPre*obspyr, tdist - 1)
#
#         # the observations that represent the post-disturbed state
#       if (obspyr == 1 | nPostMin == nPostMax) {
#         ts_post <-  seq(tdist + (nPostMin*obspyr), tdist + (nPostMax*obspyr))
#       } else {
#         ts_post <-  seq(tdist + (nPostMin*obspyr), tdist + (nPostMax*obspyr) - 1)
#       }
#         # the observations that define the post-disturbance state for YrYr metric (relative to start disturbance)
#       # the post-disturbance period for the YrYr metric is defined as the last year of the regular post-disturbance period
#       if(obspyr==1){
#         deltat <- nPostMax
#       }else{
#         deltat <- ts_post[(length(ts_post) - obspyr+1):length(ts_post)]-tdist
#       }
#         # Derive recovery indicators
#       RRI <- rri(ts, ys, tpert, ts_pre, ts_post)
#       R80P <- r80p(ts, ys, r = 0.8, ts_pre, ts_post)
#       YrYr <- yryr(ts, ys, tpert, deltat)
#       # make list of recovery indicators as output of the function
#       lst <- list(RRI, R80P, YrYr)
#       names(lst) <- c('RRI', 'R80P', 'YrYr')
#       # give NA as output if not able to calculate the recovery indicators
#     } else {
#       lst <- list(NA, NA, NA)
#       names(lst) <- c('RRI', 'R80P', 'YrYr')
#     }
#   lst
# }

#' Post-disturbance slope and recovery metrics derived from BFAST0n trend segments. The calcBFASTrec function derives a set of recovery indicators after fitting a segmented trend in the time series. Using the breakpoints function of the strucchange package, a segmented trend is fitted (hereafter called BFAST0n trend segments). The detected break showing the largest change (in absolute values) is assumed to represent the disturbance. Using the segmented trend and detected disturbance date, the RRI, R80p, YrYr and the slope of the post-disturbance trend segment are derived as recovery indicators.
#'
#' @param tsio vector of observations (time series)
#' @param obspyr number of observations in one year
#' @param h This parameter defines the minimal segment size either given as fraction relative to the sample size or as an integer giving the minimal number of observations in each segment.
#' @param nPre  number of years prior to the disturbance used to calculate the pre-disturbance value
#' @param nDist  number of months used to quantify the time series value during the disturbance
#' @param nPostMin  min number of years after the disturbance used to quantify the recovery
#' @param nPostMax  max number of years after the disturbance used to quantify the recovery
#' @param seas TRUE or FALSE, include seasonal term when detecting breaks?
#' @param breaks 'BIC' or 'LWZ': criteria used to define the optimal number of breaks (desactivated)
#'
#' @return a list containing  the RRI, R80p, YrYr recovery indicator derived from the BFAST0n trend segments and slope of the trend segment after the disturbance (sl).
#' @export
#' @import strucchange
#' @import stats
calcBFASTrec <- function(tsio, obspyr, h, nPre, nDist, nPost, nPostStart, nDelta, nDeltaStart, seas = F) {
  # Create time series object, needed as input for the piecewise regression
  tsi <- ts(tsio, frequency = obspyr)
  # Convert the time series object into a dataframe (with response, trend, harmonic), needed for the breakpoints function (dependent and independent variables in regression)
  if( obspyr>1 ) {# if temporal frequency is higher than 1 observation per year, provide the response, linear trend, and harmonic
    datapp <- bfastpp(tsi, order = 1, lag = NULL, slag = NULL,
                      na.action = na.omit, stl = 'none')
  } else if(!seas) {# if annual time series, generate a dataframe with response and linear trend
    datapp <- data.frame(response = tsio, trend = seq(1:length(tsio)))
  } else {
    stop('No seasonal term allowed for time series with one observation per year or less.')
  }

  # Test if enough observations are available to fit piecewise model
  nreg <- switch(seas+1, 2, 5)
  if(floor(length(tsio[is.na(tsio)==F]) * h) > nreg) {
    # set_fast_options()
    # Apply BFAST0n on time series: find breaks in the regression (ie fit piecewise model)
    if (seas) {
      bp <- breakpoints(response ~ trend + harmon, data = datapp, h = h)#, breaks = breaks
    } else {
      bp <- breakpoints(response ~ trend, data = datapp, h = h)##, breaks = breaks
    }
    # Check if BFAST0n found breakpoints
    if( is.na(bp$breakpoints[1]) ){
        # if no breakpoint was found, no recovery indicators can be derived (because no disturbance date could be detected by the break)
      frz <- list(NA, NA, NA, NA)
      names(frz) <- c('RRI', 'R80P', 'YrYr', 'Sl')
    }else{# at least one breakpoint found
      # Extract trend component and breaks
      cf <- coef(bp)#, breaks = breaks
      # Extract observation number of break
      tbp <- bp$breakpoints #observation number of break (for a time series where the NA values are removed)
      #tr <- rep(NA,length(tsi))
      # correct observation number of breaks for missing values
      indna <- which(is.na(tsi)==F)
      tbp <- indna[tbp]
      #tr[is.na(tsi)==F] <- fitted(bptst, length(tbptst))
      #Derive trend component
      bpf <- c(0, tbp, length(tsi))# boundaries (in terms of observation number) of each trend component
      trf <- rep(NA,length(tsi))# initialise the vector with the trend component
      for(ti in 1:(length(bpf)-1)){# iterate over the segments and derive its trend component
        trf[(bpf[ti]+1):bpf[ti+1]] <- cf[ti,1] + ((cf[ti,2]*((bpf[ti]+1):bpf[ti+1])))
      }
      # Find the major break (is assumed to represent the simulated disturbance date)
      dbr <- trf[tbp+1]-trf[tbp]
      tbp <- tbp[which(abs(dbr) == max(abs(dbr)))]
      # Calculate Frazier recovery metrics on trend component
      frz <- calcFrazier(as.numeric(trf), (tbp+1), floor(obspyr), nPre, nDist, nPost, nPostStart, nDelta, nDeltaStart)
      # Calculate the post-disturbance slope of the trend component (first segment after break)
      sl <- (trf[tbp+3] - trf[tbp+2])
      frz <- c(frz, sl)
      names(frz) <- c('RRI', 'R80P', 'YrYr', 'Sl')
    }
  } else {
    frz <- list(NA, NA, NA, NA)
    names(frz) <- c('RRI', 'R80P', 'YrYr', 'Sl')
  }

  frz
}

#
#' convert to matrix of performance indicators to a dataframe
#'
#' @param mat matrix of recovery indicators (rows represent the values of the tested time series characteristic, colums represent the different recovery indicators)
#' @param setvr time series characteristic being evaluated
#' @param metric recovery metric being evaluated
#' @param funSet settings list for the recovery analysis
#' @param input vector of preprocessing techniques per recovery indicator
#' @param nDist vector of the time span during the disturbance used to measure recovery
#'
#' @return dataframe
#' @export
toDF <- function(mat, setvr, metric, funSet){#freq, input, nPostMin, seas
  tst <- as.data.frame(t(mat))
  names(tst) <- setvr
  tst$Metric <- factor(metric)
  tst$Dense <- factor(funSet$freq)
  tst$Smooth <- factor(funSet$input)
  tst$h <- factor(funSet$h)
  tst$nDist <- factor(funSet$nDist)
  tst$nPre <- factor(funSet$nPre)
  tst$nPostStart <- factor(funSet$nPostStart)
  tst$nDeltaStart <- factor(funSet$nDeltaStart)
  tst$nPost <- factor(funSet$nPost)
  tst$nDelta <- factor(funSet$nDelta)
  # tst$nPostMin <- factor(funSet$nPostMin)#revalue(factor(nPostMin), c("1"="Short", "4"="Long"))#factor(recSttngs$nDist)#
  # tst$nPostMax <- factor(funSet$nPostMax)#revalue(factor(nPostMin), c("1"="Short", "4"="Long"))#factor(recSttngs$nDist)#
  tst$Seas <- factor(funSet$seas)
  tst
}

#' Run sensitivity analysis for particular parameter of interest
#'
#' @param evr char: name of parameter that will be evaluated in the simulation (for instance: "distT")
#' @param sttngs list of settings
#' @param funSet list of recovery indicator settings
#' @param basename basename for the output files
#' @param ofolder (optional) folder where the output files should be stored. If nothing is provided, no file is saved
#'
#' @return performance indicators (R2, RMSE, MAPE) and indicators of the relation between the measured and simulated recovery indicators (slope, intercept, p value of normality test of residuals). For each performace indicator a matrix is saved where the rows refer to the evaluated values of the paramter of interest and the columns to the various recovery indicator settings.
#' @export
#'
evalParam <- function(evr, sttngs, funSet, basename, ofolder = '') {

  # TODO:
  # The main to do is to pass `R80p`, `YrYr` and so on as parameters.
  #
  # Currently calcBFASTrec and calcFrazier calculate all the cases by default.
  # window size for smoothing
  winsize <- c(365, 4, 1)
  names(winsize) <- c('dense', 'quarterly', 'annual')

  # Extract case parameters
  pars <- setParamValues(sttngs)
  parvr <- pars[[evr]]

  # Initialize data containers
  empty <- matrix(NA, length(parvr), length(funSet[[1]]))
  RRI_rmse <- empty
  RRI_mape <- empty
  RRI_rsq <- empty
  RRI_nTS <- empty

  R80p_rmse <- empty
  R80p_mape <- empty
  R80p_rsq <- empty
  R80p_nTS <- empty

  YrYr_rmse <- empty
  YrYr_mape <- empty
  YrYr_rsq <- empty
  YrYr_nTS <- empty

  # iterate over values of evaluated parameter and simulate nrep time series per combination of all other variables
  for (i in 1:length(parvr)) {

    empty2 <- matrix(NA,nrow = length(funSet[[1]]), ncol = length(parvr[[1]][[1]]))
      # the recovery indicators extracted/measured from the time series (m stands for measured)
    m_RRIi <- empty2
    m_R80pi <-  empty2
    m_YrYri <-  empty2
    # m_SLi <-  empty2
    # the true value for the recovery indicators (s stands for simulated)
    s_RRIi <- empty2
    s_R80pi <-  empty2
    s_YrYri <-  empty2
    # s_SLi <-  empty2

    # iterate over the parameter settings and simulate each time a time series
    for (pari in 1: length(parvr[[1]][[1]])){ # TODO: what is this length? --> number of simulated time series per evaluated parameter value
      # simulate the pari'th time series of evaluated parameter i using the parameter settings given by the parameters in parvr
      sc <- simulCase(parvr[[i]]$nrep[pari], parvr[[i]]$nyr[pari], parvr[[i]]$nobsYr[pari], parvr[[i]]$nDr[pari], parvr[[i]]$seasAv[[1]], parvr[[i]]$seasAmp[pari],parvr[[i]]$trAv[pari], parvr[[i]]$remSd[pari], c(parvr[[i]]$distMag[pari],parvr[[i]]$distMag[pari]), parvr[[i]]$distT[pari], c(parvr[[i]]$distRec[pari],parvr[[i]]$distRec[pari]), parvr[[i]]$missVal[pari], parvr[[i]]$DistMissVal[pari], parvr[[i]]$distType[pari])

      # iterate over the recovery indicator settings
      for (rset in 1:length(funSet[[1]])){ # TODO: what is this length?

        # Extract simulation's key parameters
        # Note: Yes, this is the correct place for this block. See issue #46 (https://github.com/RETURN-project/BenchmarkRecovery/issues/46)
        tsi <- sc[[1]][1,] # simulated time series = seasonality + trend + noise component (contains missing values)
        tsseas <-sc[[2]][1,] # simulated seasonality component of time series  (contains no missing values)
        obspyr <- sc[[5]][1,]$obs_per_year # number of observations per year
        tdist <- sc[[5]][1,]$dist_time# timing of the disturbance (observation number)
        tsref <- sc[[3]][1,]# the simulated trend component = reference time series to measure the recovery indicators for validation (contains no missing values)
        nobs <- (sc[[5]][1,]$number_yrs)* obspyr# total number of observations = number of years * number of observations per year
        tm <- 1:nobs# the time is represented by the observation number

        # Most of the descriptions of this variables are in vignettes/sensitivity_analysis.Rmd
        # ============================================
        inp <- funSet$input[rset]# 'smoothed', 'raw', 'segmented'. Defines the type of time series that is used for the recovery indicators. For 'raw', the simulated time series are directly used to calculate recovery, for 'smooth' a time series smoothing algorithm is used before recovery calculation, for 'BFAST' trend segmentation (BFAST0n) is used.
        frq <- funSet[['freq']][rset]# 'dense', 'quarterly', or 'annual'. Defines the observation frequency. For 'dense' the original frequency is used. For 'annual' or 'quarterly', the time series are converted to annual or quarterly frequency, respectively.
        nPre <- funSet[['nPre']][rset]# the number of years before the disturbance used to derive the pre-disturbance values
        nDist <- funSet[['nDist']][rset]# the number of years after the disturbance used to derive the value during the disturbance
        # nPostMin <- funSet[['nPostMin']][rset]# the post-disturbance values are derived between nPostMin and nPostMax years after the disturbance
        # nPostMax <- funSet[['nPostMax']][rset]#  the post-disturbance values are derived between nPostMin and nPostMax years after the disturbance
        nPostStart = funSet[['nPostStart']][rset]
        nDeltaStart = funSet[['nDeltaStart']][rset]
        nPost = funSet[['nPost']][rset]
        nDelta = funSet[['nDelta']][rset]
        h <- funSet[['h']][rset]# only relevant if inp equals 'segmented', the h value is used in the piecewise regression to define the minimal segment size either given as fraction relative to the sample size or as an integer giving the minimal number of observations in each segment
        seas <- funSet[['seas']][rset]# only relevant if inp equals 'segmented', seas denotes whether a seasonal term needs to be used in the piecewise regression
        breaks <- funSet[['breaks']][rset]# only relevant if inp equals 'segmented', the criterium given by breaks is used in the piecewise regression to define the optimal number of segments. Can be set to 'BIC' or 'LWZ'. This option has been deactivated
        # ============================================

        # change temporal resolution
        if (frq == 'annual'){
          #convert time series to annual values by selecting date closest to seasonal max
          tsi <- toAnnualTS(tsseas, tsi, obspyr, dtmax = 2/12)
          tsref <- toAnnualTS(tsseas, tsref, obspyr, dtmax = 2/12)
          tdist <- which(tsref == min(tsref, na.rm = T))#ceiling(tdist/obspyr)
          obspyr <- 1
        } else if (frq == 'quarterly'){
          #convert time series to quarterly resolution
          dts <- seq(as.Date('2000-01-01'), by = '1 days', length = nobs)
          tsi <- toRegularTS(tsi, dts, 'mean', 'quart')
          tsref <- toRegularTS(tsref, dts, 'mean', 'quart')
          tdist <- which(tsref == min(tsref, na.rm = T))#ceiling(tdist/obspyr)
          obspyr <- 4
        }

        # If the noise removal technique is 'smoothed', use rolling mean to smooth the simulated time series
        if (inp == 'smoothed'){
            # smooth the time series using rolling mean with window size given by winsize
          temp.zoo<-zoo(tsi,(1:length(tsi)))# generate a zoo time series object
          m.av<-rollapply(temp.zoo, as.numeric(winsize[frq]), mean, na.rm = T, fill = NA)# smooth time series using rolling mean
          tsi <- as.numeric(m.av)
        }
        # calculate recovery indicators from simulated time series
        if((inp == 'smoothed') | (inp == 'raw')){
          outp <- calcFrazier(tsi, tdist, obspyr, nPre, nDist, nPost, nPostStart, nDelta, nDeltaStart)
          m_RRIi[rset,pari] <- outp$RRI# measured RRI
          m_R80pi[rset,pari] <- outp$R80P# measured R80p
          m_YrYri[rset,pari] <- outp$YrYr# measured YrYR
          rm(outp)
        } else if(funSet$input[rset] == 'segmented'){
          outp <- calcBFASTrec(tsi, obspyr, h, nPre, nDist, nPost, nPostStart, nDelta, nDeltaStart, seas)
          m_RRIi[rset,pari] <- outp$RRI# measured RRI
          m_R80pi[rset,pari] <- outp$R80P# measured R80p
          m_YrYri[rset,pari] <- outp$YrYr# measured YrYR
          rm(outp)
        }

        # calculate reference indicators (assumed truth for validation) based on simulated reference component
        outp <- calcFrazier(tsref, tdist, obspyr, nPre, nDist, nPost, nPostStart, nDelta, nDeltaStart)
        s_RRIi[rset,pari] <- outp$RRI#simulated (true) RRI
        s_R80pi[rset,pari] <- outp$R80P#simulated (true) R80p
        s_YrYri[rset,pari] <- outp$YrYr
      }
      rm(sc)
    }

    # evaluate recovery indicators
    # R2
    RRI_rsq[i,] <- sapply(1:dim(s_RRIi)[1], function(it) rsq(s_RRIi[it,], m_RRIi[it,]))
    R80p_rsq[i,] <- sapply(1:dim(s_R80pi)[1], function(it) rsq(s_R80pi[it,], m_R80pi[it,]))
    YrYr_rsq[i,] <- sapply(1:dim(s_YrYri)[1], function(it) rsq(s_YrYri[it,], m_YrYri[it,]))

    # MAPE
    RRI_mape[i,] <- sapply(1:dim(s_RRIi)[1], function(it) mape(s_RRIi[it,], m_RRIi[it,]))
    R80p_mape[i,] <- sapply(1:dim(s_R80pi)[1], function(it) mape(s_R80pi[it,], m_R80pi[it,]))
    YrYr_mape[i,] <- sapply(1:dim(s_YrYri)[1], function(it) mape(s_YrYri[it,], m_YrYri[it,]))

    # RMSE
    RRI_rmse[i,] <- sapply(1:dim(s_RRIi)[1], function(it) rmse(s_RRIi[it,], m_RRIi[it,]))
    R80p_rmse[i,] <- sapply(1:dim(s_R80pi)[1], function(it) rmse(s_R80pi[it,], m_R80pi[it,]))
    YrYr_rmse[i,] <- sapply(1:dim(s_YrYri)[1], function(it) rmse(s_YrYri[it,], m_YrYri[it,]))

    # nTS - the fraction of time series for which a recovery indicator could be derived (thus not a NA value)
    RRI_nTS[i,] <- apply(m_RRIi, 1, function(x){sum(is.na(x)==F)/length(x)})
    R80p_nTS[i,] <- apply(m_R80pi, 1, function(x){sum(is.na(x)==F)/length(x)})
    YrYr_nTS[i,] <- apply(m_YrYri, 1, function(x){sum(is.na(x)==F)/length(x)})

  }

  RRI_rsqDF <- toDF(RRI_rsq, names(parvr), 'RRI', funSet)
  R80p_rsqDF <- toDF(R80p_rsq, names(parvr), 'R80p', funSet)
  YrYr_rsqDF <- toDF(YrYr_rsq, names(parvr), 'YrYr', funSet)

  RRI_mapeDF <- toDF(RRI_mape, names(parvr), 'RRI', funSet)
  R80p_mapeDF <- toDF(R80p_mape, names(parvr), 'R80p', funSet)
  YrYr_mapeDF <- toDF(YrYr_mape, names(parvr), 'YrYr', funSet)

  RRI_rmseDF <- toDF(RRI_rmse, names(parvr), 'RRI', funSet)
  R80p_rmseDF <- toDF(R80p_rmse, names(parvr), 'R80p', funSet)
  YrYr_rmseDF <- toDF(YrYr_rmse, names(parvr), 'YrYr', funSet)

  RRI_nTSDF <- toDF(RRI_nTS, names(parvr), 'RRI', funSet)
  R80p_nTSDF <- toDF(R80p_nTS, names(parvr), 'R80p', funSet)
  YrYr_nTSDF <- toDF(YrYr_nTS, names(parvr), 'YrYr', funSet)

  # Save the performance indicators (if desired)
  save_results = (ofolder != '')
  if(save_results) {
    save(RRI_rsqDF, file = file.path(ofolder, paste0(basename, '_RRI_R2_' , evr, '.rda')))
    save(R80p_rsqDF, file = file.path(ofolder, paste0(basename, '_R80p_R2_' , evr, '.rda')))
    save(YrYr_rsqDF, file = file.path(ofolder, paste0(basename, '_YrYr_R2_' , evr, '.rda')))

    save(RRI_rmseDF, file = file.path(ofolder, paste0(basename, '_RRI_RMSE_' , evr, '.rda')))
    save(R80p_rmseDF, file = file.path(ofolder, paste0(basename, '_R80p_RMSE_' , evr, '.rda')))
    save(YrYr_rmseDF, file = file.path(ofolder, paste0(basename, '_YrYr_RMSE_' , evr, '.rda')))

    save(RRI_mapeDF, file = file.path(ofolder, paste0(basename, '_RRI_MAPE_' , evr, '.rda')))
    save(R80p_mapeDF, file = file.path(ofolder, paste0(basename, '_R80p_MAPE_' , evr, '.rda')))
    save(YrYr_mapeDF, file = file.path(ofolder, paste0(basename, '_YrYr_MAPE_' , evr, '.rda')))

    save(RRI_nTSDF, file = file.path(ofolder, paste0(basename, '_RRI_nTS_' , evr, '.rda')))
    save(R80p_nTSDF, file = file.path(ofolder, paste0(basename, '_R80p_nTS_' , evr, '.rda')))
    save(YrYr_nTSDF, file = file.path(ofolder, paste0(basename, '_YrYr_nTS_' , evr, '.rda')))
  }

  # Output the performance indicators
  return(
    list(RRI_rsq = RRI_rsqDF, R80p_rsq = R80p_rsqDF, YrYr_rsq = YrYr_rsqDF,
         RRI_rmse = RRI_rmseDF, R80p_rmse = R80p_rmseDF, YrYr_rmse = YrYr_rmseDF,
         RRI_mape = RRI_mapeDF, R80p_mape = R80p_mapeDF, YrYr_mape = YrYr_mapeDF,
         RRI_nTS = RRI_nTSDF, R80p_nTS = R80p_nTSDF, YrYr_nTS = YrYr_nTSDF)
  )
}

