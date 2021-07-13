context("Recovery indicators")

test_that("Frazier - annual - too short time series", {

  tsio <- c(rep(0,1), seq(-1, 0), rep(0,1))
  tdist <- 2
  obspyr <- 1
  shortDenseTS <- FALSE
  nPre <- 2
  nDist <- 1
  nPostStart <- 4
  nPost <- 2
  nDeltaStart <- 5
  nDelta <- 1

  metrics <- calcFrazier(tsio, tdist, obspyr, nPre, nDist, nPost, nPostStart, nDelta, nDeltaStart)

  expect_equal(metrics$RRI, NA)
  expect_equal(metrics$R80P, NA)
  expect_equal(metrics$YrYr, NA)
})

test_that("Frazier - annual", {

  tsio <- c(rep(1,2), seq(-5, -1), rep(-2,1))
  tdist <- 3
  obspyr <- 1
  shortDenseTS <- FALSE
  nPre <- 2
  nDist <- 1
  nPostStart <- 4
  nPost <- 2
  nDeltaStart <- 5
  nDelta <- 1

  metrics <- calcFrazier(tsio, tdist, obspyr, nPre, nDist, nPost, nPostStart, nDelta, nDeltaStart)
  pre <- 1
  dnbr <- 6
  ari <- 4

  rrim <- ari/dnbr
  r80pm <- -1/(0.8*pre)
  yryrm <- (-2 + 5)/5

  expect_equal(metrics$RRI, rrim, tolerance = 1e-4)
  expect_equal(metrics$R80P, r80pm, tolerance = 1e-4)
  expect_equal(metrics$YrYr, yryrm, tolerance = 1e-4)
})

test_that("Frazier - dense", {

  tsio <- c(rep(1,24), seq(-5, -1, length.out=72), rep(-2,12))
  tdist <- 25
  obspyr <- 12
  shortDenseTS <- TRUE
  nPre <- 2
  nDist <- 1
  nPostStart <- 4
  nPost <- 2
  nDeltaStart <- 5
  nDelta <- 1

  metrics <- calcFrazier(tsio, tdist, obspyr, nPre, nDist, nPost, nPostStart, nDelta, nDeltaStart)
  pre <- 1
  dist <- mean(tsio[25:36])
  post <- max(tsio[73:96])
  dnbr <- pre-dist
  ari <- post-dist

  rrim <- ari/dnbr
  r80pm <- post/(0.8*pre)
  yryrm <- (mean(tsio[85:96]) - dist)/(5*12)

  expect_equal(metrics$RRI, rrim, tolerance = 1e-4)
  expect_equal(metrics$R80P, r80pm, tolerance = 1e-4)
  expect_equal(metrics$YrYr, yryrm, tolerance = 1e-4)
})

test_that("Frazier - segmented", {

  tsio <- c(rep(1,24), seq(-5, -1, length.out=60), rep(-2,12))
  tdist <- 25
  obspyr <- 12
  nPre <- 2
  nDist <- 1
  nPost <- 1
  nPostStart <- 4
  nDelta <- 1
  nDeltaStart <- 4
  # nPostMin <- 4
  # nPostMax <- 5
  h <- 0.1

  metrics <- calcBFASTrec(tsio, obspyr, h, nPre, nDist, nPost, nPostStart, nDelta, nDeltaStart)
  pre <- 1
  dist <- mean(tsio[25:36])
  post <- max(tsio[73:84])
  dnbr <- pre-dist
  ari <- post-dist

  rri <- ari/dnbr
  r80p <- post/(0.8*pre)
  yryr <- (mean(tsio[73:84]) - dist)/(4*12)

  expect_equal(metrics$RRI, rri, tolerance = 1e-4)
  expect_equal(metrics$R80P, r80p, tolerance = 1e-4)
  expect_equal(metrics$YrYr, yryr, tolerance = 1e-4)
})


test_that("Frazier - segmented annual - long", {

  tsio <- c(rep(1,8), seq(-5, 0, by = 0.5), rep(0,8))
  tdist <- 9
  obspyr <- 1
  nPre <- 2
  nDist <- 1
  nPost <- 2
  nPostStart <- 4
  nDelta <- 1
  nDeltaStart <- 5
  h <- 0.2
  seas <- F

  metrics <- calcBFASTrec(tsio, obspyr, h, nPre, nDist, nPost, nPostStart, nDelta, nDeltaStart, seas = F)
  pre <- 1
  dnbr <- 6
  ari <- 2.5

  rrim <- ari/dnbr
  r80pm <- tsio[14]/(0.8*pre)
  yryrm <- (tsio[14] - tsio[9])/5

  expect_equal(metrics$RRI, rrim, tolerance = 1e-4)
  expect_equal(metrics$R80P, r80pm, tolerance = 1e-4)
  expect_equal(metrics$YrYr, yryrm, tolerance = 1e-4)
})

test_that("Frazier - segmented annual - short", {

  tsio <- c(rep(1,8), seq(-5, 0, by = 0.5), rep(0,8))
  tdist <- 9
  obspyr <- 1
  nPre <- 2
  nDist <- 1
  nPost <- 1
  nPostStart <- 1
  nDelta <- 1
  nDeltaStart <- 1
  # nPostMin <- 1
  # nPostMax <- 1
  h <- 0.2
  seas <- F

  metrics <- calcBFASTrec(tsio, obspyr, h, nPre, nDist, nPost, nPostStart, nDelta, nDeltaStart, seas = F)
  pre <- 1
  dnbr <- 6
  ari <- 0.5

  rrim <- ari/dnbr
  r80pm <- tsio[10]/(0.8*pre)
  yryrm <- (tsio[10] - tsio[9])/1

  expect_equal(metrics$RRI, rrim, tolerance = 1e-4)
  expect_equal(metrics$R80P, r80pm, tolerance = 1e-4)
  expect_equal(metrics$YrYr, yryrm, tolerance = 1e-4)
})

test_that("Evaluate recovery indicators - perfect fit", {
  #-------------------------------------------------
  #  settings  simulation

  basename <- 'test'
  # generate some time series characteristics - no noise, no missing values and small seasonal amplitude
  set.seed(197)
  STnobsYr <- 365# observations per year
  Vqntl <- c(0.75,.95)# set of quantiles used to derive realistic values
  seasVImax <- rep(0,100) #runif(100, 0.1, 0.2)# seasonal amplitude
  tsVIMissVal <- rep(0,100)#runif(100,1/12,6/12)# missing values
  Rem_VIsd <- rep(0,100)#runif(100,0.00001,0.00002)# sd remainder
  TrVImean <- 0.7# offset
  seasVImean <- 0.01* rep(sin(seq(0,pi*2,pi/183))[1:365],6)# seasonal pattern

  # generate a settings list to simulate time series
  sttngs <- list()
  sttngs$'seasAmp' <- list(type = 'dist', vals = c(0,0,quantile(seasVImax, Vqntl)), obs = seasVImax, fix = quantile(seasVImax, c(.4, .6)))
  sttngs$'missVal' <- list(type = 'dist',  vals = c(1-1/16,1-1/16,quantile(tsVIMissVal, Vqntl)), obs = tsVIMissVal, fix = quantile(tsVIMissVal, c(.4, .6)))
  sttngs$'remSd' <- list(type = 'dist',  vals = c(0,0,quantile(Rem_VIsd, Vqntl)), obs = Rem_VIsd, fix = quantile(Rem_VIsd, c(.4, .6)))
  sttngs$'nyr' <- list(type = 'range', vals = c(25), fix = 25)#seq(6,36, by = 6)
  sttngs$'distMag' <- list(type = 'range', vals = -c(0.5, 0.5), fix = -c(0.25,0.35))
  sttngs$'distT' <- list(type = 'range', vals = c(15,16), fix = c(15,16))#seq(3,33, by = 6)
  sttngs$'distRec' <- list(type = 'range', vals = c(0.5,2,4,6.5)*STnobsYr, fix = c(2.25,3.75)*STnobsYr) #seq(0.5,6.5,by=0.5)
  sttngs$'nDr' <- list(type = 'range',  vals = c(0), fix = 0)
  sttngs$'distType' <- list(type = 'cat',  vals = c('piecewise'), fix = 'piecewise')#piecewise
  sttngs$'DistMissVal' <- list(type = 'cat', vals = 'random', fix = 'random')
  sttngs$'trAv' <- list(type = 'range', vals = TrVImean, fix = TrVImean)
  sttngs$'general' <- list(
    eval = c('distMag'),#parameters to be evaluated, can be  'distT','distRec', 'missVal' 'distMag', 'seasAmp', 'remSd'
    nTS = 10,
    nobsYr = STnobsYr,
    seasAv = seasVImean,
    parSetUp = 'int')#parameter set-up: can be avg dist, comb, or int

  # settings of the recovery indicators
  funSet <- list('freq' = 'dense',#rep('dense',3),
                 'input' = 'raw',#c('raw', 'smoothed','segmented'),# settings for the recovery indicators
                 'nPre' = 2,# rep(2,3),
                 'nDist' = 1,#rep(1,3),
                 'nPost' = 2,
                 'nPostStart' = 4,
                 'nDelta' = 1,
                 'nDeltaStart' = 5,
                 'h' = 0.15,#rep(0.15,3),
                 'seas' = T)#rep(T,3))
  # evaluate recovery indicators
  perf <- evalParam('distMag', sttngs, funSet, basename, ofolder = '')

  # we expect a nearly perfect performance
  expect_equal(length(perf), 12, tolerance = 1e-4)
  expect_equal(perf$RRI_nTS$`-0.5`, 1, tolerance = 1e-4)
  expect_equal(perf$R80p_nTS$`-0.5`, 1, tolerance = 1e-4)
  expect_equal(perf$YrYr_nTS$`-0.5`, 1, tolerance = 1e-4)
  expect_equal(perf$RRI_rsq$`-0.5`, 1, tolerance = 1e-4)
  expect_equal(perf$R80p_rsq$`-0.5`, 1, tolerance = 1e-4)
  expect_equal(perf$YrYr_rsq$`-0.5`, 1, tolerance = 1e-4)
  expect_equal(perf$RRI_rmse$`-0.5`, 0, tolerance = 1e-4)
  expect_equal(perf$R80p_rmse$`-0.5`, 0, tolerance = 1e-4)
  expect_equal(perf$YrYr_rmse$`-0.5`, 0, tolerance = 1e-4)
  expect_equal(perf$RRI_mape$`-0.5`, 0, tolerance = 1e-4)
  expect_equal(perf$R80p_mape$`-0.5`, 0, tolerance = 1e-4)
  expect_equal(perf$YrYr_mape$`-0.5`, 0, tolerance = 1e-4)
})

test_that("Evaluate recovery indicators - perfect fit", {
  #-------------------------------------------------
  #  settings  simulation
  basename <- 'test'
  # generate some time series characteristics - no noise, no missing values and small seasonal amplitude
  set.seed(197)
  STnobsYr <- 365# observations per year
  Vqntl <- c(0.75,.95)# set of quantiles used to derive realistic values
  seasVImax <- rep(0,100) #runif(100, 0.1, 0.2)# seasonal amplitude
  tsVIMissVal <- rep(0,100)#runif(100,1/12,6/12)# missing values
  Rem_VIsd <- rep(0,100)#runif(100,0.00001,0.00002)# sd remainder
  TrVImean <- 0.7# offset
  seasVImean <- 0.01* rep(sin(seq(0,pi*2,pi/183))[1:365],6)# seasonal pattern

  # generate a settings list to simulate time series
  sttngs <- list()
  sttngs$'seasAmp' <- list(type = 'dist', vals = c(0,0,quantile(seasVImax, Vqntl)), obs = seasVImax, fix = quantile(seasVImax, c(.4, .6)))
  sttngs$'missVal' <- list(type = 'dist',  vals = c(1-1/16,1-1/16,quantile(tsVIMissVal, Vqntl)), obs = tsVIMissVal, fix = quantile(tsVIMissVal, c(.4, .6)))
  sttngs$'remSd' <- list(type = 'dist',  vals = c(0,0,quantile(Rem_VIsd, Vqntl)), obs = Rem_VIsd, fix = quantile(Rem_VIsd, c(.4, .6)))
  sttngs$'nyr' <- list(type = 'range', vals = c(25), fix = 25)#seq(6,36, by = 6)
  sttngs$'distMag' <- list(type = 'range', vals = -c(0.5, 0.5), fix = -c(0.25,0.35))
  sttngs$'distT' <- list(type = 'range', vals = c(15,16), fix = c(15,16))#seq(3,33, by = 6)
  sttngs$'distRec' <- list(type = 'range', vals = c(0.5,2,4,6.5)*STnobsYr, fix = c(2.25,3.75)*STnobsYr) #seq(0.5,6.5,by=0.5)
  sttngs$'nDr' <- list(type = 'range',  vals = c(0), fix = 0)
  sttngs$'distType' <- list(type = 'cat',  vals = c('piecewise'), fix = 'piecewise')#piecewise
  sttngs$'DistMissVal' <- list(type = 'cat', vals = 'random', fix = 'random')
  sttngs$'trAv' <- list(type = 'range', vals = TrVImean, fix = TrVImean)
  sttngs$'general' <- list(
    eval = c('distMag'),#parameters to be evaluated, can be  'distT','distRec', 'missVal' 'distMag', 'seasAmp', 'remSd'
    nTS = 10,
    nobsYr = STnobsYr,
    seasAv = seasVImean,
    parSetUp = 'int')#parameter set-up: can be avg dist, comb, or int

  # settings of the recovery indicators
  funSet <- list('freq' = 'dense',#rep('dense',3),
                 'input' = 'raw',#c('raw', 'smoothed','segmented'),# settings for the recovery indicators
                 'shortDenseTS' = T,# rep(TRUE,3),
                 'nPre' = 2,# rep(2,3),
                 'nDist' = 1,#rep(1,3),
                 'nPost' = 2,
                 'nPostStart' = 4,
                 'nDelta' = 1,
                 'nDeltaStart' = 5,
                 'nPostMin' = 4,#c(4,4,4),
                 'nPostMax' = 6,#rep(6,3),
                 'h' = 0.15,#rep(0.15,3),
                 'seas' = T)#rep(T,3))
  # evaluate recovery indicators
  perf <- evalParam('distMag', sttngs, funSet, basename, ofolder = '')

  # we expect a nearly perfect performance
  expect_equal(length(perf), 12, tolerance = 1e-4)
  expect_equal(perf$RRI_nTS$`-0.5`, 1, tolerance = 1e-4)
  expect_equal(perf$R80p_nTS$`-0.5`, 1, tolerance = 1e-4)
  expect_equal(perf$YrYr_nTS$`-0.5`, 1, tolerance = 1e-4)
  expect_equal(perf$RRI_rsq$`-0.5`, 1, tolerance = 1e-4)
  expect_equal(perf$R80p_rsq$`-0.5`, 1, tolerance = 1e-4)
  expect_equal(perf$YrYr_rsq$`-0.5`, 1, tolerance = 1e-4)
  expect_equal(perf$RRI_rmse$`-0.5`, 0, tolerance = 1e-4)
  expect_equal(perf$R80p_rmse$`-0.5`, 0, tolerance = 1e-4)
  expect_equal(perf$YrYr_rmse$`-0.5`, 0, tolerance = 1e-4)
  expect_equal(perf$RRI_mape$`-0.5`, 0, tolerance = 1e-4)
  expect_equal(perf$R80p_mape$`-0.5`, 0, tolerance = 1e-4)
  expect_equal(perf$YrYr_mape$`-0.5`, 0, tolerance = 1e-4)
})

test_that("Evaluate recovery indicators - temporal aggregation to quarterly time series", {
  #-------------------------------------------------
  #  settings  simulation
  basename <- 'test'
  # generate some time series characteristics - no noise, no missing values and small seasonal amplitude
  set.seed(197)
  STnobsYr <- 365# observations per year
  Vqntl <- c(0.75,.95)# set of quantiles used to derive realistic values
  seasVImax <- runif(100, 0.1, 0.2)# seasonal amplitude
  tsVIMissVal <- runif(100,1/12,6/12)# missing values
  Rem_VIsd <- runif(100,0.01,0.02)# sd remainder
  TrVImean <- 0.7# offset
  seasVImean <- 0.01* rep(sin(seq(0,pi*2,pi/183))[1:365],6)# seasonal pattern

  # generate a settings list to simulate time series
  sttngs <- list()
  sttngs$'seasAmp' <- list(type = 'dist', vals = c(0,0,quantile(seasVImax, Vqntl)), obs = seasVImax, fix = quantile(seasVImax, c(.4, .6)))
  sttngs$'missVal' <- list(type = 'dist',  vals = c(1-1/16,1-1/16,quantile(tsVIMissVal, Vqntl)), obs = tsVIMissVal, fix = quantile(tsVIMissVal, c(.4, .6)))
  sttngs$'remSd' <- list(type = 'dist',  vals = c(0,0,quantile(Rem_VIsd, Vqntl)), obs = Rem_VIsd, fix = quantile(Rem_VIsd, c(.4, .6)))
  sttngs$'nyr' <- list(type = 'range', vals = c(25), fix = 25)#seq(6,36, by = 6)
  sttngs$'distMag' <- list(type = 'range', vals = -c(0.5, 0.5), fix = -c(0.25,0.35))
  sttngs$'distT' <- list(type = 'range', vals = c(15,16), fix = c(15,16))#seq(3,33, by = 6)
  sttngs$'distRec' <- list(type = 'range', vals = c(0.5,2,4,6.5)*STnobsYr, fix = c(2.25,3.75)*STnobsYr) #seq(0.5,6.5,by=0.5)
  sttngs$'nDr' <- list(type = 'range',  vals = c(0), fix = 0)
  sttngs$'distType' <- list(type = 'cat',  vals = c('piecewise'), fix = 'piecewise')#piecewise
  sttngs$'DistMissVal' <- list(type = 'cat', vals = 'random', fix = 'random')
  sttngs$'trAv' <- list(type = 'range', vals = TrVImean, fix = TrVImean)
  sttngs$'general' <- list(
    eval = c('distMag'),#parameters to be evaluated, can be  'distT','distRec', 'missVal' 'distMag', 'seasAmp', 'remSd'
    nTS = 1000,
    nobsYr = STnobsYr,
    seasAv = seasVImean,
    parSetUp = 'int')#parameter set-up: can be avg dist, comb, or int

  # settings of the recovery indicators
  funSet <- list('freq' = 'quarterly',#rep('dense',3),
                 'input' = 'raw',#c('raw', 'smoothed','segmented'),# settings for the recovery indicators
                 'shortDenseTS' = T,# rep(TRUE,3),
                 'nPre' = 2,# rep(2,3),
                 'nDist' = 1,#rep(1,3),
                 'nPost' = 2,
                 'nPostStart' = 4,
                 'nDelta' = 1,
                 'nDeltaStart' = 5,
                 'h' = 0.15,#rep(0.15,3),
                 'seas' = T)#rep(T,3))
  # evaluate recovery indicators
  perf <- evalParam('distMag', sttngs, funSet, basename, ofolder = '')

  # pars <- setParamValues(sttngs)
  # m_RRIi <- rep(NA,1000)
  # m_R80pi <- rep(NA,1000)
  # m_YrYri <- rep(NA,1000)
  # s_RRIi <- rep(NA,1000)
  # s_R80pi <- rep(NA,1000)
  # s_YrYri <- rep(NA,1000)
  # for (pari in 1: 1000){
  #   sc <- simulCase(pars[[1]][[1]]$nrep[pari], pars[[1]][[1]]$nyr[pari], pars[[1]][[1]]$nobsYr[pari], pars[[1]][[1]]$nDr[pari], pars[[1]][[1]]$seasAv[[1]], pars[[1]][[1]]$seasAmp[pari],pars[[1]][[1]]$trAv[pari], pars[[1]][[1]]$remSd[pari], c(pars[[1]][[1]]$distMag[pari],pars[[1]][[1]]$distMag[pari]), pars[[1]][[1]]$distT[pari], c(pars[[1]][[1]]$distRec[pari],pars[[1]][[1]]$distRec[pari]), pars[[1]][[1]]$missVal[pari], pars[[1]][[1]]$DistMissVal[pari], pars[[1]][[1]]$distType[pari])
  #   tsi <- sc[[1]][1,] # simulated time series = seasonality + trend + noise component (contains missing values)
  #   obspyr <- sc[[5]][1,]$obs_per_year # number of observations per year
  #   tdist <- sc[[5]][1,]$dist_time# timing of the disturbance (observation number)
  #   tsref <- sc[[3]][1,]# the simulated trend componend = reference time series to measure the recovery indicators for validation (contains no missing values)
  #   nobs <- (sc[[5]][1,]$number_yrs)* obspyr# total number of observations = number of years * number of observations per year
  #
  #   # temporal aggregation
  #   dts <- seq(as.Date('2000-01-01'), by = '1 days', length = nobs)
  #   tsi <- toRegularTS(tsi, dts, 'mean', 'quart')
  #   tsref <- toRegularTS(tsref, dts, 'mean', 'quart')
  #   tdist <- which(tsref == min(tsref, na.rm = T))#ceiling(tdist/obspyr)
  #   obspyr <- 4
  #
  #   outp <- calcFrazier(tsi, tdist, obspyr, funSet[['nPre']], funSet[['nDist']], funSet[['nPost']], funSet[['nPostStart']], funSet[['nDelta']], funSet[['nDeltaStart']])
  #   m_RRIi[pari] <- outp$RRI# measured RRI
  #   m_R80pi[pari] <- outp$R80P# measured R80p
  #   m_YrYri[pari] <- outp$YrYr# measured YrYR
  #   outp <- calcFrazier(tsref, tdist, obspyr, funSet[['nPre']], funSet[['nDist']], funSet[['nPost']], funSet[['nPostStart']], funSet[['nDelta']], funSet[['nDeltaStart']])
  #   s_RRIi[pari] <- outp$RRI#simulated (true) RRI
  #   s_R80pi[pari] <- outp$R80P#simulated (true) R80p
  #   s_YrYri[pari] <- outp$YrYr
  # }


# check if the results agree with expectations
expect_equal(length(perf), 12, tolerance = 1e-4)
expect_equal(perf$RRI_nTS$`-0.5`, 1, tolerance = 1e-2)
expect_equal(perf$R80p_nTS$`-0.5`, 1, tolerance = 1e-2)
expect_equal(perf$YrYr_nTS$`-0.5`, 1, tolerance = 1e-2)

expect_equal(perf$RRI_rsq$`-0.5`, 0.883, tolerance = 1e-2)
expect_equal(perf$R80p_rsq$`-0.5`, 0.875, tolerance = 1e-2)
expect_equal(perf$YrYr_rsq$`-0.5`, 0.997, tolerance = 1e-2)
expect_equal(perf$RRI_rmse$`-0.5`, 0.216, tolerance = 1e-2)
expect_equal(perf$R80p_rmse$`-0.5`, 0.167, tolerance = 1e-2)
expect_equal(perf$YrYr_rmse$`-0.5`, 0.0001, tolerance = 1e-2)
expect_equal(perf$RRI_mape$`-0.5`, 0.208, tolerance = 1e-2)
expect_equal(perf$R80p_mape$`-0.5`, 0.129, tolerance = 1e-2)
expect_equal(perf$YrYr_mape$`-0.5`, 0.004, tolerance = 1e-2)
})

# test_that("Evaluate recovery indicators - temporal aggregation to annual time series", {
  #-------------------------------------------------
  #  settings  simulation
  basename <- 'test'
  # generate some time series characteristics - no noise, no missing values and small seasonal amplitude
  set.seed(197)
  STnobsYr <- 365# observations per year
  Vqntl <- c(0.75,.95)# set of quantiles used to derive realistic values
  seasVImax <- runif(100, 0.1, 0.2)# seasonal amplitude
  tsVIMissVal <- runif(100,1/12,6/12)# missing values
  Rem_VIsd <- runif(100,0.01,0.02)# sd remainder
  TrVImean <- 0.7# offset
  seasVImean <- 0.01* rep(sin(seq(0,pi*2,pi/183))[1:365],6)# seasonal pattern

  # generate a settings list to simulate time series
  sttngs <- list()
  sttngs$'seasAmp' <- list(type = 'dist', vals = c(0,0,quantile(seasVImax, Vqntl)), obs = seasVImax, fix = quantile(seasVImax, c(.4, .6)))
  sttngs$'missVal' <- list(type = 'dist',  vals = c(1-1/16,1-1/16,quantile(tsVIMissVal, Vqntl)), obs = tsVIMissVal, fix = quantile(tsVIMissVal, c(.4, .6)))
  sttngs$'remSd' <- list(type = 'dist',  vals = c(0,0,quantile(Rem_VIsd, Vqntl)), obs = Rem_VIsd, fix = quantile(Rem_VIsd, c(.4, .6)))
  sttngs$'nyr' <- list(type = 'range', vals = c(25), fix = 25)#seq(6,36, by = 6)
  sttngs$'distMag' <- list(type = 'range', vals = -c(0.5, 0.5), fix = -c(0.25,0.35))
  sttngs$'distT' <- list(type = 'range', vals = c(15,16), fix = c(15,16))#seq(3,33, by = 6)
  sttngs$'distRec' <- list(type = 'range', vals = c(0.5,2,4,6.5)*STnobsYr, fix = c(2.25,3.75)*STnobsYr) #seq(0.5,6.5,by=0.5)
  sttngs$'nDr' <- list(type = 'range',  vals = c(0), fix = 0)
  sttngs$'distType' <- list(type = 'cat',  vals = c('piecewise'), fix = 'piecewise')#piecewise
  sttngs$'DistMissVal' <- list(type = 'cat', vals = 'random', fix = 'random')
  sttngs$'trAv' <- list(type = 'range', vals = TrVImean, fix = TrVImean)
  sttngs$'general' <- list(
    eval = c('distMag'),#parameters to be evaluated, can be  'distT','distRec', 'missVal' 'distMag', 'seasAmp', 'remSd'
    nTS = 1000,
    nobsYr = STnobsYr,
    seasAv = seasVImean,
    parSetUp = 'int')#parameter set-up: can be avg dist, comb, or int

  # settings of the recovery indicators
  funSet <- list('freq' = 'annual',#rep('dense',3),
                 'input' = 'raw',#c('raw', 'smoothed','segmented'),# settings for the recovery indicators
                 'shortDenseTS' = T,# rep(TRUE,3),
                 'nPre' = 2,# rep(2,3),
                 'nDist' = 1,#rep(1,3),
                 'nPost' = 2,
                 'nPostStart' =4,
                 'nDelta' = 1,
                 'nDeltaStart'=5,
                 # 'nPostMin' = 4,#c(4,4,4),
                 # 'nPostMax' = 6,#rep(6,3),
                 'h' = 0.15,#rep(0.15,3),
                 'seas' = T)#rep(T,3))
  # evaluate recovery indicators
  perf <- evalParam('distMag', sttngs, funSet, basename, ofolder = '')

  pars <- setParamValues(sttngs)
  m_RRIi <- rep(NA,1000)
  m_R80pi <- rep(NA,1000)
  m_YrYri <- rep(NA,1000)
  s_RRIi <- rep(NA,1000)
  s_R80pi <- rep(NA,1000)
  s_YrYri <- rep(NA,1000)
  for (pari in 1: 1000){
    sc <- simulCase(pars[[1]][[1]]$nrep[pari], pars[[1]][[1]]$nyr[pari], pars[[1]][[1]]$nobsYr[pari], pars[[1]][[1]]$nDr[pari], pars[[1]][[1]]$seasAv[[1]], pars[[1]][[1]]$seasAmp[pari],pars[[1]][[1]]$trAv[pari], pars[[1]][[1]]$remSd[pari], c(pars[[1]][[1]]$distMag[pari],pars[[1]][[1]]$distMag[pari]), pars[[1]][[1]]$distT[pari], c(pars[[1]][[1]]$distRec[pari],pars[[1]][[1]]$distRec[pari]), pars[[1]][[1]]$missVal[pari], pars[[1]][[1]]$DistMissVal[pari], pars[[1]][[1]]$distType[pari])
    tsi <- sc[[1]][1,] # simulated time series = seasonality + trend + noise component (contains missing values)
    tsseas <-sc[[2]][1,] # simulated seasonality component of time series  (contains no missing values)
    obspyr <- sc[[5]][1,]$obs_per_year # number of observations per year
    tdist <- sc[[5]][1,]$dist_time# timing of the disturbance (observation number)
    tsref <- sc[[3]][1,]# the simulated trend componend = reference time series to measure the recovery indicators for validation (contains no missing values)
    nobs <- (sc[[5]][1,]$number_yrs)* obspyr# total number of observations = number of years * number of observations per year

    # temporal aggregation
    tsi <- toAnnualTS(tsseas, tsi, obspyr, dtmax = 2/12)
    tsref <- toAnnualTS(tsseas, tsref, obspyr, dtmax = 2/12)
    tdist <- which(tsref == min(tsref, na.rm = T))#ceiling(tdist/obspyr)
    obspyr <- 1

    outp <- calcFrazier(tsi, tdist, obspyr, funSet[['nPre']], funSet[['nDist']], funSet[['nPost']], funSet[['nPostStart']], funSet[['nDelta']], funSet[['nDeltaStart']])
    m_RRIi[pari] <- outp$RRI# measured RRI
    m_R80pi[pari] <- outp$R80P# measured R80p
    m_YrYri[pari] <- outp$YrYr# measured YrYR
    outp <- calcFrazier(tsref, tdist, obspyr,  funSet[['nPre']], funSet[['nDist']], funSet[['nPost']], funSet[['nPostStart']], funSet[['nDelta']], funSet[['nDeltaStart']])
    s_RRIi[pari] <- outp$RRI#simulated (true) RRI
    s_R80pi[pari] <- outp$R80P#simulated (true) R80p
    s_YrYri[pari] <- outp$YrYr
  }
  rsq(s_RRIi, m_RRIi)
  rsq(s_R80pi, m_R80pi)

  # check if the results agree with expectations
  expect_equal(length(perf), 12, tolerance = 1e-4)
  expect_equal(perf$RRI_nTS$`-0.5`, 1, tolerance = 1e-2)
  expect_equal(perf$R80p_nTS$`-0.5`, 1, tolerance = 1e-2)
  expect_equal(perf$YrYr_nTS$`-0.5`, 1, tolerance = 1e-2)

  # expect_equal(perf$RRI_rsq$`-0.5`, 0.624, tolerance = 1e-2)
  # rexpect_equal(perf$R80p_rsq$`-0.5`, 0.654, tolerance = 1e-2)
  # expect_equal(perf$YrYr_rsq$`-0.5`, 0.794, tolerance = 5e-2)
  expect_equal(perf$RRI_rmse$`-0.5`, 0.0458, tolerance = 1e-2)
  expect_equal(perf$R80p_rmse$`-0.5`, 0.028, tolerance = 1e-2)
  expect_equal(perf$YrYr_rmse$`-0.5`, 0.0034, tolerance = 1e-2)
  expect_equal(perf$RRI_mape$`-0.5`, 0.037, tolerance = 1e-2)
  expect_equal(perf$R80p_mape$`-0.5`, 0.018, tolerance = 1e-2)
  expect_equal(perf$YrYr_mape$`-0.5`, 0.0314, tolerance = 1e-2)
})

test_that("Evaluate recovery indicators - smoothed time series", {
  #-------------------------------------------------

  #  settings  simulation
  basename <- 'test'
  # generate some time series characteristics - no noise, no missing values and small seasonal amplitude
  set.seed(197)
  STnobsYr <- 365# observations per year
  Vqntl <- c(0.75,.95)# set of quantiles used to derive realistic values
  seasVImax <- runif(100, 0.1, 0.2)# seasonal amplitude
  tsVIMissVal <- runif(100,1/12,6/12)# missing values
  Rem_VIsd <- runif(100,0.01,0.02)# sd remainder
  TrVImean <- 0.7# offset
  seasVImean <- 0.01* rep(sin(seq(0,pi*2,pi/183))[1:365],6)# seasonal pattern

  # generate a settings list to simulate time series
  sttngs <- list()
  sttngs$'seasAmp' <- list(type = 'dist', vals = c(0,0,quantile(seasVImax, Vqntl)), obs = seasVImax, fix = quantile(seasVImax, c(.4, .6)))
  sttngs$'missVal' <- list(type = 'dist',  vals = c(1-1/16,1-1/16,quantile(tsVIMissVal, Vqntl)), obs = tsVIMissVal, fix = quantile(tsVIMissVal, c(.4, .6)))
  sttngs$'remSd' <- list(type = 'dist',  vals = c(0,0,quantile(Rem_VIsd, Vqntl)), obs = Rem_VIsd, fix = quantile(Rem_VIsd, c(.4, .6)))
  sttngs$'nyr' <- list(type = 'range', vals = c(25), fix = 25)#seq(6,36, by = 6)
  sttngs$'distMag' <- list(type = 'range', vals = -c(0.5, 0.5), fix = -c(0.25,0.35))
  sttngs$'distT' <- list(type = 'range', vals = c(15,16), fix = c(15,16))#seq(3,33, by = 6)
  sttngs$'distRec' <- list(type = 'range', vals = c(0.5,2,4,6.5)*STnobsYr, fix = c(2.25,3.75)*STnobsYr) #seq(0.5,6.5,by=0.5)
  sttngs$'nDr' <- list(type = 'range',  vals = c(0), fix = 0)
  sttngs$'distType' <- list(type = 'cat',  vals = c('piecewise'), fix = 'piecewise')#piecewise
  sttngs$'DistMissVal' <- list(type = 'cat', vals = 'random', fix = 'random')
  sttngs$'trAv' <- list(type = 'range', vals = TrVImean, fix = TrVImean)
  sttngs$'general' <- list(
    eval = c('distMag'),#parameters to be evaluated, can be  'distT','distRec', 'missVal' 'distMag', 'seasAmp', 'remSd'
    nTS = 1000,
    nobsYr = STnobsYr,
    seasAv = seasVImean,
    parSetUp = 'int')#parameter set-up: can be avg dist, comb, or int

  # settings of the recovery indicators
  funSet <- list('freq' = 'dense',#rep('dense',3),
                 'input' = 'smoothed',#c('raw', 'smoothed','segmented'),# settings for the recovery indicators
                 'nPre' = 2,# rep(2,3),
                 'nDist' = 1,#rep(1,3),
                 'nPost' = 2,
                 'nPostStart' = 4,
                 'nDelta' = 1,
                 'nDeltaStart' = 5,
                 # 'nPostMin' = 4,#c(4,4,4),
                 # 'nPostMax' = 6,#rep(6,3),
                 'h' = 0.15,#rep(0.15,3),
                 'seas' = T)#rep(T,3))
  # evaluate recovery indicators
  perf <- evalParam('distMag', sttngs, funSet, basename, ofolder = '')

  # pars <- setParamValues(sttngs)
  # winsize <- c(365, 4, 1)
  # names(winsize) <- c('dense', 'quarterly', 'annual')
  # m_RRIi <- rep(NA,1000)
  # m_R80pi <- rep(NA,1000)
  # m_YrYri <- rep(NA,1000)
  # s_RRIi <- rep(NA,1000)
  # s_R80pi <- rep(NA,1000)
  # s_YrYri <- rep(NA,1000)
  # frq <- funSet[['freq']][1]
  # for (pari in 1: 1000){
  #   sc <- simulCase(pars[[1]][[1]]$nrep[pari], pars[[1]][[1]]$nyr[pari], pars[[1]][[1]]$nobsYr[pari], pars[[1]][[1]]$nDr[pari], pars[[1]][[1]]$seasAv[[1]], pars[[1]][[1]]$seasAmp[pari],pars[[1]][[1]]$trAv[pari], pars[[1]][[1]]$remSd[pari], c(pars[[1]][[1]]$distMag[pari],pars[[1]][[1]]$distMag[pari]), pars[[1]][[1]]$distT[pari], c(pars[[1]][[1]]$distRec[pari],pars[[1]][[1]]$distRec[pari]), pars[[1]][[1]]$missVal[pari], pars[[1]][[1]]$DistMissVal[pari], pars[[1]][[1]]$distType[pari])
  #   tsi <- sc[[1]][1,] # simulated time series = seasonality + trend + noise component (contains missing values)
  #   tsseas <-sc[[2]][1,] # simulated seasonality component of time series  (contains no missing values)
  #   obspyr <- sc[[5]][1,]$obs_per_year # number of observations per year
  #   tdist <- sc[[5]][1,]$dist_time# timing of the disturbance (observation number)
  #   tsref <- sc[[3]][1,]# the simulated trend componend = reference time series to measure the recovery indicators for validation (contains no missing values)
  #   nobs <- (sc[[5]][1,]$number_yrs)* obspyr# total number of observations = number of years * number of observations per year
  #
  #   # smoothing
  #   temp.zoo<-zoo(tsi,(1:length(tsi)))# generate a zoo time series object
  #   m.av<-rollapply(temp.zoo, as.numeric(winsize[frq]), mean, na.rm = T, fill = NA)# smooth time series using rolling mean
  #   tsi <- as.numeric(m.av)
  #
  #   outp <- calcFrazier(tsi, tdist, obspyr, funSet[['shortDenseTS']], funSet[['nPre']], funSet[['nDist']], funSet[['nPostMin']], funSet[['nPostMax']])
  #   m_RRIi[pari] <- outp$RRI# measured RRI
  #   m_R80pi[pari] <- outp$R80P# measured R80p
  #   m_YrYri[pari] <- outp$YrYr# measured YrYR
  #   outp <- calcFrazier(tsref, tdist, obspyr, funSet[['shortDenseTS']], funSet[['nPre']], funSet[['nDist']], funSet[['nPostMin']], funSet[['nPostMax']])
  #   s_RRIi[pari] <- outp$RRI#simulated (true) RRI
  #   s_R80pi[pari] <- outp$R80P#simulated (true) R80p
  #   s_YrYri[pari] <- outp$YrYr
  # }

  # check if the results agree with expectations
  expect_equal(length(perf), 12, tolerance = 1e-4)
  expect_equal(perf$RRI_nTS$`-0.5`, 1, tolerance = 1e-2)
  expect_equal(perf$R80p_nTS$`-0.5`, 1, tolerance = 1e-2)
  expect_equal(perf$YrYr_nTS$`-0.5`, 1, tolerance = 1e-2)

  expect_equal(perf$RRI_rsq$`-0.5`, 0.973, tolerance = 1e-2)
  expect_equal(perf$R80p_rsq$`-0.5`, 0.975, tolerance = 1e-2)
  expect_equal(perf$YrYr_rsq$`-0.5`, 0.989, tolerance = 1e-2)
  expect_equal(perf$RRI_rmse$`-0.5`, 0.073, tolerance = 1e-2)
  expect_equal(perf$R80p_rmse$`-0.5`, 0.055, tolerance = 1e-2)
  expect_equal(perf$YrYr_rmse$`-0.5`, 0.000, tolerance = 1e-2)
  expect_equal(perf$RRI_mape$`-0.5`, 0.066, tolerance = 1e-2)
  expect_equal(perf$R80p_mape$`-0.5`, 0.043, tolerance = 1e-2)
  expect_equal(perf$YrYr_mape$`-0.5`, 0.215, tolerance = 1e-1)
})
