# Synthetic Data - TimeChange

# Sampling Scheme: 48h every 2h
# Replicates: 2,3
# Noise Levels: 0.1, 0.2, 0.4 as percent of Amplitude
# 1000 of each waveform shape
# 1) Sin
# 2) Sin 2x Amp
# 3) Saw Tooth
# 4) Linear Trend
# 5) Damped
# 6) Contractile

## ---------- Parameters
set.seed(123)
interval <- 2
len <- 48
xVals <- seq(0, len, interval)
reps <- c(2, 3)
period <- 24
nLevel <- c(0.1, 0.2, 0.4)
phaseShift <- 0
nRuns <- 1000

## ---------- Functions

# Sin Function
sinwave <- function(xVals, amp, phaseShift, period, nL) {
  return(amp / 2 * sin((xVals + phaseShift) * (2 * pi / period)) + nL * amp / 2 * rnorm(n = length(xVals), mean = 0, sd = 1))
}
# Sawtooth Function
saw <- function(xVals, amp, period, phaseShift, asym, nL) {
  up <- seq(from = 0, to = amp, length.out = round(asym * period))
  down <- rev(seq(from = 0, to = amp, length.out = (period - round(asym * period) + 2)))
  output <- c(rep(c(up, down[-c(1, length(down))]), length(xVals) / period), down[length(down)])
  phaseShift <- phaseShift %% 24

  if (phaseShift == 1) {
    output <- c(output[2:length(output)], output[2])
  } else if (phaseShift > 0) {
    output <- c(output[-c(1:phaseShift - 1)], output[2:length(1:phaseShift + 1)])
  }

  # add noise
  output <- output + nL * amp / 2 * rnorm(n = length(xVals), mean = 0, sd = 1)
  return(output)
}

# LinearTrend Function
linearTrendSinwave <- function(xVals, amp, period, phaseShift, lTrend, nL) {
  return(amp / 2 * sin((xVals + phaseShift) * (2 * pi / period)) + (lTrend * (xVals + phaseShift)) + nL * amp / 2 * rnorm(n = length(xVals), mean = 0, sd = 1))
}

# Damped Function
dampedSinwave <- function(xVals, amp, period, phaseShift, d, nL) {
  return(amp / 2 * sin((xVals + phaseShift) * (2 * pi / period)) * exp(-d * xVals) + nL * amp / 2 * rnorm(n = length(xVals), mean = 0, sd = 1))
}

# Contractile Function
contract <- function(xVals, amp, period, phaseShift, k, nL) {
  return(amp / 2 * sin((xVals + phaseShift) * (2 * pi / period) * ((xVals)^k)) + nL * amp / 2 * rnorm(n = length(xVals), mean = 0, sd = 1))
}

## -----------  Condition Data Set (1000 of each type of waveform)
# 1) Sin - Control
# 2) Sin 2x Amp
# 3) Saw Tooth
# 4) Linear Trend
# 5) Damped
# 6) Contractile

# random parameters
asym <- runif(n = nRuns, min = 0.1, max = 1)
lTrend <- sample(c(runif(n = nRuns, min = 0.01, max = 0.1), runif(n = nRuns, min = -0.1, max = -0.01)), size = nRuns)
d <- runif(n = nRuns, min = 0.01, max = 0.1)
k <- runif(n = nRuns, min = 0.1, max = 0.2)

getRandomizedWaveforms <- function(noiseLevel = noiseLevel, replicates = replicates) {
  # 1) Sin - control
  sinWaves <- t(replicate(n = nRuns, expr = {
    c(t(replicate(n = replicates, expr = {
      sinwave(xVals = xVals, amp = 1, phaseShift = 0, period = period, nL = noiseLevel)
    })))
  }))
  sinWaves <- as.data.frame(sinWaves)


  # 2) Sin 2x Amp
  sinWaves2Amp <- t(replicate(n = nRuns, expr = {
    c(t(replicate(n = replicates, expr = {
      sinwave(xVals = xVals, amp = 2, phaseShift = 0, period = period, nL = noiseLevel)
    })))
  }))
  sinWaves2Amp <- as.data.frame(sinWaves2Amp)



  # 3) Saw Tooth
  sawtooth <- t(sapply(1:nRuns, FUN = function(i) {
    c(t(replicate(n = replicates, expr = {
      saw(xVals = xVals, amp = 1, period = period / 2, phaseShift = 0, asym = asym[i], nL = noiseLevel)
    })))
  }))
  sawtooth <- as.data.frame(sawtooth)



  # 4) Linear Trend
  linearTrend <- t(sapply(1:nRuns, FUN = function(i) {
    c(t(replicate(n = replicates, expr = {
      linearTrendSinwave(xVals = xVals, amp = 1, period = period, phaseShift = 0, lTrend = lTrend[i], nL = noiseLevel)
    })))
  }))
  linearTrend <- as.data.frame(linearTrend)



  # 5) Damped
  damped <- t(sapply(1:nRuns, FUN = function(i) {
    c(t(replicate(n = replicates, expr = {
      dampedSinwave(xVals = xVals, amp = 1, period = period, phaseShift = 0, d = d[i], nL = noiseLevel)
    })))
  }))
  damped <- as.data.frame(damped)



  # 6) Contractile
  contractile <- t(sapply(1:nRuns, FUN = function(i) {
    c(t(replicate(n = replicates, expr = {
      contract(xVals = xVals, amp = 1, period = period, phaseShift = 0, k = k[i], nL = noiseLevel)
    })))
  }))
  contractile <- as.data.frame(contractile)


  # bind data into single dataframe
  condition <- rbind(sinWaves, sinWaves2Amp, sawtooth, linearTrend, damped, contractile)
}

## ----------- Format data for export

gridOfparams <- expand.grid(nLevel, reps)

apply(gridOfparams, 1, function(params){
  control <- getRandomizedWaveforms(noiseLevel = params[1], replicates = params[2])
  condition <- getRandomizedWaveforms(noiseLevel = params[1], replicates = params[2])

  colNamesZTreps <- c(t(sapply(1:params[2], function(rep) {
    paste0("ZT_", xVals, "_Rep_", rep)
  })))

  colnames(control) <- colNamesZTreps
  colnames(condition) <- colNamesZTreps

  rownames(control) <- c(
    paste0("sin_", 1:nRuns),
    paste0("sin2X_", 1:nRuns),
    paste0("saw_", 1:nRuns),
    paste0("LT_", 1:nRuns),
    paste0("Damped_", 1:nRuns),
    paste0("Contract_", 1:nRuns)
  )
  rownames(condition) <- c(
    paste0("sin_", 1:nRuns),
    paste0("sin2X_", 1:nRuns),
    paste0("saw_", 1:nRuns),
    paste0("LT_", 1:nRuns),
    paste0("Damped_", 1:nRuns),
    paste0("Contract_", 1:nRuns)
  )

  syntheticData <- list(control, condition)


  fileName <- paste("syntheticData", interval, len, "rep", params[2], "nl", params[1], sep = "_")

  save(syntheticData, file = paste0("~/Desktop/TimeChange-data/Data/syntheticData_general/", fileName, ".RData"))
})
