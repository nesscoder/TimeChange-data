# Synthetic Data - TimeChange

# Sampling Scheme: 48h every 2h
# Replicates: 2,3
# Noise Levels: 0.1, 0.2, 0.4 as percent of Amplitude
# 500 of each waveform shape
# 1) Sin - phaseShifted every 2 hours

## ---------- Parameters
set.seed(123)
interval <- 2
len <- 48
xVals <- seq(0, len, interval)
reps <- c(2, 3)
period <- 24
nLevel <- c(0.1, 0.2, 0.4)
phaseShift <- 0:23
nRuns <- 500

## ---------- Functions

# Sin Function
sinwave <- function(xVals, amp, phaseShift, period, nL) {
  return(amp / 2 * sin((xVals + phaseShift) * (2 * pi / period)) + nL * amp / 2 * rnorm(n = length(xVals), mean = 0, sd = 1))
}

getRandomizedWaveforms <- function(phaseShift = phaseShift, noiseLevel = noiseLevel, replicates = replicates) {
  # 1) Sin - control
  sinWaves <- t(replicate(n = nRuns, expr = {
    c(t(replicate(n = replicates, expr = {
      sinwave(xVals = xVals, amp = 1, phaseShift = phaseShift, period = period, nL = noiseLevel)
    })))
  }))
  sinWaves <- as.data.frame(sinWaves)
}


## ----------- Format data for export

gridOfparams <- expand.grid(nLevel, reps)

apply(gridOfparams, 1, function(params) {
  control <- lapply(phaseShift, function(phase) {
    getRandomizedWaveforms(phaseShift = phase, noiseLevel = params[1], replicates = params[2])
  })
  control <- do.call(what = rbind, args = control)

  condition <- lapply(phaseShift, function(phase) {
    getRandomizedWaveforms(phaseShift = phase, noiseLevel = params[1], replicates = params[2])
  })
  condition <- do.call(what = rbind, args = condition)

  colNamesZTreps <- c(t(sapply(1:params[2], function(rep) {
    paste0("ZT_", xVals, "_Rep_", rep)
  })))

  colnames(control) <- colNamesZTreps
  colnames(condition) <- colNamesZTreps

  rownames(control) <- as.vector(sapply(paste0("phase", phaseShift), function(i) {
    paste0(i, "_", 1:nRuns)
  }))

  rownames(condition) <- as.vector(sapply(paste0("phase", phaseShift), function(i) {
    paste0(i, "_", 1:nRuns)
  }))

  syntheticData <- list(control, condition)

  fileName <- paste("syntheticData_dphase", interval, len, "rep", params[2], "nl", params[1], sep = "_")

  save(syntheticData, file = paste0("~/Desktop/TimeChange-data/Data/syntheticData_dphase/", fileName, ".RData"))
})
