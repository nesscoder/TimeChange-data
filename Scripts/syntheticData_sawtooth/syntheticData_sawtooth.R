# Synthetic Data - TimeChange

# Sampling Scheme: 48h every 2h
# Replicates: 3,2
# Noise Levels: 0, 0.2, 0.4 as percent of Amplitude
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
nRuns <- 500

## ---------- Functions

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

getRandomizedWaveforms <- function(asym = asym, noiseLevel = noiseLevel, replicates = replicates) {

  # 3) Saw Tooth
  sawtooth <- t(sapply(1:nRuns, FUN = function(i) {
    c(t(replicate(n = replicates, expr = {
      saw(xVals = xVals, amp = 1, period = period / 2, phaseShift = 0, asym = asym, nL = noiseLevel)
    })))
  }))
  sawtooth <- as.data.frame(sawtooth)
}

# Half  parameters
nPeriods <- (len / interval) / (len / period)
asym <- seq(2 / nPeriods, 1, 1 / nPeriods)

## ----------- Format data for export

gridOfparams <- expand.grid(nLevel, reps)

apply(gridOfparams, 1, function(params) {
  control <- lapply(asym, function(peaks) {
    getRandomizedWaveforms(asym = peaks, noiseLevel = params[1], replicates = params[2])
  })
  control <- do.call(what = rbind, args = control)

  condition <- lapply(asym, function(peaks) {
    getRandomizedWaveforms(asym = peaks, noiseLevel = params[1], replicates = params[2])
  })
  condition <- do.call(what = rbind, args = condition)

  colNamesZTreps <- c(t(sapply(1:params[2], function(rep) {
    paste0("ZT_", xVals, "_Rep_", rep)
  })))

  colnames(control) <- colNamesZTreps
  colnames(condition) <- colNamesZTreps

  rownames(control) <- as.vector(sapply(paste0("peak", 2:nPeriods), function(i) {
    paste0(i, "_", 1:nRuns)
  }))

  rownames(condition) <- as.vector(sapply(paste0("peak", 2:nPeriods), function(i) {
    paste0(i, "_", 1:nRuns)
  }))

  syntheticData <- list(control, condition)
  fileName <- paste("syntheticData_sawtooth", interval, len, "rep", params[2], "nl", params[1], sep = "_")
  rm(list = setdiff(ls(), c("fileName", "syntheticData")))

  save(syntheticData, file = paste0("~/Desktop/TimeChange-data/Data/syntheticData_sawtooth/", fileName, ".RData"))
})
