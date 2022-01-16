#get arguments from command line
args <- commandArgs(trailingOnly = TRUE)
print(args)

library('annotate')
library('tidyverse')
library('data.table')
library('foreach')
library('ggplot2')
library('knitr')
library('limma')
library('limorhyde')
library('readr')
library('stringr')
library('pROC')


source(system.file('extdata', 'vignette_functions.R', package = 'limorhyde'))

# load Data
file <- args[1]
load(args[2])

contExpression <- syntheticData[[1]]
mutExpression <- syntheticData[[2]]

reps <- parse_number(str_extract(string = file, pattern = "rep_._"))
phaseShift <- paste0('phase', 0:23)

period <- 24
startZT <- 0
len <- 48
interval <- 2

ZTtp <- seq(startZT,len, interval)
baseTPs <- length(ZTtp)
repLabel <- rep(reps, ncol(contExpression) / reps)
idxNum <- table(sort(ZTtp%%period))
tpReps <- split(colnames(contExpression), rep(1:length(repLabel), repLabel))
df <- data.frame(idx = 1:baseTPs, rep = (ZTtp%%period/interval+1)) %>%
  arrange(rep)
tpIdx <- as.vector(sapply(df$idx, FUN = function(i){
  tpReps[[i]]
}))
repLabel <- as.vector(unlist(lapply(split(repLabel, rep(names(idxNum),unname(idxNum))),sum)))
contExpression <- contExpression[,tpIdx]
mutExpression <- mutExpression[,tpIdx]

# use each waveform in turn as control
output <- lapply(X = 1:length(phaseShift), FUN = function(phase) {
  contIdx <- (500 * (phase - 1) + 1):(500 * phase)
  pairwiseComp <- replicate(n = length(phaseShift), contExpression[contIdx, ], simplify = FALSE)
  contExpression <- do.call(rbind.data.frame, pairwiseComp)
  rownames(contExpression) <- paste0(phaseShift[phase], "_", 1:nrow(contExpression))


  ## --- Run LimoRhyde
  start <- Sys.time()
  emat <- cbind(contExpression, mutExpression)
  timeZT <- parse_number(colnames(emat))
  sm <- data.frame(cond = c(rep("wt", 25*reps), rep("mut", 25*reps)), time = timeZT, limorhyde(timeZT, 'time_'))
  design <- model.matrix(~ cond * (time_cos + time_sin), data = sm)
  fit <- lmFit(emat, design)
  fit <- eBayes(fit, trend = TRUE)
  results <- data.table(topTable(fit, coef = 5:6, number = Inf,sort.by = "none"), keep.rownames = TRUE)
  end <- Sys.time()
  totalTime <- end - start

  ## --- ROC & AUC Check
  expected <- rep(1, 12000)
  expected[contIdx] <- 0
  pred <- as.numeric(as.vector(unlist(results$P.Value)))
  ROC <- roc(expected, pred)

  fileName <- paste0(phaseShift[phase], "results_", file)
  return(list(fileName = fileName,
              inputData = list(control = contExpression, mutant = mutExpression),
              results = results,
              totalTime = totalTime,
              ROC = ROC))
})

names(output) <- phaseShift
rm(list = setdiff(ls(), c("args", "output", "syntheticData")))

save.image(file = args[3])
