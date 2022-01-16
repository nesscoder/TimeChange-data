#' Resampled Time-Series
#'
#' Generates a \code{data.frame} of resampled time-series conditioned on transcription and degradation rates.
#' Each row represents a new resampling, with the number of resampled time-series defined by the \code{resamplings} parameter.
#' See \code{\link{conditionedTPsim},\link{replicateTPsim}} for more details.
#'
#' @param timeSeries a \code{vector} of numeric gene expression values.
#' @param repLabel a \code{vector} defining the number of replicates at each time-point.
#' @param nbFit a fitted model \code{object} of class \code{negbin} inheriting from \code{glm} and \code{lm}. See \code{MASS} package for details.
#' @param resamplings a \code{numeric} specifying the number of resampled time-series to produce. The Default is \code{5}.
#'
#' @seealso \code{\link{conditionedTPsim},\link{replicateTPsim}}.
#'
#' @return a \code{data.frame} of resampled numeric gene expression time-series conditioned on transcription and degradation rates.
#'
getSampledTimeSeries <- function(timeSeries, repLabel, nbFit, resamplings = 5) {
  timePointsSplitByRepsPerTimePoint <- split(timeSeries, rep(1:length(repLabel), repLabel))
  out <- t(replicate(n = resamplings, conditionedTPsim(timePointsSplitByRepsPerTimePoint, nbFit)))
  return(out)
}





#' Biologically Constrained Resampling
#'
#' Creates a vector of resampled time-series via Monte Carlo simulation.
#' Based on the median expression of a given time-point, the accompany standard deviation is computed using the negative binomial fit.
#' Simulated Gene expression is then drawn from a normal distribution with median (\eqn{\mu}) and standard deviation (\eqn{\sigma}) for each time-point.
#' Resampled time-points are conditioned on transcription and degradation rates (e.g. changes in expression cannot be larger or smaller than the difference in expression between time-points in the raw data)
#' If a gene is outside the conditioned transcription and degradation rates, a new point is simulated until it fits inside the bounds.
#'
#' @param timePointsSplitByRepsPerTimePoint a \code{list} of numeric gene expression values split by time-point.
#' @param nbFit a fitted model \code{object} of class \code{negbin} inheriting from \code{glm} and \code{lm}. See \code{MASS} package for details.
#'
#' @return a \code{vector} of resampled numeric gene expression values conditioned on transcription and degradation rates.
#'
conditionedTPsim <- function(timePointsSplitByRepsPerTimePoint, nbFit) {
  tsDiff <- diff(unname(unlist(lapply(timePointsSplitByRepsPerTimePoint, median))))

  oldTP <- NA
  simulatedTimeSeries <- lapply(timePointsSplitByRepsPerTimePoint, function(tp) {
    if (is.na(oldTP)) {
      oldTP <- replicateTPsim(replicateTPs = tp, nbFit)
      tempTP <- oldTP
    } else {
      repeat{
        newTP <- replicateTPsim(replicateTPs = tp, nbFit)
        tpChange <- diff(c(tempTP, newTP))
        # transcription and degradation rates constraint
        if (tpChange > min(tsDiffs) & tpChange < max(tsDiffs)) {
          tempTP <- newTP
          return(newTP)
        }
      }
    }
  })

  return(unname(unlist(simulatedTimeSeries)))
}





#' Simulate Replicate Time-Point
#'
#' Helper functions used by \code{conditionedTPsim()}, for simulating replicate time-points from the negative binomial fit.
#'
#' @param replicateTPs a \code{vector} defining the expression of replicate at a single time-point.
#' @param nbFit a fitted model \code{object} of class \code{negbin} inheriting from \code{glm} and \code{lm}. See \code{MASS} package for details.
#'
#' @seealso \code{\link{conditionedTPsim}}
#'
#' @return a \code{numeric} defining a single simulated time-serie based on the negative binomial fit.
#'
replicateTPsim <- function(replicateTPs, nbFit) {
  meanRep <- mean(replicateTPs)
  sdRep <- predict(object = nbFit, newdata = data.frame(mean = meanRep), type = "response")
  randomTP <- rnorm(n = 1, mean = meanRep, sd = sdRep)
  return(randomTP)
}





#' Median Curvature Across Replicates
#'
#' Computes the median curvature across embedded time-series replicates. Curvature values are computed on a per lag basis.
#'
#' @param TimeSeriesDF a \code{data.frame} of numeric gene expression over time (row = genes \emph{x} col = ZT times).
#' @param minLag a \code{numeric} specifying the min lag to check in the 2-D embedding. Default is \code{2}.
#' @param maxLag a \code{numeric} specifying the max lag to check in the 2-D embedding. Default is  \code{5}.
#' @param dimension a \code{numeric} specifying the dimension to use for in the time-delayed embedding (i.e. 2-D or 3-D). Default is \code{2}.
#'
#' @seealso \code{\link{buildTakens_ndim}, \link{getCurvature}}.
#'
#' @return a \code{data.frame} of median curvature across all replicates on a per lag basis.
#'
getMedianCurvatureFromReplicateResamplings <- function(TimeSeriesDF, minLag = 2, maxLag = 5, dimension = 2) {
  meanCurvaturePerLag <- apply(TimeSeriesDF, 1, function(ts) {
    m <- sapply(minLag:maxLag, function(lag) {
      embedding <- buildTakens_ndim(ts, dim = dimension, delay = lag)
      curviture <- sapply(1:(dim(embedding)[1] - 2), function(i) {
        getCurvature(embedding[i:(i + 2), ])
      })
      mean(curviture)
    })
  })

  # take the median over all replicates to compute the average curvature on a per lags basis
  medianCurvatureOverReplicates <- apply(meanCurvaturePerLag, 1, median)
  return(medianCurvatureOverReplicates)
}





#' Embed Time-Series
#'
#' Takes in a list of resampled time-series and Embeds them as defined by the embedding \code{dimension} and lag constant \code{tau}.
#'
#' @param data a \code{data.frame} of numeric gene expression over time (row = genes \emph{x} col = ZT times).
#' @param label a \code{character} defining a name for the output embedding.
#' @param tau a \code{numeric} specifying the lag to use for in the n-dimensional time delayed embedding specified by \code{dim}.
#' @param dimension a \code{numeric} specifying the dimension to use for in the time-delayed embedding (i.e. 2-D or 3-D).
#'
#' @seealso \code{\link{buildTakens_ndim}}.
#'
#' @return a \code{list} of embedded time-series.
#'
getEmbedding <- function(data, label, tau, dimension = 2) {
  emData <- lapply(1:nrow(data), function(i) {
    embedding <- buildTakens_ndim(data[i, ],
      dim = dimension,
      delay = tau
    )

    embedding <- data.frame(embedding, condition = label)
    colnames(embedding) <- c("x", "y", "condition")
    return(embedding)
  })

  emDataOutput <- do.call(what = rbind, args = emData)

  return(emDataOutput)
}





#' Generates Takens' Embedding For a Time-Series
#'
#' Generates the Takens' embedding for a time-series given a specified delay (i.e. lag) and dimension (i.e. 2-D or 3-D).
#'
#' @param x a \code{vector} of \code{numeric} time-series expression values.
#' @param dim a \code{numeric} specifying the dimension to use for in the time-delayed embedding (i.e. 2-D or 3-D).
#' @param delay a \code{numeric} specifying the lag to use for in the n-dimensional time delayed embedding specified by \code{dim}.
#'
#' @return a \code{data.frame} of the n-dimensional Takens' embedding. Columns defined from (1-D to n-D).
#'
buildTakens_ndim <- function(x, dim = 2, delay) {
  n <- length(x) - (dim - 1) * delay
  X <- seq_along(x)
  if (n <= 0) {
    stop("Insufficient observations for the requested embedding")
  }
  out <- matrix(rep(X[seq_len(n)], dim), ncol = dim)
  out[, -1] <- out[, -1, drop = FALSE] +
    rep(seq_len(dim - 1) * delay, each = nrow(out))

  out <- matrix(x[out], ncol = dim)

  return(out)
}





#' Takens' Embedding Curvature
#'
#' Computes the Curvature along the Takens' Embedding For a Time-Series
#'
#' @param df a \code{data.frame} with 3 rows and 2 columns. The Columns define the x and y coordinates of a point. Each row defines a point the define a triangle.
#'
#' @references{
#' \itemize{
#'    \item{Varad Deshmukh, Elizabeth Bradley, Joshua Garland, and James D. Meiss , "Using curvature to select the time lag for delay reconstruction", Chaos 30, 063143 (2020) \url{https://doi.org/10.1063/5.0005890}}
#'    }
#'}
#'
#' @return a \code{numeric} defining the curvature of the circumscribed triangle.
#'
getCurvature <- function(df) {
  # calculate area of circle
  tri_sides <- as.vector(dist(df))
  # heron formula
  s <- sum(tri_sides) / 2
  Area <- sqrt(s * (s - tri_sides[1]) * (s - tri_sides[2]) * (s - tri_sides[3]))
  curvature <- 4 * Area / prod(tri_sides)
  return(curvature)
}





#' Negative Binomial Fit
#'
#' Gets the Negative Binomial Fit of the Data
#'
#' @param data a \code{data.frame} of numeric gene expression over time (row = genes \emph{x} col = ZT times).
#' @param repLabel a \code{vector} defining the number of replicates at each time-point.
#'
#' @seealso \code{\link{getMeanVsSd}}.
#'
#' @return a fitted model \code{object} of class \code{negbin} inheriting from \code{glm} and \code{lm}. See \code{MASS} package for details.
#'
getNBfit <- function(data, repLabel) {
  sd_all <- apply(data, 1, function(i) {
    getMeanVsSd(timeSeries = i, repLabel = repLabel)
  })
  sd_all <- do.call(rbind, sd_all)
  nbfit <- MASS::glm.nb(sd ~ mean, data = sd_all)
  return(nbfit)
}




#' Time-Series Mean and Standard Deviation
#'
#' Computes the mean and variance at each time-point in a time-series.
#'
#' @param timeSeries a \code{vector} of numeric gene expression values.
#' @param repLabel a \code{vector} defining the number of replicates at each time-point.
#'
#' @return a \code{list} of all means and standard deviation of expression at each time-point.
#'
getMeanVsSd <- function(timeSeries, repLabel) {
  timePointsSplitByRepsPerTimePoint <- split(timeSeries, rep(1:length(repLabel), repLabel))
  out <- lapply(timePointsSplitByRepsPerTimePoint, function(tp) {
    data.frame(
      mean = mean(tp),
      sd = sd(tp)
    )
  })
  return(do.call(rbind, out))
}
