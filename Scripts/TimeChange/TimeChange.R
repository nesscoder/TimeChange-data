#' Main TimeChange Function
#'
#' Main function of the \pkg{TimeChange} package used for detecting differential rhythmic signals in time-series gene expression sets.
#' For additional help with parameter selection, see TimeChange's vignette:
#' \code{vignette("TimeChange")}.
#'
#' @param sample1 a \code{data.frame} of numeric gene expression over time (row = genes \emph{x} col = ZT times).
#' @param sample2 a \code{data.frame} of numeric gene expression over time (row = genes \emph{x} col = ZT times).
#' @param repLabelS1 a \code{vector} defining the number of replicates at each time-point in sample 1.
#' @param repLabelS2 a \code{vector} defining the number of replicates at each time-point in sample 2. If not specified, \code{repLabelS2} is assumed to have the same sampling scheme as \code{repLabelS1}.
#' @param minLag a \code{numeric} specifying the min lag to check in the 2-D embedding. Default is \code{2}.
#' @param maxLag a \code{numeric} specifying the max lag to check in the 2-D embedding. Default is  \code{5}.
#' @param dimension a \code{numeric} specifying the dimensions in which to embed the time-series. Default is \code{2}. TimeChange has not been adequately tested with values other than 2.
#' @param nPermute a \code{numeric} defining the number of permuted samples for the Fasano-Franceschini test to generate when computing the empirical p-value (note this procedure is slow and computationally expensive on the order of \code{nPermute}*O(n^2). Default is set to 0.
#' If nPermute is 0, the Fasano Franceschini distributional approximation defined by is used for defining the p-value. See Fasano and Franceschini test (1987) for details.
#' @param plotEmbedding a \code{boolean} variable indicating where to visualize the embedded data. Default is  \code{FALSE}. This should not be used when analyzing many samples, and is mainly meant for low-throughput data and debugging code.
#' @param sample1Name a \code{character} defining a name for the Sample 1 \code{data.frame}. Only used in embedding plots when \code{plotEmbedding} is \code{TRUE}.
#' @param sample2Name a \code{character} defining a name for the Sample 2 \code{data.frame}. Only used in embedding plots when \code{plotEmbedding} is \code{TRUE}
#'
#' @references{
#' \itemize{
#'    \item{TimeChange's Manuscript is available from BioRxiv at \url{INSERT LINK TO MS}}
#'    }
#'    \subsection{Takens' Theorem and Embedding Parameter Selection}{
#'    \itemize{
#'    \item{Floris Takens. Detecting strange attractors in turbulence. pages 366–381. (1981) \url{http://www.crcv.ucf.edu/gauss/info/Takens.pdf}.}
#'    \item{Joshua Garland, Ryan G. James, and Elizabeth Bradley. Leveraging information storage to select forecast-optimal parameters for delay-coordinate reconstructions. Physical Review E, 93(2). (2016) \doi{10.1103/PhysRevE.93.022221}.}
#'    \item{Varad Deshmukh, Elizabeth Bradley, Joshua Garland, and James D. Meiss. Using curvature to select the time lag for delay reconstruction. Chaos 30, 063143. (2020) \doi{10.1063/5.0005890}. }
#'    }
#'    }
#'    \subsection{Fasano-Franceschini Test}{
#'    \itemize{
#'    \item{Fasano-Franceschini Test R implementation Manuscript is available from BioRxiv at \url{https://arxiv.org/abs/2106.10539}}
#'    \item{Fasano, G., Franceschini, A. A multidimensional version of the Kolmogorov-Smirnov test. Monthly Notices of the Royal Astronomical Society 225:155-170. (1987) \doi{10.1093/mnras/225.1.155}.}
#'    \item{Peacock J.A. Two-dimensional goodness-of-fit testing in astronomy. Monthly Notices of the Royal Astronomical Society 202:615-627. (1983) \doi{10.1093/mnras/202.3.615}.}
#'    \item{Press, W. H., Teukolsky, S. A., Vetterling, W. T., Flannery, B. P. (2007). Numerical Recipes 3rd Edition: The Art of Scientific Computing. Cambridge University Press. ISBN: 0521880688}
#' }
#'
#'
#'    }
#' }
#'
#'
#'
#' @return a tidy \code{data.frame} of processed results for each gene:
#' \tabular{cccc}{
#'  \strong{sample1} \tab  \strong{sample2} \tab \strong{D.stat} \tab \strong{p.value}\cr
#'  Sample 1 gene ID \tab  Sample 2 gene ID  \tab Fasano-Franceschini Test D-stat \tab \emph{p}-value \cr
#' }
#'
#' @examples
#' # > library(TimeChange)
#'
#' #define seed for reproducible example
#' # > set.seed(2020)
#'
#' #define time-points and replicates per time-point
#' # > xVal <- seq(0,48,2)
#' # > reps <- 2
#'
#' #generate replicate sine waves with 10% noise
#' # > sine <- replicate(n = reps, expr = {
#' # >   sin(xVal*(2*pi/24)) + rnorm(length(xVal), sd = .1)
#' # > })
#' # > sine <- c(cbind(sine[,1],sine[,2]))
#'
#' #generate replicate 2x sine waves with 10% noise
#' # > sine2X <- replicate(n = reps, expr = {
#' # >    2*sin(xVal*(2*pi/24)) + rnorm(length(xVal), sd = .1)
#' # > })
#' # > sine2X <- c(cbind(sine2X[,1],sine2X[,2]))
#'
#' # pairwise compare
#' # sin vs sin
#' # sin 2x vs sin 2x
#' # sin vs sin 2x
#' # >   sample1 <- t(data.frame(
#' # >      "gene_sine" = sine,
#' # >      "gene_2Xsine" = sine2X,
#' # >      "gene_sine" = sine
#' # >   ))
#' # >  sample2 <- t(data.frame(
#' # >      "gene_sine" = sine,
#' # >      "gene_2Xsine" =  sine2X,
#' # >      "gene_2Xsine" = sine2X
#' # >   ))
#'
#' # run TimeChange on the data
#' # > TimeChange(sample1, sample2, repLabelS1 = rep(2, length(xVal)))
#'
#' @export
#'
TimeChange <- function(sample1, sample2, repLabelS1, repLabelS2 = NULL, sample1Name = "Sample 1", sample2Name = "Sample 2", minLag = 2, maxLag = 5, dimension = 2, nPermute = 0, plotEmbedding = F) {
  if (nrow(sample1) != nrow(sample2)) {
    stop("There are an uneven number of samples, make sure they the number of rows in both samples data frames are equal")
  }

  sample1_ID <- rownames(sample1)
  sample2_ID <- rownames(sample2)

  if (is.null(repLabelS2)) {
    repLabelS2 <- repLabelS1
  }

  # get Negative Binomial fit to data to approximate variance
  sample1_nbFit <- getNBfit(sample1, repLabelS1)
  sample2_nbFit <- getNBfit(sample2, repLabelS2)

  output <- sapply(1:nrow(sample1), function(i) {
    S1 <- unlist(unname(sample1[i, ]))
    S2 <- unlist(unname(sample2[i, ]))

    contExpression <- S1 - mean(S1)
    mutExpression <- S2 - mean(S2)

    # Use Replicates to resample a distribution of trajectories for each gene
    control <- getSampledTimeSeries(timeSeries = contExpression, repLabel = repLabelS1, nbFit = sample1_nbFit)
    condition <- getSampledTimeSeries(timeSeries = mutExpression, repLabel = repLabelS2, nbFit = sample2_nbFit)

    # select embedding with min curvature for each gene
    controlCurvatures <- getMedianCurvatureFromReplicateResamplings(TimeSeriesDF = control, minLag = minLag, maxLag = maxLag, dimension = dimension)
    conditionCurvatures <- getMedianCurvatureFromReplicateResamplings(TimeSeriesDF = condition, minLag = minLag, maxLag = maxLag, dimension = dimension)
    sel <- which.min(c(min(controlCurvatures, na.rm = T), min(conditionCurvatures, na.rm = T)))
    tau <- c(minLag:maxLag)[which.min(data.frame(controlCurvatures, conditionCurvatures)[, sel])]

    pDataPlot <- rbind(
      getEmbedding(data = control, label = sample1Name, tau = tau),
      getEmbedding(data = condition, label = sample2Name, tau = tau)
    )

    if (plotEmbedding) {
      plot(
        x = pDataPlot[, 1],
        y = pDataPlot[, 2],
        ylim = range(c(pDataPlot[, 1], pDataPlot[, 2])),
        pch = 21,
        cex = 2,
        bg = c("#537FAE", "#F09838")[as.factor(pDataPlot[, 3])],
        main = paste0(sample1Name, ": ", sample1_ID[i],"\n versus \n", sample2Name, ": ", sample2_ID[i]),
        xlab = "f(t)",
        ylab = paste0("f(t+", tau, ")"),
        asp = 1
      )
      legend("topleft",
        legend = c(sample1Name, sample2Name),
        fill = c("#537FAE", "#F09838"), box.lwd = 0, cex = 0.8, inset = 0.01
      )
    }

    # 2-D KS Test
    contData <- pDataPlot %>%
      filter(condition == sample1Name)

    mutData <- pDataPlot %>%
      filter(condition == sample2Name)

    results <- fasano.franceschini.test::fasano.franceschini.test(S1 = contData[, 1:2], S2 = mutData[, 1:2], nPermute = nPermute)
    return(c(results$statistic, results$p.value))
  })

  timeChangeOutput <- cbind(sample1_ID, sample2_ID, t(output[1:2, ]))
  return(data.frame(timeChangeOutput))
}
