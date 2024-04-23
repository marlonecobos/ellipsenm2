#' Helper function to sub-sample occurrence data
#'
#' @param data data.frame of occurrence records. Columns must be longitude and
#' latitude.
#' @param replicates (numeric) number of subsamples to be produced. If
#' \code{replicate_type} = "jackknife" the maximum number of subsamples produced
#' equal the number of records in \code{data}. Default = 1 returns the complete
#' set of data as the only element of a list.
#' @param replicate_type (character) name of the sampling method to be used.
#' Options are "bootstrap", "jackknife", or "subsample". Default = "subsample".
#' @param percentage (numeric) percentage of data to be bootstrapped
#' for each replicate. Default = 50. Valid if \code{replicates} > 1 and
#' \code{replicate_type} = "bootstrap".
#'
#' @usage
#' data_subsample(data, replicates = 1, replicate_type = "subsample",
#'                percentage = 50)
#'
#' @export
#'
#' @examples
#' # reading data
#' occurrences <- read.csv(system.file("extdata", "occurrences.csv",
#'                                     package = "ellipsenm"))
#'
#' # subsampling by bootstrap
#' subsamples <- data_subsample(data = occurrences[, -1], replicates = 5,
#'                              replicate_type = "bootstrap", percentage = 70)
#'
#' lapply(subsamples, head)
#' lapply(subsamples, dim)
#'
#' # subsampling by "jackknife" one record out at the time
#' subsamples1 <- data_subsample(data = occurrences[, -1], replicates = 5,
#'                               replicate_type = "jackknife")
#'
#' lapply(subsamples1, head)
#' lapply(subsamples1, dim)

data_subsample <- function(data, replicates = 1, replicate_type = "subsample",
                           percentage = 75) {
  # -----------
  # detecting potential errors
  if (missing(data)) {
    stop("Argument 'data' is necessary to perform the analysis")
  }

  # -----------
  # running
  if (replicates == 1) {
    samples <- list(data)

  } else {
    ndata <- nrow(data)

    if (replicate_type %in% c("bootstrap", "jackknife", "subsample")) {

      if (replicate_type == "bootstrap") {
        part <- round(ndata * percentage / 100)
        samples <- lapply(1:replicates, function(x) {
          set.seed(x)
          data[sample(ndata, part, replace = TRUE), ]
        })
      }

      if (replicate_type == "jackknife") {
        if (replicates > ndata) {
          warning(paste0("Maximum number of replicates under \"jackknife\" replicate_type",
                         "\ncannot be > ", ndata, ", only ", ndata, " replicates will be produced."))
          replicates <- ndata
        }

        samples <- lapply(1:replicates, function(x) {
          data[-x, ]
        })
      }

      if (replicate_type == "subsample") {
        part <- round(ndata * percentage / 100)
        samples <- lapply(1:replicates, function(x) {
          set.seed(x)
          data[sample(ndata, part, replace = FALSE), ]
        })
      }

    } else {
      stop("'replicate_type' must be: 'bootstrap', 'jackknife', or 'subsample'.")
    }
  }

  return(samples)
}
