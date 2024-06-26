% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{occ_blocksplit}
\alias{occ_blocksplit}
\title{Split occurrences in training and testing data using blocks}
\usage{
occ_blocksplit(
  data,
  longitude,
  latitude,
  train_proportion = 0.75,
  raster_layer = NULL,
  background_n = 10000
)
}
\arguments{
\item{data}{matrix or data.frame with the occurrences to be split. Columns
may vary but species, longitude, and latitude are recommended.}

\item{longitude}{(character) name of the column with longitude data.}

\item{latitude}{(character) name of the column with latitude data.}

\item{train_proportion}{(numeric) proportion of data to be used as training
occurrences. Available options are 0.25, 0.5, and 0.75. The remaining data
will be used for testing. Default = 0.75.}

\item{raster_layer}{optional SpatRaster to prepare background data.}

\item{background_n}{(numeric) optional number of coordinates to be extracted
using the \code{raster_layer}. Default = 10000.}
}
\value{
List with all occurrences (all), training occurrences (train), testing
(test) occurrences, and all occurrences with assigned blocks (block).

If a raster layer is given in \code{raster_layer}, background coordinates
will be returned as part of this list. Data will be named as bg_all, bg_train,
bg_test, and bg_block, for all, training, testing, and all background with
assigned blocks, respectively.
}
\description{
occ_blocksplit splits a set of occurrences to obtain training
and testing data using blocks.
}
