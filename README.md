ellipsenm: Ecological Niche Characterizations Using Ellipsoids
================
Marlon E. Cobos, Luis Osorio-Olvera, Jorge Soberón, A. Townsend Peterson

- [Package description](#package-description)
- [Installing the package](#installing-the-package)
- [Exploring ellipsenm](#exploring-ellipsenm)
  - [Setting a directory](#setting-a-directory)
  - [Ecological niches in ellipsenm](#ecological-niches-in-ellipsenm)
  - [Functions in ellipsenm](#functions-in-ellipsenm)
  - [Model calibration in ellipsenm](#model-calibration-in-ellipsenm)
  - [Modeling ecological niches using
    ellipsoids](#modeling-ecological-niches-using-ellipsoids)
  - [Niche overlap using ellipsenm](#niche-overlap-using-ellipsenm)
- [Project description](#project-description)
  - [Status of the project](#status-of-the-project)

<br>

## Package description

The **ellipsenm** R package implements multiple tools to help in using
ellipsoid envelopes to model ecological niches of species. Handly
options for calibrating and selecting models, producing models with
replicates and projections, and assessing niche overlap are included as
part of this package. Other functions implemented here are useful to
perform a series of pre- and post-modeling analyses.

<br>

## Installing the package

**ellipsenm** is in a GitHub repository and can be installed and/or
loaded using the code below (make sure to have Internet connection). One
of the functions to evaluate model performance in this package needs
compilation. That is why you must install a compilation tools before
installing the package, **Rtools** for Windows or other tools in other
Operative Systems. A guide for downloading and installing Rtools can be
found
<a href="http://jtleek.com/modules/01_DataScientistToolbox/02_10_rtools/#1" target="_blank">here</a>.
IMPORTANT note: Add Rtools to the **PATH** during its installation.

Try the code below first… If you have any problem during the
installation, restart your R session, close other sessions you may have
open, and try again. If during the installation you are asked to update
packages, please do it (select the option that says All). If any of the
packages gives an error, please install it alone using
install.packages(), then try re-installing **ellipsenm** again. Also, it
may be a good idea to update R and RStudio (if you are using it).

``` r
# Installing and loading packages
if(!require(remotes)){
    install.packages("remotes")
}
if(!require(ellipsenm)){
    remotes::install_github("marlonecobos/ellipsenm2")
}
library(ellipsenm)
```

<br>

## Exploring ellipsenm

### Setting a directory

The main functions of the **ellipsenm** package produce results that
need to be written in a directory in your computer. Writing the results
outside the R environment helps to avoid problems related to RAM
limitations. That is why, setting a working directory is recommended
before starting. You ca do that using the code below:

``` r
setwd("Your/Directory") # change the characters accordingly
```

<br>

### Ecological niches in ellipsenm

An ecological niche, from a Grinnellian perspective, is the set of
environmental conditions that allow a species to maintain populations
for long periods of time, without immigration events. Models created
with the ellipsenm package are ellipsoid envelope models and assume that
a species ecological niche is convex, has an only optimum, and the
species response to each variable covaries with the response to other
variables. Mahalanobis distances are used to represent how far is each
combination of environmental conditions from the optimum (ellipsoid
centroid). Suitability values result by default from a multivariate
normal transformation of Mahalanobis distances. Therefore, maximum
values of suitability will be close to the centroid and minimum values
will be close to the border of the ellipsoid.

<br>

### Functions in ellipsenm

A complete list of the main functions in the **ellipsenm** package can
be found in the package documentation. Use the following code to see the
list.

``` r
help(ellipsenm)
```

<br>

### Model calibration in ellipsenm

Model calibration is a process in which various candidate models,
created with distinct (coherent) parameter settings, are tested based on
metrics that reflect their performance. After that models the best
models are selected, based on user defined selection criteria, to
represent the phenomenon of interest. The function
*ellipsoid_calibration* will help to create candidate models, evaluate
them, and select the best according to the selected criteria. The
performance of models in this process is assessed based on statistical
significance (partial ROC), predictive power (omission rates), and
prevalence. Where good models are the ones statistically significant,
with low omission rates, and with low prevalence.

We encourage the users to check the function’s help before using it.
This is possible using the code below:

``` r
help(ellipsoid_calibration)
```

The code below helps users to perform a small exercise using various
functions to prepare the data and to perform the model calibration
process.

``` r
# reading data
occurrences <- read.csv(system.file("extdata", "occurrences.csv",
                                    package = "ellipsenm"))
colnames(occurrences)

# raster layers of environmental data (this ones are masked to the accessible area)
# users must prepare their layers accordingly if using other data
vars <- terra::rast(list.files(system.file("extdata", package = "ellipsenm"),
                               pattern = "bio", full.names = TRUE))

# preparing training and testing data
data_split <- split_data(occurrences, method = "random", longitude = "longitude",
                         latitude = "latitude", train_proportion = 0.75)

# sets of variables (example)
sets <- list(set_1 = c("bio_1", "bio_7", "bio_15"),
             set_2 = c("bio_1", "bio_12", "bio_15"),
             set3 = c("bio_1", "bio_7", "bio_12", "bio_15")) # change as needed

variable_sets <- prepare_sets(vars, sets)

# methods to create ellipsoids
methods <- c("covmat")

# model calibration process (Make sure to define your working directory first)
calib <- ellipsoid_calibration(data = data_split, species = "species",
                               longitude = "longitude", latitude = "latitude",
                               variables = variable_sets, methods = methods,
                               level = 99, selection_criteria = "S_OR_P",
                               error = 5, iterations = 500, percentage = 50,
                               output_directory = "calibration_results")

class(calib)
# check your directory, folder "calibration_results"
```

<br>

### Modeling ecological niches using ellipsoids

Once you have decided what are the best parameter settings for your
models either using a model calibration process or other approach, your
final models can be produced. This models will help to represent the
ecological niche of a species in environmental and geographic space. The
function *ellipsoid_model* will help to perform all necessary analyses.

Please check the function’s help with the code below:

``` r
help(ellipsoid_model)
```

Now the examples. Three distinct ways to create ecological niche models
using **ellipsenm** are presented below:

<br>

**Creating a simple ecological niche model**

When models are created this way, the whole dataset is used for fitting
the ellipsoids.

``` r
# reading data
occurrences <- read.csv(system.file("extdata", "occurrences.csv",
                                    package = "ellipsenm"))

# raster layers of environmental data
vars <- terra::rast(list.files(system.file("extdata", package = "ellipsenm"),
                               pattern = "bio", full.names = TRUE))

# creating the model with no replicates
ell_model <- ellipsoid_model(data = occurrences, species = "species",
                             longitude = "longitude", latitude = "latitude",
                             raster_layers = vars, method = "covmat", level = 99,
                             replicates = 1, prediction = "suitability",
                             return_numeric = TRUE, format = "GTiff",
                             overwrite = FALSE,
                             output_directory = "ellipsenm_model")

class(ell_model)
# check your directory, folder "ellipsenm_model"
```

<br>

**Modeling an ecological niche with replicates**

When models are replicated, sub-samples of the data are used to create
each replicate. Mean, minimum, and maximum ellipsoid models are also
calculated. See more details about sub-sampling in the function’s help.

``` r
# creating the model with replicates
ell_model1 <- ellipsoid_model(data = occurrences, species = "species",
                              longitude = "longitude", latitude = "latitude",
                              raster_layers = vars, method = "covmat", level = 99,
                              replicates = 5, prediction = "suitability",
                              return_numeric = TRUE, format = "GTiff",
                              overwrite = FALSE,
                              output_directory = "ellipsenm_model1")

class(ell_model1)
# check your directory, folder "ellipsenm_model1"
```

<br>

**Example of an ecological niche model with projections**

Ellipsoid envelope models can also be projected to other scenarios. This
is one example using a RasterStack of layers. Projections can be done to
multiple scenarios in same modeling process by changing the type of
*projection_variables*. Check the function’s help for more details.

``` r
# creating the model with projections
pr_vars <- terra::rast(system.file("extdata", "proj_variables.tif",
                                   package = "ellipsenm"))

ell_model2 <- ellipsoid_model(data = occurrences, species = "species",
                              longitude = "longitude", latitude = "latitude",
                              raster_layers = vars, method = "covmat", level = 99,
                              replicates = 3, replicate_type = "subsample",
                              percentage = 75, projection_variables = pr_vars,
                              prediction = "suitability", return_numeric = TRUE,
                              format = "GTiff", overwrite = FALSE,
                              output_directory = "ellipsenm_model2")

class(ell_model2)
# check your directory, folder "ellipsenm_model2"
```

<br>

### Niche overlap using ellipsenm

Niche overlap analyses that can be performed in **ellipsenm** allow to
compare the degree of overlap between ecological niches in environmental
space. The function *ellipsoid_overlap* helps is performing needed
analyses. To plot results use the function *plot_overlap*.

Please check the function’s help with the code below:

``` r
help(ellipsoid_overlap)
```

The code below helps to prepare data and perform an example of overlap
analysis.

``` r
# preparing data
vext <- terra::ext(vars)
ext1 <- terra::ext(vext[1], (mean(vext[1:2]) + 0.2), vext[3:4])
ext2 <- terra::ext((mean(vext[1:2]) + 0.2), vext[2], vext[3:4])

# croping variables and splitting occurrences
vars1 <- terra::crop(vars, ext1)
vars2 <- terra::crop(vars, ext2)

occurrences1 <- occurrences[occurrences$longitude < (mean(vext[1:2]) + 0.2), ]
occurrences2 <- occurrences[!occurrences$longitude %in% occurrences1$longitude, ]

# preparing overlap objects to perform analyses
niche1 <- overlap_object(occurrences1, species =  "species", longitude = "longitude",
                         latitude = "latitude", method = "covmat", level = 95,
                         variables = vars1)

niche2 <- overlap_object(occurrences2, species =  "species", longitude = "longitude",
                         latitude = "latitude", method = "covmat", level = 95,
                         variables = vars2)

# niche overlap analysis
overlap <- ellipsoid_overlap(niche1, niche2)

# niche overlap analysis with test of significance
overlap_st <- ellipsoid_overlap(niche1, niche2, overlap_type = "back_union",
                                significance_test = TRUE, replicates = 100)
```

Use the following lines of code for plotting results using distinct
options.

``` r
# plotting only ellipsoids
plot_overlap(overlap)

# plotting ellispods and background for full overlap
plot_overlap(overlap, background = TRUE, proportion = 0.6, background_type = "full")

# plotting ellispods and background for overlap based on accessible environments
plot_overlap(overlap, background = TRUE,  proportion = 1, background_type = "back_union")
```

<br>

**Detailed examples for other functions of the package can be found in
their respective documentation.**

<br>

**This repository is for the project “Grinnellian ecological niches and
ellipsoids in R” developed during the program GSoC 2019.**

<br>

## Project description

Student: *Marlon E. Cobos*

GSoC Mentors: *Luis Osorio-Olvera, Vijay Barve, and Narayani Barve*

Motivation:

Ecological niche modeling represents a set of methodological tools that
allow researchers to characterize and analyze species ecological niches.
Currently, these methods are used widely and their applications include
disease risk mapping, climate change risk predictions, conservation
biology, among others. Physiological data suggests that Grinnellian
ecological niches are convex in nature and they may probably have an
ellipsoidal form when multiple dimensions are considered. However, among
the diverse available software to characterize ecological niches,
algorithms to model these niches as ellipsoids in the environmental
space are scarce. Given the need for more solutions, the **ellipsenm**
package aimed for developing specialized tools to perform multiple
analyses related to ecological niches of species.

### Status of the project

At the moment we have completed the three phases of the project. We have
made few modifications to the list of original products that have helped
us to improve the package functionality. The package is fully functional
but improvements can be done and in the future, they will be included.
Future improvements will include: options to allow reproducibility,
production of HTML reports of other main analyses, and generic methods
for summarizing and plotting resultant objects. Future plans also
include submitting the package to CRAN.

All commits made can be seen at the
<a href="https://github.com/marlonecobos/ellipsenm/commits/master" target="_blank">complete
list of commits</a>.

Following you can find a brief description of this R package, as well as
some examples of its use (still in development).
