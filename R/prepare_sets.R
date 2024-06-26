#' Prepare sets of variables for calibration
#'
#' @param variables (character or SpatRaster) if character, name of the folder
#' containing only the variables to be divided into sets; if SpatRaster, stack of
#' variables to be divided into sets.
#' @param sets named list of character vectors with the names of the variables
#'  per each subset. Default = NULL.
#' @param all_combinations (logical) whether or not to prepare sets based on all
#' potential combinations of two or more \code{variables}. Ignored if sets is
#' defined and default = FALSE changes to TRUE if sets = NULL.
#' @param minvar_perset (integer) if \code{all_combinations} = TRUE and sets is
#' not defined, minimum number of variables per combination. This number must
#' be > 1. Default = 2.
#' @param format_in (character) if \code{variables} is character, format of the
#' variables found in the folder. Default = NULL.
#' @param save (logical) whether or not to write sets of variables in a subfolder
#' in the working directory. Default = FALSE.
#' @param format_out (character) format of the variables to be written in
#' \code{output_directory} if \code{save} = TRUE. Default = "GTiff".
#' @param overwrite (logical) whether or not to overwrite existent results in
#' \code{output_directory} if \code{save} = TRUE. Default = FALSE.
#' @param output_directory name of the folder were subsets will be written if
#' \code{save} = TRUE.
#'
#' @export
#'
#' @return
#' A list of character vectors containing the names of the variables to be used
#' in each set. A folder with subfolders (sets) and variables divided into sets.
#'
#' @usage
#' prepare_sets(variables, sets = NULL, all_combinations = FALSE,
#'              minvar_perset = 2, format_in = NULL, save = FALSE,
#'              format_out = "GTiff", overwrite = FALSE, output_directory)

prepare_sets <- function(variables, sets = NULL, all_combinations = FALSE,
                         minvar_perset = 2, format_in = NULL, save = FALSE,
                         format_out = "GTiff", overwrite = FALSE,
                         output_directory) {
  # -----------
  # detecting potential errors
  if (missing(variables)) {
    stop("Argument 'variables' is necessary to perform the analysis.")
  }
  if (missing(sets) & all_combinations == FALSE) {
    stop("Argument 'sets' must be provided if all_cambinations = FALSE. See function's help.")
  }
  clvar <- class(variables)[1]
  if (clvar == "character" | clvar == "SpatRaster") {
    if (clvar == "character") {
      if (is.null(format_in)) {
        stop("Argument 'fomat_in' cannot be NULL if variables is of class character.")
      }
      patt <- paste0(rformat_type(format_in), "$")
      rtype <- rformat_type(format_in)
      var_dir <- variables
      vars <- list.files(variables, pattern = patt, full.names = TRUE)
      variables <- terra::rast(vars)
    }
    var_names <- names(variables)
  } else {
    stop("'variables' must be either character or SpatRaster. See function's help.")
  }

  # -----------
  # Preparing sets
  message("Preparing sets of variables...")
  if (!missing(sets)) {
    if (!missing(sets) & all_combinations == TRUE) {
      message("Argument 'sets' was provided, all_cambinations = TRUE will be ignored.")
    }
    names(sets) <- paste0("set_", 1:length(sets))
  } else {
    var_comb <- lapply(minvar_perset:length(var_names), function(x) {
      comb <- combn(x = var_names, m = x)
      split(comb, col(comb))
    })

    sets <- do.call(c, var_comb)
    names(sets) <- paste0("set_", 1:length(sets))
  }

  message(paste0("\t", length(sets)), " sets of variables prepared")

  # -----------
  # writing sets in output directory
  if (save == TRUE) {
    if (missing(output_directory)) {
      stop("'output_directory' must be defined if 'save' = TRUE.")
    }

    message("Writing sets of variables in output directory:")
    dir.create(output_directory)
    sub_paths <- paste0(output_directory, "/", names(sets))

    lapply(1:length(sub_paths), function(x) {
      dir.create(sub_paths[x])
      out_type <- rformat_type(format_out)
      vars_set <- paste0(sub_paths[x], "/", sets[[x]], out_type)

      if (clvar == "character") {
        vars_comb <- paste0(var_dir, "/", sets[[x]], rtype)
        if (format_in == format_out) {
          suppressMessages(file.copy(from = vars_comb, to = vars_set,
                                     overwrite = overwrite))
        } else {
          lapply(1:length(vars_set), function(y) {
            var_name <- sets[[x]][y]
            terra::writeRaster(variables[[var_name]], filename = vars_set[y],
                               overwrite = overwrite)
          })
        }
      } else {
        lapply(1:length(vars_set), function(y) {
          var_name <- sets[[x]][y]
          terra::writeRaster(variables[[var_name]], filename = vars_set[y],
                             overwrite = overwrite)
        })
      }

      message("\tset ", x, " of ", length(sub_paths), " created")
    })
  }

  # -----------
  # returning results
  results <- list(raster_layers = variables,
                  variable_sets = sets)

  return(results)
}
