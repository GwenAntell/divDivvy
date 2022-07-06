#' Convert point environment data to a raster of majority-environment classes
#'
#' Given point occurrences of environmental categories, `classRast` generates
#' a raster grid with cell values specifying the majority environment therein.
#'
#' The `cutoff` threshold is an inclusive bound: environmental incidence
#' proportions greater than or equal to the `cutoff` will assign cell values
#' to the majority environmental class. If no single category attains a
#' sufficient majority in a cell, the cell value is indeterminate (`indet.`).
#' Cells lacking environmental occurrences altogether have `NA` values.
#'
#' Missing environment values in the point data should be coded as `NA`,
#' not e.g. `'unknown'`. `classRast` ignores `NA` occurrences when tallying
#' environmental occurrences against the `cutoff`. However, `NA` occurrences
#' still count when determining `NA` status of cells in the raster: a cell
#' containing occurrences of only `NA` value is classified as `indet.`, not `NA`.
#' That is, any grid cell encompassing original point data is non-`NA`.
#'
#' Antell and others (2020) set a `cutoff` of 0.8, based on the same threshold
#' Nürnberg and Aberhan (2013) used to classify environmental preferences for taxa.
#'
#' @param grid A `SpatRaster` to use as a template for the
#' resolution, extent, and coordinate reference system of the returned object.
#' Values can be empty.
#' @param dat Either a `data.frame` or `matrix` for which `xy` and `env` are
#' column names, or an empty argument.
#' @param xy A vector specifying the name or numeric position of columns
#' in `dat` containing coordinates, if `dat` is supplied, or a 2-column
#' `data.frame` or `matrix` of coordinate values.
#' @param env The name or numeric position of the column in `dat` containing a
#' categorical environmental variable, if `dat` is supplied, or a vector of
#' environmental values.
#' @param cutoff The (decimal) proportion of incidences of an environmental
#' category above which a cell will be assigned as that category.
#' `cutoff` must be greater than 0.5.
#'
#' @return A raster of class `SpatRaster` defined by the `terra` package
#'
#' @references
#' \insertRef{Antell2020}{divvy}
#'
#' \insertRef{Nurnberg2013}{divvy}
#'
#' @export

classRast <- function(grid, dat = NULL, xy, env, cutoff){
  if (! is.null(dat)){
    xy  <- dat[, xy]
    env <- dat[, env]
  }
  naEnv <- is.na(env)
  env <- env[!naEnv]
  xy  <-  xy[!naEnv,]

  # tally occurrences with each enviro within raster cells
  lvls <- unique(env)
  nLvl <- length(lvls)
  lyrL <- lapply(lvls, function(lvl){
    xyMat <- data.matrix( xy[env==lvl,] )
    terra::rasterize(xyMat, grid, fun = length)
  }  )
  lyrs <- terra::rast(lyrL)
  names(lyrs) <- lvls

  if (cutoff <= 0.5){ stop('cutoff must be greater than 0.5') }
  # set NA values to 0 so when summed across layers, result is not NA
  reclass <- matrix(c(NA, 0), ncol = 2)
  lyrsNoNA <- terra::classify(lyrs, reclass)
  prop <- lyrs / sum(lyrsNoNA) # convert to proportion
  cuts <- rbind(c(0, cutoff, 0), # binarise
                c(cutoff, 1, 1))
  # values = cutoff classified as above the cutoff (left-closed interval)
  mjtys <- terra::classify(prop, cuts,
                    include.lowest = TRUE,
                    right = FALSE)
  terra::values(grid) <- NA # ensure blank starting raster
  for (l in 1:nLvl){
    grid[ mjtys[[l]]==1 ] <- l
  }

  # change NA cells to indet. if they do have occurrence data
  # but no single environmental class exceeds cutoff
  lvls <- c(lvls, 'indet.')
  nLvl <- nLvl + 1
  reclass <- cbind(1, nLvl)
  anydat <- terra::rasterize( data.matrix(xy), grid)
  anydat <- terra::classify(anydat, reclass)
  grid <- terra::cover(grid, anydat)
  attrs <- data.frame('num' = 1:nLvl)
  attrs$mainClass <- lvls
  levels(grid) <- attrs
  return(grid)
}
