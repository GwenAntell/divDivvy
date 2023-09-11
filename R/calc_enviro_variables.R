#' Convert point environment data to a raster of majority-environment classes
#'
#' Given point occurrences of environmental categories, `classRast` generates
#' a raster grid with cell values specifying the majority environment therein.
#'
#' The `cutoff` threshold is an inclusive bound: environmental incidence
#' proportions greater than or equal to the `cutoff` will assign cell values
#' to the majority environmental class. For instance, if category A represents
#' 65% of occurrences in a cell and `cutoff = 0.65`, the returned value for the
#' cell will be A. If no single category in a cell meets or exceeds the
#' representation necessary to reach the given `cutoff`, the value returned
#' for the cell is `indet.`, indeterminate.
#' Cells lacking environmental occurrences altogether return `NA` values.
#'
#' The `env` object can contain more than two classes, but in many cases it will
#' be less likely for any individual class to attain an absolute majority the
#' more finely divided classes are. For example, if there are three classes,
#' A, B, and C, with relative proportions of 20%, 31%, and 49%, the cell value
#' will be returned as `indet.` because no single class can attain a `cutoff`
#' above 50%, despite class C having the largest relative representation.
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
#' The coordinates associated with points should be given with respect to the
#' same coordinate reference system (CRS) of the target raster grid, e.g. both
#' given in latitude-longitude, Equal Earth projected coordinates, or other CRS.
#' The CRS of a `SpatRaster` object can be retrieved with `terra` function `crs`
#' (with the optional but helpful argument `describe = TRUE`).
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
#' @examples
#' library(terra)
#' # work in Equal Earth projected coordinates
#' prj <- 'EPSG:8857'
#' # generate point occurrences in a small area of Northern Africa
#' n <- 100
#' set.seed(5)
#' x <- runif(n,  0, 30)
#' y <- runif(n, 10, 30)
#' # generate an environmental variable with a latitudinal gradient
#' # more habitat type 0 (e.g. rock) near equator, more 1 (e.g. grassland) to north
#' env <- rbinom(n, 1, prob = (y-10)/20)
#' env[env == 0] <- 'rock'
#' env[env == 1] <- 'grass'
#' # units for Equal Earth are meters, so if we consider x and y as given in km,
#' x <- x * 1000
#' y <- y * 1000
#' ptsDf <- data.frame(x, y, env)
#' # raster for study area at 5-km resolution
#' r <- rast(resolution = 5*1000, crs = prj,
#'           xmin = 0, xmax = 30000, ymin = 10000, ymax = 30000)
#'
#' binRast <- classRast(grid = r, dat = ptsDf, xy = c('x', 'y'),
#'                      env = 'env', cutoff = 0.6)
#' binRast
#'
#' # plot environment classification vs. original points
#' plot(binRast, col = c('lightgreen', 'grey60', 'white'))
#' points(ptsDf[env=='rock', ], pch = 16, cex = 1.2) # occurrences of given habitat
#' points(ptsDf[env=='grass',], pch =  1, cex = 1.2)
#'
#' # classRast can also accept more than 2 environmental classes:
#'
#' # add a 3rd environmental class with maximum occurrence in bottom-left grid cell
#' newEnv <- data.frame('x' = rep(0,       10),
#'                      'y' = rep(10000,   10),
#'                      'env' = rep('new', 10))
#' ptsDf <- rbind(ptsDf, newEnv)
#' binRast <- classRast(grid = r, dat = ptsDf, xy = c('x', 'y'),
#'                      env = 'env', cutoff = 0.6)
#' plot(binRast, col = c('lightgreen', 'grey60', 'purple', 'white'))
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
