#' Paleobiology Database occurrences of Pliocene fossil bivalves
#'
#' A dataset containing the (paleao)coordinates and genus identifications of
#' nearly 10,000 marine bivalves from the Pliocene (ca. 5.3-2.6 Ma).
#'
#' @format A data frame with 9269 rows and 8 variables:
#' \describe{
#'   \item{genus}{Latin genus name of a fossil, synonymised according to the
#'   latest taxonomy of the database at time of download}
#'   \item{paleolng, paleolat}{Coordinates of an occurrence, rotated to
#'   its palaeogeographic location}
#'   \item{geoplate}{Identifier of the geologic plate on which the occurrence
#'   falls, specified by the paleogeographic model of
#'    \link[=https://www.earthbyte.org/]{GPlates}}
#'   \item{max_ma, min_ma}{Bounds of the age estimate for an occurrence}
#'   \item{collection_no, reference_no}{Unique identifiers for the collection and published
#'   reference containing the occurrence}
#' }
#' @source \url{https://paleobiodb.org/}
'bivalves'
