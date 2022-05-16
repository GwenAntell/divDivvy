#' Paleobiology Database occurrences of Pliocene fossil bivalves
#'
#' A dataset containing the (palaeo)coordinates and genus identifications of
#' 8,000 marine bivalves from the Pliocene (ca. 5.3-2.6 Ma). Records with
#' uncertain or unaccepted taxonomic names, non-marine palaeo-environments,
#' or missing coordinates are excluded from the original download
#' (13 May 2022).
#'
#' @format A data frame with 8072 rows and 10 variables:
#' \describe{
#'   \item{genus}{Latin genus identification. Subgenera are not elevated.}
#'   \item{paleolng, paleolat}{Coordinates of an occurrence, rotated to
#'   its palaeogeographic location with the tectonic plate model of
#'    \link[=https://www.earthbyte.org/]{GPlates}}
#'   \item{collection_no, reference_no}{Unique identifiers for the collection and published
#'   reference containing the occurrence}
#'   \item{environment}{One of 23 marine environment categories}
#'   \item{geology_comments}{Verbatim geological description notes}
#'   \item{max_ma, min_ma}{Bounds of the age estimate for an occurrence}
#'   \item{accepted_name}{Original identification, including subgenus and
#'   species epithet if applicable, according to the latest PBDB
#'   accepted taxonomy at time of download}
#' }
#' @source \url{https://paleobiodb.org/}
'bivalves'
