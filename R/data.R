#' Paleobiology Database occurrences of Pliocene fossil bivalves
#'
#' A dataset containing the (palaeo)coordinates and genus identifications of
#' 8,000 marine bivalves from the Pliocene (ca. 5.3-2.6 Ma). Records with
#' uncertain or unaccepted taxonomic names, non-marine palaeo-environments,
#' or missing coordinates are excluded from the original download
#' (24 June 2022).
#'
#' @format A data frame with 8095 rows and 9 variables:
#' \describe{
#'   \item{genus}{Latin genus identification. Subgenera are not elevated.}
#'   \item{paleolng, paleolat}{Coordinates of an occurrence, rotated to
#'   its palaeogeographic location with the tectonic plate model of
#'    \href{https://www.earthbyte.org/}{GPlates}}
#'   \item{collection_no, reference_no}{Unique identifiers for the collection and published
#'   reference containing the occurrence}
#'   \item{environment}{One of 23 marine environment categories}
#'   \item{max_ma, min_ma}{Bounds of the age estimate for an occurrence}
#'   \item{accepted_name}{Original identification, including subgenus and
#'   species epithet if applicable, according to the latest PBDB
#'   accepted taxonomy at time of download}
#' }
#' @source \url{https://paleobiodb.org/}
'bivalves'

#' Paleobiology Database collections of Silurian marine fossils
#'
#' A dataset containing the (palaeo)coordinates and recorded marine environment
#' of 8,000 PBDB fossil collections from the Silurian, formatted and downloaded
#' from the Paleobiology Database on 24 June 2022.
#'
#' @format A data frame with 8345 rows and 7 variables:
#' \describe{
#'   \item{paleolng, paleolat}{Coordinates of a collection, rotated to
#'   its palaeogeographic location with the tectonic plate model of
#'    \href{https://www.earthbyte.org/}{GPlates}}
#'   \item{collection_no, reference_no}{Unique identifier for the collection
#'   and its published reference}
#'   \item{environment}{One of 23 marine environment categories}
#'   \item{max_ma, min_ma}{Bounds of the age estimate for a collection}
#' }
#' @source \url{https://paleobiodb.org/}
'collSilur'

#' Paleobiology Database occurrences of Silurian fossil brachiopods
#'
#' A dataset containing the (palaeo)coordinates and genus identifications of
#' 16,000 marine brachiopods from the Silurian (443.1-419 Ma). Records with
#' uncertain or unaccepted taxonomic names, non-marine palaeo-environments,
#' or missing coordinates are excluded from the original download (24 June 2022).
#' Taxonomic synonymisation and removal of stratigraphic outliers
#' follows the `fossilbrush` vignette example of cross-correlation with the
#' Sepkoski range-through database (`fossilbrush::sepkoski`).
#'
#' @format A data frame with 15909 rows and 11 variables:
#' \describe{
#'   \item{order, family, genus}{Latin order, family, and genus name,
#'   as synonymised against Sepkoski database}
#'   \item{paleolng, paleolat}{Coordinates of an occurrence, rotated to
#'   its palaeogeographic location with the tectonic plate model of
#'    \href{https://www.earthbyte.org/}{GPlates}}
#'   \item{collection_no, reference_no}{Unique identifiers for the collection
#'   and published reference containing the occurrence}
#'   \item{environment}{One of 23 marine environment categories}
#'   \item{max_ma, min_ma}{Bounds of the age estimate for an occurrence,
#'   according to the ICS 2013 geologic time scale.}
#'   \item{accepted_name}{Original identification, including subgenus and
#'   species epithet if applicable, according to the latest PBDB
#'   accepted taxonomy at time of download}
#' }
#' @source \url{https://paleobiodb.org/}
'occSilur'
