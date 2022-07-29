library(usethis)
library(fossilbrush)

# URL downloads -----------------------------------------------------------

# For csv files previously downloaded into local /raw-data folder:
  # datFls <- list.files('./data-raw') |>
  #    (\(.) grep('csv', ., value = TRUE))()
  # bivNm <- grep('Bivalvia', datFls, value = TRUE) |>
  #   (\(.) paste0('data-raw/', .))()
  # bivRaw <- read.csv(bivNm, skip = 21, header = TRUE)
# if metadata included in download file:
  # read.csv(header = TRUE, skip = 19)

# Pliocene bivalve occurrence data
bivPth <- paste0(
  'https://paleobiodb.org/data1.2/occs/list.csv?datainfo',
  '&rowcount&base_name=Bivalvia&taxon_reso=genus&idqual=certain',
  '&interval=Pliocene,Pliocene&time_rule=contain&envtype=marine',
  '&show=coll,paleoloc,env,geo,ref,acconly'
)
bivRaw <- read.csv( url(bivPth), skip = 21, header = TRUE)

# Silurian brachiopod occurrence data
brachPath <- paste0(
  'https://paleobiodb.org/data1.2/occs/list.csv?datainfo',
  '&rowcount&base_name=Brachiopoda&taxon_reso=genus&idqual=certain',
  '&interval=Silurian,Silurian&time_rule=contain&envtype=marine',
  '&show=class,coll,paleoloc,geo,acconly' # 'class' gives all higher taxon units
)
brachRaw <- read.csv( url(brachPath), skip = 21, header = TRUE)

# Silurian marine collections data
collPath <- paste0(
  'https://paleobiodb.org/data1.2/colls/list.csv?datainfo',
  '&rowcount&interval=Silurian,Silurian&time_rule=contain&envtype=marine',
  '&show=paleoloc,geo' # includes environment and geology notes fields
)
collRaw <- read.csv( url(collPath), skip = 18, header = TRUE)
keepRows <- complete.cases( collRaw[, c('paleolng','paleolat')] )
collDat <- collRaw[keepRows,]

# Bivalve cleaning and export ---------------------------------------------

cleanr <- function(dat, keepCols){
  # manually extract genus name (remove subgenera) - replace existing genus col
  dat$genus <- sapply(dat$accepted_name, function(a){
    strsplit(a, ' ') |> unlist() |> head(1)
  })
  # test <- divDyn::cleansp(dat$accepted_name, debug = TRUE, # debug = F returns vector
  #    subgenera = FALSE, misspells = FALSE, stems = FALSE)

  keepRows <- complete.cases(dat[, c('paleolng','paleolat')])
  dat <- dat[keepRows, keepCols]
  rownames(dat) <- 1:nrow(dat)
  return(dat)
}

keepCols <- c('genus', 'paleolng', 'paleolat',
              'collection_no','reference_no',
              'environment',
              # 'geology_comments', # for manual vetting of enviro
              'max_ma', 'min_ma', 'accepted_name'
)
bivalves <- cleanr(bivRaw, keepCols)

# save individual object as rda in /data
use_data(bivalves, internal = FALSE, overwrite = TRUE)

# Brachiopod cleaning and export ------------------------------------------

# brachs   <- cleanr(brachRaw)

# substitute new brach download into verbatim code from fossilbrush vignette:

data(sepkoski)
# trim the Sepkoski Compendium to the relevant entries
sepkoski <- sepkoski[which(sepkoski$PHYLUM == "Brachiopoda"),]

# convert ICS2013 to GTS2020 timescale
brachios <- chrono_scale(brachRaw, srt = "early_interval", end = "late_interval",
                         max_ma = "max_ma", min_ma = "min_ma", verbose = FALSE)
# chrono_scale adds newfad and newlad columns, but doesn't adjust original ages
  # identical(brachios$max_ma, brachRaw$max_ma)
  # identical(brachios$min_ma, brachRaw$min_ma)
  # > [1] TRUE

# combine the datasets
sepkNA <- rep(NA, nrow(sepkoski))
occs <- cbind.data.frame(phylum = c(brachios$phylum, sepkoski$PHYLUM),
                         class = c(brachios$class, sepkoski$CLASS),
                         order = c(brachios$order, sepkoski$ORDER),
                         family = c(brachios$family, sepkNA),
                         genus = c(brachios$genus, sepkoski$GENUS),
                         max_ma = c(brachios$newFAD, sepkoski$RANGE_BASE),
                         min_ma = c(brachios$newLAD, sepkoski$RANGE_TOP),
                         collection_no = c(brachios$collection_no, sepkNA),
                         reference_no = c(brachios$reference_no, sepkNA),
                         paleolng = c(brachios$paleolng, sepkNA),
                         paleolat = c(brachios$paleolat, sepkNA),
                         environment = c(brachios$environment, sepkNA),
                         accepted_name = c(brachios$accepted_name, sepkNA)
                       # geology_comments = c(brachios$geology_comments, sepkNA)
                         )
# define the taxonomic ranks used in the dataset (re-used elsewhere)
b_ranks <- c("phylum", "class", "order", "family", "genus")
# define a list of suffixes to be used at each taxonomic level when scanning for synonyms
b_suff = list(NULL, NULL, NULL, NULL, c("ina", "ella", "etta"))

# manually resolve the clear synonymous spellings
occs$family[which(occs$family == "Disciniidae")] <- "Discinidae"
# occs$genus[which(occs$genus == "Ptychomalotoechia")] <- "Ptychomaletoechia" # none in Silurian
occs$genus[which(occs$genus == "Hipparionix")] <- "Hipparionyx"
# occs$genus[which(occs$genus == "Sphaenospira")] <- "Sphenospira" # none in Silurian
# strip out the PBDB 'missing taxon' format
for(i in c("phylum", "class", "order", "family", "genus")) {
  occs[grep("^NO_", occs[,i]),i] <- NA
}

# clean the data, resolving classifications
occs_c <- check_taxonomy(occs, suff_set = b_suff, ranks = b_ranks, verbose = FALSE,
                         clean_name = TRUE, resolve_duplicates = TRUE, jump = 5)
# extract Sepkoski
sepkoski_c <- occs_c$data[(nrow(brachios) + 1):nrow(occs_c$data),]
# extract PBDB
brachios_c <- occs_c$data[1:nrow(brachios),]
# drop occurrences with older LADs than FADs
brachios_c <- brachios_c[brachios_c$max_ma > brachios_c$min_ma,]

# flag and resolve against the Sepkoski Compendium, collection-wise
revrng <- flag_ranges(x = brachios_c, y = sepkoski_c,
                      xcols = c("genus", "max_ma", "min_ma"),
                      verbose = F
                      )
# Error flag codes from revise_ranges documentation:
# "000" - unchecked, "R1R" - valid,
# "0R0" - both FAD and LAD exceeded,
# "00R" - totally older than range,
# "R00" - totally younger than range,
# "01R" - FAD exceeded, "1R0" - LAD exceeded
brachios_c$errcode <- revrng$occurrence$status
noFlag <- brachios_c$errcode %in% c('000', 'R1R')
brachios_c <- brachios_c[ noFlag, ]

# there are weirdly some PBDB accepted genera converted to NA -
# Flannery-Sutherland attributes to inconsistencies in those occurrences'
# taxonomy fields in original PBDB entries
noName <- is.na(brachios_c$genus)
brachios_c <- brachios_c[ !noName, ]

# geology comments only necessary if manually correct enviro classes
occCols <- c('order', 'family', keepCols)
occSilur <- brachios_c[, occCols]
collCols <- setdiff(keepCols, c('genus','accepted_name'))
collSilur <- collDat[, collCols]
use_data(occSilur, collSilur, internal = FALSE, overwrite = TRUE)
