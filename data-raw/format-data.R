library(magrittr)
library(dplyr)
library(usethis)

# data examples: 1. bivalves, 2. dinos, 3. trilos

# For csv files previously downloaded into local /raw-data folder:
  # datFls <- list.files('./data-raw') %>%
  #   grep('csv', ., value = TRUE)
  # bivNm <- paste0('data-raw/', datFls[3])
  # bivRaw <- read.csv(bivNm, skip = 21, header = TRUE)
  # # read.csv(header = TRUE, skip = 19) # if metadata included in download file

bivPth <- paste0(
  'https://paleobiodb.org/data1.2/occs/list.csv?datainfo',
  '&rowcount&base_name=Bivalvia&taxon_reso=genus&idqual=certain',
  '&interval=Pliocene,Pliocene&time_rule=contain&envtype=marine',
  '&show=coll,paleoloc,env,geo,ref,acconly'
)
bivRaw <- url(bivPth) %>% read.csv( ., skip = 21, header = TRUE)


# Clean and save data -----------------------------------------------------

cleanr <- function(dat){
  # manually extract genus name (remove subgenera) - replace existing genus col
  dat$genus <- sapply(dat$accepted_name, function(a){
    strsplit(a, ' ') %>% unlist %>% .[1]
  })
  # test <- divDyn::cleansp(dat$accepted_name, debug = TRUE, # debug = F returns vector
  #    subgenera = FALSE, misspells = FALSE, stems = FALSE)

  keepRows <- complete.cases(dat[, c('paleolng','paleolat')])
  keepCols <- c('genus', 'paleolng', 'paleolat',
                'collection_no','reference_no',
                'environment', 'geology_comments',
                'max_ma', 'min_ma', 'accepted_name' # 'geoplate'
  )
  dat <- dat[keepRows, keepCols]
  rownames(dat) <- 1:nrow(dat)
  return(dat)
}

bivalves <- cleanr(bivRaw)

# save individual object as rda in /data
use_data(bivalves, internal = FALSE, overwrite = TRUE)
