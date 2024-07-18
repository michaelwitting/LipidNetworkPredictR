## first, load the metaboliteIdmapping package into your R session, 
## when the package is loaded, the data will be available as tibble
library(metaboliteIDmapping)
metabolitesMapping
system.file(package = "metaboliteIDmapping", "/scripts/make-data.R")

## second, search for the mapping table in the AnnotationHub resource interface
library(AnnotationHub)
ah <- AnnotationHub()
datasets <- query(ah, "metaboliteIDmapping")

## Currently, there are two versions of the mapping table.
## AH79817 represents the original ID mapping containing 9 different ID formats
## AH83115 mapping table which also includes common names for each compound
## AH91792 current version of the mapping table that also accounts for tautomers
## For implementing this data in your code, it is recommended to use the 
## AHid for retrieval:
annotation_set <- ah[["AH91792"]]

## RHEA

## load the files
## rhea_reaction_smiles
file <- file.path(
    path.package("LipidNetworkPredictR"), "extdata", "rhea-reaction-smiles.tsv",
    fsep = .Platform$file.sep)
if (!file.exists(file)) {
    file <- "https://ftp.expasy.org/databases/rhea/tsv/rhea-reaction-smiles.tsv"
}
rhea_reaction_smiles <- read.csv(file, sep = "\t", header = FALSE)

## rhea_directions
file <- file.path(
    path.package("LipidNetworkPredictR"), "extdata", "rhea-directions.tsv",
    fsep = .Platform$file.sep)
if (!file.exists(file)) {
    file <- "https://ftp.expasy.org/databases/rhea/tsv/rhea-directions.tsv"
}
rhea_directions <- read.csv(file, sep = "\t", header = TRUE)


## rhea_chebi_smiles
file <- file.path(
    path.package("LipidNetworkPredictR"), "data", "rhea-chebi-smiles.tsv",
    fsep = .Platform$file.sep)
if (!file.exists(file)) {
    file <- "https://ftp.expasy.org/databases/rhea/tsv/rhea-chebi-smiles.tsv"
}
rhea_chebi_smiles <- read.csv(file, sep = "\t", header = FALSE)
