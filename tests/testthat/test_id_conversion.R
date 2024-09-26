## define the mapping table
library(metaboliteIDmapping)
library(AnnotationHub)
ah <- AnnotationHub()
datasets <- query(ah, "metaboliteIDmapping")
annotation_set <- ah[["AH91792"]]

## function translate_id_to_ChEBI
test_that("translate_id_to_ChEBI works", {
    
    ids <- c("HMDB0004947", "HMDB0004949", "HMDB0011763", "HMDB0004974", "HMDB0011594")
    id_type <- "HMDB"
    
    ## run the function
    translations <- translate_id_to_ChEBI(ids = ids, id_type = "HMDB",
        annotation_set = annotation_set)
    expect_equal(translations[["HMDB0004947"]], "72956")
    expect_equal(translations[["HMDB0004949"]], "72956")
    expect_equal(translations[["HMDB0011763"]], "74100")
    expect_equal(translations[["HMDB0004974"]], "62107")
    expect_equal(translations[["HMDB0011594"]], "62109")
    translations <- translate_id_to_ChEBI(ids = ids, id_type = "CAS",
        annotation_set = annotation_set)
    expect_equal(translations[["HMDB0004947"]], "character(0)")
    expect_equal(translations[["HMDB0004949"]], "character(0)")
    expect_equal(translations[["HMDB0011763"]], "character(0)")
    expect_equal(translations[["HMDB0004974"]], "character(0)")
    expect_equal(translations[["HMDB0011594"]], "character(0)")
    
    expect_error(translate_id_to_ChEBI(ids = ids, id_type = "foo",
        annotation_set = annotation_set), "should be one of ")
    expect_error(translate_id_to_ChEBI(ids = ids, id_type = "HMDB",
        annotation_set = NULL), "should be one of ")
    expect_error(translate_id_to_ChEBI(ids = ids, id_type = NULL,
        annotation_set = NULL), "applied to an object of class")
})

## function translate_ChEBI_to_SMILES
test_that("translate_ChEBI_to_SMILES works", {
    
    ids <- list("HMDB0004947" = "72956", "HMDB0004949" = "72959", 
        "HMDB0011763" = "74100", "HMDB0004974" = "62107", 
        "HMDB0011594" = "62109")
    
    ## load rhea_chebi_smiles
    file <- file.path(
        path.package("LipidNetworkPredictR"), "extdata", "rhea-chebi-smiles.tsv",
        fsep = .Platform$file.sep)
    rhea_chebi_smiles <- read.csv(file, sep = "\t", header = FALSE)
    
    ## run the function
    translations <- translate_ChEBI_to_SMILES(ids = ids, rhea_chebi_smiles = rhea_chebi_smiles)
    expect_equal(translations[["HMDB0004947"]], "CCCCCCCCCCCCC/C=C/[C@@H](O)[C@H](CO)NC(=O)CCCCCCCCCCC")
    expect_equal(translations[["HMDB0004949"]], "CCCCCCCCCCCCC/C=C/[C@@H](O)[C@H](CO)NC(=O)CCCCCCCCCCCCCCC")
    expect_equal(translations[["HMDB0011763"]], "CCCCCCCC/C=C\\CCCCCCCC(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCCCC")
    expect_equal(translations[["HMDB0004974"]], NA)
    expect_equal(translations[["HMDB0011594"]], NA)
    
    translations <- translate_ChEBI_to_SMILES(ids = ids, rhea_chebi_smiles = NULL)
    expect_equal(translations[["HMDB0004947"]], NA)
    expect_equal(translations[["HMDB0004949"]], NA)
    expect_equal(translations[["HMDB0011763"]], NA)
    expect_equal(translations[["HMDB0004974"]], NA)
    expect_equal(translations[["HMDB0011594"]], NA)
    expect_error(
        translate_ChEBI_to_SMILES(ids = NULL, rhea_chebi_smiles = rhea_chebi_smiles), 
        "is not a list ")
    expect_equal(
        translate_ChEBI_to_SMILES(ids = list(), rhea_chebi_smiles = rhea_chebi_smiles), 
        list())
})

## function select_substrates_or_products
test_that("select_substrates_or_products works", {
    
    ## given a rection that describes the change of (3S)-3-Hydroxy-2-butanone
    ## to its stereoisomer (3R)-3-Hydroxy-2-butanone
    reaction <- list(c("CC(=O)[C@H](C)O",  "CC(=O)[C@@H](C)O"))
    
    ## cols = 2
    cols <- 2
    res <- LipidNetworkPredictR:::select_substrates_or_products(reaction = reaction, 
        cols = cols, type = "substrates")
    expect_equal(length(res), 1)
    expect_equal(res[[1]], "CC(=O)[C@@H](C)O")
    res <- LipidNetworkPredictR:::select_substrates_or_products(reaction = reaction, 
        cols = cols, type = "products")
    expect_equal(length(res), 1)
    expect_equal(res[[1]], "CC(=O)[C@H](C)O")
    
    ## cols = 3
    cols <- 3
    res <- LipidNetworkPredictR:::select_substrates_or_products(reaction = reaction, 
        cols = cols, type = "substrates")
    expect_equal(length(res), 1)
    expect_equal(res[[1]], "CC(=O)[C@@H](C)O")
    res <- LipidNetworkPredictR:::select_substrates_or_products(reaction = reaction, 
        cols = cols, type = "products")
    expect_equal(length(res), 1)
    expect_equal(res[[1]], "CC(=O)[C@H](C)O")
})

## function find_RHEA_ids_from_SMILES
test_that("find_RHEA_ids_from_SMILES works", {
    
    ids <- list(
        "HMDB0004947" = "CCCCCCCCCCCCC/C=C/[C@@H](O)[C@H](CO)NC(=O)CCCCCCCCCCC",
        "HMDB0004949" = "CCCCCCCCCCCCC/C=C/[C@@H](O)[C@H](CO)NC(=O)CCCCCCCCCCCCCCC",
        "HMDB0011763" = "CCCCCCCC/C=C\\CCCCCCCC(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCCCC",
        "HMDB0004974" = NA, 
        "HMDB0011594" = NA)
    
    ## rhea_reaction_smiles
    file <- file.path(
        path.package("LipidNetworkPredictR"), "extdata", "rhea-reaction-smiles.tsv",
        fsep = .Platform$file.sep)
    rhea_reaction_smiles <- read.csv(file, sep = "\t", header = FALSE)
 
    ## rhea_directions
    file <- file.path(
        path.package("LipidNetworkPredictR"), "extdata", "rhea-directions.tsv",
        fsep = .Platform$file.sep)
    rhea_directions <- read.csv(file, sep = "\t", header = TRUE)
    
    ## run the function
    ## type = "both"
    res <- find_RHEA_ids_from_SMILES(ids = ids, 
        rhea_reaction_smiles = rhea_reaction_smiles,
        rhea_directions = rhea_directions, type = "both")
    expect_equal(res[["HMDB0004947"]], 
        c(41292, 41293, 70256, 70257, 70308, 70309))
    expect_equal(res[["HMDB0004949"]], 
        c(38892, 38893, 45721, 45722, 46341, 46342, 36688, 36689, 45645, 45646,
            58297, 58298, 58317, 58318, 69700, 69701, 69704, 69705))
    expect_equal(res[["HMDB0011763"]], 
        c(45373, 45374, 36576, 36577))
    expect_equal(res[["HMDB0004974"]], NA)
    expect_equal(res[["HMDB0011594"]], NA)
    
    ## type = "substrates"
    res <- find_RHEA_ids_from_SMILES(ids = ids, 
        rhea_reaction_smiles = rhea_reaction_smiles,
        rhea_directions = rhea_directions, type = "substrates")
    expect_equal(res[["HMDB0004947"]], 
        c(41292, 41293))
    expect_equal(res[["HMDB0004949"]], 
        c(38892, 38893, 45721, 45722, 46341, 46342))
    expect_equal(res[["HMDB0011763"]], 
        c(45373, 45374))
    expect_equal(res[["HMDB0004974"]], NA)
    expect_equal(res[["HMDB0011594"]], NA)
    
    ## type = "products"
    res <- find_RHEA_ids_from_SMILES(ids = ids, 
        rhea_reaction_smiles = rhea_reaction_smiles,
        rhea_directions = rhea_directions, type = "products")
    expect_equal(res[["HMDB0004947"]], 
        c(70256, 70257, 70308, 70309))
    expect_equal(res[["HMDB0004949"]], 
        c(36688, 36689, 45645, 45646, 58297, 58298, 58317, 58318, 
            69700, 69701, 69704, 69705))
    expect_equal(res[["HMDB0011763"]], 
                 c(36576, 36577))
    expect_equal(res[["HMDB0004974"]], NA)
    expect_equal(res[["HMDB0011594"]], NA)
})

## function find_RHEA_ids_from_ids
test_that("find_RHEA_ids_from_ids works", {
    
    ids <- c("HMDB0004947", "HMDB0004949", "HMDB0011763", "HMDB0004974", "HMDB0011594")
    id_type <- "HMDB"
    
    ## load rhea_chebi_smiles
    file <- file.path(
        path.package("LipidNetworkPredictR"), "extdata", "rhea-chebi-smiles.tsv",
        fsep = .Platform$file.sep)
    rhea_chebi_smiles <- read.csv(file, sep = "\t", header = FALSE)
    
    ## rhea_reaction_smiles
    file <- file.path(
        path.package("LipidNetworkPredictR"), "extdata", "rhea-reaction-smiles.tsv",
        fsep = .Platform$file.sep)
    rhea_reaction_smiles <- read.csv(file, sep = "\t", header = FALSE)
    
    ## rhea_directions
    file <- file.path(
        path.package("LipidNetworkPredictR"), "extdata", "rhea-directions.tsv",
        fsep = .Platform$file.sep)
    rhea_directions <- read.csv(file, sep = "\t", header = TRUE)
    
    ## run the function
    ## type = "both"
    res <- find_RHEA_ids_from_ids(ids = ids, id_type = id_type, 
        annotation_set = annotation_set, rhea_chebi_smiles = rhea_chebi_smiles,
        rhea_reaction_smiles = rhea_reaction_smiles, 
        rhea_directions = rhea_directions, type = "both")
    expect_equal(res[["HMDB0004947"]], 
        c(41292, 41293, 70256, 70257, 70308, 70309))
    expect_equal(res[["HMDB0004949"]], 
        c(38892, 38893, 45721, 45722, 46341, 46342, 36688, 36689, 45645, 45646,
            58297, 58298, 58317, 58318, 69700, 69701, 69704, 69705))
    expect_equal(res[["HMDB0011763"]], 
        c(45373, 45374, 36576, 36577))
    expect_equal(res[["HMDB0004974"]], NA)
    expect_equal(res[["HMDB0011594"]], NA)
    
    ## type = "substrates"
    res <- find_RHEA_ids_from_ids(ids = ids, id_type = id_type, 
        annotation_set = annotation_set, rhea_chebi_smiles = rhea_chebi_smiles,
        rhea_reaction_smiles = rhea_reaction_smiles, 
        rhea_directions = rhea_directions, type = "substrates")
    expect_equal(res[["HMDB0004947"]], 
        c(41292, 41293))
    expect_equal(res[["HMDB0004949"]], 
        c(38892, 38893, 45721, 45722, 46341, 46342))
    expect_equal(res[["HMDB0011763"]], 
        c(45373, 45374))
    expect_equal(res[["HMDB0004974"]], NA)
    expect_equal(res[["HMDB0011594"]], NA)
    
    ## type = "products"
    res <- find_RHEA_ids_from_ids(ids = ids, id_type = id_type, 
        annotation_set = annotation_set, rhea_chebi_smiles = rhea_chebi_smiles,
        rhea_reaction_smiles = rhea_reaction_smiles, 
        rhea_directions = rhea_directions, type = "products")
    expect_equal(res[["HMDB0004947"]], 
        c(70256, 70257, 70308, 70309))
    expect_equal(res[["HMDB0004949"]], 
        c(36688, 36689, 45645, 45646, 58297, 58298, 58317, 58318, 
            69700, 69701, 69704, 69705))
    expect_equal(res[["HMDB0011763"]], 
        c(36576, 36577))
    expect_equal(res[["HMDB0004974"]], NA)
    expect_equal(res[["HMDB0011594"]], NA)
})
