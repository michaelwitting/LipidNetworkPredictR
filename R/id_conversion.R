
## a user would typically have a matrix with samples in columns and metabolites
## in the rows, the rownames would be of kind 'id', e.g. HMDB ids or Name

## define a function that translates the rownames/ids into ChEBI ids using annotation_set


#' @name translate_id_to_ChEBI
#' 
#' @title Translate ID to ChEBI
#' 
#' @description
#' The function \code{translate_id_to_ChEBI} will translate identifiers 
#' (\code{ids}) of type \code{id_type} to ChEBI identifiers. 
#' 
#' @details
#' The function will return a named list. Each entry will contain the
#' corresponding ChEBI ids. In case there are no corresponding ChEBI ids
#' to \code{ids}, the list will contain the entry \code{NA}.
#' 
#' @param ids character, identifiers to be translated to ChEBI ids
#' @param id_type character(1), one of \code{colnames(annotation_set)}
#' @param annotation_set data.frame containing the mappings between different
#' ids
#' 
#' @author Thomas Naake
#' 
#' @importFrom dplyr filter
#' 
#' @return list
#' 
#' @export
#' 
#' @examples
#' ## define the identifiers and the id type 
#' ids <- c("HMDB0004947", "HMDB0004949", "HMDB0011763", "HMDB0004974", "HMDB0011594")
#' id_type <- "HMDB"
#' 
#' ## define the mapping table
#' library(AnnotationHub)
#' ah <- AnnotationHub()
#' datasets <- query(ah, "metaboliteIDmapping")
#' annotation_set <- ah[["AH91792"]]
#' 
#' ## run the function
#' translate_id_to_ChEBI(ids = ids, id_type = id_type,
#'     annotation_set = annotation_set)
translate_id_to_ChEBI <- function(ids, id_type = colnames(annotation_set),
    annotation_set = annotation_set) {

    id_type <- match.arg(id_type)

    ## to reduce the size of annotation_set, filter the rows that are not
    ## NA in the column id_type, get all ids of id_type
    annotation_set <- annotation_set |>
        dplyr::filter( !is.na(get(id_type)) )
    ids_all <- annotation_set[[id_type]]

    ## get all the mappings from ids in annotation_set, then select the ChEBI 
    ## ids
    mappings <- lapply(ids, function(id_i) {
        annotation_set |>
            dplyr::filter(ids_all == id_i)
    })
    chebi_ids <- lapply(mappings, function(mappings_i) {
        mappings_i |>
            dplyr::select(ChEBI) |>
            as.character()
    })
    names(chebi_ids) <- ids

    ## return the chebi_ids
    chebi_ids
}


## from the ChEBI ids find the SMILES using rhea_chebi_smiles

#' @name translate_ChEBI_to_SMILES
#' 
#' @title Translate ChEBI to SMILES
#' 
#' @description
#' The function \code{translate_ChEBI_to_SMILES} takes ChEBI identifiers
#' (\code{ids}) and translate these via the dictionary \code{rhea_chebi_smiles}
#' to SMILES.
#'  
#' @details
#' The function will return a named list. Each entry will contain the
#' corresponding SMILES ids. In case there are no corresponding SMILES ids
#' to \code{ids}, the list will contain the entry \code{NA}.
#' 
#' @param ids list, entries contain the 
#' @param rhea_chebi_smiles data.frame containing the mappings between ChEBI 
#' identifiers and SMILES identifiers
#' 
#' @author Thomas Naake
#' 
#' @return list
#' 
#' @export
#' 
#' @examples
#' ids <- list("HMDB0004947" = "72956", "HMDB0004949" = "72959", 
#'     "HMDB0011763" = "74100", "HMDB0004974" = "62107", 
#'     "HMDB0011594" = "62109")
#'     
#' ## load rhea_chebi_smiles
#' file <- file.path(
#'     path.package("LipidNetworkPredictR"), "extdata", "rhea-chebi-smiles.tsv",
#'     fsep = .Platform$file.sep)
#' rhea_chebi_smiles <- read.csv(file, sep = "\t", header = FALSE)
#' 
#' ## run the function
#' translate_ChEBI_to_SMILES(ids = ids, rhea_chebi_smiles = rhea_chebi_smiles)
translate_ChEBI_to_SMILES <- function(ids, rhea_chebi_smiles) {
    
    if (!is.list(ids)) 
        stop("'ids' is not a list")
    
    ids <- lapply(ids, function(ids_i) paste0("CHEBI:", ids_i))
    
    smiles <- lapply(ids, function(ids_i) {
        smiles <- rhea_chebi_smiles[rhea_chebi_smiles[, 1] == ids_i, 2]
        if (length(smiles) == 0) 
            smiles <- NA
        smiles
    })
    
    smiles
}

## from the SMILES find the RHEA reaction id using rhea_reaction_smiles

#' @name select_substrates_or_products
#' 
#' @title Select substrates or products
#' 
#' @description
#' The function \code{select_substrates_or_products} selects either 
#' substrates or products of (a) reaction(s) stored in a \code{list} object
#' depending on the \code{cols} argument.
#' 
#' Depending on the \code{type} argument the function will either return 
#' a list of substrates (\code{type = "substrates"}) or products
#' (\code{type = "products"}).
#' 
#' @details
#' Each list entry stores the substrates/products or products/substrates of
#' a reaction. The direction of this reaction is codified in the 
#' \code{cols} argument. 
#' If \code{cols} = 1: first element is equal to substrates, 
#' if \code{cols} = 2: first entry is equal to substrates,
#' if \code{cols} = 3: first entry is equal to products, 
#' if \code{cols} = 4: first entry is equal to substrates.
#'  
#' @param reaction list, each entry contains a character vector of length 2
#' @param cols numeric of length \code{length(list)}
#' @param type character(1) either substrates or products
#' 
#' @author Thomas Naake
#' 
#' @return list
#' 
#' @examples
#' ## given a rection that describes the change of (3S)-3-Hydroxy-2-butanone
#' ## to its stereoisomer (3R)-3-Hydroxy-2-butanone
#' reaction <- list(c("CC(=O)[C@H](C)O",  "CC(=O)[C@@H](C)O"))
#'
#' ## cols = 2
#' cols <- 2
#' LipidNetworkPredictR:::select_substrates_or_products(reaction = reaction, 
#'     cols = cols, type = "substrates")
#' LipidNetworkPredictR:::select_substrates_or_products(reaction = reaction, 
#'     cols = cols, type = "products")
#' 
#' ## cols = 3
#' cols <- 3
#' LipidNetworkPredictR:::select_substrates_or_products(reaction = reaction, 
#'     cols = cols, type = "substrates")
#' LipidNetworkPredictR:::select_substrates_or_products(reaction = reaction, 
#'     cols = cols, type = "products")
select_substrates_or_products <- function(reaction, cols, 
    type = c("substrates", "products")) {
    
    if (!is.list(reaction))
        stop("'reaction' has to be a list")
    
    if (length(reaction) != length(cols))
        stop("'length(reaction)' has to be identical to 'length(cols)'")
    
    if (!all(cols %in% c(1, 2, 3, 4))) 
        stop("'cols' has to be either 1, 2, 3, or 4")
    
    type <- match.arg(type)
    
    ## check the direction of the reaction
    ## if column 1: first entry = substrates, column 2: first entry = substrates
    ## column 3: first entry = products, column 4: first entry = substrates
    lapply(seq_along(reaction), 
        function(entry_i) {
            if (type == "substrates")
                cols_i <- ifelse(cols[entry_i] %in% c(1, 2, 4), 1, 2)
            if (type == "products")
                cols_i <- ifelse(cols[entry_i] %in% c(1, 2, 4), 2, 1)
            reaction[[entry_i]][cols_i]
        }) |>
        unlist() |>
        strsplit(split = "[.]")
}


#' @name obtain_substrates_products
#' 
#' @title Obtain substrates and products
#' 
#' @description
#' The function \code{obtain_substrates_products} obtains from the 
#' reactions stored in the \code{rhea_reaction_smiles} object the 
#' substrates and products.
#'  
#' @details
#' The function \code{obtain_substrates_products} uses the information on the 
#' directionality of the reactions stored in the \code{rhea_directions}
#' object. This object stores in the first column reactions that have no
#' information on directionality, in the second column left-to-right reactions,
#' in the third column right-to-left reactions, and in the fourth column 
#' bidirectional reactions. To comply with downstream functions 
#' \code{obtain_substrates_products} will interpret the entries of the reactions
#' in \code{rhea_reaction_smiles} as in the following:
#' 
#' \itemize{
#'     \item if reaction id in column 1 of \code{rhea_directions}: 
#'         first entry = substrates, 
#'     \item if reaction id in column 2 of \code{rhea_directions}:
#'         first entry = substrates,
#'     \item if reaction id in column 3 of \code{rhea_directions}:
#'         first entry = products, 
#'     \item if reaction id in column 4 of \code{rhea_directions}:
#'         first entry = substrates.
#' }
#' 
#' @param rhea_reaction_smiles data.frame containing the mappings between RHEA
#' reaction ids and SMILES identifiers
#' @param rhea_directions data.frame containing information on the 
#' directionality of RHEA directions
#' 
#' @author Thomas Naake
#' 
#' @return list
#' 
#' @examples
#' ## rhea_reaction_smiles
#' file <- file.path(
#'     path.package("LipidNetworkPredictR"), "extdata", "rhea-reaction-smiles.tsv",
#'     fsep = .Platform$file.sep)
#' rhea_reaction_smiles <- read.csv(file, sep = "\t", header = FALSE)
#'  
#' ## rhea_directions
#' file <- file.path(
#'     path.package("LipidNetworkPredictR"), "extdata", "rhea-directions.tsv",
#'     fsep = .Platform$file.sep)
#' rhea_directions <- read.csv(file, sep = "\t", header = TRUE)
#' 
#' ## run the function
#' LipidNetworkPredictR:::obtain_substrates_products(
#'     rhea_reaction_smiles = rhea_reaction_smiles,
#'     rhea_directions = rhea_directions)
obtain_substrates_products <- function(rhea_reaction_smiles, rhea_directions) {
    
    ## split each reaction into a first and second part, from this part it 
    ## is not directly derivable which part refers to products/substrates
    reaction <- rhea_reaction_smiles[, 2] |>
        strsplit(split = ">>")
    
    ## check now the first column of rhea_reaction_smiles and check the direction
    ## of the reaction
    ## if column 1: first entry = substrates, column 2: first entry = substrates
    ## column 3: first entry = products, column 4: first entry = substrates
    ## get the column
    cols <- lapply(rhea_reaction_smiles[, 1], 
        function(entry_i) {
            as.numeric(which(rhea_directions == entry_i, arr.ind = TRUE)[, "col"])
    }) |>
        unlist()
    
    ## choose the substrate/product based on the column
    substrates <- select_substrates_or_products(reaction = reaction, 
        cols = cols, type = "substrates")
    products <- select_substrates_or_products(reaction = reaction, 
        cols = cols, type = "products")
    
    list(substrates = substrates, products = products)
}


#' @name find_RHEA_ids_from_SMILES
#' 
#' @title Find RHEA ids from SMILES  
#' 
#' @description
#' The function \code{find_RHEA_ids_from_SMILES} returns the RHEA ids
#' that are involved in reactions of the given SMILES (\code{ids}). 
#'  
#' @details
#' The returned result might differ depending on the \code{type} object:
#' 
#' \itemize{
#'     \item if \code{type = "both"}, the function will return the RHEA 
#'     identifiers for the SMILES found both in the substrates or products of 
#'     the reactions,
#'     \item if \code{type = "substrates"}, the function will return the RHEA 
#'     identifiers for the SMILES found in the substrates of the reactions,
#'     \item if \code{type = "products"}, the function will return the RHEA 
#'     identifiers for the SMILES found in the products of the reactions.
#' }
#' 
#' @param ids list, entries contain the SMILES identifiers.
#' @param rhea_reaction_smiles data.frame containing the mappings between ChEBI 
#' identifiers and SMILES identifiers
#' @param rhea_directions data.frame containing information on the 
#' directionality of RHEA directions
#' @param type character(1), one of "both", "substrates", or "products"
#' 
#' @author Thomas Naake
#' 
#' @return list
#' 
#' @export
#' 
#' @examples
#' ids <- list(
#'     "HMDB0004947" = "CCCCCCCCCCCCC/C=C/[C@@H](O)[C@H](CO)NC(=O)CCCCCCCCCCC",
#'     "HMDB0004949" = "CCCCCCCCCCCCC/C=C/[C@@H](O)[C@H](CO)NC(=O)CCCCCCCCCCCCCCC",
#'     "HMDB0011763" = "CCCCCCCC/C=C\\CCCCCCCC(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCCCC",
#'     "HMDB0004974" = NA, 
#'     "HMDB0011594" = NA)
#'     
#' ## rhea_reaction_smiles
#' file <- file.path(
#'     path.package("LipidNetworkPredictR"), "extdata", "rhea-reaction-smiles.tsv",
#'     fsep = .Platform$file.sep)
#' rhea_reaction_smiles <- read.csv(file, sep = "\t", header = FALSE)
#' 
#' ## rhea_directions
#' file <- file.path(
#'     path.package("LipidNetworkPredictR"), "extdata", "rhea-directions.tsv",
#'     fsep = .Platform$file.sep)
#' rhea_directions <- read.csv(file, sep = "\t", header = TRUE)
#'
#' ## run the function
#' ## type = "both"
#' find_RHEA_ids_from_SMILES(ids = ids, 
#'     rhea_reaction_smiles = rhea_reaction_smiles,
#'     rhea_directions = rhea_directions, type = "both")
#' ## type = "substrates"
#' find_RHEA_ids_from_SMILES(ids = ids, 
#'     rhea_reaction_smiles = rhea_reaction_smiles,
#'     rhea_directions = rhea_directions, type = "substrates")
#' ## type = "products"
#' find_RHEA_ids_from_SMILES(ids = ids, 
#'     rhea_reaction_smiles = rhea_reaction_smiles,
#'     rhea_directions = rhea_directions, type = "products")
find_RHEA_ids_from_SMILES <- function(ids, rhea_reaction_smiles, 
        rhea_directions, type = c("both", "substrates", "products")) {
    
    if (!is.list(ids)) 
        stop("'ids' is not a list")
    type <- match.arg(type)
    
    ## obtain the substrates and products from the reactions, the entries
    ## wll be of type SMILES
    substrates_products <- obtain_substrates_products(
        rhea_reaction_smiles = rhea_reaction_smiles, 
        rhea_directions = rhea_directions)
    substrates <- substrates_products[[1]]
    products <- substrates_products[[2]]
    
    substrates_inds <- lapply(ids, function(ids_i) {
        if (length(ids_i) > 0) {
            lapply(substrates, function(substrates_i)
                any(substrates_i %in% ids_i)) |>
                unlist() |>
                which()
    }})
    products_inds <- lapply(ids, function(ids_i) {
        if (length(ids_i) > 0) {
            lapply(products, function(products_i)
                any(products_i %in% ids_i)) |>
                unlist() |>
                which()
    }})
    
    ## finally, obtain the RHEA ids and return the object
    lapply(seq_along(ids), function(i) {
        if (type == "substrates")
            inds <- substrates_inds[[i]]
        if (type == "products")
            inds <- products_inds[[i]]
        if (type == "both")
            inds <- c(substrates_inds[[i]], products_inds[[i]])
        
        inds <- unique(inds)
        
        ## get the RHEA ids from inds
        ids <- rhea_reaction_smiles[inds, 1]
        
        ## set ids to NA if no corresponding RHEA ids are found and return
        if (length(ids) == 0)
            ids <- NA
        ids
    })
}


#' @name find_RHEA_ids_from_ids
#' 
#' @title Find RHEA ids from ids
#' 
#' @description
#' The function \code{find_RHEA_ids_from_ids} returns the RHEA ids
#' that are involved in reactions of the given identifiers (\code{ids})
#' of type \code{id_type}. 
#'  
#' @details
#' If one or several arguments (\code{annotation_set}, 
#' \code{rhea_chebi_smiles}, \code{rhea_reaction_smiles}, 
#' \code{rhea_directions}) are set to NULL, the objects will be loaded
#' as by default (the default values are identical to the values in the
#' Example section below).
#' 
#' The returned result might differ depending on the \code{type} object:
#' 
#' \itemize{
#'      \item if \code{type = "both"}, the function will return the RHEA 
#'     identifiers for the SMILES found both in the substrates or products of 
#'     the reactions,
#'     \item if \code{type = "substrates"}, the function will return the RHEA 
#'     identifiers for the SMILES found in the substrates of the reactions,
#'     \item if \code{type = "products"}, the function will return the RHEA 
#'     identifiers for the SMILES found in the products of the reactions.
#' }
#' 
#' @param ids character, identifiers to be translated to ChEBI ids
#' @param id_type character(1), one of \code{colnames(annotation_set)}
#' @param annotation_set data.frame containing the mappings between different
#' ids
#' @param rhea_chebi_smiles data.frame containing the mappings between ChEBI 
#' identifiers and SMILES identifiers 
#' @param rhea_reaction_smiles data.frame containing the mappings between ChEBI 
#' identifiers and SMILES identifiers
#' @param rhea_directions data.frame containing information on the 
#' directionality of RHEA directions
#' @param type character(1), one of "both", "substrates", or "products"
#' 
#' @author Thomas Naake
#' 
#' @export
#' 
#' @return list
#' 
#' @import metaboliteIDmapping
#' @importFrom AnnotationHub AnnotationHub query
#' 
#' @examples
#' ids <- c("HMDB0004947", "HMDB0004949", "HMDB0011763", "HMDB0004974", "HMDB0011594")
#' id_type <- "HMDB"
#' 
#' ## load annotation_set
#' library(metaboliteIDmapping)
#' library(AnnotationHub)
#' ah <- AnnotationHub()
#' datasets <- query(ah, "metaboliteIDmapping")
#' 
#' ## Currently, there are two versions of the mapping table.
#' ## AH79817 represents the original ID mapping containing 9 different ID formats
#' ## AH83115 mapping table which also includes common names for each compound
#' ## AH91792 current version of the mapping table that also accounts for tautomers
#' ## For implementing this data in your code, it is recommended to use the 
#' ## AHid for retrieval:
#' annotation_set <- ah[["AH91792"]]
#' 
#' ## load rhea_chebi_smiles
#' file <- file.path(
#'     path.package("LipidNetworkPredictR"), "extdata", "rhea-chebi-smiles.tsv",
#'     fsep = .Platform$file.sep)
#' rhea_chebi_smiles <- read.csv(file, sep = "\t", header = FALSE)
#' 
#' ## load rhea_reaction_smiles
#' file <- file.path(
#'     path.package("LipidNetworkPredictR"), "extdata", "rhea-reaction-smiles.tsv",
#'     fsep = .Platform$file.sep)
#' rhea_reaction_smiles <- read.csv(file, sep = "\t", header = FALSE)
#' 
#' ## load rhea_directions
#' file <- file.path(
#'     path.package("LipidNetworkPredictR"), "extdata", "rhea-directions.tsv",
#'     fsep = .Platform$file.sep)
#' rhea_directions <- read.csv(file, sep = "\t", header = TRUE)
#' 
#' ## run the function
#' find_RHEA_ids_from_ids(ids = ids, id_type = id_type, 
#'     annotation_set = annotation_set, rhea_chebi_smiles = rhea_chebi_smiles,
#'     rhea_reaction_smiles = rhea_reaction_smiles, 
#'     rhea_directions = rhea_directions, type = "substrates")
find_RHEA_ids_from_ids <- function(ids, id_type, annotation_set = NULL, 
    rhea_chebi_smiles = NULL, rhea_reaction_smiles = NULL, 
    rhea_directions = NULL, type = c("substrates", "products", "both")) {

    ## load the databases in case the objects are not given
    if (is.null(annotation_set)) {
   
        ## search for the mapping table in the AnnotationHub resource interface,
        ## metaboliteIDmapping comes from the package metaboliteIDmapping
        ah <- AnnotationHub()
        datasets <- query(ah, "metaboliteIDmapping")
        
        ## Currently, there are two versions of the mapping table.
        ## AH79817 represents the original ID mapping containing 9 different ID formats
        ## AH83115 mapping table which also includes common names for each compound
        ## AH91792 current version of the mapping table that also accounts for tautomers
        ## For implementing this data in your code, it is recommended to use the 
        ## AHid for retrieval:
        annotation_set <- ah[["AH91792"]]
    }

    if (is.null(rhea_chebi_smiles)) {
        file <- file.path(
            path.package("LipidNetworkPredictR"), "extdata", "rhea-chebi-smiles.tsv",
            fsep = .Platform$file.sep)
        rhea_chebi_smiles <- read.csv(file, sep = "\t", header = FALSE)
    }

    if (is.null(rhea_reaction_smiles)) {
        file <- file.path(
            path.package("LipidNetworkPredictR"), "extdata", "rhea-reaction-smiles.tsv",
            fsep = .Platform$file.sep)
        rhea_reaction_smiles <- read.csv(file, sep = "\t", header = FALSE)
    }

    if (is.null(rhea_directions)) {
        file <- file.path(
            path.package("LipidNetworkPredictR"), "extdata", "rhea-directions.tsv",
            fsep = .Platform$file.sep)
        rhea_directions <- read.csv(file, sep = "\t", header = TRUE)
    }

    ## check the type argument
    type <- match.arg(type)

    ## translate the ids to ChEBI ids
    chebi_ids <- translate_id_to_ChEBI(ids = ids, id_type = id_type, 
        annotation_set = annotation_set)

    ## translate the ChEBI ids to SMILES ids
    smiles_ids <- translate_ChEBI_to_SMILES(ids = chebi_ids, 
        rhea_chebi_smiles = rhea_chebi_smiles)

    ## find the RHEA ids from the SMILES ids
    rhea_ids <- find_RHEA_ids_from_SMILES(ids = smiles_ids, 
        rhea_reaction_smiles = rhea_reaction_smiles, 
        rhea_directions = rhea_directions, type = type) 

    ## return the rhea_ids object
    rhea_ids
}
