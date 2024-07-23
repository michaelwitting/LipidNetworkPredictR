## function create_via_chebi
test_that("create_via_chebi works", {
    
    file <- file.path(
        path.package("LipidNetworkPredictR"), "extdata", "rhea_reactions.tsv",
        fsep = .Platform$file.sep)
    rhea_reactions_order <- read.table(file, sep = "\t", header = TRUE)
    rhea_reactions_order <- rhea_reactions_order[c(1, 5, 9), ]

    ## load rhea_reactions
    file <- file.path(
        path.package("LipidNetworkPredictR"), "extdata", "rhea-reactions.txt",
        fsep = .Platform$file.sep)
    rhea_reactions <- read.table(file, sep = "\t", header = FALSE, quote = "")

    ## run function
    ##df <- create_via_chebi(rhea_reactions_order = rhea_reactions_order, 
    ##    rhea_reactions = rhea_reactions)
    
})
