#' @name create_via_chebi 
#' 
#' @title create_via_chebi 
#' 
#' @description
#' Creates ...
#' 
#' @details
#' Uses ...
#' 
#' @param rhea_reactions_order data.frame
#' @param rhea_reactions data.frame
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom rlang sym
#'    
#' @examples
#' ## load rhea_reactions_order
#' file <- file.path(
#'     path.package("LipidNetworkPredictR"), "extdata", "rhea_reactions.tsv",
#'     fsep = .Platform$file.sep)
#' rhea_reactions_order <- read.table(file, sep = "\t", header = TRUE)
#' rhea_reactions_order <- rhea_reactions_order[c(1, 5, 9), ]
#' 
#' ## load rhea_reactions
#' file <- file.path(
#'     path.package("LipidNetworkPredictR"), "extdata", "rhea-reactions.txt",
#'     fsep = .Platform$file.sep)
#' rhea_reactions <- read.table(file, sep = "\t", header = FALSE, quote = "")
#'
#' ## run function
#' create_via_chebi(rhea_reactions_order = rhea_reactions_order, 
#'     rhea_reactions = rhea_reactions)
create_via_chebi <- function(rhea_reactions_order, rhea_reactions) {
    
    rhea_pasted <- paste(rhea_reactions_order[, "rhea"], collapse = "|")
    
    ## only keep the entries of interest
    rhea_reactions <- rhea_reactions |>
        dplyr::filter(grepl(x = !!rlang::sym("V1"), pattern = "ENTRY|DEFINITION|EQUATION"))
    
    ## split the dataset into columns
    rhea_reactions_l <- strsplit(rhea_reactions[, 1], split = "       |    |  ")
    rhea_reactions_cut <- do.call("rbind", rhea_reactions_l)
    colnames(rhea_reactions_cut) <- c("category", "value")
    
    ## 
    inds_rhea <- match(
        rhea_reactions_order[["rhea"]], rhea_reactions_cut[, "value"])
    rhea_reactions_cut[sort(c(inds_rhea, inds_rhea + 1, inds_rhea + 2)), ]
    
    ## create matching table: CHEBI id --> name
    inds_rhea <- grep(x = rhea_reactions_cut[, "value"], pattern = "RHEA:")
    split_pattern <- " = | <=> | => | <= | [+] "
    reactants_l <- strsplit(rhea_reactions_cut[inds_rhea + 1, "value"], 
        split = split_pattern)
    chebi_l <- strsplit(rhea_reactions_cut[inds_rhea + 2, "value"], 
        split = split_pattern)
    
    reactants_l <- lapply(seq_along(chebi_l), function(i) {
        reactants_i_l <- reactants_l[[i]]
        names(reactants_i_l) <- chebi_l[[i]]
        reactants_i_l
    })
    reactants_l <- unlist(reactants_l)
    reactants_l <- reactants_l[!duplicated(reactants_l)]
    
    chebi_name_df <- data.frame(chebi = names(reactants_l), 
        metabolite = as.character(reactants_l))
}

# inds <- grep(x = rhea_reactions_cut[, "value"], "RHEA:44604") ## -> create_template.R
# rhea_reactions_cut[c(inds, inds + 1, inds + 2), ]
# 
# ## matching table: CHEBI --> shorthand
# 
# matching_table <- rbind(
#     c("CHEBI:15850", "1-MG-O", "1-O-alkyl-sn-glycerol"),
#     c("CHEBI:17135", "FAO", "a long chain fatty alcohol"),
#     c("CHEBI:17389", "2-MG", "a 2-acylglycerol"),
#     c("CHEBI:17636", "SM", "a sphingomyelin"),
#     c("CHEBI:17815", "1,2-DG", "a 1,2-diacyl-sn-glycerol"),
#     c("CHEBI:22801", "GlcCer", "a beta-D-glucosyl-(1<->1')-N-acylsphing-4-enine"),
#     c("CHEBI:28868", "FA", "a fatty acid"),
#     c("CHEBI:29067", "FA", "a carboxylate"),
#     c("CHEBI:30909", "1-LPC-O", "1-O-alkyl-sn-glycero-3-phosphocholine"),
#     c("CHEBI:35759", "1-MG", "a 1-acylglycerol"),
#     c("CHEBI:36702", "PC-O", "1-O-alkyl-2-acyl-sn-glycero-3-phosphocholine"),
#     c("CHEBI:49172", "1,2-DG", "a 1,2-diacylglycerol"),
#     c("CHEBI:52595", "DG-O", "1-O-alkyl-2-acyl-sn-glycerol"),
#     c("CHEBI:52639", "Cer", "an N-acylsphing-4-enine"),
#     c("CHEBI:57262", "PS", "a 1,2-diacyl-sn-glycero-3-phospho-L-serine"),
#     c("CHEBI:57534", "AcylDHAP", "a 1-acylglycerone 3-phosphate"),
#     c("CHEBI:57560", "FA", "a long-chain fatty acid"),
#     c("CHEBI:57597", "AcylCoA", "sn-glycerol 3-phosphate"),
#     c("CHEBI:57643", "PC", "a 1,2-diacyl-sn-glycero-3-phosphocholine"),
#     c("CHEBI:57674", "CerP", "an N-acylsphing-4-enine 1-phosphate"),
#     c("CHEBI:57875", "2-LPC", "a 2-acyl-sn-glycero-3-phosphocholine"),
#     c("CHEBI:57880", "PI", "a 1,2-diacyl-sn-glycero-3-phospho-(1D-myo-inositol)"),
#     c("CHEBI:57970", "LPA", "a 1-acyl-sn-glycero-3-phosphate"),
#     c("CHEBI:58014", "LPA-O", "1-O-alkyl-sn-glycero-3-phosphate"),
#     c("CHEBI:58168", "1-LPC", "a 1-acyl-sn-glycero-3-phosphocholine"),
#     c("CHEBI:58332", "CDP-DG", "a CDP-1,2-diacyl-sn-glycerol"),
#     c("CHEBI:58342", "AcylCoA", "an acyl-CoA"),
#     c("CHEBI:58608", "PA", "a 1,2-diacyl-sn-glycero-3-phosphate"),
#     c("CHEBI:60520", "PE-O", "1-O-alkyl-2-acyl-sn-glycero-3-phosphoethanolamine"),
#     c("CHEBI:60110", "PGP", "1,2-diacyl-sn-glycero-3-phospho-(1'-sn-glycero-3'-phosphate"), 
#     c("CHEBI:62237", "CL", "a cardiolipin"),
#     c("CHEBI:64381", "1-LPE", "a 1-acyl-sn-glycero-3-phosphoethanolamine"),
#     c("CHEBI:64583", "DhSM", "an N-(acyl)-sphingosylphosphocholine"),
#     c("CHEBI:64612", "PE", "a 1,2-diacyl-sn-glycero-3-phosphoethanolamine"),
#     c("CHEBI:64615", "TG", "a triacyl-sn-glycerol"),
#     c("CHEBI:64683", "1-MG", "a 1-acyl-sn-glycerol"),
#     c("CHEBI:64716", "PG", "1,2-diacyl-sn-glycero-3-phospho-(1'-sn-glycerol)"),
#     c("CHEBI:64743", "1,2,4-LCL", "1'-[1,2-diacyl-sn-glycero-3-phospho],3'-[1-acyl-sn-glycero-3-phospho]-glycerol"),
#     c("CHEBI:64771", "1-LPI", "a 1-acyl-sn-glycero-3-phospho-(1D-myo-inositol)"),
#     c("CHEBI:73315", "AlkylDHAP", "1-O-alkylglycerone 3-phosphate"),
#     c("CHEBI:73332", "PA-O", "1-O-alkyl-2-acyl-sn-glycero-3-phosphate"),
#     c("CHEBI:73359", "FAL", "a 2,3-saturated aldehyde"),
#     c("CHEBI:75028", "PE-O", "1-(1,2-saturated alkyl)-2-acyl-sn-glycero-3-phosphoethanolamine"),
#     c("CHEBI:77283", "LPA-P", "1-O-(1Z-alkenyl)-sn-glycero-3-phosphate"),
#     c("CHEBI:77288", "LPE-P", "1-O-(1Z-alkenyl)-sn-glycero-3-phosphoethanolamine"),
#     c("CHEBI:77290", "PE-P", "1-O-(1Z-alkenyl)-2-acyl-sn-glycero-3-phosphoethanolamine"),
#     c("CHEBI:77297", "1-MG-P", "1-O-(1Z-alkenyl)-sn-glycerol"),
#     c("CHEBI:77396", "FAO", "a long-chain primary fatty alcohol"),
#     c("CHEBI:77636", "AcylCoA", "a fatty acyl-CoA"),
#     c("CHEBI:83139", "AcylCoA", "a long-chain fatty acyl-CoA"),
#     c("CHEBI:83273", "DhCer", "an N-acyl-sphingoid base"),
#     c("CHEBI:85216", "LNAPE", "N,1-diacyl-sn-glycero-3-phosphoethanolamine"),
#     c("CHEBI:85225", "GPNAE", "N-acyl-sn-glycero-3-phosphoethanolamine")
# )
