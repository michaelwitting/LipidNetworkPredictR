#'
#'
#' @export
pe_to_nape_sn1 <- function(PEs, PCs, template = NA, constraints = c("", ""), negate = c(FALSE, FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_pe + M_pc <=> M_nape + M_2agpc",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 2) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  PEs <- PEs[stringr::str_detect(PEs, constraints[1], negate = negate[1])]
  PCs <- PCs[stringr::str_detect(PCs, constraints[2], negate = negate[2])]
  
  # make combinatorics
  NAPEs <- data.frame(expand.grid(PEs = PEs,
                                  PCs = PCs, stringsAsFactors = FALSE),
                      stringsAsFactors = FALSE)
  
  NAPEs$LPCs <- unlist(lapply(NAPEs$PCs,
                              function(x) {
                                paste0("PC(0:0/", lipidomicsUtils::isolate_fatty_acyls(x)[2], ")")
                              }))
  
  
  NAPEs$NAPEs <- stringr::str_replace_all(NAPEs$PEs,
                                          "\\)$",
                                          unlist(lapply(NAPEs$PCs,
                                                        function(x) {
                                                          paste0("/",lipidomicsUtils::isolate_fatty_acyls(x)[1], ")")
                                                        }
                                          )))
  
  
  NAPEs$NAPEs <- stringr::str_replace_all(NAPEs$NAPEs, "PE", "NAPE")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(nrow(NAPEs)),
                            ReactionFormula = character(nrow(NAPEs)),
                            IsReversible = character(nrow(NAPEs)),
                            GeneAssociation = character(nrow(NAPEs)),
                            Pathway = character(nrow(NAPEs)))
  
  # fill with data
  reaction_df$Name <- template$reaction_name
  reaction_df$ReactionFormula <- template$reaction_formula
  reaction_df$IsReversible <- template$reaction_isReversible
  reaction_df$GeneAssociation <- template$reaction_geneAssociation
  reaction_df$Pathway <- template$reaction_pathway
  
  # adjust fixed substrates and products
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_pc",
                                                          NAPEs$PCs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_pe",
                                                          NAPEs$PEs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_nape",
                                                          NAPEs$NAPEs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_2agpc",
                                                          NAPEs$LPCs)
  
  # return results
  return(list(unique(NAPEs$NAPEs), unique(NAPEs$LPCs), reaction_df))
}


#'
#'
#' @export
pe_to_nape_sn2 <- function(PEs, PCs, template = NA, constraints = c("", ""), negate = c(FALSE, FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_pe + M_pc <=> M_nape + M_ag3pc",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 2) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  PEs <- PEs[stringr::str_detect(PEs, constraints[1], negate = negate[1])]
  PCs <- PCs[stringr::str_detect(PCs, constraints[2], negate = negate[2])]
  
  # make combinatorics
  NAPEs <- data.frame(expand.grid(PEs = PEs,
                                PCs = PCs, stringsAsFactors = FALSE),
                    stringsAsFactors = FALSE)
  
  NAPEs$LPCs <- unlist(lapply(NAPEs$PCs,
                       function(x) {
                         paste0("PC(", lipidomicsUtils::isolate_fatty_acyls(x)[1], "/0:0)")
                       }))
  
  
  NAPEs$NAPEs <- stringr::str_replace_all(NAPEs$PEs,
                                      "\\)$",
                                      unlist(lapply(NAPEs$PCs,
                                                    function(x) {
                                                      paste0("/",lipidomicsUtils::isolate_fatty_acyls(x)[2], ")")
                                                    }
                                      )))
  
  NAPEs$NAPEs <- stringr::str_replace_all(NAPEs$NAPEs, "PE", "NAPE")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(nrow(NAPEs)),
                            ReactionFormula = character(nrow(NAPEs)),
                            IsReversible = character(nrow(NAPEs)),
                            GeneAssociation = character(nrow(NAPEs)),
                            Pathway = character(nrow(NAPEs)))
  
  # fill with data
  reaction_df$Name <- template$reaction_name
  reaction_df$ReactionFormula <- template$reaction_formula
  reaction_df$IsReversible <- template$reaction_isReversible
  reaction_df$GeneAssociation <- template$reaction_geneAssociation
  reaction_df$Pathway <- template$reaction_pathway
  
  # adjust fixed substrates and products

  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_pc",
                                                          NAPEs$PCs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_pe",
                                                          NAPEs$PEs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_nape",
                                                          NAPEs$NAPEs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_ag3pc",
                                                          NAPEs$LPCs)
  
  # return results
  return(list(unique(NAPEs$NAPEs), unique(NAPEs$LPCs), reaction_df))
}
