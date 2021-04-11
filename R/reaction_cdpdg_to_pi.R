#'
#'
#' @export
cdpdg_to_pi <- function(CDPDGs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_inost + M_cdpdag <=> M_h + M_cmp + M_pail",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  CDPDGs <- CDPDGs[stringr::str_detect(CDPDGs, constraints, negate = negate)]
  
  # make new list for LPAs
  PIs <- stringr::str_replace_all(CDPDGs, "CDP-DG", "PI")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(PIs)),
                            ReactionFormula = character(length(PIs)),
                            IsReversible = character(length(PIs)),
                            GeneAssociation = character(length(PIs)),
                            Pathway = character(length(PIs)))
  
  # fill with data
  reaction_df$Name <- template$reaction_name
  reaction_df$ReactionFormula <- template$reaction_formula
  reaction_df$IsReversible <- template$reaction_isReversible
  reaction_df$GeneAssociation <- template$reaction_geneAssociation
  reaction_df$Pathway <- template$reaction_pathway
  
  # adjust fixed substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_inost",
                                                          "Inositol")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_h",
                                                          "H+")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_cmp",
                                                          "CMP")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_cdpdag",
                                                          CDPDGs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_pail",
                                                          PIs)
  
  # return results
  return(list(unique(PIs), reaction_df))
  
}