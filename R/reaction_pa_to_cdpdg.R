#'
#'
#' @export
pa_to_cdpdg <- function(PAs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_ctp + M_pa_pl <=> M_ppi + M_cdpdag",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  PAs <- PAs[stringr::str_detect(PAs, constraints, negate = negate)]
  
  # make new list for LPAs
  CDPDGs <- stringr::str_replace_all(PAs, "PA", "CDP-DG")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(CDPDGs)),
                            ReactionFormula = character(length(CDPDGs)),
                            IsReversible = character(length(CDPDGs)),
                            GeneAssociation = character(length(CDPDGs)),
                            Pathway = character(length(CDPDGs)))
  
  # fill with data
  reaction_df$Name <- template$reaction_name
  reaction_df$ReactionFormula <- template$reaction_formula
  reaction_df$IsReversible <- template$reaction_isReversible
  reaction_df$GeneAssociation <- template$reaction_geneAssociation
  reaction_df$Pathway <- template$reaction_pathway
  
  # adjust fixed substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_ctp",
                                                          "CTP")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_ppi",
                                                          "PPi")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_pa_pl",
                                                          PAs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_cdpdag",
                                                          CDPDGs)
  
  # return results
  return(list(unique(CDPDGs), reaction_df))
  
}