#'
#'
#' @export
cer_to_c1p <- function(CERs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_atp + M_c17isocrm <=> M_h + M_adp + M_c17isocrmp",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  CERs <- CERs[stringr::str_detect(CERs, constraints, negate = negate)]
  
  # make new list for LCERs
  C1Ps <- stringr::str_replace_all(CERs, "Cer", "C1P")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(C1Ps)),
                            ReactionFormula = character(length(C1Ps)),
                            IsReversible = character(length(C1Ps)),
                            GeneAssociation = character(length(C1Ps)),
                            Pathway = character(length(C1Ps)))
  
  # fill with data
  reaction_df$Name <- template$reaction_name
  reaction_df$ReactionFormula <- template$reaction_formula
  reaction_df$IsReversible <- template$reaction_isReversible
  reaction_df$GeneAssociation <- template$reaction_geneAssociation
  reaction_df$Pathway <- template$reaction_pathway
  
  # adjust fixed substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_atp",
                                                          "ATP")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_h",
                                                          "H+")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_adp",
                                                          "ADP")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_c17isocrm",
                                                          CERs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_c17isocrmp",
                                                          C1Ps)
  
  # return results
  return(list(unique(C1Ps), reaction_df))
  
}