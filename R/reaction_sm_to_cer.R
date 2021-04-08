#'
#'
#' @export
sm_to_cer <- function(SMs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h2o + M_c17isosphmyln <=> M_cholp + M_h + M_c17isocrm",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  SMs <- SMs[stringr::str_detect(SMs, constraints, negate = negate)]
  
  # make new list for LCERs
  CERs <- stringr::str_replace_all(SMs, "SM", "Cer")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(SMs)),
                            ReactionFormula = character(length(SMs)),
                            IsReversible = character(length(SMs)),
                            GeneAssociation = character(length(SMs)),
                            Pathway = character(length(SMs)))
  
  # fill with data
  reaction_df$Name <- template$reaction_name
  reaction_df$ReactionFormula <- template$reaction_formula
  reaction_df$IsReversible <- template$reaction_isReversible
  reaction_df$GeneAssociation <- template$reaction_geneAssociation
  reaction_df$Pathway <- template$reaction_pathway
  
  # adjust fixed substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_h2o",
                                                          "H2O")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_cholp",
                                                          "P-Choline")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_h",
                                                          "H+")
  
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_c17isocrm",
                                                          CERs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_c17isosphmyln",
                                                          SMs)
  
  # return results
  return(list(unique(CERs), reaction_df))
  
}