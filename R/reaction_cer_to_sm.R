#'
#'
#' @export
cer_to_sm <- function(CERs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_pchol + M_c17isocrm <=> M_12dag + M_c17isosphmyln",
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
  SMs <- stringr::str_replace_all(CERs, "Cer", "SM")
  
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
                                                          "M_pchol",
                                                          "PC")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_12dag",
                                                          "DG")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_c17isocrm",
                                                          CERs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_c17isosphmyln",
                                                          SMs)
  
  # return results
  return(list(unique(SMs), reaction_df))
  
}