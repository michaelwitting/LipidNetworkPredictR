#'
#'
#' @export
dhsm_to_dhcer <- function(DHSMs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h2o + M_c17isodhsphmyln <=> M_cholp + M_h + M_c17isodhcrm",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  DHSMs <- DHSMs[stringr::str_detect(DHSMs, constraints, negate = negate)]
  
  # make new list for LDHCERs
  DHCERs <- stringr::str_replace_all(DHSMs, "SM", "Cer")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(DHSMs)),
                            ReactionFormula = character(length(DHSMs)),
                            IsReversible = character(length(DHSMs)),
                            GeneAssociation = character(length(DHSMs)),
                            Pathway = character(length(DHSMs)))
  
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
                                                          "M_c17isodhcrm",
                                                          DHCERs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_c17isodhsphmyln",
                                                          DHSMs)
  
  # return results
  return(list(unique(DHCERs), reaction_df))
  
}