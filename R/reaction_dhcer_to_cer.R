#'
#'
#' @export
dhcer_to_cer <- function(DHCERs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h + M_nadh + M_o2 + M_c17isodhcrm <=> 2.0 M_h2o + M_nad + M_c17isocrm",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  DHCERs <- DHCERs[stringr::str_detect(DHCERs, constraints, negate = negate)]
  
  # make new list for LDHCERs
  CERs <- stringr::str_replace_all(DHCERs, "Cer\\(d16:0\\(3OH,4OH\\)\\(15Me\\)\\/", "Cer(d16:1(4E)(3OH,4OH)(15Me)/")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(CERs)),
                            ReactionFormula = character(length(CERs)),
                            IsReversible = character(length(CERs)),
                            GeneAssociation = character(length(CERs)),
                            Pathway = character(length(CERs)))
  
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
                                                          "M_h",
                                                          "H+")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_nadh",
                                                          "NADH")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_o2",
                                                          "O2")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_nad",
                                                          "NAD+")
  
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_c17isodhcrm",
                                                          DHCERs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_c17isocrm",
                                                          CERs)
  
  # return results
  return(list(unique(CERs), reaction_df))
  
}