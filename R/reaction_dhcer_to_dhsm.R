#'
#'
#' @export
dhcer_to_dhsm <- function(DHCERs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_pchol + M_c17isodhcrm <=> M_12dag + M_c17isodhsphmyln",
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
  DHSMs <- stringr::str_replace_all(DHCERs, "Cer", "SM")
  
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
                                                          "M_pchol",
                                                          "PC")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_12dag",
                                                          "DG")

  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_c17isodhcrm",
                                                          DHCERs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_c17isodhsphmyln",
                                                          DHSMs)
  
  # return results
  return(list(unique(DHSMs), reaction_df))
  
}