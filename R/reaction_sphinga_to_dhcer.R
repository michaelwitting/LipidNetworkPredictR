#'
#'
#' @export
sphinga_to_dhcer <- function(CoAs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_fataccoa + M_c17isosphgn <=> M_coa + M_c17isodhcrm",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  CoAs <- CoAs[stringr::str_detect(CoAs, constraints, negate = negate)]
  
  # make new list for DHCERs
  DHCERs <- stringr::str_replace_all(CoAs, "CoA\\(", "Cer(d16:0(3OH,4OH)(15Me)/")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(CoAs)),
                            ReactionFormula = character(length(CoAs)),
                            IsReversible = character(length(CoAs)),
                            GeneAssociation = character(length(CoAs)),
                            Pathway = character(length(CoAs)))
  
  # fill with data
  reaction_df$Name <- template$reaction_name
  reaction_df$ReactionFormula <- template$reaction_formula
  reaction_df$IsReversible <- template$reaction_isReversible
  reaction_df$GeneAssociation <- template$reaction_geneAssociation
  reaction_df$Pathway <- template$reaction_pathway
  
  # adjust fixed substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_c17isosphgn",
                                                          "SPH(d16:0(1OH,3OH)(15Me))")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_coa",
                                                          "CoA")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fataccoa",
                                                          CoAs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_alpa_pl",
                                                          DHCERs)
  
  # return results
  return(list(unique(DHCERs), reaction_df))
  
}