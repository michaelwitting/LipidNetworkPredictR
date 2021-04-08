#'
#'
#' @export
ps_to_pe <- function(PSs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h + M_ps <=> M_co2 + M_pe",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  PSs <- PSs[stringr::str_detect(PSs, constraints, negate = negate)]
  
  # make new list for LPAs
  PEs <- stringr::str_replace_all(PSs, "PS", "PE")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(PSs)),
                            ReactionFormula = character(length(PSs)),
                            IsReversible = character(length(PSs)),
                            GeneAssociation = character(length(PSs)),
                            Pathway = character(length(PSs)))
  
  # fill with data
  reaction_df$Name <- template$reaction_name
  reaction_df$ReactionFormula <- template$reaction_formula
  reaction_df$IsReversible <- template$reaction_isReversible
  reaction_df$GeneAssociation <- template$reaction_geneAssociation
  reaction_df$Pathway <- template$reaction_pathway
  
  # adjust fixed substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_h",
                                                          "H+")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_co2",
                                                          "CO2")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_ps",
                                                          PSs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_pe",
                                                          PEs)
  
  # return results
  return(list(unique(PEs), reaction_df))
  
}