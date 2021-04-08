#'
#'
#' @export
pe_to_dg <- function(PEs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h2o + M_pe <=> M_ethamp + M_12dag",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  PEs <- PEs[stringr::str_detect(PEs, constraints, negate = negate)]
  
  # make new list for LPAs
  DGs <- stringr::str_replace_all(PEs, "PE", "DG")
  DGs <- stringr::str_replace(DGs, "\\)$", "/0:0\\)")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(PEs)),
                            ReactionFormula = character(length(PEs)),
                            IsReversible = character(length(PEs)),
                            GeneAssociation = character(length(PEs)),
                            Pathway = character(length(PEs)))
  
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
                                                          "M_ethamp",
                                                          "P-Ethanolamine")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_pe",
                                                          PEs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_12dag",
                                                          DGs)
  
  # return results
  return(list(reaction_df))
  
}