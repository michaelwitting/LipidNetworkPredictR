#'
#'
#' @export
pc_to_dg <- function(PCs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h2o + M_pchol <=> M_cholp + M_12dag",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  PCs <- PCs[stringr::str_detect(PCs, constraints, negate = negate)]
  
  # make new list for LPAs
  DGs <- stringr::str_replace_all(PCs, "PC", "DG")
  DGs <- stringr::str_replace(DGs, "\\)$", "/0:0\\)")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(PCs)),
                            ReactionFormula = character(length(PCs)),
                            IsReversible = character(length(PCs)),
                            GeneAssociation = character(length(PCs)),
                            Pathway = character(length(PCs)))
  
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
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_pchol",
                                                          PCs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_12dag",
                                                          DGs)
  
  # return results
  return(list(reaction_df))
  
}