#'
#'
#' @export
pao_to_dgo <- function(PAOs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h2o + M_akac2gp <=> M_pi + M_akac2g",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  PAOs <- PAOs[stringr::str_detect(PAOs, constraints, negate = negate)]
  
  # make new list for LPAOs
  DGOs <- stringr::str_replace_all(PAOs, "PA", "DG")
  DGOs <- stringr::str_replace_all(DGOs, "\\)$", "/0:0\\)")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(DGOs)),
                            ReactionFormula = character(length(DGOs)),
                            IsReversible = character(length(DGOs)),
                            GeneAssociation = character(length(DGOs)),
                            Pathway = character(length(DGOs)))
  
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
                                                          "M_pi",
                                                          "Pi")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_akac2gp",
                                                          PAOs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_akac2g",
                                                          DGOs)
  
  # return results
  return(list(unique(DGOs), reaction_df))
  
}