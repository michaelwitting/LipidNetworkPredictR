#'
#'
#' @export
dgo_to_peo <- function(DGOs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_cdpea + M_akac2g <=> M_h + M_cmp + M_akac2gpe",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  DGOs <- DGOs[stringr::str_detect(DGOs, constraints, negate = negate)]
  
  # make new list for LPAs
  PEOs <- stringr::str_replace_all(DGOs, "DG", "PE")
  PEOs <- stringr::str_replace_all(PEOs, "/0:0\\)$", "\\)")
  
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
                                                          "M_cdpea",
                                                          "CDP-Ethn")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_h",
                                                          "H+")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_cmp",
                                                          "CMP")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_akac2g",
                                                          DGOs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_akac2gpe",
                                                          PEOs)
  
  # return results
  return(list(unique(PEOs), reaction_df))
  
}