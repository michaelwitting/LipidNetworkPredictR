#'
#'
#' @export
cdpdg_to_pgp <- function(CDPDGs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_glyc3p + M_cdpdag <=> M_h + M_cmp + M_pgp",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  CDPDGs <- CDPDGs[stringr::str_detect(CDPDGs, constraints, negate = negate)]
  
  # make new list for LPAs
  PGPs <- stringr::str_replace_all(CDPDGs, "CDPDG", "PGP")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(PGPs)),
                            ReactionFormula = character(length(PGPs)),
                            IsReversible = character(length(PGPs)),
                            GeneAssociation = character(length(PGPs)),
                            Pathway = character(length(PGPs)))
  
  # fill with data
  reaction_df$Name <- template$reaction_name
  reaction_df$ReactionFormula <- template$reaction_formula
  reaction_df$IsReversible <- template$reaction_isReversible
  reaction_df$GeneAssociation <- template$reaction_geneAssociation
  reaction_df$Pathway <- template$reaction_pathway
  
  # adjust fixed substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_glyc3p",
                                                          "Glycerol-3-P")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_h",
                                                          "H+")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_cmp",
                                                          "CMP")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_cdpdag",
                                                          CDPDGs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_pail",
                                                          PGPs)
  
  # return results
  return(list(unique(PGPs), reaction_df))
  
}