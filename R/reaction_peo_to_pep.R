#'
#'
#' @export
peo_to_pep <- function(PEOs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_akac2gpe + M_nadph + M_h + M_o2 <=> M_alkenac2gpe + M_nadp + 2 M_h2o",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  PEOs <- PEOs[stringr::str_detect(PEOs, constraints, negate = negate)]
  
  # make new list for LPAs
  PEPs <- stringr::str_replace_all(PEOs, "PE\\(O-", "PE(P-")

  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(PEOs)),
                            ReactionFormula = character(length(PEOs)),
                            IsReversible = character(length(PEOs)),
                            GeneAssociation = character(length(PEOs)),
                            Pathway = character(length(PEOs)))
  
  # fill with data
  reaction_df$Name <- template$reaction_name
  reaction_df$ReactionFormula <- template$reaction_formula
  reaction_df$IsReversible <- template$reaction_isReversible
  reaction_df$GeneAssociation <- template$reaction_geneAssociation
  reaction_df$Pathway <- template$reaction_pathway
  
  # adjust fixed substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_nadph",
                                                          "NADPH")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_h",
                                                          "H+")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_o2",
                                                          "O2")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_nadp",
                                                          "NADP+")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_h2o",
                                                          "H2O")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_akac2gpe",
                                                          PEOs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_alkenac2gpe_c",
                                                          PEPs)
  
  # return results
  return(list(unique(PEPs), reaction_df))
  
}