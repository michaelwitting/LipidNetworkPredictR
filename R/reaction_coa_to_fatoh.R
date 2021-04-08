#'
#'
#' @export
coa_to_fatoh <- function(CoAs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_fataccoa + 2 M_napdh + 2 M_h <=> M_alkylR1oh + 2 M_nadp + M_coa",
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
  
  # make new list for FATOHs
  FATOHs <- stringr::str_replace_all(CoAs, "CoA", "FATOH")
  
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
                                                          "M_napdh",
                                                          "NADPH")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_h",
                                                          "H+")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_nadp",
                                                          "NADP+")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_coa",
                                                          "CoA")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fataccoa",
                                                          CoAs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_alkylR1oh",
                                                          FATOHs)
  
  # return results
  return(list(unique(FATOHs), reaction_df))
  
}