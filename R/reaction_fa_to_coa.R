#'
#'
#' @export
fa_to_coa <- function(FAs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_atp + M_coa + M_fatacid <=> M_ppi + M_amp + M_fataccoa",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  

  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  FAs <- FAs[stringr::str_detect(FAs, constraints[1], negate = negate[1])]
  
  # make new list for coas
  CoAs <- stringr::str_replace_all(FAs, "FA", "CoA")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(FAs)),
                            ReactionFormula = character(length(FAs)),
                            IsReversible = character(length(FAs)),
                            GeneAssociation = character(length(FAs)),
                            Pathway = character(length(FAs)))
  
  # fill with data
  reaction_df$Name <- template$reaction_name
  reaction_df$ReactionFormula <- template$reaction_formula
  reaction_df$IsReversible <- template$reaction_isReversible
  reaction_df$GeneAssociation <- template$reaction_geneAssociation
  reaction_df$Pathway <- template$reaction_pathway
  
  # adjust fixed substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_atp",
                                                          "ATP")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_coa",
                                                          "CoA")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_ppi",
                                                          "PPi")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_amp",
                                                          "AMP")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fatacid",
                                                          FAs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fataccoa",
                                                          CoAs)
  
  # return results
  return(list(unique(CoAs), reaction_df))
}
