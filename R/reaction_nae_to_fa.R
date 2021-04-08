#'
#'
#' @export
nae_to_fa <- function(NAEs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h2o + M_nae <=> M_etha + M_h + M_fatacid",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  NAEs <- NAEs[stringr::str_detect(NAEs, constraints[1], negate = negate[1])]
  
  # "sn2 loss"
  FAs <- stringr::str_replace_all(NAEs, "NAE", "FA")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(NAEs)),
                            ReactionFormula = character(length(NAEs)),
                            IsReversible = character(length(NAEs)),
                            GeneAssociation = character(length(NAEs)),
                            Pathway = character(length(NAEs)))
  
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
                                                          "M_h",
                                                          "H+")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_etha",
                                                          "Ethanolamine")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_nae",
                                                          NAEs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fatacid",
                                                          FAs)
  
  # return results
  return(list(unique(FAs), reaction_df))
}