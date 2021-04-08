#'
#'
#' @export
sn1mg_to_lpa <- function(sn1MGs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_atp + M_1magol <=> M_h + M_adp + M_alpa_pl",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  sn1MGs <- sn1MGs[stringr::str_detect(sn1MGs, constraints, negate = negate)]
  
  # make new list for LPAs
  LPAs <- stringr::str_replace_all(sn1MGs, "MG", "PA")
  LPAs <- stringr::str_replace_all(LPAs, "/0:0\\)$", "\\)")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(sn1MGs)),
                            ReactionFormula = character(length(sn1MGs)),
                            IsReversible = character(length(sn1MGs)),
                            GeneAssociation = character(length(sn1MGs)),
                            Pathway = character(length(sn1MGs)))
  
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
                                                          "M_h",
                                                          "H+")
  
    reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_adp",
                                                          "ADP")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_1magol",
                                                          sn1MGs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_alpa_pl",
                                                          LPAs)
  
  # return results
  return(list(unique(LPAs), reaction_df))
  
}