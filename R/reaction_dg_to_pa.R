#'
#'
#' @export
dg_to_pa <- function(DGs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_atp + M_12dag <=> M_h + M_adp + M_pa_pl",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  DGs <- DGs[stringr::str_detect(DGs, constraints, negate = negate)]
  
  # make new list for LDGs
  PAs <- stringr::str_replace_all(DGs, "DG", "PA")
  PAs <- stringr::str_replace_all(PAs, "/0:0\\)$", "\\)")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(DGs)),
                            ReactionFormula = character(length(DGs)),
                            IsReversible = character(length(DGs)),
                            GeneAssociation = character(length(DGs)),
                            Pathway = character(length(DGs)))
  
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
                                                          "M_pa_pl",
                                                          PAs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_12dag",
                                                          DGs)
  
  # return results
  return(list(unique(PAs), reaction_df))
  
}