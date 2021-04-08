#'
#'
#' @export
alkyldhap_to_lpao <- function(ALKYLDHAPs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h + M_nadph + M_akdhap <=> M_alkgp + M_nadp",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  ALKYLDHAPs <- ALKYLDHAPs[stringr::str_detect(ALKYLDHAPs, constraints, negate = negate)]
  
  # make new list for LPAOs
  LPAOs <- stringr::str_replace_all(ALKYLDHAPs, "DHAP", "PA")
  LPAOs <- stringr::str_replace_all(LPAOs, "\\)$", "/0:0\\)")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(ALKYLDHAPs)),
                            ReactionFormula = character(length(ALKYLDHAPs)),
                            IsReversible = character(length(ALKYLDHAPs)),
                            GeneAssociation = character(length(ALKYLDHAPs)),
                            Pathway = character(length(ALKYLDHAPs)))
  
  # fill with data
  reaction_df$Name <- template$reaction_name
  reaction_df$ReactionFormula <- template$reaction_formula
  reaction_df$IsReversible <- template$reaction_isReversible
  reaction_df$GeneAssociation <- template$reaction_geneAssociation
  reaction_df$Pathway <- template$reaction_pathway
  
  # adjust fixed substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_H",
                                                          "H+")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_nadph",
                                                          "NADPH")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_napd",
                                                          "NADP+")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_akdhap",
                                                          ALKYLDHAPs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_alpa_pl",
                                                          LPAOs)
  
  # return results
  return(list(unique(LPAOs), reaction_df))
  
}