#'
#'
#' @export
pco_to_lpco <- function(PCOs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h2o + M_akac2gchol <=> M_h + M_ak2lgchol + M_fatacid",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  PCOs <- PCOs[stringr::str_detect(PCOs, constraints[1], negate = negate[1])]
  
  # make DF
  LPCOs <- data.frame(PC = character(length(PCOs)),
                        LPCO = character(length(PCOs)),
                        FA = character(length(PCOs)))
  
  LPCOs$PCO <- PCOs
  
  # "sn2 loss"
  LPCOs$LPCO <- unlist(lapply(LPCOs$PCO, function(pco) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(pco)
    lpco <- paste0("PC(", fatty_acyls[1], "/0:0)")
    return(lpco)
  }))
  
  LPCOs$FA <- unlist(lapply(LPCOs$PCO, function(pco) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(pco)
    fa <- paste0("FA(", fatty_acyls[2], ")")
    return(fa)
  }))
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(PCOs)),
                            ReactionFormula = character(length(PCOs)),
                            IsReversible = character(length(PCOs)),
                            GeneAssociation = character(length(PCOs)),
                            Pathway = character(length(PCOs)))
  
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
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_akac2gchol",
                                                          LPCOs$PCO)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_ak2lgchol",
                                                          LPCOs$LPCO)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fatacid",
                                                          LPCOs$FA)
  
  # return results
  return(list(unique(LPCOs$LPCO), reaction_df))
}