#'
#'
#' @export
pep_to_lpep <- function(PEPs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h2o + M_alkenac2gpe <=> M_h + M_alken2gpe + M_fatacid",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  PEPs <- PEPs[stringr::str_detect(PEPs, constraints[1], negate = negate[1])]
  
  # make DF
  LPEPs <- data.frame(PEP = character(length(PEPs)),
                      LPEP = character(length(PEPs)),
                      FA = character(length(PEPs)))
  
  LPEPs$PEP <- PEPs
  
  # "sn2 loss"
  LPEPs$LPEP <- unlist(lapply(LPEPs$PEP, function(pep) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(pep)
    lPEP <- paste0("PE(", fatty_acyls[1], "/0:0)")
    return(lPEP)
  }))
  
  LPEPs$FA <- unlist(lapply(LPEPs$PEP, function(pep) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(pep)
    fa <- paste0("FA(", fatty_acyls[2], ")")
    return(fa)
  }))
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(PEPs)),
                            ReactionFormula = character(length(PEPs)),
                            IsReversible = character(length(PEPs)),
                            GeneAssociation = character(length(PEPs)),
                            Pathway = character(length(PEPs)))
  
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
                                                          "M_alkenac2gpe",
                                                          LPEPs$PEP)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_alken2gpe",
                                                          LPEPs$LPEP)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fatacid",
                                                          LPEPs$FA)
  
  # return results
  return(list(unique(LPEPs$LPEP), reaction_df))
}