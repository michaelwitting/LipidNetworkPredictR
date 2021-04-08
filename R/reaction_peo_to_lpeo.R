#'
#'
#' @export
peo_to_lpeo <- function(PEOs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h2o + M_akac2gpe <=> M_h + M_ak2lgpe + M_fatacid",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  PEOs <- PEOs[stringr::str_detect(PEOs, constraints[1], negate = negate[1])]
  
  # make DF
  LPEOs <- data.frame(PEO = character(length(PEOs)),
                      LPEO = character(length(PEOs)),
                      FA = character(length(PEOs)))
  
  LPEOs$PEO <- PEOs
  
  # "sn2 loss"
  LPEOs$LPEO <- unlist(lapply(LPEOs$PEO, function(peo) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(peo)
    lpeo <- paste0("PE(", fatty_acyls[1], "/0:0)")
    return(lpeo)
  }))
  
  LPEOs$FA <- unlist(lapply(LPEOs$PEO, function(peo) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(peo)
    fa <- paste0("FA(", fatty_acyls[2], ")")
    return(fa)
  }))
  
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
                                                          "M_h2o",
                                                          "H2O")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_h",
                                                          "H+")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_akac2gpe",
                                                          LPEOs$PEO)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_ak2lgpe",
                                                          LPEOs$LPEO)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fatacid",
                                                          LPEOs$FA)
  
  # return results
  return(list(unique(LPEOs$LPEO), reaction_df))
}