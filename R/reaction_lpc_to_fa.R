#'
#'
#' @export
sn1lpc_to_fa <- function(sn1LPCs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h2o + M_ag3pc <=> M_g3pc + M_h + M_fatacid",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  sn1LPCs <- sn1LPCs[stringr::str_detect(sn1LPCs, constraints[1], negate = negate[1])]

  # "sn2 loss"
  FAs <- unlist(lapply(sn1LPCs, function(sn1lpc) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(sn1lpc)
    fa <- paste0("FA(", fatty_acyls[1], ")")
    return(fa)
  }))
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(sn1LPCs)),
                            ReactionFormula = character(length(sn1LPCs)),
                            IsReversible = character(length(sn1LPCs)),
                            GeneAssociation = character(length(sn1LPCs)),
                            Pathway = character(length(sn1LPCs)))
  
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
                                                          "M_g3pc",
                                                          "Glycerophosphocholine")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_ag3pc",
                                                          sn1LPCs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fatacid",
                                                          FAs)
  
  # return results
  return(list(unique(FAs), reaction_df))
}

#'
#'
#' @export
sn2lpc_to_fa <- function(sn2LPCs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h2o +  M_2agpc <=> M_g3pc + M_h + M_fatacid",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  sn2LPCs <- sn2LPCs[stringr::str_detect(sn2LPCs, constraints[1], negate = negate[1])]
  
  # "sn2 loss"
  FAs <- unlist(lapply(sn2LPCs, function(sn2lpc) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(sn2lpc)
    fa <- paste0("FA(", fatty_acyls[2], ")")
    return(fa)
  }))
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(sn2LPCs)),
                            ReactionFormula = character(length(sn2LPCs)),
                            IsReversible = character(length(sn2LPCs)),
                            GeneAssociation = character(length(sn2LPCs)),
                            Pathway = character(length(sn2LPCs)))
  
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
                                                          "M_g3pc",
                                                          "Glycerophosphocholine")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_2agpc",
                                                          sn2LPCs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fatacid",
                                                          FAs)
  
  # return results
  return(list(unique(FAs), reaction_df))
}