#'
#'
#' @export
pc_to_sn1lpc <- function(PCs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h2o +  M_pchol <=> M_ag3pc + M_fatacid",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  PCs <- PCs[stringr::str_detect(PCs, constraints[1], negate = negate[1])]
  
  # make DF
  sn1LPCs <- data.frame(PC = character(length(PCs)),
                        sn1LPC = character(length(PCs)),
                        FA = character(length(PCs)))
  
  sn1LPCs$PC <- PCs
  
  # "sn2 loss"
  sn1LPCs$sn1LPC <- unlist(lapply(sn1LPCs$PC, function(pc) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(pc)
    sn1lpc <- paste0("PC(", fatty_acyls[1], "/0:0)")
    return(sn1lpc)
  }))
  
  sn1LPCs$FA <- unlist(lapply(sn1LPCs$PC, function(pc) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(pc)
    fa <- paste0("FA(", fatty_acyls[2], ")")
    return(fa)
  }))
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(PCs)),
                            ReactionFormula = character(length(PCs)),
                            IsReversible = character(length(PCs)),
                            GeneAssociation = character(length(PCs)),
                            Pathway = character(length(PCs)))
  
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

  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_pchol",
                                                          sn1LPCs$PC)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_ag3pc",
                                                          sn1LPCs$sn1LPC)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fatacid",
                                                          sn1LPCs$FA)
  
  # return results
  return(list(unique(sn1LPCs$sn1LPC), reaction_df))
}

#'
#'
#' @export
pc_to_sn2lpc <- function(PCs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h2o + M_pchol <=> M_2agpc + M_fatacid",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  PCs <- PCs[stringr::str_detect(PCs, constraints[1], negate = negate[1])]
  
  # make DF
  sn2LPCs <- data.frame(PC = character(length(PCs)),
                        sn2LPC = character(length(PCs)),
                        FA = character(length(PCs)))
  
  sn2LPCs$PC <- PCs
  
  # "sn2 loss"
  sn2LPCs$sn2LPC <- unlist(lapply(sn2LPCs$PC, function(pc) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(pc)
    sn2lpc <- paste0("PC(0:0/", fatty_acyls[2], ")")
    return(sn2lpc)
  }))
  
  sn2LPCs$FA <- unlist(lapply(sn2LPCs$PC, function(pc) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(pc)
    fa <- paste0("FA(", fatty_acyls[1], ")")
    return(fa)
  }))
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(PCs)),
                            ReactionFormula = character(length(PCs)),
                            IsReversible = character(length(PCs)),
                            GeneAssociation = character(length(PCs)),
                            Pathway = character(length(PCs)))
  
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
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_pchol",
                                                          sn2LPCs$PC)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_2agpc",
                                                          sn2LPCs$sn2LPC)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fatacid",
                                                          sn2LPCs$FA)
  
  # return results
  return(list(unique(sn2LPCs$sn2LPC), reaction_df))
}
