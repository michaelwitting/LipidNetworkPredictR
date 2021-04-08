#'
#'
#' @export
pe_to_sn1lpe <- function(PEs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h2o + M_pe <=> M_h + M_fatacid + M_ag3pe",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  PEs <- PEs[stringr::str_detect(PEs, constraints[1], negate = negate[1])]
  
  # make DF
  sn1LPEs <- data.frame(PE = character(length(PEs)),
                        sn1LPE = character(length(PEs)),
                        FA = character(length(PEs)))
  
  sn1LPEs$PE <- PEs
  
  # "sn2 loss"
  sn1LPEs$sn1LPE <- unlist(lapply(sn1LPEs$PE, function(PE) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(PE)
    sn1LPE <- paste0("PE(", fatty_acyls[1], "/0:0)")
    return(sn1LPE)
  }))
  
  sn1LPEs$FA <- unlist(lapply(sn1LPEs$PE, function(PE) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(PE)
    fa <- paste0("FA(", fatty_acyls[2], ")")
    return(fa)
  }))
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(PEs)),
                            ReactionFormula = character(length(PEs)),
                            IsReversible = character(length(PEs)),
                            GeneAssociation = character(length(PEs)),
                            Pathway = character(length(PEs)))
  
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
                                                          "M_pe",
                                                          sn1LPEs$PE)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_ag3pe",
                                                          sn1LPEs$sn1LPE)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fatacid",
                                                          sn1LPEs$FA)
  
  # return results
  return(list(unique(sn1LPEs$sn1LPE), reaction_df))
}

#'
#'
#' @export
pe_to_sn2lpe <- function(PEs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h2o + M_pe <=> M_h + M_fatacid + M_2agpe",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  PEs <- PEs[stringr::str_detect(PEs, constraints[1], negate = negate[1])]
  
  # make DF
  sn2LPEs <- data.frame(PE = character(length(PEs)),
                        sn1LPE = character(length(PEs)),
                        FA = character(length(PEs)))
  
  sn2LPEs$PE <- PEs
  
  # "sn2 loss"
  sn2LPEs$sn2LPE <- unlist(lapply(sn2LPEs$PE, function(PE) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(PE)
    sn2LPE <- paste0("PE(0:0/", fatty_acyls[1], ")")
    return(sn2LPE)
  }))
  
  sn2LPEs$FA <- unlist(lapply(sn2LPEs$PE, function(PE) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(PE)
    fa <- paste0("FA(", fatty_acyls[2], ")")
    return(fa)
  }))
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(PEs)),
                            ReactionFormula = character(length(PEs)),
                            IsReversible = character(length(PEs)),
                            GeneAssociation = character(length(PEs)),
                            Pathway = character(length(PEs)))
  
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
                                                          "M_pe",
                                                          sn2LPEs$PE)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_2agpe",
                                                          sn2LPEs$sn2LPE)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fatacid",
                                                          sn2LPEs$FA)
  
  # return results
  return(list(unique(sn2LPEs$sn2LPE), reaction_df))
}