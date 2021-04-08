#'
#'
#' @export
sn1lpe_to_fa <- function(sn1LPEs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h2o + M_ag3pe <=> M_h + M_fatacid + M_g3pe",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  sn1LPEs <- sn1LPEs[stringr::str_detect(sn1LPEs, constraints[1], negate = negate[1])]
  
  # "sn2 loss"
  FAs <- unlist(lapply(sn1LPEs, function(sn1lpe) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(sn1lpe)
    fa <- paste0("FA(", fatty_acyls[1], ")")
    return(fa)
  }))
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(sn1LPEs)),
                            ReactionFormula = character(length(sn1LPEs)),
                            IsReversible = character(length(sn1LPEs)),
                            GeneAssociation = character(length(sn1LPEs)),
                            Pathway = character(length(sn1LPEs)))
  
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
                                                          "Glycerophosphoethanolamine")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_ag3pe",
                                                          sn1LPEs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fatacid",
                                                          FAs)
  
  # return results
  return(list(unique(FAs), reaction_df))
}

#'
#'
#' @export
sn2lpe_to_fa <- function(sn2LPEs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h2o +  M_2agpe <=> M_h + M_fatacid + M_g3pe",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  sn2LPEs <- sn2LPEs[stringr::str_detect(sn2LPEs, constraints[1], negate = negate[1])]
  
  # "sn2 loss"
  FAs <- unlist(lapply(sn2LPEs, function(sn2lpe) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(sn2lpe)
    fa <- paste0("FA(", fatty_acyls[2], ")")
    return(fa)
  }))
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(sn2LPEs)),
                            ReactionFormula = character(length(sn2LPEs)),
                            IsReversible = character(length(sn2LPEs)),
                            GeneAssociation = character(length(sn2LPEs)),
                            Pathway = character(length(sn2LPEs)))
  
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
                                                          "Glycerophosphoethanolamine")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_2agpe",
                                                          sn2LPEs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fatacid",
                                                          FAs)
  
  # return results
  return(list(unique(FAs), reaction_df))
}