#'
#'
#' @export
dg_to_sn1mg <- function(DGs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h2o + M_12dag <=> M_h + M_1magol + M_fatacid",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  DGs <- DGs[stringr::str_detect(DGs, constraints[1], negate = negate[1])]
  
  # make DF
  MGs <- data.frame(DG = character(length(DGs)),
                    MG = character(length(DGs)),
                    FA = character(length(DGs)))
  
  MGs$DG <- DGs

  # "sn2 loss"
  MGs$MG <- unlist(lapply(MGs$DG, function(dg) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(dg)
    mg <- paste0("MG(", fatty_acyls[1], "/0:0/0:0)")
    return(mg)
  }))
  
  MGs$FA <- unlist(lapply(MGs$DG, function(dg) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(dg)
    mg <- paste0("FA(", fatty_acyls[2], ")")
    return(mg)
  }))
  
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
                                                          "M_h2o",
                                                          "H2O")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_h",
                                                          "H+")
  
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_12dag",
                                                          MGs$DG)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_1magol",
                                                          MGs$MG)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fatacid",
                                                          MGs$FA)
  
  # return results
  return(list(unique(MGs$MG), reaction_df))
}



#'
#'
#' @export
dg_to_sn2mg <- function(DGs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h2o + M_12dag <=> M_h + M_mag + M_fatacid",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  DGs <- DGs[stringr::str_detect(DGs, constraints[1], negate = negate[1])]
  
  # make DF
  MGs <- data.frame(DG = character(length(DGs)),
                    MG = character(length(DGs)),
                    FA = character(length(DGs)))
  
  MGs$DG <- DGs
  
  # "sn2 loss"
  MGs$MG <- unlist(lapply(MGs$DG, function(dg) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(dg)
    mg <- paste0("MG(0:0/", fatty_acyls[2], "/0:0)")
    return(mg)
  }))
  
  MGs$FA <- unlist(lapply(MGs$DG, function(dg) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(dg)
    mg <- paste0("FA(", fatty_acyls[1], ")")
    return(mg)
  }))
  
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
                                                          "M_h2o",
                                                          "H2O")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_h",
                                                          "H+")
  
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_12dag",
                                                          MGs$DG)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_mag",
                                                          MGs$MG)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fatacid",
                                                          MGs$FA)
  
  # return results
  return(list(unique(MGs$MG), reaction_df))
}
