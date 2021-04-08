#'
#'
#' @export
tg_to_dg <- function(TGs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h2o + M_tag <=> M_h + M_fatacid + M_12dag",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  TGs <- TGs[stringr::str_detect(TGs, constraints[1], negate = negate[1])]
  
  # make DF
  DGs <- data.frame(TGs = character(length(TGs)),
                    sn1Loss = character(length(TGs)),
                    sn2Loss = character(length(TGs)))
  
  DGs$TGs <- TGs
  
  # "sn1 loss"
  DGs$sn1Loss_dg <- unlist(lapply(DGs$TGs, .sn1Loss_dg))
  DGs$sn1Loss_fa <- unlist(lapply(DGs$TGs, .sn1Loss_fa))
  
  # "sn3 loss"
  DGs$sn3Loss_dg <- unlist(lapply(DGs$TGs, .sn3Loss_dg))
  DGs$sn3Loss_fa <- unlist(lapply(DGs$TGs, .sn3Loss_fa))
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(2 *length(TGs)),
                            ReactionFormula = character(2 * length(TGs)),
                            IsReversible = character(2 * length(TGs)),
                            GeneAssociation = character(2 * length(TGs)),
                            Pathway = character(2 * length(TGs)))
  
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
                                                          "M_tag",
                                                          c(DGs$TGs, DGs$TGs))
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_12dag",
                                                          c(DGs$sn1Loss_dg, DGs$sn3Loss_dg))
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fatacid",
                                                          c(DGs$sn1Loss_fa, DGs$sn3Loss_fa))
  
  # return results
  return(list(unique(DGs$sn1Loss_dg, DGs$sn1Loss_fa), reaction_df))
}

.sn1Loss_dg <- function(TG) {
  
  fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(TG)
  dg <- paste0("DG(", fatty_acyls[3], "/", fatty_acyls[2], "/0:0)")
  
  return(dg)
  
}

.sn3Loss_dg <- function(TG) {
  
  fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(TG)
  dg <- paste0("DG(", fatty_acyls[1], "/", fatty_acyls[2], "/0:0)")
  
  return(dg)
  
}

.sn1Loss_fa <- function(TG) {
  
  fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(TG)
  fa <- paste0("FA(", fatty_acyls[1], ")")
  
  return(fa)
  
}

.sn3Loss_fa <- function(TG) {
  
  fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(TG)
  fa <- paste0("FA(", fatty_acyls[3], ")")
  
  return(fa)
  
}