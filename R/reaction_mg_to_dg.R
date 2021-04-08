#'
#'
#' @export
sn1mg_to_dg <- function(sn1MGs, CoAs, template = NA, constraints = c("", ""), negate = c(FALSE, FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_1magol + M_fataccoa <=> M_coa + M_12dag",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 2) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  sn1MGs <- sn1MGs[stringr::str_detect(sn1MGs, constraints[1], negate = negate[1])]
  CoAs <- CoAs[stringr::str_detect(CoAs, constraints[2], negate = negate[2])]
  
  # make combinatorics
  DGs <- data.frame(expand.grid(sn1MGs = sn1MGs,
                                CoAs = CoAs, stringsAsFactors = FALSE),
                    stringsAsFactors = FALSE)
  
  
  DGs$DGs <- stringr::str_replace_all(DGs$sn1MGs,
                                      "/0:0/0:0",
                                      unlist(lapply(DGs$CoAs,
                                                    function(x) {
                                                      paste0("/",lipidomicsUtils::isolate_fatty_acyls(x), "/0:0")
                                                    }
                                      )))
  
  DGs$DGs <- stringr::str_replace_all(DGs$DGs, "MG", "DG")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(nrow(DGs)),
                            ReactionFormula = character(nrow(DGs)),
                            IsReversible = character(nrow(DGs)),
                            GeneAssociation = character(nrow(DGs)),
                            Pathway = character(nrow(DGs)))
  
  # fill with data
  reaction_df$Name <- template$reaction_name
  reaction_df$ReactionFormula <- template$reaction_formula
  reaction_df$IsReversible <- template$reaction_isReversible
  reaction_df$GeneAssociation <- template$reaction_geneAssociation
  reaction_df$Pathway <- template$reaction_pathway
  
  # adjust fixed substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_coa",
                                                          "CoA")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fataccoa",
                                                          DGs$CoAs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_1magol",
                                                          DGs$sn1MGs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_12dag",
                                                          DGs$DGs)
  
  # return results
  return(list(unique(DGs$DGs), reaction_df))
}

#'
#'
#' @export
sn2mg_to_dg <- function(sn2MGs, CoAs, template = NA, constraints = c("", ""), negate = c(FALSE, FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_fataccoa + M_mag <=> M_coa + M_12dag",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 2) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  sn2MGs <- sn2MGs[stringr::str_detect(sn2MGs, constraints[1], negate = negate[1])]
  CoAs <- CoAs[stringr::str_detect(CoAs, constraints[2], negate = negate[2])]
  
  # make combinatorics
  DGs <- data.frame(expand.grid(sn2MGs = sn2MGs,
                                CoAs = CoAs, stringsAsFactors = FALSE),
                    stringsAsFactors = FALSE)
  
  DGs$DGs <- stringr::str_replace_all(DGs$sn2MGs,
                                      "\\(0:0/",
                                      unlist(lapply(DGs$CoAs,
                                                    function(x) {
                                                      paste0("(", lipidomicsUtils::isolate_fatty_acyls(x), "/")
                                                    }
                                      )))
  
  DGs$DGs <- stringr::str_replace_all(DGs$DGs, "MG", "DG")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(nrow(DGs)),
                            ReactionFormula = character(nrow(DGs)),
                            IsReversible = character(nrow(DGs)),
                            GeneAssociation = character(nrow(DGs)),
                            Pathway = character(nrow(DGs)))
  
  # fill with data
  reaction_df$Name <- template$reaction_name
  reaction_df$ReactionFormula <- template$reaction_formula
  reaction_df$IsReversible <- template$reaction_isReversible
  reaction_df$GeneAssociation <- template$reaction_geneAssociation
  reaction_df$Pathway <- template$reaction_pathway
  
  # adjust fixed substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_coa",
                                                          "CoA")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fataccoa",
                                                          DGs$CoAs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_mag",
                                                          DGs$sn2MGs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_12dag",
                                                          DGs$DGs)
  
  # return results
  return(list(unique(DGs$DGs), reaction_df))
}
