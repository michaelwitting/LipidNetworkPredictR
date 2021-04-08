#'
#'
#' @export
lpao_to_pao <- function(LPAOs, CoAs, template = NA, constraints = c("", ""), negate = c(FALSE, FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_alkgp + M_fataccoa <=> M_akac2gp + M_coa",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 2) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  LPAOs <- LPAOs[stringr::str_detect(LPAOs, constraints[1], negate = negate[1])]
  CoAs <- CoAs[stringr::str_detect(CoAs, constraints[2], negate = negate[2])]
  
  # make combinatorics
  PAOs <- data.frame(expand.grid(LPAOs = LPAOs,
                                CoAs = CoAs, stringsAsFactors = FALSE),
                    stringsAsFactors = FALSE)
  
  
  PAOs$PAOs <- stringr::str_replace_all(PAOs$LPAOs,
                                      "/0:0",
                                      unlist(lapply(PAOs$CoAs,
                                                    function(x) {
                                                      paste0("/",lipidomicsUtils::isolate_fatty_acyls(x))
                                                    }
                                      )))
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(nrow(PAOs)),
                            ReactionFormula = character(nrow(PAOs)),
                            IsReversible = character(nrow(PAOs)),
                            GeneAssociation = character(nrow(PAOs)),
                            Pathway = character(nrow(PAOs)))
  
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
                                                          PAOs$CoAs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_alpa_pl",
                                                          PAOs$LPAOs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_pa_pl",
                                                          PAOs$PAOs)
  
  # return results
  return(list(unique(PAOs$PAOs), reaction_df))
}