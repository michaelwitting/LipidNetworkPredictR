#'
#'
#' @export
lpa_to_pa <- function(LPAs, CoAs, template = NA, constraints = c("", ""), negate = c(FALSE, FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_fataccoa + M_alpa_pl <=> M_coa + M_pa_pl",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 2) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  LPAs <- LPAs[stringr::str_detect(LPAs, constraints[1], negate = negate[1])]
  CoAs <- CoAs[stringr::str_detect(CoAs, constraints[2], negate = negate[2])]

  # make combinatorics
  PAs <- data.frame(expand.grid(LPAs = LPAs,
                                CoAs = CoAs, stringsAsFactors = FALSE),
                    stringsAsFactors = FALSE)

  
  PAs$PAs <- stringr::str_replace_all(PAs$LPAs,
                                      "/0:0",
                                      unlist(lapply(PAs$CoAs,
                                                    function(x) {
                                                      paste0("/",lipidomicsUtils::isolate_fatty_acyls(x))
                                                    }
                                                    )))

  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(nrow(PAs)),
                            ReactionFormula = character(nrow(PAs)),
                            IsReversible = character(nrow(PAs)),
                            GeneAssociation = character(nrow(PAs)),
                            Pathway = character(nrow(PAs)))
  
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
                                                          PAs$CoAs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_alpa_pl",
                                                          PAs$LPAs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_pa_pl",
                                                          PAs$PAs)
  
  # return results
  return(list(unique(PAs$PAs), reaction_df))
}