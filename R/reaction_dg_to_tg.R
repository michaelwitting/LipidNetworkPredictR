#'
#'
#' @export
dg_to_tg <- function(DGs, CoAs, template = NA, constraints = c("", ""), negate = c(FALSE, FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_fataccoa + M_12dag <=> M_coa + M_tag",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 2) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  DGs <- DGs[stringr::str_detect(DGs, constraints[1], negate = negate[1])]
  CoAs <- CoAs[stringr::str_detect(CoAs, constraints[2], negate = negate[2])]
  
  # make combinatorics
  TGs <- data.frame(expand.grid(DGs = DGs,
                                CoAs = CoAs, stringsAsFactors = FALSE),
                    stringsAsFactors = FALSE)
  
  
  TGs$TGs <- stringr::str_replace_all(TGs$DGs,
                                      "/0:0",
                                      paste0("/", unlist(lapply(TGs$CoAs,
                                                    lipidomicsUtils::isolate_fatty_acyls))))
  
  TGs$TGs <- stringr::str_replace_all(TGs$TGs, "^DG", "TG")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(nrow(TGs)),
                            ReactionFormula = character(nrow(TGs)),
                            IsReversible = character(nrow(TGs)),
                            GeneAssociation = character(nrow(TGs)),
                            Pathway = character(nrow(TGs)))
  
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
                                                          TGs$CoAs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_12dag",
                                                          TGs$DGs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_tag",
                                                          TGs$TGs)

  # return results
  return(list(unique(TGs$TGs), reaction_df))
}