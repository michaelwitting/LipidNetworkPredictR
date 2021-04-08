#'
#'
#' @export
lpep_to_pep <- function(LPEPs, CoAs, template = NA, constraints = c("", ""), negate = c(FALSE, FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_fataccoa + M_alken2gpe <=> M_coa + M_alkenac2gpe",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 2) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  LPEPs <- LPEPs[stringr::str_detect(LPEPs, constraints[1], negate = negate[1])]
  CoAs <- CoAs[stringr::str_detect(CoAs, constraints[2], negate = negate[2])]
  
  # make combinatorics
  PEPs <- data.frame(expand.grid(LPEPs = LPEPs,
                                 CoAs = CoAs, stringsAsFactors = FALSE),
                     stringsAsFactors = FALSE)
  
  
  PEPs$PEPs <- stringr::str_replace_all(PEPs$LPEPs,
                                        "/0:0",
                                        unlist(lapply(PEPs$CoAs,
                                                      function(x) {
                                                        paste0("/",lipidomicsUtils::isolate_fatty_acyls(x))
                                                      }
                                        )))
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(nrow(PEPs)),
                            ReactionFormula = character(nrow(PEPs)),
                            IsReversible = character(nrow(PEPs)),
                            GeneAssociation = character(nrow(PEPs)),
                            Pathway = character(nrow(PEPs)))
  
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
                                                          PEPs$CoAs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_alken2gpe",
                                                          PEPs$LPEPs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_alkenac2gpe",
                                                          PEPs$PEPs)
  
  # return results
  return(list(unique(PEPs$PEPs), reaction_df))
}
