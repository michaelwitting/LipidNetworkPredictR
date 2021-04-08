#'
#'
#' @export
lpeo_to_peo <- function(LPEOs, CoAs, template = NA, constraints = c("", ""), negate = c(FALSE, FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_fataccoa + M_ak2lgpe <=> M_coa + M_akac2gpe",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 2) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  LPEOs <- LPEOs[stringr::str_detect(LPEOs, constraints[1], negate = negate[1])]
  CoAs <- CoAs[stringr::str_detect(CoAs, constraints[2], negate = negate[2])]
  
  # make combinatorics
  PEOs <- data.frame(expand.grid(LPEOs = LPEOs,
                                CoAs = CoAs, stringsAsFactors = FALSE),
                    stringsAsFactors = FALSE)
  
  
  PEOs$PEOs <- stringr::str_replace_all(PEOs$LPEOs,
                                      "/0:0",
                                      unlist(lapply(PEOs$CoAs,
                                                    function(x) {
                                                      paste0("/",lipidomicsUtils::isolate_fatty_acyls(x))
                                                    }
                                      )))
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(nrow(PEOs)),
                            ReactionFormula = character(nrow(PEOs)),
                            IsReversible = character(nrow(PEOs)),
                            GeneAssociation = character(nrow(PEOs)),
                            Pathway = character(nrow(PEOs)))
  
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
                                                          PEOs$CoAs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_ak2lgpe",
                                                          PEOs$LPEOs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_akac2gpe",
                                                          PEOs$PEOs)
  
  # return results
  return(list(unique(PEOs$PEOs), reaction_df))
}
