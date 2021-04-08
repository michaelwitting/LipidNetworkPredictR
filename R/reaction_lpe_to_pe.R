#'
#'
#' @export
sn1lpe_to_pe <- function(sn1LPEs, CoAs, template = NA, constraints = c("", ""), negate = c(FALSE, FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_fataccoa + M_acg3pe <=> M_coa + M_pe",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 2) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  sn1LPEs <- sn1LPEs[stringr::str_detect(sn1LPEs, constraints[1], negate = negate[1])]
  CoAs <- CoAs[stringr::str_detect(CoAs, constraints[2], negate = negate[2])]
  
  # make combinatorics
  PEs <- data.frame(expand.grid(sn1LPEs = sn1LPEs,
                                CoAs = CoAs, stringsAsFactors = FALSE),
                    stringsAsFactors = FALSE)
  
  
  PEs$PEs <- stringr::str_replace_all(PEs$sn1LPEs,
                                      "/0:0",
                                      unlist(lapply(PEs$CoAs,
                                                    function(x) {
                                                      paste0("/",lipidomicsUtils::isolate_fatty_acyls(x))
                                                    }
                                      )))
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(nrow(PEs)),
                            ReactionFormula = character(nrow(PEs)),
                            IsReversible = character(nrow(PEs)),
                            GeneAssociation = character(nrow(PEs)),
                            Pathway = character(nrow(PEs)))
  
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
                                                          PEs$CoAs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_acg3pe",
                                                          PEs$sn1LPEs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_pe",
                                                          PEs$PEs)
  
  # return results
  return(list(unique(PEs$PEs), reaction_df))
}
