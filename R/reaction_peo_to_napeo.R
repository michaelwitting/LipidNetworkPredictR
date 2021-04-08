#'
#'
#' @export
peo_to_napeo_sn1 <- function(PEOs, PCs, template = NA, constraints = c("", ""), negate = c(FALSE, FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_akac2gpe + M_pc <=> M_akac2nape + M_2agpc",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 2) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  PEOs <- PEOs[stringr::str_detect(PEOs, constraints[1], negate = negate[1])]
  PCs <- PCs[stringr::str_detect(PCs, constraints[2], negate = negate[2])]
  
  # make combinatorics
  NAPEOs <- data.frame(expand.grid(PEOs = PEOs,
                                  PCs = PCs, stringsAsFactors = FALSE),
                      stringsAsFactors = FALSE)
  
  NAPEOs$LPCs <- unlist(lapply(NAPEOs$PCs,
                              function(x) {
                                paste0("PC(0:0/", lipidomicsUtils::isolate_fatty_acyls(x)[2], ")")
                              }))
  
  
  NAPEOs$NAPEOs <- stringr::str_replace_all(NAPEOs$PEOs,
                                          "\\)$",
                                          unlist(lapply(NAPEOs$PCs,
                                                        function(x) {
                                                          paste0("/",lipidomicsUtils::isolate_fatty_acyls(x)[1], ")")
                                                        }
                                          )))
  
  
  NAPEOs$NAPEOs <- stringr::str_replace_all(NAPEOs$NAPEOs, "PE", "NAPE")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(nrow(NAPEOs)),
                            ReactionFormula = character(nrow(NAPEOs)),
                            IsReversible = character(nrow(NAPEOs)),
                            GeneAssociation = character(nrow(NAPEOs)),
                            Pathway = character(nrow(NAPEOs)))
  
  # fill with data
  reaction_df$Name <- template$reaction_name
  reaction_df$ReactionFormula <- template$reaction_formula
  reaction_df$IsReversible <- template$reaction_isReversible
  reaction_df$GeneAssociation <- template$reaction_geneAssociation
  reaction_df$Pathway <- template$reaction_pathway
  
  # adjust fixed substrates and products
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_pc",
                                                          NAPEOs$PCs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_akac2gpe",
                                                          NAPEOs$PEOs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_akac2nape",
                                                          NAPEOs$NAPEOs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_2agpc",
                                                          NAPEOs$LPCs)
  
  # return results
  return(list(unique(NAPEOs$NAPEOs), unique(NAPEOs$LPCs), reaction_df))
}


#'
#'
#' @export
peo_to_napeo_sn2 <- function(PEOs, PCs, template = NA, constraints = c("", ""), negate = c(FALSE, FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_akac2gpe + M_pc <=> M_akac2nape + M_ag3pc",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 2) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  PEOs <- PEOs[stringr::str_detect(PEOs, constraints[1], negate = negate[1])]
  PCs <- PCs[stringr::str_detect(PCs, constraints[2], negate = negate[2])]
  
  # make combinatorics
  NAPEOs <- data.frame(expand.grid(PEOs = PEOs,
                                  PCs = PCs, stringsAsFactors = FALSE),
                      stringsAsFactors = FALSE)
  
  NAPEOs$LPCs <- unlist(lapply(NAPEOs$PCs,
                              function(x) {
                                paste0("PC(", lipidomicsUtils::isolate_fatty_acyls(x)[1], "/0:0)")
                              }))
  
  
  NAPEOs$NAPEOs <- stringr::str_replace_all(NAPEOs$PEOs,
                                          "\\)$",
                                          unlist(lapply(NAPEOs$PCs,
                                                        function(x) {
                                                          paste0("/",lipidomicsUtils::isolate_fatty_acyls(x)[2], ")")
                                                        }
                                          )))
  
  NAPEOs$NAPEOs <- stringr::str_replace_all(NAPEOs$NAPEOs, "PE", "NAPE")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(nrow(NAPEOs)),
                            ReactionFormula = character(nrow(NAPEOs)),
                            IsReversible = character(nrow(NAPEOs)),
                            GeneAssociation = character(nrow(NAPEOs)),
                            Pathway = character(nrow(NAPEOs)))
  
  # fill with data
  reaction_df$Name <- template$reaction_name
  reaction_df$ReactionFormula <- template$reaction_formula
  reaction_df$IsReversible <- template$reaction_isReversible
  reaction_df$GeneAssociation <- template$reaction_geneAssociation
  reaction_df$Pathway <- template$reaction_pathway
  
  # adjust fixed substrates and products
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_pc",
                                                          NAPEOs$PCs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_akac2gpe",
                                                          NAPEOs$PEOs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_akac2nape",
                                                          NAPEOs$NAPEOs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_ag3pc",
                                                          NAPEOs$LPCs)
  
  # return results
  return(list(unique(NAPEOs$NAPEOs), unique(NAPEOs$LPCs), reaction_df))
}
