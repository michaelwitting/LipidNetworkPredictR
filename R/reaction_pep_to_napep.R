#'
#'
#' @export
pep_to_napep_sn1 <- function(PEPs, PCs, template = NA, constraints = c("", ""), negate = c(FALSE, FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_alkenac2gpe + M_pc <=> M_alkenac2nape + M_2agpc",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 2) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  PEPs <- PEPs[stringr::str_detect(PEPs, constraints[1], negate = negate[1])]
  PCs <- PCs[stringr::str_detect(PCs, constraints[2], negate = negate[2])]
  
  # make combinatorics
  NAPEPs <- data.frame(expand.grid(PEPs = PEPs,
                                   PCs = PCs, stringsAsFactors = FALSE),
                       stringsAsFactors = FALSE)
  
  NAPEPs$LPCs <- unlist(lapply(NAPEPs$PCs,
                               function(x) {
                                 paste0("PC(0:0/", lipidomicsUtils::isolate_fatty_acyls(x)[2], ")")
                               }))
  
  
  NAPEPs$NAPEPs <- stringr::str_replace_all(NAPEPs$PEPs,
                                            "\\)$",
                                            unlist(lapply(NAPEPs$PCs,
                                                          function(x) {
                                                            paste0("/",lipidomicsUtils::isolate_fatty_acyls(x)[1], ")")
                                                          }
                                            )))
  
  
  NAPEPs$NAPEPs <- stringr::str_replace_all(NAPEPs$NAPEPs, "PE", "NAPE")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(nrow(NAPEPs)),
                            ReactionFormula = character(nrow(NAPEPs)),
                            IsReversible = character(nrow(NAPEPs)),
                            GeneAssociation = character(nrow(NAPEPs)),
                            Pathway = character(nrow(NAPEPs)))
  
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
                                                          NAPEPs$PCs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_alkenac2gpe",
                                                          NAPEPs$PEPs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_alkenac2nape",
                                                          NAPEPs$NAPEPs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_2agpc",
                                                          NAPEPs$LPCs)
  
  # return results
  return(list(unique(NAPEPs$NAPEPs), unique(NAPEPs$LPCs), reaction_df))
}


#'
#'
#' @export
pep_to_napep_sn2 <- function(PEPs, PCs, template = NA, constraints = c("", ""), negate = c(FALSE, FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_alkenac2gpe + M_pc <=> M_alkenac2nape + M_ag3pc",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 2) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  PEPs <- PEPs[stringr::str_detect(PEPs, constraints[1], negate = negate[1])]
  PCs <- PCs[stringr::str_detect(PCs, constraints[2], negate = negate[2])]
  
  # make combinatorics
  NAPEPs <- data.frame(expand.grid(PEPs = PEPs,
                                   PCs = PCs, stringsAsFactors = FALSE),
                       stringsAsFactors = FALSE)
  
  NAPEPs$LPCs <- unlist(lapply(NAPEPs$PCs,
                               function(x) {
                                 paste0("PC(", lipidomicsUtils::isolate_fatty_acyls(x)[1], "/0:0)")
                               }))
  
  
  NAPEPs$NAPEPs <- stringr::str_replace_all(NAPEPs$PEPs,
                                            "\\)$",
                                            unlist(lapply(NAPEPs$PCs,
                                                          function(x) {
                                                            paste0("/",lipidomicsUtils::isolate_fatty_acyls(x)[2], ")")
                                                          }
                                            )))
  
  NAPEPs$NAPEPs <- stringr::str_replace_all(NAPEPs$NAPEPs, "PE", "NAPE")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(nrow(NAPEPs)),
                            ReactionFormula = character(nrow(NAPEPs)),
                            IsReversible = character(nrow(NAPEPs)),
                            GeneAssociation = character(nrow(NAPEPs)),
                            Pathway = character(nrow(NAPEPs)))
  
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
                                                          NAPEPs$PCs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_alkenac2gpe",
                                                          NAPEPs$PEPs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_alkenac2nape",
                                                          NAPEPs$NAPEPs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_ag3pc",
                                                          NAPEPs$LPCs)
  
  # return results
  return(list(unique(NAPEPs$NAPEPs), unique(NAPEPs$LPCs), reaction_df))
}
