#' @name .create_template
#' 
#' @title Create template with mininum information or return a given template
#' 
#' @description 
#' Helper function for \code{create_reaction}.
#' 
#' The function \code{.create_template} returns a data.frame of the template
#' of a given \code{reaction}. The returned data.frame is either 
#' the pass-through of the given \code{template} or will be created and 
#' filled with minimum information using the information from
#' \code{reaction}.
#' 
#' @details 
#' The function \code{.create_template} will check the \code{template} object
#' for the following information:
#' 
#' \itemize{
#'     \item column "reaction_name" with character entry,
#'     \item column "reaction_formula" with character entry,
#'     \item column "reaction_isReversible" with character entry, 
#'     \item column "reaction_geneAssociation" with character entry,
#'     \item column "reaction_pathway" with character entry.
#' }    
#' 
#' If \code{template} is \code{NA} or \code{NULL}, \code{template} will be set 
#' to an empty list.
#' 
#' @param template \code{list}, \code{NA} or \code{NULL}
#' @param reaction \code{character(1)}
#' 
#' @return data.frame
#' 
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#'    and Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @importFrom stringi stri_replace_all_fixed stri_replace_all_regex stri_split
#' 
#' @examples 
#' LipidNetworkPredictR:::.create_template(template = list(), reaction = "RHEA:15421")
.create_template <- function(template = list(), reaction = "RHEA:15421") {
    
    ## checks on template, convert to empty list if template is NA or NULL
    if (length(template) == 1) {
        if (is.na(template)) 
            template <- list()
    }
    
    if (length(template) == 0) {
        if (is.null(template))
            template <- list()
    }

    if (!is.list(template)) 
        stop("'template' has to be a list")
    
    ## convert to list if template is a data.frame
    if (is.data.frame(template)) {
        template <- as.list(template)
    }

    ## write reaction to reation_RHEA entry
    template[["reaction_RHEA"]] <- reaction
        
    ## check integrity (all names are present in template)
    if (!"reaction_name" %in% names(template)) 
        template[["reaction_name"]] <- ""
    if (!"reaction_isReversible" %in% names(template)) 
        template[["reaction_isReversible"]] <- ""
    if (!"reaction_geneAssociation" %in% names(template)) 
        template[["reaction_geneAssociation"]] <- ""
    if (!"reaction_pathway" %in% names(template))
        template[["reaction_pathway"]] <- ""
    if (!"reaction_constraints" %in% names(template))
        template[["reaction_constraints"]] <- ""
    if (!"reaction_constraints_negate" %in% names(template))
        template[["reaction_constraints_negate"]] <- FALSE
    
    ## all entries are of mode character
    if (!is.vector(template$reaction_name, "character"))
        stop("entry 'reaction_name' is not of mode 'character'")
    if (!is.vector(template$reaction_isReversible, "character"))
        stop("entry 'reaction_isReversible' is not of mode 'character'")
    if (!is.vector(template$reaction_geneAssociation, "character"))
        stop("entry 'reaction_geneAssociation' is not of mode 'character'")
    if (!is.vector(template$reaction_pathway, "character"))
        stop("entry 'reaction_pathway' is not of mode 'character'")

    if (is.null(template$reaction_formula)) {
        ## this event happens if the param template is not a list
        
        ## check if reaction is in the parameter space
        .check_reaction(reaction)
        
        ## acyldhap_to_alkyldhap
        if (reaction == "RHEA:36171") 
            .formula <- "M_AcylDHAP + M_FAO = M_AlkylDHAP + M_FA + M_H+"
        if (reaction == "RHEA:36172") 
            .formula <- "M_AcylDHAP + M_FAO => M_AlkylDHAP + M_FA + M_H+"
        if (reaction == "RHEA:36173") 
            .formula <- "M_AcylDHAP + M_FAO <= M_AlkylDHAP + M_FA + M_H+"
        if (reaction == "RHEA:36174") 
            .formula <- "M_AcylDHAP + M_FAO <=> M_AlkylDHAP + M_FA + M_H+"
        
        ## alkyldhap_to_lpao
        if (reaction == "RHEA:36175") 
            .formula <- "M_AlkylDHAP + M_H+ + M_NADPH = M_LPA-O + M_NADP"
        if (reaction == "RHEA:36176") 
            .formula <- "M_AlkylDHAP + M_H+ + M_NADPH => M_LPA-O + M_NADP"
        if (reaction == "RHEA:36177") 
            .formula <- "M_AlkylDHAP + M_H+ + M_NADPH <= M_LPA-O + M_NADP"
        if (reaction == "RHEA:36178") 
            .formula <- "M_AlkylDHAP + M_H+ + M_NADPH <=> M_LPA-O + M_NADP"
        
        ## cerp_to_cer
        if (reaction == "RHEA:33743")
            .formula <- "M_H2O + M_CerP = M_Pi + M_Cer"
        if (reaction == "RHEA:33744")
            .formula <- "M_H2O + M_CerP => M_Pi + M_Cer"
        if (reaction == "RHEA:33745")
            .formula <- "M_H2O + M_CerP <= M_Pi + M_Cer"
        if (reaction == "RHEA:33746")
            .formula <- "M_H2O + M_CerP <=> M_Pi + M_Cer"
        
        ## cdpdg_to_pgp
        if (reaction == "RHEA:12593") 
            .formula <- "M_CDP-DG + M_Glycerol-3-P = M_PGP + M_CMP + M_H+"
        if (reaction == "RHEA:12594") 
            .formula <- "M_CDP-DG + M_Glycerol-3-P => M_PGP + M_CMP + M_H+"
        if (reaction == "RHEA:12595") 
            .formula <- "M_CDP-DG + M_Glycerol-3-P <= M_PGP + M_CMP + M_H+"
        if (reaction == "RHEA:12596") 
            .formula <- "M_CDP-DG + M_Glycerol-3-P <=> M_PGP + M_CMP + M_H+"
        
        ## cdpdg_to_pi
        if (reaction == "RHEA:11580") 
            .formula <- "M_CDP-DG + M_myo-Inositol = M_PI + M_CMP + M_H+"
        if (reaction == "RHEA:11581") 
            .formula <- "M_CDP-DG + M_myo-Inositol => M_PI + M_CMP + M_H+"
        if (reaction == "RHEA:11582") 
            .formula <- "M_CDP-DG + M_myo-Inositol <= M_PI + M_CMP + M_H+"
        if (reaction == "RHEA:11583") 
            .formula <- "M_CDP-DG + M_myo-Inositol <=> M_PI + M_CMP + M_H+"
        
        ## cer_to_cerp
        if (reaction == "RHEA:17929") 
            .formula <- "M_Cer + M_ATP = M_ADP + M_CerP + M_H+"
        if (reaction == "RHEA:17930") 
            .formula <- "M_Cer + M_ATP => M_ADP + M_CerP + M_H+"
        if (reaction == "RHEA:17931") 
            .formula <- "M_Cer + M_ATP <= M_ADP + M_CerP + M_H+"
        if (reaction == "RHEA:17932") 
            .formula <- "M_Cer + M_ATP <=> M_ADP + M_CerP + M_H+"
        
        ## cer_to_glccer
        if (reaction == "RHEA:12088")
            .formula <- "M_Cer + M_UDP-Glucose = M_GlcCer + M_H+ + M_UDP"
        if (reaction == "RHEA:12089")
            .formula <- "M_Cer + M_UDP-Glucose => M_GlcCer + M_H+ + M_UDP"
        if (reaction == "RHEA:12090")
            .formula <- "M_Cer + M_UDP-Glucose <= M_GlcCer + M_H+ + M_UDP"
        if (reaction == "RHEA:12091")
            .formula <- "M_Cer + M_UDP-Glucose <=> M_GlcCer + M_H+ + M_UDP"
        
        ## cer_to_sm
        if (reaction == "RHEA:18765") 
            .formula <- "M_PC + M_Cer = M_1,2-DG + M_SM"
        if (reaction == "RHEA:18766") 
            .formula <- "M_PC + M_Cer => M_1,2-DG + M_SM"
        if (reaction == "RHEA:18767") 
            .formula <- "M_PC + M_Cer <= M_1,2-DG + M_SM"
        if (reaction == "RHEA:18768") 
            .formula <- "M_PC + M_Cer <=> M_1,2-DG + M_SM"
        
        ## cl_to_lcl
        if (reaction == "RHEA:32935")
            .formula <- "M_CL + M_H2O = M_1,2,4-LCL + M_FA + M_H+"
        if (reaction == "RHEA:32936")
            .formula <- "M_CL + M_H2O => M_1,2,4-LCL + M_FA + M_H+"
        if (reaction == "RHEA:32937")
            .formula <- "M_CL + M_H2O <= M_1,2,4-LCL + M_FA + M_H+"
        if (reaction == "RHEA:32938")
            .formula <- "M_CL + M_H2O <=> M_1,2,4-LCL + M_FA + M_H+"
        
        ## coa_to_acdhap
        if (reaction == "RHEA:17657")
            .formula <- "M_AcylCoA + M_Dihydroxyacetone-P = M_AcylDHAP + M_CoA"
        if (reaction == "RHEA:17658")
            .formula <- "M_AcylCoA + M_Dihydroxyacetone-P => M_AcylDHAP + M_CoA"
        if (reaction == "RHEA:17659")
            .formula <- "M_AcylCoA + M_Dihydroxyacetone-P <= M_AcylDHAP + M_CoA"
        if (reaction == "RHEA:17660")
            .formula <- "M_AcylCoA + M_Dihydroxyacetone-P <=> M_AcylDHAP + M_CoA"
    
        ## coa_to_ce
        if (reaction == "RHEA:17729")
            .formula <- "M_Cholesterol + M_AcylCoA = M_CE + M_CoA"
        if (reaction == "RHEA:17730")
            .formula <- "M_Cholesterol + M_AcylCoA => M_CE + M_CoA"
        if (reaction == "RHEA:17731")
            .formula <- "M_Cholesterol + M_AcylCoA <= M_CE + M_CoA"
        if (reaction == "RHEA:17732")
            .formula <- "M_Cholesterol + M_AcylCoA <=> M_CE + M_CoA"
        
        ## coa_to_fao
        if (reaction == "RHEA:52716") 
            .formula <- "M_AcylCoA + 2 M_H+ + 2 M_NADPH = M_FAO + M_CoA + 2 M_NADP"
        if (reaction == "RHEA:52717") 
            .formula <- "M_AcylCoA + 2 M_H+ + 2 M_NADPH => M_FAO + M_CoA + 2 M_NADP"
        if (reaction == "RHEA:52718") 
            .formula <- "M_AcylCoA + 2 M_H+ + 2 M_NADPH <= M_FAO + M_CoA + 2 M_NADP"
        if (reaction == "RHEA:52719") 
            .formula <- "M_AcylCoA + 2 M_H+ + 2 M_NADPH <=> M_FAO + M_CoA + 2 M_NADP"
        
        ## coa_to_lpa
        if (reaction == "RHEA:15325") 
            .formula <- "M_AcylCoA + M_Glycerol-3-P = M_LPA + M_CoA"
        if (reaction == "RHEA:15326") 
            .formula <- "M_AcylCoA + M_Glycerol-3-P => M_LPA + M_CoA"
        if (reaction == "RHEA:15327") 
            .formula <- "M_AcylCoA + M_Glycerol-3-P <= M_LPA + M_CoA"
        if (reaction == "RHEA:15328")
            .formula <- "M_AcylCoA + M_Glycerol-3-P <=> M_LPA + M_CoA"

        ## dg_to_sn1mg
        if (reaction %in% c("RHEA:44712", "RHEA:35663"))
            .formula <- "M_1,2-DG + M_H2O = M_1-MG + M_FA + M_H+"
        if (reaction %in% c("RHEA:44713", "RHEA:35664"))
            .formula <- "M_1,2-DG + M_H2O => M_1-MG + M_FA + M_H+"
        if (reaction %in% c("RHEA:44714", "RHEA:35665"))
            .formula <- "M_1,2-DG + M_H2O <= M_1-MG + M_FA + M_H+"
        if (reaction %in% c("RHEA:44715", "RHEA:35666"))
            .formula <- "M_1,2-DG + M_H2O <=> M_1-MG + M_FA + M_H+"

        ## dg_to_sn2mg
        if (reaction == "RHEA:33275")
            .formula <- "M_1,2-DG + M_H2O = M_2-MG + M_FA + M_H+"
        if (reaction == "RHEA:33276")
            .formula <- "M_1,2-DG + M_H2O => M_2-MG + M_FA + M_H+"
        if (reaction == "RHEA:33277")
            .formula <- "M_1,2-DG + M_H2O <= M_2-MG + M_FA + M_H+"
        if (reaction == "RHEA:33278")
            .formula <- "M_1,2-DG + M_H2O <=> M_2-MG + M_FA + M_H+"
        
        ## dg_to_pa
        if (reaction == "RHEA:10272")
            .formula <- "M_1,2-DG + M_ATP = M_PA + M_ADP + M_H+"
        if (reaction == "RHEA:10273")
            .formula <- "M_1,2-DG + M_ATP => M_PA + M_ADP + M_H+"
        if (reaction == "RHEA:10274")
            .formula <- "M_1,2-DG + M_ATP <= M_PA + M_ADP + M_H+"
        if (reaction == "RHEA:10275")
            .formula <- "M_1,2-DG + M_ATP <=> M_PA + M_ADP + M_H+"
        
        ## dg_to_pc
        if (reaction == "RHEA:32939") 
            .formula <- "M_1,2-DG + M_CDP-Choline = M_PC + M_CMP + M_H+"
        if (reaction == "RHEA:32940") 
            .formula <- "M_1,2-DG + M_CDP-Choline => M_PC + M_CMP + M_H+"
        if (reaction == "RHEA:32941") 
            .formula <- "M_1,2-DG + M_CDP-Choline <= M_PC + M_CMP + M_H+"
        if (reaction == "RHEA:32942") 
            .formula <- "M_1,2-DG + M_CDP-Choline <=> M_PC + M_CMP + M_H+"
        
        ## dg_to_pe
        if (reaction == "RHEA:32943") 
            .formula <- "M_1,2-DG + M_CDP-Ethanolamine = M_PE + M_CMP + M_H+"
        if (reaction == "RHEA:32944") 
            .formula <- "M_1,2-DG + M_CDP-Ethanolamine => M_PE + M_CMP + M_H+"
        if (reaction == "RHEA:32945") 
            .formula <- "M_1,2-DG + M_CDP-Ethanolamine <= M_PE + M_CMP + M_H+"
        if (reaction == "RHEA:32946") 
            .formula <- "M_1,2-DG + M_CDP-Ethanolamine <=> M_PE + M_CMP + M_H+"
        
        ## dg_to_tg
        if (reaction == "RHEA:10868") 
            .formula <- "M_1,2-DG + M_AcylCoA = M_TG + M_CoA"
        if (reaction == "RHEA:10869") 
            .formula <- "M_1,2-DG + M_AcylCoA => M_TG + M_CoA"
        if (reaction == "RHEA:10870") 
            .formula <- "M_1,2-DG + M_AcylCoA <= M_TG + M_CoA"
        if (reaction == "RHEA:10871") 
            .formula <- "M_1,2-DG + M_AcylCoA <=> M_TG + M_CoA"
        
        ## dgo_to_pco
        if (reaction == "RHEA:36179")
            .formula <- "M_DG-O + M_CDP-Choline = M_PC-O + M_CMP + M_H+"
        if (reaction == "RHEA:36180")
            .formula <- "M_DG-O + M_CDP-Choline => M_PC-O + M_CMP + M_H+"
        if (reaction == "RHEA:36181")
            .formula <- "M_DG-O + M_CDP-Choline <= M_PC-O + M_CMP + M_H+"
        if (reaction == "RHEA:36182")
            .formula <- "M_DG-O + M_CDP-Choline <=> M_PC-O + M_CMP + M_H+"
        
        ## dgo_to_peo
        if (reaction == "RHEA:36187")
            .formula <- "M_DG-O + M_CDP-Ethanolamine = M_PE-O + M_CMP + M_H+"
        if (reaction == "RHEA:36188")
            .formula <- "M_DG-O + M_CDP-Ethanolamine => M_PE-O + M_CMP + M_H+"
        if (reaction == "RHEA:36189")
            .formula <- "M_DG-O + M_CDP-Ethanolamine <= M_PE-O + M_CMP + M_H+"
        if (reaction == "RHEA:36190")
            .formula <- "M_DG-O + M_CDP-Ethanolamine <=> M_PE-O + M_CMP + M_H+"
        
        ## dhcer_to_cer
        if (reaction == "RHEA:46544")
            .formula <- "M_DhCer + 2 M_Fe2+-cytochrome_b5 + 2 M_H+ + M_O2 = 2 M_Fe3+-cytochrome_b5 + M_Cer + 2 M_H2O"
        if (reaction == "RHEA:46545")
            .formula <- "M_DhCer + 2 M_Fe2+-cytochrome_b5 + 2 M_H+ + M_O2 => 2 M_Fe3+-cytochrome_b5 + M_Cer + 2 M_H2O"
        if (reaction == "RHEA:46546")
            .formula <- "M_DhCer + 2 M_Fe2+-cytochrome_b5 + 2 M_H+ + M_O2 <= 2 M_Fe3+-cytochrome_b5 + M_Cer + 2 M_H2O"
        if (reaction == "RHEA:46547")
            .formula <- "M_DhCer + 2 M_Fe2+-cytochrome_b5 + 2 M_H+ + M_O2 <=> 2 M_Fe3+-cytochrome_b5 + M_Cer + 2 M_H2O"
        
        ## dhcer_to_dhsm
        if (reaction == "RHEA:44620")
            .formula <- "M_PC + M_DhCer = M_1,2-DG + M_DhSM"
        if (reaction == "RHEA:44621")
            .formula <- "M_PC + M_DhCer => M_1,2-DG + M_DhSM"
        if (reaction == "RHEA:44622")
            .formula <- "M_PC + M_DhCer <= M_1,2-DG + M_DhSM"
        if (reaction == "RHEA:44623")
            .formula <- "M_PC + M_DhCer <=> M_1,2-DG + M_DhSM"
        
        ## dhsm_to_dhcer
        if (reaction == "RHEA:45300")
            .formula <- "M_DhSM + M_H2O = M_DhCer + M_H+ + M_Phosphocholine"
        if (reaction == "RHEA:45301")
            .formula <- "M_DhSM + M_H2O => M_DhCer + M_H+ + M_Phosphocholine"
        if (reaction == "RHEA:45302")
            .formula <- "M_DhSM + M_H2O <= M_DhCer + M_H+ + M_Phosphocholine"
        if (reaction == "RHEA:45303")
            .formula <- "M_DhSM + M_H2O <=> M_DhCer + M_H+ + M_Phosphocholine"
        
        ## fa_to_coa
        if (reaction %in% c("RHEA:15421", "RHEA:38883")) 
            .formula <- "M_FA + M_ATP + M_CoA = M_AcylCoA + M_AMP + M_PPi"
        if (reaction %in% c("RHEA:15422", "RHEA:38884"))
            .formula <- "M_FA + M_ATP + M_CoA => M_AcylCoA + M_AMP + M_PPi"
        if (reaction %in% c("RHEA:15423", "RHEA:38885"))
            .formula <- "M_FA + M_ATP + M_CoA <= M_AcylCoA + M_AMP + M_PPi"
        if (reaction %in% c("RHEA:15424", "RHEA:38886"))
            .formula <- "M_FA + M_ATP + M_CoA <=> M_AcylCoA + M_AMP + M_PPi"
        
        ## lcl_to_cl
        if (reaction == "RHEA:35839")
            .formula <- "M_1,2,4-LCL + M_AcylCoA = M_CL + M_CoA"
        if (reaction == "RHEA:35840")
            .formula <- "M_1,2,4-LCL + M_AcylCoA => M_CL + M_CoA"
        if (reaction == "RHEA:35841")
            .formula <- "M_1,2,4-LCL + M_AcylCoA <= M_CL + M_CoA"
        if (reaction == "RHEA:35842")
            .formula <- "M_1,2,4-LCL + M_AcylCoA <=> M_CL + M_CoA"
        
        ## lnape_to_gpnae
        if (reaction == "RHEA:45420")
            .formula <- "M_H2O + M_LNAPE = M_FA + M_H+ + M_GPNAE"
        if (reaction == "RHEA:45421")
            .formula <- "M_H2O + M_LNAPE => M_FA + M_H+ + M_GPNAE"
        if (reaction == "RHEA:45422")
            .formula <- "M_H2O + M_LNAPE <= M_FA + M_H+ + M_GPNAE"
        if (reaction == "RHEA:45423")
            .formula <- "M_H2O + M_LNAPE <=> M_FA + M_H+ + M_GPNAE"
        
        ## lpa_to_pa
        if (reaction == "RHEA:19709")
            .formula <- "M_LPA + M_AcylCoA = M_PA + M_CoA"
        if (reaction == "RHEA:19710")
            .formula <- "M_LPA + M_AcylCoA => M_PA + M_CoA"
        if (reaction == "RHEA:19711") 
            .formula <- "M_LPA + M_AcylCoA <= M_PA + M_CoA"
        if (reaction == "RHEA:19712") 
            .formula <- "M_LPA + M_AcylCoA <=> M_PA + M_CoA"
        
        ## lpao_to_pao
        if (reaction == "RHEA:36235") 
             .formula <- "M_LPA-O + M_AcylCoA = M_PA-O + M_CoA"
        if (reaction == "RHEA:36236") 
            .formula <- "M_LPA-O + M_AcylCoA => M_PA-O + M_CoA"
        if (reaction == "RHEA:36237") 
            .formula <- "M_LPA-O + M_AcylCoA <= M_PA-O + M_CoA"
        if (reaction == "RHEA:36238") 
            .formula <- "M_LPA-O + M_AcylCoA <=> M_PA-O + M_CoA"
        
        ## sn1lpc_to_fa
        if (reaction == "RHEA:15177")
            .formula <- "M_1-LPC + M_H2O = M_FA + M_H+ + M_Glycerophosphocholine"
        if (reaction == "RHEA:15178")
            .formula <- "M_1-LPC + M_H2O => M_FA + M_H+ + M_Glycerophosphocholine"
        if (reaction == "RHEA:15179")
            .formula <- "M_1-LPC + M_H2O <= M_FA + M_H+ + M_Glycerophosphocholine"
        if (reaction == "RHEA:15180")
            .formula <- "M_1-LPC + M_H2O <=> M_FA + M_H+ + M_Glycerophosphocholine"
        
        ## sn2lpc_to_fa
        if (reaction == "RHEA:44696") 
            .formula <- "M_2-LPC + M_H2O = M_FA + M_H+ + M_Glycerophosphocholine"
        if (reaction == "RHEA:44697") 
            .formula <- "M_2-LPC + M_H2O => M_FA + M_H+ + M_Glycerophosphocholine"
        if (reaction == "RHEA:44698") 
            .formula <- "M_2-LPC + M_H2O <= M_FA + M_H+ + M_Glycerophosphocholine"
        if (reaction == "RHEA:44699") 
            .formula <- "M_2-LPC + M_H2O <=> M_FA + M_H+ + M_Glycerophosphocholine"
        
        ## sn1lpc_to_pc
        if (reaction == "RHEA:12937") 
            .formula <- "M_1-LPC + M_AcylCoA = M_PC + M_CoA"
        if (reaction == "RHEA:12938") 
            .formula <- "M_1-LPC + M_AcylCoA => M_PC + M_CoA"
        if (reaction == "RHEA:12939") 
            .formula <- "M_1-LPC + M_AcylCoA <= M_PC + M_CoA"
        if (reaction == "RHEA:12940") 
            .formula <- "M_1-LPC + M_AcylCoA <=> M_PC + M_CoA"
        
        ## lpco_to_pco
        if (reaction == "RHEA:23992") 
            .formula <- "M_1-LPC-O + M_AcylCoA = M_PC-O + M_CoA"
        if (reaction == "RHEA:23993") 
            .formula <- "M_1-LPC-O + M_AcylCoA => M_PC-O + M_CoA"
        if (reaction == "RHEA:23994") 
            .formula <- "M_1-LPC-O + M_AcylCoA <= M_PC-O + M_CoA"
        if (reaction == "RHEA:23995") 
            .formula <- "M_1-LPC-O + M_AcylCoA <=> M_PC-O + M_CoA"
        
        ## sn1lpe_to_fa
        if (reaction == "RHEA:32967") 
            .formula <- "M_1-LPE + M_H2O = M_FA + M_H+ + M_Glycerophosphoethanolamine"
        if (reaction == "RHEA:32968") 
            .formula <- "M_1-LPE + M_H2O => M_FA + M_H+ + M_Glycerophosphoethanolamine"
        if (reaction == "RHEA:32969") 
            .formula <- "M_1-LPE + M_H2O <= M_FA + M_H+ + M_Glycerophosphoethanolamine"
        if (reaction == "RHEA:32970") 
            .formula <- "M_1-LPE + M_H2O <=> M_FA + M_H+ + M_Glycerophosphoethanolamine"
        
        ## sn2lpe_to_fa
        if (reaction == "sn2lpe_to_fa")
            .formula <- "M_2-LPE + M_H2O <=> M_FA + M_H+ + M_Glycerophosphoethanolamine"
        
        ## sn1lpe_to_pe
        if (reaction == "RHEA:32995") 
            .formula <- "M_1-LPE + M_AcylCoA = M_PE + M_CoA"
        if (reaction == "RHEA:32996") 
            .formula <- "M_1-LPE + M_AcylCoA => M_PE + M_CoA"
        if (reaction == "RHEA:32997") 
            .formula <- "M_1-LPE + M_AcylCoA <= M_PE + M_CoA"
        if (reaction == "RHEA:32998") 
            .formula <- "M_1-LPE + M_AcylCoA <=> M_PE + M_CoA"
        
        ## sn1lpi_to_pi
        if (reaction == "RHEA:33195")
            .formula <- "M_1-LPI + M_AcylCoA = M_PI + M_CoA"
        if (reaction == "RHEA:33196")
            .formula <- "M_1-LPI + M_AcylCoA => M_PI + M_CoA"
        if (reaction == "RHEA:33197")
            .formula <- "M_1-LPI + M_AcylCoA <= M_PI + M_CoA"
        if (reaction == "RHEA:33198")
            .formula <- "M_1-LPI + M_AcylCoA <=> M_PI + M_CoA"
        
        ## lpeo_to_peo
        if (reaction == "lpeo_to_peo")
            .formula <- "M_AcylCoA + M_LPE-O <=> M_CoA + M_PE-O"
        
        ## lpep_to_pep
        if (reaction == "RHEA:16245")
            .formula <- "M_LPE-P + M_AcylCoA = M_PE-P + M_CoA"
        if (reaction == "RHEA:16246")    
            .formula <- "M_LPE-P + M_AcylCoA => M_PE-P + M_CoA"
        if (reaction == "RHEA:16247")
            .formula <- "M_LPE-P + M_AcylCoA <= M_PE-P + M_CoA"
        if (reaction == "RHEA:16248")
            .formula <- "M_LPE-P + M_AcylCoA <=> M_PE-P + M_CoA"
        
        ## sn1mg_to_dg
        if (reaction %in% c("RHEA:38463", "RHEA:39943"))
            .formula <- "M_1-MG + M_AcylCoA = M_1,2-DG + M_CoA"
        if (reaction %in% c("RHEA:38464", "RHEA:39944"))
            .formula <- "M_1-MG + M_AcylCoA => M_1,2-DG + M_CoA"
        if (reaction %in% c("RHEA:38465", "RHEA:39945"))
            .formula <- "M_1-MG + M_AcylCoA <= M_1,2-DG + M_CoA"
        if (reaction %in% c("RHEA:38466", "RHEA:39946"))
            .formula <- "M_1-MG + M_AcylCoA <=> M_1,2-DG + M_CoA"
        
        ## sn2mg_to_dg
        if (reaction %in% c("RHEA:32947", "RHEA:16741"))
            .formula <- "M_2-MG + M_AcylCoA = M_1,2-DG + M_CoA"
        if (reaction %in% c("RHEA:32948", "RHEA:16742"))
            .formula <- "M_2-MG + M_AcylCoA => M_1,2-DG + M_CoA"
        if (reaction %in% c("RHEA:32949", "RHEA:16743"))
            .formula <- "M_2-MG + M_AcylCoA <= M_1,2-DG + M_CoA"
        if (reaction %in% c("RHEA:32950", "RHEA:16744"))
            .formula <- "M_2-MG + M_AcylCoA <=> M_1,2-DG + M_CoA"
        
        ## sn1mg_to_fa
        if (reaction == "RHEA:34019")
            .formula <- "M_1-MG + M_H2O = M_FA + M_Glycerol + M_H+"
        if (reaction == "RHEA:34020")
            .formula <- "M_1-MG + M_H2O => M_FA + M_Glycerol + M_H+"
        if (reaction == "RHEA:34021")
            .formula <- "M_1-MG + M_H2O <= M_FA + M_Glycerol + M_H+"
        if (reaction == "RHEA:34022")
            .formula <- "M_1-MG + M_H2O <=> M_FA + M_Glycerol + M_H+"

        ## sn2mg_to_fa
        if (reaction %in% c("RHEA:32871"))
            .formula <-  "M_H2O + M_2-MG = M_FA + M_Glycerol + M_H+"
        if (reaction %in% c("RHEA:32872"))
            .formula <-  "M_H2O + M_2-MG => M_FA + M_Glycerol + M_H+"
        if (reaction %in% c("RHEA:32873"))
            .formula <-  "M_H2O + M_2-MG <= M_FA + M_Glycerol + M_H+"
        if (reaction %in% c("RHEA:32874"))
            .formula <-  "M_H2O + M_2-MG <=> M_FA + M_Glycerol + M_H+"
        
        ## sn1mg_to_lpa
        if (reaction == "RHEA:33747")
            .formula <- "M_1-MG + M_ATP = M_LPA + M_ADP + M_H+"
        if (reaction == "RHEA:33748")
            .formula <- "M_1-MG + M_ATP => M_LPA + M_ADP + M_H+"
        if (reaction == "RHEA:33749")
            .formula <- "M_1-MG + M_ATP <= M_LPA + M_ADP + M_H+"
        if (reaction == "RHEA:33750")
            .formula <- "M_1-MG + M_ATP <=> M_LPA + M_ADP + M_H+"
        
        ## sn2mg_to_sn1mg
        if (reaction == "sn2mg_to_sn1mg")
            .formula <- "M_2-MG <=> M_1-MG"
        
        ## nae_to_fa
        if (reaction %in% c("RHEA:17505", "RHEA:39995"))
            .formula <- "M_H2O + M_NAE = M_FA + M_Ethanolamine"
        if (reaction %in% c("RHEA:17506", "RHEA:39996"))
            .formula <- "M_H2O + M_NAE => M_FA + M_Ethanolamine"
        if (reaction %in% c("RHEA:17507", "RHEA:39997"))
            .formula <- "M_H2O + M_NAE <= M_FA + M_Ethanolamine"
        if (reaction %in% c("RHEA:17508", "RHEA:39998"))
            .formula <- "M_H2O + M_NAE <=> M_FA + M_Ethanolamine"
        
        ## nape_to_lnape
        if (reaction == "RHEA:45460")
            .formula <- "M_H2O + M_NAPE = M_FA + M_H+ + M_LNAPE"
        if (reaction == "RHEA:45461")
            .formula <- "M_H2O + M_NAPE => M_FA + M_H+ + M_LNAPE"
        if (reaction == "RHEA:45462")
            .formula <- "M_H2O + M_NAPE <= M_FA + M_H+ + M_LNAPE"
        if (reaction == "RHEA:45463")
            .formula <- "M_H2O + M_NAPE <=> M_FA + M_H+ + M_LNAPE"
        
        ## nape_to_nae
        if (reaction == "RHEA:33159")
            .formula <- "M_H2O + M_NAPE = M_PA + M_NAE + M_H+"
        if (reaction == "RHEA:33160")
            .formula <- "M_H2O + M_NAPE => M_PA + M_NAE + M_H+"
        if (reaction == "RHEA:33161")
            .formula <- "M_H2O + M_NAPE <= M_PA + M_NAE + M_H+"
        if (reaction == "RHEA:33162")
            .formula <- "M_H2O + M_NAPE <=> M_PA + M_NAE + M_H+"
        
        ## nape_to_pnae
        if (reaction == "nape_to_pnae")
            .formula <- "M_NAPE + M_H2O <=> M_PNAE + M_1,2-DG"
        
        ## napeo_to_nae
        if (reaction == "napeo_to_nae")
            .formula <- "M_NAPEO + M_H2O <=> M_NAE + M_PA-O"
        
        ## pa_to_cdpdg
        if (reaction == "RHEA:16229")
            .formula <- "M_PA + M_CTP + M_H+ = M_CDP-DG + M_PPi"
        if (reaction == "RHEA:16230")
            .formula <- "M_PA + M_CTP + M_H+ => M_CDP-DG + M_PPi"
        if (reaction == "RHEA:16231")
            .formula <- "M_PA + M_CTP + M_H+ <= M_CDP-DG + M_PPi"
        if (reaction == "RHEA:16232")
            .formula <- "M_PA + M_CTP + M_H+ <=> M_CDP-DG + M_PPi"
        
        ## pa_to_dg
        if (reaction == "RHEA:27429")
            .formula <- "M_PA + M_H2O = M_1,2-DG + M_Pi"
        if (reaction == "RHEA:27430")
            .formula <- "M_PA + M_H2O => M_1,2-DG + M_Pi"
        if (reaction == "RHEA:27431")
            .formula <- "M_PA + M_H2O <= M_1,2-DG + M_Pi"
        if (reaction == "RHEA:27432")
            .formula <- "M_PA + M_H2O <=> M_1,2-DG + M_Pi"
        
        ## pao_to_dgo
        if (reaction == "RHEA:36239")
            .formula <- "M_PA-O + M_H2O = M_DG-O + M_Pi"
        if (reaction == "RHEA:36240")
            .formula <- "M_PA-O + M_H2O => M_DG-O + M_Pi"
        if (reaction == "RHEA:36241")
            .formula <- "M_PA-O + M_H2O <= M_DG-O + M_Pi"
        if (reaction == "RHEA:36242")
            .formula <- "M_PA-O + M_H2O <=> M_DG-O + M_Pi"
        
        ## pc_to_dg
        if (reaction == "RHEA:10604")
            .formula <- "M_PC + M_H2O = M_1,2-DG + M_H+ + M_Phosphocholine"
        if (reaction == "RHEA:10605")
            .formula <- "M_PC + M_H2O => M_1,2-DG + M_H+ + M_Phosphocholine"
        if (reaction == "RHEA:10606")
            .formula <- "M_PC + M_H2O <= M_1,2-DG + M_H+ + M_Phosphocholine"
        if (reaction == "RHEA:10607")
            .formula <- "M_PC + M_H2O <=> M_1,2-DG + M_H+ + M_Phosphocholine"
        
        ## pc_to_sn1lpc
        if (reaction == "RHEA:15801")
            .formula <- "M_PC + M_H2O = M_1-LPC + M_FA + M_H+"
        if (reaction == "RHEA:15802")
            .formula <- "M_PC + M_H2O => M_1-LPC + M_FA + M_H+"
        if (reaction == "RHEA:15803")
            .formula <- "M_PC + M_H2O <= M_1-LPC + M_FA + M_H+"
        if (reaction == "RHEA:15804")
            .formula <- "M_PC + M_H2O <=> M_1-LPC + M_FA + M_H+"
        
        ## pc_to_sn2lpc
        if (reaction == "RHEA:18689") 
            .formula <- "M_PC + M_H2O = M_2-LPC + M_FA + M_H+"
        if (reaction == "RHEA:18690") 
            .formula <- "M_PC + M_H2O => M_2-LPC + M_FA + M_H+"
        if (reaction == "RHEA:18691") 
            .formula <- "M_PC + M_H2O <= M_2-LPC + M_FA + M_H+"
        if (reaction == "RHEA:18692") 
            .formula <- "M_PC + M_H2O <=> M_2-LPC + M_FA + M_H+"
        
        ## pc_to_pa
        if (reaction == "RHEA:14445")
            .formula <- "M_PC + M_H2O = M_PA + M_Choline + M_H+"
        if (reaction == "RHEA:14446")
            .formula <- "M_PC + M_H2O => M_PA + M_Choline + M_H+"
        if (reaction == "RHEA:14447")
            .formula <- "M_PC + M_H2O <= M_PA + M_Choline + M_H+"
        if (reaction == "RHEA:14448")
            .formula <- "M_PC + M_H2O <=> M_PA + M_Choline + M_H+"
        
        ## pc_to_ps
        if (reaction == "RHEA:45088") 
            .formula <- "M_PC + M_L-Serine = M_PS + M_Choline"
        if (reaction == "RHEA:45089") 
            .formula <- "M_PC + M_L-Serine => M_PS + M_Choline"
        if (reaction == "RHEA:45090") 
            .formula <- "M_PC + M_L-Serine <= M_PS + M_Choline"
        if (reaction == "RHEA:45091")
            .formula <- "M_PC + M_L-Serine <=> M_PS + M_Choline"

        ## pco_to_lpco
        if (reaction == "RHEA:36231")
            .formula <- "M_PC-O + M_H2O = M_LPC-O + M_FA + M_H+"
        if (reaction == "RHEA:36232")
            .formula <- "M_PC-O + M_H2O => M_LPC-O + M_FA + M_H+"
        if (reaction == "RHEA:36233")
            .formula <- "M_PC-O + M_H2O <= M_LPC-O + M_FA + M_H+"
        if (reaction == "RHEA:36234")
            .formula <- "M_PC-O + M_H2O <=> M_LPC-O + M_FA + M_H+"
        
        ## lpco_to_lpao
        if (reaction == "RHEA:39927")
            .formula <- "M_1-LPC-O + M_H2O = M_1-LPA-O + M_Choline + M_H+"
        if (reaction == "RHEA:39928")
            .formula <- "M_1-LPC-O + M_H2O => M_1-LPA-O + M_Choline + M_H+"
        if (reaction == "RHEA:39929")
            .formula <- "M_1-LPC-O + M_H2O <= M_1-LPA-O + M_Choline + M_H+"
        if (reaction == "RHEA:39930")
            .formula <- "M_1-LPC-O + M_H2O <=> M_1-LPA-O + M_Choline + M_H+"
        
        ## lpco_to_mgo
        if (reaction == "RHEA:36083")
            .formula <- "M_1-LPC-O + M_H2O = M_1-MG-O + M_H+ + M_Phosphocholine"
        if (reaction == "RHEA:36084")
            .formula <- "M_1-LPC-O + M_H2O => M_1-MG-O + M_H+ + M_Phosphocholine"
        if (reaction == "RHEA:36085")
            .formula <- "M_1-LPC-O + M_H2O <= M_1-MG-O + M_H+ + M_Phosphocholine"
        if (reaction == "RHEA:36086")
            .formula <- "M_1-LPC-O + M_H2O <=> M_1-MG-O + M_H+ + M_Phosphocholine"
        
        ## lpco_to_pco
        if (reaction == "RHEA:23992")
            .formula <- "M_1-LPC-O + M_AcylCoA = M_PC-O + M_CoA"
        if (reaction == "RHEA:23993")
            .formula <- "M_1-LPC-O + M_AcylCoA => M_PC-O + M_CoA"
        if (reaction == "RHEA:23994")
            .formula <- "M_1-LPC-O + M_AcylCoA <= M_PC-O + M_CoA"
        if (reaction == "RHEA:23995")
            .formula <- "M_1-LPC-O + M_AcylCoA <=> M_PC-O + M_CoA"
        
        ## pe_to_dg
        if (reaction == "RHEA:78951")
            .formula <- "M_H2O + M_PE = M_P-Ethanolamine + M_1,2-DG + M_H+"
        if (reaction == "RHEA:78952")
            .formula <- "M_H2O + M_PE => M_P-Ethanolamine + M_1,2-DG + M_H+"
        if (reaction == "RHEA:78953")
            .formula <- "M_H2O + M_PE <= M_P-Ethanolamine + M_1,2-DG + M_H+"
        if (reaction == "RHEA:78954")
            .formula <- "M_H2O + M_PE <=> M_P-Ethanolamine + M_1,2-DG + M_H+"
        
        ## pe_to_sn1lpe
        if (reaction == "RHEA:44604")
            .formula <- "M_PE + M_H2O = M_1-LPE + M_FA + M_H+"
        if (reaction == "RHEA:44605")
            .formula <- "M_PE + M_H2O => M_1-LPE + M_FA + M_H+"
        if (reaction == "RHEA:44606")
            .formula <- "M_PE + M_H2O <= M_1-LPE + M_FA + M_H+"
        if (reaction == "RHEA:44607")
            .formula <- "M_PE + M_H2O <=> M_1-LPE + M_FA + M_H+"
        
        ## pe_to_sn2lpe
        if (reaction == "RHEA:44408") 
            .formula <- "M_PE + M_H2O = M_2-LPE + M_FA + M_H+"
        if (reaction == "RHEA:44409") 
            .formula <- "M_PE + M_H2O => M_2-LPE + M_FA + M_H+"
        if (reaction == "RHEA:44410") 
            .formula <- "M_PE + M_H2O <= M_2-LPE + M_FA + M_H+"
        if (reaction == "RHEA:44411") 
            .formula <- "M_PE + M_H2O <=> M_2-LPE + M_FA + M_H+"
        
        ## pe_to_nape_sn1
        if (reaction == "RHEA:45188")
            .formula <- "M_PC + M_PE = M_2-LPC + M_H+ + M_NAPE"
        if (reaction == "RHEA:45189")
            .formula <- "M_PC + M_PE => M_2-LPC + M_H+ + M_NAPE"
        if (reaction == "RHEA:45190")
            .formula <- "M_PC + M_PE <= M_2-LPC + M_H+ + M_NAPE"
        if (reaction == "RHEA:45191")
            .formula <- "M_PC + M_PE <=> M_2-LPC + M_H+ + M_NAPE"
        
        ## pe_to_nape_sn2
        if (reaction == "RHEA:45192")
            .formula <- "M_PC + M_PE = M_1-LPC + M_H+ + M_NAPE"
        if (reaction == "RHEA:45193")
            .formula <- "M_PC + M_PE => M_1-LPC + M_H+ + M_NAPE"
        if (reaction == "RHEA:45194")
            .formula <- "M_PC + M_PE <= M_1-LPC + M_H+ + M_NAPE"
        if (reaction == "RHEA:45195")
            .formula <- "M_PC + M_PE <=> M_1-LPC + M_H+ + M_NAPE"
        
        ## pe_to_pa
        if (reaction == "pe_to_pa")
            .formula <- "M_PE + M_H2O <=> M_PA + M_Ethanolamine + M_H+"
        
        ## pe_to_ps
        if (reaction == "RHEA:27606")
            .formula <- "M_PE + M_L-Serine = M_PS + M_Ethanolamine"
        if (reaction == "RHEA:27607")
            .formula <- "M_PE + M_L-Serine => M_PS + M_Ethanolamine"
        if (reaction == "RHEA:27608")
            .formula <- "M_PE + M_L-Serine <= M_PS + M_Ethanolamine"
        if (reaction == "RHEA:27609")
            .formula <- "M_PE + M_L-Serine <=> M_PS + M_Ethanolamine"
        
        ## peo_to_lpeo
        if (reaction == "peo_to_lpeo")
            .formula <- "M_PE-O + M_H2O <=> M_LPE-O + M_FA + M_H+"
        
        ## peo_to_napeo_sn1
        if (reaction == "peo_to_napeo_sn1")
            .formula <- "M_PE-O + M_PC <=> M_NAPEO + M_2-LPC"
        
        ## peo_to_napeo_sn2
        if (reaction == "peo_to_napeo_sn2")
            .formula <- "M_PE-O + M_PC <=> M_NAPEO + M_1-LPC"
        
        ## peo_to_pep
        if (reaction == "RHEA:22956")
            .formula <- "M_PE-O + M_Fe2+-cytochrome_b5 + 2 M_H+ + M_O2 = M_PE-P + M_Fe3+-cytochrome_b5 + 2 M_H2O"
        if (reaction == "RHEA:22957")
            .formula <- "M_PE-O + M_Fe2+-cytochrome_b5 + 2 M_H+ + M_O2 => M_PE-P + M_Fe3+-cytochrome_b5 + 2 M_H2O"
        if (reaction == "RHEA:22958")
            .formula <- "M_PE-O + M_Fe2+-cytochrome_b5 + 2 M_H+ + M_O2 <= M_PE-P + M_Fe3+-cytochrome_b5 + 2 M_H2O"
        if (reaction == "RHEA:22959")
            .formula <- "M_PE-O + M_Fe2+-cytochrome_b5 + 2 M_H+ + M_O2 <=> M_PE-P + M_Fe3+-cytochrome_b5 + 2 M_H2O"
        
        ## pep_to_lpep
        if (reaction == "RHEA:36195")
            .formula <- "M_PE-P + M_H2O = M_LPE-P + M_FA + M_H+"
        if (reaction == "RHEA:36196")
            .formula <- "M_PE-P + M_H2O => M_LPE-P + M_FA + M_H+"
        if (reaction == "RHEA:36197")
            .formula <- "M_PE-P + M_H2O <= M_LPE-P + M_FA + M_H+"
        if (reaction == "RHEA:36198")
            .formula <- "M_PE-P + M_H2O <=> M_LPE-P + M_FA + M_H+"
        
        ## lpep_to_fal
        if (reaction == "RHEA:16905")
            .formula <- "M_1-LPE-P + M_H2O = M_FAL + M_Glycerophosphoethanolamine"
        if (reaction == "RHEA:16906")
            .formula <- "M_1-LPE-P + M_H2O => M_FAL + M_Glycerophosphoethanolamine"
        if (reaction == "RHEA:16907")
            .formula <- "M_1-LPE-P + M_H2O <= M_FAL + M_Glycerophosphoethanolamine"
        if (reaction == "RHEA:16908")
            .formula <- "M_1-LPE-P + M_H2O <=> M_FAL + M_Glycerophosphoethanolamine"
        
        ## lpep_to_lpap
        if (reaction == "RHEA:36203")
            .formula <- "M_1-LPE-P + M_H2O = M_LPA-P + M_Ethanolamine + M_H+"
        if (reaction == "RHEA:36204")
            .formula <- "M_1-LPE-P + M_H2O => M_LPA-P + M_Ethanolamine + M_H+"
        if (reaction == "RHEA:36205")
            .formula <- "M_1-LPE-P + M_H2O <= M_LPA-P + M_Ethanolamine + M_H+"
        if (reaction == "RHEA:36206")
            .formula <- "M_1-LPE-P + M_H2O <=> M_LPA-P + M_Ethanolamine + M_H+"
        
        ## lpep_to_mgp
        if (reaction == "RHEA:36199")
            .formula <- "M_1-LPE-P + M_H2O = M_1-MG-P + M_H+ + M_Phosphoethanolamine"
        if (reaction == "RHEA:36200")
            .formula <- "M_1-LPE-P + M_H2O => M_1-MG-P + M_H+ + M_Phosphoethanolamine"
        if (reaction == "RHEA:36201")
            .formula <- "M_1-LPE-P + M_H2O <= M_1-MG-P + M_H+ + M_Phosphoethanolamine"
        if (reaction == "RHEA:36202")
            .formula <- "M_1-LPE-P + M_H2O <=> M_1-MG-P + M_H+ + M_Phosphoethanolamine"
        
        ## pep_to_napep_sn1
        if (reaction == "RHEA:63596")
            .formula <- "M_PE-P + M_PC = M_2-LPC + M_H+ + M_NAPEP"
        if (reaction == "RHEA:63597")
            .formula <- "M_PE-P + M_PC => M_2-LPC + M_H+ + M_NAPEP"
        if (reaction == "RHEA:63598")
            .formula <- "M_PE-P + M_PC <= M_2-LPC + M_H+ + M_NAPEP"
        if (reaction == "RHEA:63599")
            .formula <- "M_PE-P + M_PC <=> M_2-LPC + M_H+ + M_NAPEP"
        
        ## pep_to_napep_sn2
        if (reaction == "pep_to_napep_sn2")
            .formula <- "M_PE-P + M_PC <=> M_NAPEP + M_1-LPC"
        
        ## pg_to_cl
        if (reaction == "RHEA:32931")
            .formula <- "M_PG + M_CDP-DG = M_CL + M_CMP + M_H+"
        if (reaction == "RHEA:32932")
            .formula <- "M_PG + M_CDP-DG => M_CL + M_CMP + M_H+"
        if (reaction == "RHEA:32933")
            .formula <- "M_PG + M_CDP-DG <= M_CL + M_CMP + M_H+"
        if (reaction == "RHEA:32934")
            .formula <- "M_PG + M_CDP-DG <=> M_CL + M_CMP + M_H+"
        
        ## pgp_to_pg
        if (reaction == "RHEA:33751")
            .formula <- "M_PGP + M_H2O = M_PG + M_Pi"
        if (reaction == "RHEA:33752")
            .formula <- "M_PGP + M_H2O => M_PG + M_Pi"
        if (reaction == "RHEA:33753")
            .formula <- "M_PGP + M_H2O <= M_PG + M_Pi"
        if (reaction == "RHEA:33754")
            .formula <- "M_PGP + M_H2O <=> M_PG + M_Pi"
        
        ## pi_to_dg
        if (reaction == "RHEA:43484")
            .formula <- "M_PI + M_H2O = M_myo-Inositol-1-P + M_1,2-DG + M_H+"
        if (reaction == "RHEA:43485")
            .formula <- "M_PI + M_H2O => M_myo-Inositol-1-P + M_1,2-DG + M_H+"
        if (reaction == "RHEA:43486")
            .formula <- "M_PI + M_H2O <= M_myo-Inositol-1-P + M_1,2-DG + M_H+"
        if (reaction == "RHEA:43487")
            .formula <- "M_PI + M_H2O <=> M_myo-Inositol-1-P + M_1,2-DG + M_H+"
        
        ## pi_to_sn1lpi
        if (reaction == "RHEA:18001")
            .formula <- "M_PI + M_H2O = M_1-LPI + M_FA + M_H+"
        if (reaction == "RHEA:18002")
            .formula <- "M_PI + M_H2O => M_1-LPI + M_FA + M_H+"
        if (reaction == "RHEA:18003")
            .formula <- "M_PI + M_H2O <= M_1-LPI + M_FA + M_H+"
        if (reaction == "RHEA:18004")
            .formula <- "M_PI + M_H2O <=> M_1-LPI + M_FA + M_H+"
            
        ## ps_to_pe
        if (reaction == "RHEA:20828")
            .formula <- "M_PS + M_H+ = M_PE + M_CO2"
        if (reaction == "RHEA:20829")
            .formula <- "M_PS + M_H+ => M_PE + M_CO2"
        if (reaction == "RHEA:20830")
            .formula <- "M_PS + M_H+ <= M_PE + M_CO2"
        if (reaction == "RHEA:20831")
            .formula <- "M_PS + M_H+ <=> M_PE + M_CO2"
        
        ## sm_to_cer
        if (reaction == "RHEA:45644")
            .formula <- "M_H2O + M_SM = M_H+ + M_Cer + M_Phosphocholine"
        if (reaction == "RHEA:45645")
            .formula <- "M_H2O + M_SM => M_H+ + M_Cer + M_Phosphocholine"
        if (reaction == "RHEA:45646")
            .formula <- "M_H2O + M_SM <= M_H+ + M_Cer + M_Phosphocholine"
        if (reaction == "RHEA:45647")
            .formula <- "M_H2O + M_SM <=> M_H+ + M_Cer + M_Phosphocholine"
        
        ## sphinga_to_dhcer
        if (reaction == "RHEA:53424")
            .formula <- "M_AcylCoA + M_Sphinganine = M_DhCer + M_CoA + M_H+"
        if (reaction == "RHEA:53425")
            .formula <- "M_AcylCoA + M_Sphinganine => M_DhCer + M_CoA + M_H+"
        if (reaction == "RHEA:53426")
            .formula <- "M_AcylCoA + M_Sphinganine <= M_DhCer + M_CoA + M_H+"
        if (reaction == "RHEA:53427")
            .formula <- "M_AcylCoA + M_Sphinganine <=> M_DhCer + M_CoA + M_H+"
        
        ## tg_to_dg
        if (reaction %in% c("RHEA:33271", "RHEA:44864"))
            .formula <- "M_TG + M_H2O = M_1,2-DG + M_FA + M_H+"
        if (reaction %in% c("RHEA:33272", "RHEA:44865"))
            .formula <- "M_TG + M_H2O => M_1,2-DG + M_FA + M_H+"
        if (reaction %in% c("RHEA:33273", "RHEA:44866"))
            .formula <- "M_TG + M_H2O <= M_1,2-DG + M_FA + M_H+"
        if (reaction %in% c("RHEA:33274", "RHEA:44867"))
            .formula <- "M_TG + M_H2O <=> M_1,2-DG + M_FA + M_H+"
            
        template$reaction_formula <- .formula
    }
    
    ## fill using reaction_formula
    ## split the template$reaction_formula into substrates and products 
    .formula_tmp <- template$reaction_formula
    .formula_tmp <- stringi::stri_replace_all_regex(str = .formula_tmp, 
        pattern = ">|<", replacement = "")
    .formula_tmp <- stringi::stri_split(str = .formula_tmp, fixed = "=")[[1]]
    .formula_substrate <- .formula_tmp[1]
    .formula_product <- .formula_tmp[2]
    
    ## remove the "+" and the trailing spaces for substrates and products
    .formula_substrate <- stringi::stri_split(str = .formula_substrate, 
        regex = "[ ][+][ ]")[[1]]
    .formula_substrate <- stringi::stri_replace_all_regex(str = .formula_substrate, 
        pattern = "^ | $", replacement = "")
    .formula_product <- stringi::stri_split(str = .formula_product, 
        regex = "[ ][+][ ]")[[1]]
    .formula_product <- stringi::stri_replace_all_regex(str = .formula_product, 
        pattern = "^ | $", replacement = "")

    ## write the substrates and products to the template
    template$reaction_substrate <- .formula_substrate
    template$reaction_product <- .formula_product
    
    ## fill using reaction_formula_chebi
    ## add ChEBI ids using the mapping table
    .mappingTable <- mappingTable()
    inds <- which(.mappingTable == reaction, arr.ind = TRUE)[, "row"]
    .mappingTable <- .mappingTable[inds, ]
    
    ## find the corresponding CHEBIs for substrates
    substrate_tmp <- stringi::stri_split(template$reaction_substrate, regex = "[ ]")
    .formula_substrate <- lapply(substrate_tmp, function(substrate_tmp_i) {
        inds_match <- match(substrate_tmp_i, .mappingTable[, 1])
        res <- .mappingTable[inds_match, 3, drop = FALSE]
        res[is.na(inds_match)] <- substrate_tmp_i[is.na(inds_match)]
        paste(res, collapse = " ")
    }) |>
        unlist()
    
    ## find the corresponding CHEBIs for products
    product_tmp <- stringi::stri_split(template$reaction_product, regex = "[ ]")
    .formula_product <- lapply(product_tmp, function(product_tmp_i) {
        inds_match <- match(product_tmp_i, .mappingTable[, 1])
        res <- .mappingTable[inds_match, 3, drop = FALSE]
        res[is.na(inds_match)] <- product_tmp_i[is.na(inds_match)]
        paste(res, collapse = " ")
    }) |>
        unlist()
    
    ## determine the reaction direction
    patterns <- c("<=>", "<=", "=>", "=")
    matches <- sapply(patterns, function(p) grepl(p, template$reaction_formula))
    reaction_type <- names(which(matches)[1])
    
    ## assemble the reaction_formula_chebi
    template$reaction_formula_chebi <- paste(
        paste(.formula_substrate, collapse = " + "), reaction_type,
        paste(.formula_product, collapse = " + "), collapse = " ")
    
    ## write the substrates and products to the template
    template$reaction_substrate_chebi <- .formula_substrate
    template$reaction_product_chebi <- .formula_product
    
    ## return the template object
    template
}
