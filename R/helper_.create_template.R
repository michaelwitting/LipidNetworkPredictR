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
#' @param template NA or \code{list}
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
#' LipidNetworkPredictR:::.create_template(template = NA, reaction = "RHEA:15421")
.create_template <- function(template = NA, reaction = "RHEA:15421") {
    
    if (is.list(template)) {
        ## check integrity (all names are present in template)
        if (!"reaction_name" %in% names(template)) 
            stop("'template' has to contain the entry 'reaction_name'")
        if (!"reaction_formula" %in% names(template)) 
            stop("'template' has to contain the entry 'reaction_formula'")
        if (!"reaction_isReversible" %in% names(template)) 
            stop("'template' has to contain the entry 'reacton_isReversible'")
        if (!"reaction_geneAssociation" %in% names(template)) 
            stop("'template' has to contain the entry 'reaction_geneAssociation'")
        if (!"reaction_pathway" %in% names(template)) 
            stop("'template' has to contain the entry 'reaction_pathway'")
        
        ## all entries are of mode character
        if (!is.vector(template$reaction_name, "character"))
            stop("entry 'reaction_name' is not of mode 'character'")
        if (!is.vector(template$reaction_formula, "character"))
            stop("entry 'reaction_formula' is not of mode 'character'")
        if (!is.vector(template$reaction_isReversible, "character"))
            stop("entry 'reaction_isReversible' is not of mode 'character'")
        if (!is.vector(template$reaction_geneAssociation, "character"))
            stop("entry 'reaction_geneAssociation' is not of mode 'character'")
        if (!is.vector(template$reaction_pathway, "character"))
            stop("entry 'reaction_pathway' is not of mode 'character'")
            
    } else {
        ## this event happens if the param template is not a list
        
        ## check if reaction is in the parameter space
        .check_reaction(reaction)
        
        ## create a template object (list)
        template <- list(reaction_name = "", reaction_formula = "",
            reaction_RHEA = reaction, reaction_isReversible = "",
            reaction_geneAssociation = "", reaction_pathway = "")
        
        ## depending on the reaction, write specific reaction .formula to the
        ## entry reaction_.formula
        
        ## acyldhap_to_alkyldhap
        if (reaction == "RHEA:36171") 
            .formula <- "M_AcylDHAP + M_FAO = M_H+ + M_FA + M_AlkylDHAP"
        if (reaction == "RHEA:36172") 
            .formula <- "M_AcylDHAP + M_FAO => M_H+ + M_FA + M_AlkylDHAP"
        if (reaction == "RHEA:36173") 
            .formula <- "M_AcylDHAP + M_FAO <= M_H+ + M_FA + M_AlkylDHAP"
        if (reaction == "RHEA:36174") 
            .formula <- "M_AcylDHAP + M_FAO <=> M_H+ + M_FA + M_AlkylDHAP"
        
        ## alkyldhap_to_lpao
        if (reaction == "RHEA:36175") 
            .formula <- "M_H+ + M_NADPH + M_AlkylDHAP = M_LPA-O + M_NADP"
        if (reaction == "RHEA:36176") 
            .formula <- "M_H+ + M_NADPH + M_AlkylDHAP => M_LPA-O + M_NADP"
        if (reaction == "RHEA:36177") 
            .formula <- "M_H+ + M_NADPH + M_AlkylDHAP <= M_LPA-O + M_NADP"
        if (reaction == "RHEA:36178") 
            .formula <- "M_H+ + M_NADPH + M_AlkylDHAP <=> M_LPA-O + M_NADP"
        
        if (reaction == "cerp_to_cer")
            .formula <- "M_H2O + M_CerP <=> M_Pi + M_Cer"
        
        ## cdpdg_to_pgp
        if (reaction == "RHEA:12593") 
            .formula <- "M_Glycerol-3-P + M_CDP-DG = M_H+ + M_CMP + M_PGP"
        if (reaction == "RHEA:12594") 
            .formula <- "M_Glycerol-3-P + M_CDP-DG => M_H+ + M_CMP + M_PGP"
        if (reaction == "RHEA:12595") 
            .formula <- "M_Glycerol-3-P + M_CDP-DG <= M_H+ + M_CMP + M_PGP"
        if (reaction == "RHEA:12596") 
            .formula <- "M_Glycerol-3-P + M_CDP-DG <=> M_H+ + M_CMP + M_PGP"
        
        ## cdpdg_to_pi
        if (reaction == "RHEA:11580") 
            .formula <- "M_myo-Inositol + M_CDP-DG = M_H+ + M_CMP + M_PI"
        if (reaction == "RHEA:11581") 
            .formula <- "M_myo-Inositol + M_CDP-DG => M_H+ + M_CMP + M_PI"
        if (reaction == "RHEA:11582") 
            .formula <- "M_myo-Inositol + M_CDP-DG <= M_H+ + M_CMP + M_PI"
        if (reaction == "RHEA:11583") 
            .formula <- "M_myo-Inositol + M_CDP-DG <=> M_H+ + M_CMP + M_PI"
        
        ## cer_to_cerp
        if (reaction == "RHEA:17929")
            .formula <- "M_ATP + M_Cer <=> M_H+ + M_ADP + M_CerP"
        
        ## cer_to_glccer
        if (reaction == "RHEA:12088")
            .formula <- "M_UDP-Glucose + M_Cer <=> M_H+ + M_UDP + M_GlcCer"
        
        ## cer_to_sm
        if (reaction == "RHEA:18765") 
            .formula <- "M_PC + M_Cer <=> M_1,2-DG + M_SM"
        
        ## cl_to_lcl
        if (reaction == "RHEA:32935")
            .formula <- "M_H2O + M_CL = M_1,2,4-LCL + M_H+ + M_FA"
        if (reaction == "RHEA:32936")
            .formula <- "M_H2O + M_CL => M_1,2,4-LCL + M_H+ + M_FA"
        if (reaction == "RHEA:32937")
            .formula <- "M_H2O + M_CL <= M_1,2,4-LCL + M_H+ + M_FA"
        if (reaction == "RHEA:32938")
            .formula <- "M_H2O + M_CL <=> M_1,2,4-LCL + M_H+ + M_FA"
        
        ## coa_to_acdhap
        if (reaction == "RHEA:17657")
            .formula <- "M_Dihydroxyacetone-P + M_AcylCoA = M_CoA + M_AcylDHAP"
        if (reaction == "RHEA:17658")
            .formula <- "M_Dihydroxyacetone-P + M_AcylCoA => M_CoA + M_AcylDHAP"
        if (reaction == "RHEA:17659")
            .formula <- "M_Dihydroxyacetone-P + M_AcylCoA <= M_CoA + M_AcylDHAP"
        if (reaction == "RHEA:17660")
            .formula <- "M_Dihydroxyacetone-P + M_AcylCoA <=> M_CoA + M_AcylDHAP"
    
        ## coa_to_fao
        if (reaction == "RHEA:52716") 
            .formula <- "M_AcylCoA + 2 M_NADPH + 2 M_H+ = M_FAO + 2 M_NADP + M_CoA"
        if (reaction == "RHEA:52717") 
            .formula <- "M_AcylCoA + 2 M_NADPH + 2 M_H+ => M_FAO + 2 M_NADP + M_CoA"
        if (reaction == "RHEA:52718") 
            .formula <- "M_AcylCoA + 2 M_NADPH + 2 M_H+ <= M_FAO + 2 M_NADP + M_CoA"
        if (reaction == "RHEA:52719") 
            .formula <- "M_AcylCoA + 2 M_NADPH + 2 M_H+ <=> M_FAO + 2 M_NADP + M_CoA"
        
        ## coa_to_lpa
        if (reaction == "RHEA:15325") 
            .formula <- "M_Glycerol-3-P + M_AcylCoA = M_CoA + M_LPA"
        if (reaction == "RHEA:15326") 
            .formula <- "M_Glycerol-3-P + M_AcylCoA => M_CoA + M_LPA"
        if (reaction == "RHEA:15327") 
            .formula <- "M_Glycerol-3-P + M_AcylCoA <= M_CoA + M_LPA"
        if (reaction == "RHEA:15328")
            .formula <- "M_Glycerol-3-P + M_AcylCoA <=> M_CoA + M_LPA"

        ## dg_to_sn1mg
        if (reaction %in% c("RHEA:44712", "RHEA:35663"))
            .formula <- "M_H2O + M_1,2-DG = M_H+ + M_1-MG + M_FA"
        if (reaction %in% c("RHEA:44713", "RHEA:35664"))
            .formula <- "M_H2O + M_1,2-DG => M_H+ + M_1-MG + M_FA"
        if (reaction %in% c("RHEA:44714", "RHEA:35665"))
            .formula <- "M_H2O + M_1,2-DG <= M_H+ + M_1-MG + M_FA"
        if (reaction %in% c("RHEA:44715", "RHEA:35666"))
            .formula <- "M_H2O + M_1,2-DG <=> M_H+ + M_1-MG + M_FA"

        ## dg_to_sn2mg
        if (reaction == "RHEA:33275")
            .formula <- "M_H2O + M_1,2-DG = M_H+ + M_2-MG + M_FA"
        if (reaction == "RHEA:33276")
            .formula <- "M_H2O + M_1,2-DG => M_H+ + M_2-MG + M_FA"
        if (reaction == "RHEA:33277")
            .formula <- "M_H2O + M_1,2-DG <= M_H+ + M_2-MG + M_FA"
        if (reaction == "RHEA:33278")
            .formula <- "M_H2O + M_1,2-DG <=> M_H+ + M_2-MG + M_FA"
        
        ## dg_to_pa
        if (reaction == "RHEA:10272")
            .formula <- "M_ATP + M_1,2-DG = M_H+ + M_ADP + M_PA"
        if (reaction == "RHEA:10273")
            .formula <- "M_ATP + M_1,2-DG => M_H+ + M_ADP + M_PA"
        if (reaction == "RHEA:10274")
            .formula <- "M_ATP + M_1,2-DG <= M_H+ + M_ADP + M_PA"
        if (reaction == "RHEA:10275")
            .formula <- "M_ATP + M_1,2-DG <=> M_H+ + M_ADP + M_PA"
        
        ## dg_to_pc
        if (reaction == "RHEA:32939") 
            .formula <- "M_CDP-Choline + M_1,2-DG = M_H+ + M_CMP + M_PC"
        if (reaction == "RHEA:32940") 
            .formula <- "M_CDP-Choline + M_1,2-DG => M_H+ + M_CMP + M_PC"
        if (reaction == "RHEA:32941") 
            .formula <- "M_CDP-Choline + M_1,2-DG <= M_H+ + M_CMP + M_PC"
        if (reaction == "RHEA:32942") 
            .formula <- "M_CDP-Choline + M_1,2-DG <=> M_H+ + M_CMP + M_PC"
        
        ## dg_to_pe
        if (reaction == "RHEA:32943") 
            .formula <- "M_CDP-Ethanolamine + M_1,2-DG = M_H+ + M_CMP + M_PE"
        if (reaction == "RHEA:32944") 
            .formula <- "M_CDP-Ethanolamine + M_1,2-DG => M_H+ + M_CMP + M_PE"
        if (reaction == "RHEA:32945") 
            .formula <- "M_CDP-Ethanolamine + M_1,2-DG <= M_H+ + M_CMP + M_PE"
        if (reaction == "RHEA:32946") 
            .formula <- "M_CDP-Ethanolamine + M_1,2-DG <=> M_H+ + M_CMP + M_PE"
        
        ## dg_to_tg
        if (reaction == "RHEA:10868") 
            .formula <- "M_AcylCoA + M_1,2-DG = M_CoA + M_TG"
        if (reaction == "RHEA:10869") 
            .formula <- "M_AcylCoA + M_1,2-DG => M_CoA + M_TG"
        if (reaction == "RHEA:10870") 
            .formula <- "M_AcylCoA + M_1,2-DG <= M_CoA + M_TG"
        if (reaction == "RHEA:10871") 
            .formula <- "M_AcylCoA + M_1,2-DG <=> M_CoA + M_TG"
        
        ## dgo_to_pco
        if (reaction == "RHEA:36179")
            .formula <- "M_CDP-Choline + M_DG-O = M_H+ + M_CMP + M_PC-O"
        if (reaction == "RHEA:36180")
            .formula <- "M_CDP-Choline + M_DG-O => M_H+ + M_CMP + M_PC-O"
        if (reaction == "RHEA:36181")
            .formula <- "M_CDP-Choline + M_DG-O <= M_H+ + M_CMP + M_PC-O"
        if (reaction == "RHEA:36182")
            .formula <- "M_CDP-Choline + M_DG-O <=> M_H+ + M_CMP + M_PC-O"
        
        ## dgo_to_peo
        if (reaction %in% c("RHEA:36187"))
            .formula <- "M_CDP-Ethanolamine + M_DG-O = M_H+ + M_CMP + M_PE-O"
        if (reaction %in% c("RHEA:36188"))
            .formula <- "M_CDP-Ethanolamine + M_DG-O => M_H+ + M_CMP + M_PE-O"
        if (reaction %in% c("RHEA:36189"))
            .formula <- "M_CDP-Ethanolamine + M_DG-O <= M_H+ + M_CMP + M_PE-O"
        if (reaction %in% c("RHEA:36190"))
            .formula <- "M_CDP-Ethanolamine + M_DG-O <=> M_H+ + M_CMP + M_PE-O"
        
        if (reaction == "dhcer_to_cer")
            .formula <- "M_H+ + M_NADH + M_O2 + M_DhCer <=> 2 M_H2O + M_NAD + M_Cer"
        
        ## dhcer_to_dhsm
        if (reaction == "RHEA:44620")
            .formula <- "M_PC + M_DhCer <=> M_1,2-DG + M_DhSM"
        
        ## dhsm_to_dhcer
        if (reaction == "RHEA:45300")
            .formula <- "M_H2O + M_DhSM = M_Phosphocholine + M_H+ + M_DhCer"
        if (reaction == "RHEA:45301")
            .formula <- "M_H2O + M_DhSM => M_Phosphocholine + M_H+ + M_DhCer"
        if (reaction == "RHEA:45302")
            .formula <- "M_H2O + M_DhSM <= M_Phosphocholine + M_H+ + M_DhCer"
        if (reaction == "RHEA:45303")
            .formula <- "M_H2O + M_DhSM <=> M_Phosphocholine + M_H+ + M_DhCer"
        
        ## fa_to_coa
        if (reaction %in% c("RHEA:15421", "RHEA:38883")) 
            .formula <- "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA"
        if (reaction %in% c("RHEA:15422", "RHEA:38884"))
            .formula <- "M_ATP + M_CoA + M_FA => M_PPi + M_AMP + M_AcylCoA"
        if (reaction %in% c("RHEA:15423", "RHEA:38885"))
            .formula <- "M_ATP + M_CoA + M_FA <= M_PPi + M_AMP + M_AcylCoA"
        if (reaction %in% c("RHEA:15424", "RHEA:38886"))
            .formula <- "M_ATP + M_CoA + M_FA <=> M_PPi + M_AMP + M_AcylCoA"
        
        ## lcl_to_cl
        if (reaction == "RHEA:35839")
            .formula <- "M_1,2,4-LCL + M_AcylCoA = M_CL + M_CoA"
        if (reaction == "RHEA:35840")
            .formula <- "M_1,2,4-LCL + M_AcylCoA => M_CL + M_CoA"
        if (reaction == "RHEA:35841")
            .formula <- "M_1,2,4-LCL + M_AcylCoA <= M_CL + M_CoA"
        if (reaction == "RHEA:35842")
            .formula <- "M_1,2,4-LCL + M_AcylCoA <=> M_CL + M_CoA"
        
        if (reaction == "RHEA:45420") ## lnape_to_gpnae
            .formula <- "M_LNAPE + M_H2O <=> M_GPNAE + M_FA"
        
        ## lpa_to_pa
        if (reaction == "RHEA:19709")
            .formula <- "M_LPA + M_AcylCoA = M_CoA + M_PA"
        if (reaction == "RHEA:19710")
            .formula <- "M_LPA + M_AcylCoA => M_CoA + M_PA"
        if (reaction == "RHEA:19711") 
            .formula <- "M_LPA + M_AcylCoA <= M_CoA + M_PA"
        if (reaction == "RHEA:19712") 
            .formula <- "M_LPA + M_AcylCoA <=> M_CoA + M_PA"
        
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
            .formula <- "M_H2O + M_1-LPC = M_Glycerophosphocholine + M_H+ + M_FA"
        if (reaction == "RHEA:15178")
            .formula <- "M_H2O + M_1-LPC => M_Glycerophosphocholine + M_H+ + M_FA"
        if (reaction == "RHEA:15179")
            .formula <- "M_H2O + M_1-LPC <= M_Glycerophosphocholine + M_H+ + M_FA"
        if (reaction == "RHEA:15180")
            .formula <- "M_H2O + M_1-LPC <=> M_Glycerophosphocholine + M_H+ + M_FA"
        
        ## sn2lpc_to_fa
        if (reaction == "RHEA:44696") 
            .formula <- "M_H2O + M_2-LPC = M_Glycerophosphocholine + M_H+ + M_FA"
        if (reaction == "RHEA:44697") 
            .formula <- "M_H2O + M_2-LPC => M_Glycerophosphocholine + M_H+ + M_FA"
        if (reaction == "RHEA:44698") 
            .formula <- "M_H2O + M_2-LPC <= M_Glycerophosphocholine + M_H+ + M_FA"
        if (reaction == "RHEA:44699") 
            .formula <- "M_H2O + M_2-LPC <=> M_Glycerophosphocholine + M_H+ + M_FA"
        
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
            .formula <- "1-LPC-O + AcylCoA = PC-O + CoA"
        if (reaction == "RHEA:23993") 
            .formula <- "1-LPC-O + AcylCoA => PC-O + CoA"
        if (reaction == "RHEA:23994") 
            .formula <- "1-LPC-O + AcylCoA <= PC-O + CoA"
        if (reaction == "RHEA:23995") 
            .formula <- "1-LPC-O + AcylCoA <=> PC-O + CoA"
        
        
        ## sn1lpe_to_fa
        if (reaction == "RHEA:32967") 
            .formula <- "M_H2O + M_1-LPE = M_H+ + M_FA + M_Glycerophosphoethanolamine"
        if (reaction == "RHEA:32968") 
            .formula <- "M_H2O + M_1-LPE => M_H+ + M_FA + M_Glycerophosphoethanolamine"
        if (reaction == "RHEA:32969") 
            .formula <- "M_H2O + M_1-LPE <= M_H+ + M_FA + M_Glycerophosphoethanolamine"
        if (reaction == "RHEA:32970") 
            .formula <- "M_H2O + M_1-LPE <=> M_H+ + M_FA + M_Glycerophosphoethanolamine"
        
        if (reaction == "sn2lpe_to_fa")
            .formula <- "M_H2O + M_2-LPE <=> M_H+ + M_FA + M_Glycerophosphoethanolamine"
        
        ## sn1lpe_to_pe
        if (reaction == "RHEA:32995") 
            .formula <- "M_AcylCoA + M_1-LPE = M_CoA + M_PE"
        if (reaction == "RHEA:32996") 
            .formula <- "M_AcylCoA + M_1-LPE => M_CoA + M_PE"
        if (reaction == "RHEA:32997") 
            .formula <- "M_AcylCoA + M_1-LPE <= M_CoA + M_PE"
        if (reaction == "RHEA:32998") 
            .formula <- "M_AcylCoA + M_1-LPE <=> M_CoA + M_PE"
        
        ## sn1lpi_to_pi
        if (reaction == "RHEA:33195")
            .formula <- "M_1-LPI + M_AcylCoA = M_PI + M_CoA"
        if (reaction == "RHEA:33196")
            .formula <- "M_1-LPI + M_AcylCoA => M_PI + M_CoA"
        if (reaction == "RHEA:33197")
            .formula <- "M_1-LPI + M_AcylCoA <= M_PI + M_CoA"
        if (reaction == "RHEA:33198")
            .formula <- "M_1-LPI + M_AcylCoA <=> M_PI + M_CoA"
        
        if (reaction == "lpeo_to_peo")
            .formula <- "M_AcylCoA + M_ak2lgpe <=> M_CoA + M_PE-O"
        
        ## lpep_to_pep
        if (reaction == "RHEA:16245")
            .formula <- "M_AcylCoA + M_LPE-P = M_CoA + M_PE-P"
        if (reaction == "RHEA:16246")    
            .formula <- "M_AcylCoA + M_LPE-P => M_CoA + M_PE-P"
        if (reaction == "RHEA:16247")
            .formula <- "M_AcylCoA + M_LPE-P <= M_CoA + M_PE-P"
        if (reaction == "RHEA:16248")
            .formula <- "M_AcylCoA + M_LPE-P <=> M_CoA + M_PE-P"
            
        
        ## sn1mg_to_dg
        if (reaction %in% c("RHEA:38463", "RHEA:39943"))
            .formula <- "M_1-MG + M_AcylCoA = M_CoA + M_1,2-DG"
        if (reaction %in% c("RHEA:38464", "RHEA:39944"))
            .formula <- "M_1-MG + M_AcylCoA => M_CoA + M_1,2-DG"
        if (reaction %in% c("RHEA:38465", "RHEA:39945"))
            .formula <- "M_1-MG + M_AcylCoA <= M_CoA + M_1,2-DG"
        if (reaction %in% c("RHEA:38466", "RHEA:39946"))
            .formula <- "M_1-MG + M_AcylCoA <=> M_CoA + M_1,2-DG"
        
        ## sn2mg_to_dg
        if (reaction %in% c("RHEA:32947", "RHEA:16741"))
            .formula <- "M_AcylCoA + M_1-MG = M_CoA + M_1,2-DG"
        if (reaction %in% c("RHEA:32948", "RHEA:16742"))
            .formula <- "M_AcylCoA + M_1-MG => M_CoA + M_1,2-DG"
        if (reaction %in% c("RHEA:32949", "RHEA:16743"))
            .formula <- "M_AcylCoA + M_1-MG <= M_CoA + M_1,2-DG"
        if (reaction %in% c("RHEA:32950", "RHEA:16744"))
            .formula <- "M_AcylCoA + M_1-MG <=> M_CoA + M_1,2-DG"
        
        ## sn1mg_to_fa
        if (reaction == "RHEA:34019")
            .formula <- "M_H2O + M_1-MG = M_Glycerol + M_H+ + M_FA"
        if (reaction == "RHEA:34020")
            .formula <- "M_H2O + M_1-MG => M_Glycerol + M_H+ + M_FA"
        if (reaction == "RHEA:34021")
            .formula <- "M_H2O + M_1-MG <= M_Glycerol + M_H+ + M_FA"
        if (reaction == "RHEA:34022")
            .formula <- "M_H2O + M_1-MG <=> M_Glycerol + M_H+ + M_FA"

        ## sn2mg_to_fa
        if (reaction %in% c("RHEA:32871"))
            .formula <-  "M_H2O + M_2-MG = M_Glycerol + M_H+ + M_FA"
        if (reaction %in% c("RHEA:32872"))
            .formula <-  "M_H2O + M_2-MG => M_Glycerol + M_H+ + M_FA"
        if (reaction %in% c("RHEA:32873"))
            .formula <-  "M_H2O + M_2-MG <= M_Glycerol + M_H+ + M_FA"
        if (reaction %in% c("RHEA:32874"))
            .formula <-  "M_H2O + M_2-MG <=> M_Glycerol + M_H+ + M_FA"
        
        ## sn1mg_to_lpa
        if (reaction == "RHEA:33747")
            .formula <- "M_ATP + M_1-MG = M_H+ + M_ADP + M_LPA"
        if (reaction == "RHEA:33748")
            .formula <- "M_ATP + M_1-MG => M_H+ + M_ADP + M_LPA"
        if (reaction == "RHEA:33749")
            .formula <- "M_ATP + M_1-MG <= M_H+ + M_ADP + M_LPA"
        if (reaction == "RHEA:33750")
            .formula <- "M_ATP + M_1-MG <=> M_H+ + M_ADP + M_LPA"
        
        if (reaction == "sn2mg_to_sn1mg")
            .formula <- "M_2-MG <=> M_1-MG"
        
        if (reaction == "nae_to_fa")
            .formula <- "M_H2O + M_NAE <=> M_Ethanolamine + M_H+ + M_FA"
        
        if (reaction == "nape_to_lnape")
            .formula <- "M_NAPE + M_H2O <=> M_LNAPE + M_FA"
        
        if (reaction == "nape_to_nae")
            .formula <- "M_NAPE + M_H2O <=> M_NAE + M_PA"
        
        if (reaction == "nape_to_pnae")
            .formula <- "M_NAPE + M_H2O <=> M_PNAE + M_1,2-DG"
        
        if (reaction == "napeo_to_nae")
            .formula <- "M_NAPEO + M_H2O <=> M_NAE + M_PA-O"
        
        ## pa_to_cdpdg
        if (reaction == "RHEA:16229")
            .formula <- "M_CTP + M_PA = M_PPi + M_CDP-DG"
        if (reaction == "RHEA:16230")
            .formula <- "M_CTP + M_PA => M_PPi + M_CDP-DG"
        if (reaction == "RHEA:16231")
            .formula <- "M_CTP + M_PA <= M_PPi + M_CDP-DG"
        if (reaction == "RHEA:16232")
            .formula <- "M_CTP + M_PA <=> M_PPi + M_CDP-DG"
        
        ## pa_to_dg
        if (reaction == "RHEA:27429")
            .formula <- "M_H2O + M_PA = M_Pi + M_1,2-DG"
        if (reaction == "RHEA:27430")
            .formula <- "M_H2O + M_PA => M_Pi + M_1,2-DG"
        if (reaction == "RHEA:27431")
            .formula <- "M_H2O + M_PA <= M_Pi + M_1,2-DG"
        if (reaction == "RHEA:27432")
            .formula <- "M_H2O + M_PA <=> M_Pi + M_1,2-DG"
        
        ## pao_to_dgo
        if (reaction == "RHEA:36239")
            .formula <- "M_H2O + M_PA-O = M_Pi + M_DG-O"
        if (reaction == "RHEA:36240")
            .formula <- "M_H2O + M_PA-O => M_Pi + M_DG-O"
        if (reaction == "RHEA:36241")
            .formula <- "M_H2O + M_PA-O <= M_Pi + M_DG-O"
        if (reaction == "RHEA:36242")
            .formula <- "M_H2O + M_PA-O <=> M_Pi + M_DG-O"
        
        ## pc_to_dg
        if (reaction == "RHEA:10604")
            .formula <- "M_H2O + M_PC = M_Phosphocholine + M_1,2-DG"
        if (reaction == "RHEA:10605")
            .formula <- "M_H2O + M_PC => M_Phosphocholine + M_1,2-DG"
        if (reaction == "RHEA:10606")
            .formula <- "M_H2O + M_PC <= M_Phosphocholine + M_1,2-DG"
        if (reaction == "RHEA:10607")
            .formula <- "M_H2O + M_PC <=> M_Phosphocholine + M_1,2-DG"
        
        ## pc_to_sn1lpc
        if (reaction == "RHEA:15801")
            .formula <- "M_H2O + M_PC = M_1-LPC + M_FA"
        if (reaction == "RHEA:15802")
            .formula <- "M_H2O + M_PC => M_1-LPC + M_FA"
        if (reaction == "RHEA:15803")
            .formula <- "M_H2O + M_PC <= M_1-LPC + M_FA"
        if (reaction == "RHEA:15804")
            .formula <- "M_H2O + M_PC <=> M_1-LPC + M_FA"
        
        ## pc_to_sn2lpc
        if (reaction == "RHEA:18689") 
            .formula <- "M_H2O + M_PC = M_2-LPC + M_FA"
        if (reaction == "RHEA:18690") 
            .formula <- "M_H2O + M_PC => M_2-LPC + M_FA"
        if (reaction == "RHEA:18691") 
            .formula <- "M_H2O + M_PC <= M_2-LPC + M_FA"
        if (reaction == "RHEA:18692") 
            .formula <- "M_H2O + M_PC <=> M_2-LPC + M_FA"
        
        ## pc_to_pa
        if (reaction == "RHEA:14445")
            .formula <- "M_H2O + M_PC = M_Choline + M_PA"
        if (reaction == "RHEA:14446")
            .formula <- "M_H2O + M_PC => M_Choline + M_PA"
        if (reaction == "RHEA:14447")
            .formula <- "M_H2O + M_PC <= M_Choline + M_PA"
        if (reaction == "RHEA:14448")
            .formula <- "M_H2O + M_PC <=> M_Choline + M_PA"
        
        ## pc_to_ps
        if (reaction == "RHEA:45088") ## 
            .formula <- "M_L-Serine + M_PC = M_Choline + M_PS"
        if (reaction == "RHEA:45089") ## 
            .formula <- "M_L-Serine + M_PC => M_Choline + M_PS"
        if (reaction == "RHEA:45090") ## 
            .formula <- "M_L-Serine + M_PC <= M_Choline + M_PS"
        if (reaction == "RHEA:45091") ## 
            .formula <- "M_L-Serine + M_PC <=> M_Choline + M_PS"

        ## pco_to_lpco
        if (reaction == "RHEA:36231")
            .formula <- "M_H2O + M_PC-O = M_H+ + M_LPC-O + M_FA"
        if (reaction == "RHEA:36232")
            .formula <- "M_H2O + M_PC-O => M_H+ + M_LPC-O + M_FA"
        if (reaction == "RHEA:36233")
            .formula <- "M_H2O + M_PC-O <= M_H+ + M_LPC-O + M_FA"
        if (reaction == "RHEA:36234")
            .formula <- "M_H2O + M_PC-O <=> M_H+ + M_LPC-O + M_FA"
        
        ## lpco_to_lpao
        if (reaction == "RHEA:39927")
            .formula <- "M_H2O + M_1-LPC-O = M_1-LPA-O + M_H+ + M_Choline"
        if (reaction == "RHEA:39928")
            .formula <- "M_H2O + M_1-LPC-O => M_1-LPA-O + M_H+ + M_Choline"
        if (reaction == "RHEA:39929")
            .formula <- "M_H2O + M_1-LPC-O <= M_1-LPA-O + M_H+ + M_Choline"
        if (reaction == "RHEA:39930")
            .formula <- "M_H2O + M_1-LPC-O <=> M_1-LPA-O + M_H+ + M_Choline"
        
        ## lpco_to_mgo
        if (reaction == "RHEA:36083")
            .formula <- "M_H2O + M_1-LPC-O = M_H+ + M_Phosphocholine + M_1-MG-O"
        if (reaction == "RHEA:36084")
            .formula <- "M_H2O + M_1-LPC-O => M_H+ + M_Phosphocholine + M_1-MG-O"
        if (reaction == "RHEA:36085")
            .formula <- "M_H2O + M_1-LPC-O <= M_H+ + M_Phosphocholine + M_1-MG-O"
        if (reaction == "RHEA:36086")
            .formula <- "M_H2O + M_1-LPC-O <=> M_H+ + M_Phosphocholine + M_1-MG-O"
        
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
        if (reaction == "pe_to_dg")
            .formula <- "M_H2O + M_PE <=> M_P-Ethanolamine + M_1,2-DG"
        
        ## pe_to_sn1lpe
        if (reaction == "RHEA:44604")
            .formula <- "M_H2O + M_PE = M_H+ + M_FA + M_1-LPE"
        if (reaction == "RHEA:44605")
            .formula <- "M_H2O + M_PE => M_H+ + M_FA + M_1-LPE"
        if (reaction == "RHEA:44606")
            .formula <- "M_H2O + M_PE <= M_H+ + M_FA + M_1-LPE"
        if (reaction == "RHEA:44607")
            .formula <- "M_H2O + M_PE <=> M_H+ + M_FA + M_1-LPE"
        
        ## pe_to_sn2lpe
        if (reaction == "RHEA:44408") 
            .formula <- "M_H2O + M_PE = M_H+ + M_FA + M_2-LPE"
        if (reaction == "RHEA:44409") 
            .formula <- "M_H2O + M_PE => M_H+ + M_FA + M_2-LPE"
        if (reaction == "RHEA:44410") 
            .formula <- "M_H2O + M_PE <= M_H+ + M_FA + M_2-LPE"
        if (reaction == "RHEA:44411") 
            .formula <- "M_H2O + M_PE <=> M_H+ + M_FA + M_2-LPE"
        
        if (reaction == "pe_to_nape_sn1")
            .formula <- "M_PE + M_PC <=> M_NAPE + M_2-LPC"
        
        if (reaction == "pe_to_nape_sn2")
            .formula <- "M_PE + M_PC <=> M_NAPE + M_1-LPC"
        
        if (reaction == "pe_to_pa")
            .formula <- "M_H2O + M_PE <=> M_Ethanolamine + M_H+ + M_PA"
        
        ## pe_to_ps
        if (reaction == "RHEA:27606")
            .formula <- "M_L-Serine + M_PE = M_Ethanolamine + M_PS"
        if (reaction == "RHEA:27607")
            .formula <- "M_L-Serine + M_PE => M_Ethanolamine + M_PS"
        if (reaction == "RHEA:27608")
            .formula <- "M_L-Serine + M_PE <= M_Ethanolamine + M_PS"
        if (reaction == "RHEA:27609")
            .formula <- "M_L-Serine + M_PE <=> M_Ethanolamine + M_PS"
        
        if (reaction == "peo_to_lpeo")
            .formula <- "M_H2O + M_PE-O <=> M_H+ + M_LPE-O + M_FA"
        
        if (reaction == "peo_to_napeo_sn1")
            .formula <- "M_PE-O + M_PC <=> M_NAPEO + M_2-LPC"
        
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
            .formula <- "M_H2O + M_PE-P = M_H+ + M_LPE-P + M_FA"
        if (reaction == "RHEA:36196")
            .formula <- "M_H2O + M_PE-P => M_H+ + M_LPE-P + M_FA"
        if (reaction == "RHEA:36197")
            .formula <- "M_H2O + M_PE-P <= M_H+ + M_LPE-P + M_FA"
        if (reaction == "RHEA:36198")
            .formula <- "M_H2O + M_PE-P <=> M_H+ + M_LPE-P + M_FA"
        
        ## lpep_to_fal
        if (reaction == "RHEA:16905")
            .formula <- "M_H2O + M_1-LPE-P = M_FAL + M_Glycerophosphoethanolamine"
        if (reaction == "RHEA:16906")
            .formula <- "M_H2O + M_1-LPE-P => M_FAL + M_Glycerophosphoethanolamine"
        if (reaction == "RHEA:16907")
            .formula <- "M_H2O + M_1-LPE-P <= M_FAL + M_Glycerophosphoethanolamine"
        if (reaction == "RHEA:16908")
            .formula <- "M_H2O + M_1-LPE-P <=> M_FAL + M_Glycerophosphoethanolamine"
        
        ## lpep_to_lpap
        if (reaction == "RHEA:36203")
            .formula <- "M_H2O + M_1-LPE-P = M_LPA-P + M_H+ + M_Ethanolamine"
        if (reaction == "RHEA:36204")
            .formula <- "M_H2O + M_1-LPE-P => M_LPA-P + M_H+ + M_Ethanolamine"
        if (reaction == "RHEA:36205")
            .formula <- "M_H2O + M_1-LPE-P <= M_LPA-P + M_H+ + M_Ethanolamine"
        if (reaction == "RHEA:36206")
            .formula <- "M_H2O + M_1-LPE-P <=> M_LPA-P + M_H+ + M_Ethanolamine"
        
        ## lpep_to_mgp
        if (reaction == "RHEA:36199")
            .formula <- "M_H2O + M_1-LPE-P = M_Phosphoethanolamine + M_H+ + M_1-MG-P"
        if (reaction == "RHEA:36200")
            .formula <- "M_H2O + M_1-LPE-P => M_Phosphoethanolamine + M_H+ + M_1-MG-P"
        if (reaction == "RHEA:36201")
            .formula <- "M_H2O + M_1-LPE-P <= M_Phosphoethanolamine + M_H+ + M_1-MG-P"
        if (reaction == "RHEA:36202")
            .formula <- "M_H2O + M_1-LPE-P <=> M_Phosphoethanolamine + M_H+ + M_1-MG-P"
        
        ## pep_to_napep_sn1
        if (reaction == "pep_to_napep_sn1")
            .formula <- "M_PE-P + M_PC <=> M_NAPEP + M_2-LPC"
        
        ## pep_to_napep_sn2
        if (reaction == "pep_to_napep_sn2")
            .formula <- "M_PE-P + M_PC <=> M_NAPEP + M_1-LPC"
        
        ## pg_to_cl
        if (reaction == "RHEA:32931")
            .formula <- "M_CDP-DG + M_PG = M_H+ + M_CMP + M_CL"
        if (reaction == "RHEA:32932")
            .formula <- "M_CDP-DG + M_PG => M_H+ + M_CMP + M_CL"
        if (reaction == "RHEA:32933")
            .formula <- "M_CDP-DG + M_PG <= M_H+ + M_CMP + M_CL"
        if (reaction == "RHEA:32934")
            .formula <- "M_CDP-DG + M_PG <=> M_H+ + M_CMP + M_CL"
        
        ## pgp_to_pg
        if (reaction == "RHEA:33751")
            .formula <- "M_H2O + M_PGP = M_Pi + M_PG"
        if (reaction == "RHEA:33752")
            .formula <- "M_H2O + M_PGP => M_Pi + M_PG"
        if (reaction == "RHEA:33753")
            .formula <- "M_H2O + M_PGP <= M_Pi + M_PG"
        if (reaction == "RHEA:33754")
            .formula <- "M_H2O + M_PGP <=> M_Pi + M_PG"
        
        ## pi_to_dg
        if (reaction == "RHEA:43484")
            .formula <- "M_H2O + M_PI = M_myo-Inositol-1-P + M_H+ + M_1,2-DG"
        if (reaction == "RHEA:43485")
            .formula <- "M_H2O + M_PI => M_myo-Inositol-1-P + M_H+ + M_1,2-DG"
        if (reaction == "RHEA:43486")
            .formula <- "M_H2O + M_PI <= M_myo-Inositol-1-P + M_H+ + M_1,2-DG"
        if (reaction == "RHEA:43487")
            .formula <- "M_H2O + M_PI <=> M_myo-Inositol-1-P + M_H+ + M_1,2-DG"
        
        ## pi_to_sn1lpi
        if (reaction == "RHEA:18001")
            .formula <- "M_H2O + M_PI = M_1-LPI + M_H+ + M_FA"
        if (reaction == "RHEA:18002")
            .formula <- "M_H2O + M_PI => M_1-LPI + M_H+ + M_FA"
        if (reaction == "RHEA:18003")
            .formula <- "M_H2O + M_PI <= M_1-LPI + M_H+ + M_FA"
        if (reaction == "RHEA:18004")
            .formula <- "M_H2O + M_PI <=> M_1-LPI + M_H+ + M_FA"
            
        ## ps_to_pe
        if (reaction == "RHEA:20828")
            .formula <- "M_H+ + M_PS = M_CO2 + M_PE"
        if (reaction == "RHEA:20829")
            .formula <- "M_H+ + M_PS => M_CO2 + M_PE"
        if (reaction == "RHEA:20830")
            .formula <- "M_H+ + M_PS <= M_CO2 + M_PE"
        if (reaction == "RHEA:20831")
            .formula <- "M_H+ + M_PS <=> M_CO2 + M_PE"
        
        if (reaction == "sm_to_cer")
            .formula <- "M_H2O + M_SM <=> M_Phosphocholine + M_H+ + M_Cer"
        
        ## sphinga_to_dhcer
        if (reaction == "RHEA:53424")
            .formula <- "M_AcylCoA + M_Sphinganine = M_CoA + M_DhCer"
        if (reaction == "RHEA:53425")
            .formula <- "M_AcylCoA + M_Sphinganine => M_CoA + M_DhCer"
        if (reaction == "RHEA:53426")
            .formula <- "M_AcylCoA + M_Sphinganine <= M_CoA + M_DhCer"
        if (reaction == "RHEA:53427")
            .formula <- "M_AcylCoA + M_Sphinganine <=> M_CoA + M_DhCer"
        
        ## tg_to_dg
        if (reaction %in% c("RHEA:33271", "RHEA:44864"))
            .formula <- "M_H2O + M_TG = M_H+ + M_FA + M_1,2-DG"
        if (reaction %in% c("RHEA:33272", "RHEA:44865"))
            .formula <- "M_H2O + M_TG => M_H+ + M_FA + M_1,2-DG"
        if (reaction %in% c("RHEA:33273", "RHEA:44866"))
            .formula <- "M_H2O + M_TG <= M_H+ + M_FA + M_1,2-DG"
        if (reaction %in% c("RHEA:33274", "RHEA:44867"))
            .formula <- "M_H2O + M_TG <=> M_H+ + M_FA + M_1,2-DG"
            
        template$reaction_formula <- .formula
    }
    
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
    
    ## return the template object
    template
}
