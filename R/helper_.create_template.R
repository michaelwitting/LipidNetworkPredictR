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
        if (reaction == "RHEA:36171") ## acdhap_to_alkyldhap
            .formula <- "M_adhap + M_alkylR1oh <=> M_H+ + M_FA + M_akdhap"
        
        if (reaction == "RHEA:36175") ## alkyldhap_to_lpao
            .formula <- "M_H+ + M_NADPH + M_akdhap <=> M_alpa_pl + M_NADP"
        
        if (reaction == "c1p_to_cer")
            .formula <- "M_H2O + M_c17isocrmp <=> M_Pi + M_c17isocrm"
        
        if (reaction == "RHEA:12593") ## cdpdg_to_pgp
            .formula <- "M_Glycerol-3-P + M_cdpdag <=> M_H+ + M_CMP + M_pgp"
        
        if (reaction == "RHEA:11580") ## cdpdg_to_pi
            .formula <- "M_Inositol + M_cdpdag <=> M_H+ + M_CMP + M_pail"
        
        if (reaction == "RHEA:17929") ## cer_to_c1p
            .formula <- "M_ATP + M_c17isocrm <=> M_H+ + M_ADP + M_c17isocrmp"
        
        if (reaction == "RHEA:12088") ## cer_to_glccer
            .formula <- "M_UDP-Glucose + M_c17isocrm <=> M_H+ + M_UDP + M_c17isogluside"
        
        if (reaction == "RHEA:18765") ## cer_to_sm
            .formula <- "M_PC + M_c17isocrm <=> M_1,2-DG + M_c17isosphmyln"
        
        if (reaction == "coa_to_acdhap")
            .formula <- "M_Dihydroxyacetone-P + M_AcylCoA <=> M_CoA + M_adhap"

        if (reaction == "RHEA:52716") ## coa_to_fatoh
            .formula <- "M_AcylCoA + 2 M_NADPH + 2 M_H+ <=> M_alkylR1oh + 2 M_NADP + M_CoA"
        
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
        if (reaction == c("RHEA:44712"))
            .formula <- "M_H2O + M_1,2-DG = M_H+ + M_1-MG + M_FA"
        if (reaction == c("RHEA:44713"))
            .formula <- "M_H2O + M_1,2-DG => M_H+ + M_1-MG + M_FA"
        if (reaction == c("RHEA:44714"))
            .formula <- "M_H2O + M_1,2-DG <= M_H+ + M_1-MG + M_FA"
        if (reaction == c("RHEA:44715"))
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
        
        
        if (reaction == "RHEA:36179") ## dgo_to_pco
            .formula <- "M_CDP-Choline + M_akac2g <=> M_H+ + M_CMP + M_akac2gchol"
        
        if (reaction == "RHEA:36187") ## dgo_to_peo
            .formula <- "M_CDP-Ethanolamine + M_akac2g <=> M_H+ + M_CMP + M_akac2gpe"
        
        if (reaction == "dhcer_to_cer")
            .formula <- "M_H+ + M_NADH + M_O2 + M_c17isodhcrm <=> 2 M_H2O + M_NAD + M_c17isocrm"
        
        if (reaction == "RHEA:44620") ## dhcer_to_dhsm
            .formula <- "M_PC + M_c17isodhcrm <=> M_1,2-DG + M_c17isodhsphmyln"
        
        if (reaction == "RHEA:19253") ## dhsm_to_dhcer
            .formula <- "M_H2O + M_c17isodhsphmyln <=> M_P-Choline + M_H+ + M_c17isodhcrm"
        
        ## fa_to_coa
        if (reaction == "RHEA:15421") 
            .formula <- "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA"
        if (reaction == "RHEA:15422")
            .formula <- "M_ATP + M_CoA + M_FA => M_PPi + M_AMP + M_AcylCoA"
        if (reaction == "RHEA:15423")
            .formula <- "M_ATP + M_CoA + M_FA <= M_PPi + M_AMP + M_AcylCoA"
        if (reaction == "RHEA:15424")
            .formula <- "M_ATP + M_CoA + M_FA <=> M_PPi + M_AMP + M_AcylCoA"
        
        if (reaction == "RHEA:45420") ## lnape_to_gpnae
            .formula <- "M_lnape + M_H2O <=> M_gpnae + M_FA"
        
        ## lpa_to_pa
        if (reaction == "RHEA:19709")
            .formula <- "M_LPA + M_AcylCoA = M_CoA + M_PA"
        if (reaction == "RHEA:19710")
            .formula <- "M_LPA + M_AcylCoA => M_CoA + M_PA"
        if (reaction == "RHEA:19711") 
            .formula <- "M_LPA + M_AcylCoA <= M_CoA + M_PA"
        if (reaction == "RHEA:19712") 
            .formula <- "M_LPA + M_AcylCoA <=> M_CoA + M_PA"
        
        if (reaction == "RHEA:36235") ## lpao_to_pao
             .formula <- "M_LPAO + M_AcylCoA <=> M_PAO + M_CoA"
        
        if (reaction == "RHEA:15177") ## sn1lpc_to_fa
            .formula <- "M_H2O + M_ag3pc <=> M_g3pc + M_H+ + M_FA"
        
        if (reaction == "RHEA:44696") ## sn2lpc_to_fa
            .formula <- "M_H2O + M_2agpc <=> M_g3pc + M_H+ + M_FA"
        
        if (reaction == "RHEA:12937") ## sn1lpc_to_pc
            .formula <- "M_ag3pc + M_AcylCoA <=> M_PC + M_CoA"
        
        if (reaction == "RHEA:32967") ## sn1lpe_to_fa
            .formula <- "M_H2O + M_ag3pe <=> M_H+ + M_FA + M_g3pe"
        
        if (reaction == "sn2lpe_to_fa")
            .formula <- "M_H2O + M_2agpe <=> M_H+ + M_FA + M_g3pe"
        
        if (reaction == "RHEA:32995") ## sn1lpe_to_pe
            .formula <- "M_AcylCoA + M_acg3pe <=> M_CoA + M_PE"
        
        if (reaction == "lpeo_to_peo")
            .formula <- "M_AcylCoA + M_ak2lgpe <=> M_CoA + M_akac2gpe"
        
        if (reaction == "lpep_to_pep")
            .formula <- "M_AcylCoA + M_alken2gpe <=> M_CoA + M_alkenac2gpe"
        
        ## sn1mg_to_dg
        if (reaction == "RHEA:38463")
            .formula <- "M_1-MG + M_AcylCoA = M_CoA + M_1,2-DG"
        if (reaction == "RHEA:38464")
            .formula <- "M_1-MG + M_AcylCoA => M_CoA + M_1,2-DG"
        if (reaction == "RHEA:38465")
            .formula <- "M_1-MG + M_AcylCoA <= M_CoA + M_1,2-DG"
        if (reaction == "RHEA:38466")
            .formula <- "M_1-MG + M_AcylCoA <=> M_CoA + M_1,2-DG"
        
        ## sn2mg_to_dg
        if (reaction == "RHEA:32947")
            .formula <- "M_AcylCoA + M_1-MG = M_CoA + M_1,2-DG"
        if (reaction == "RHEA:32948")
            .formula <- "M_AcylCoA + M_1-MG => M_CoA + M_1,2-DG"
        if (reaction == "RHEA:32949")
            .formula <- "M_AcylCoA + M_1-MG <= M_CoA + M_1,2-DG"
        if (reaction == "RHEA:32950")
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
            .formula <- "M_mag <=> M_1magol"
        
        if (reaction == "nae_to_fa")
            .formula <- "M_H2O + M_NAE <=> M_Ethanolamine + M_H+ + M_FA"
        
        if (reaction == "nape_to_lnape")
            .formula <- "M_NAPE + M_H2O <=> M_lnape + M_FA"
        
        if (reaction == "nape_to_nae")
            .formula <- "M_NAPE + M_H2O <=> M_NAE + M_PA"
        
        if (reaction == "nape_to_pnae")
            .formula <- "M_NAPE + M_H2O <=> M_PNAE + M_1,2-DG"
        
        if (reaction == "napeo_to_nae")
            .formula <- "M_akac2nape + M_H2O <=> M_NAE + M_akac2gp"
        
        if (reaction == "pa_to_cdpdg")
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
        
        if (reaction == "pao_to_dgo")
            .formula <- "M_H2O + M_akac2gp <=> M_Pi + M_akac2g"
        
        if (reaction == "RHEA:10604") ## pc_to_dg
            .formula <- "M_H2O + M_PC <=> M_P-Choline + M_1,2-DG"
        
        if (reaction == "RHEA:15801") ## pc_to_sn1lpc
            .formula <- "M_H2O + M_PC <=> M_ag3pc + M_FA"
        
        if (reaction == "RHEA:18689") ## pc_to_sn2lpc
            .formula <- "M_H2O + M_PC <=> M_2agpc + M_FA"
        
        if (reaction == "RHEA:14445") ## pc_to_pa
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
        
        if (reaction == "pco_to_lpco")
            .formula <- "M_H2O + M_akac2gchol <=> M_H+ + M_ak2lgchol + M_FA"
        
        if (reaction == "pe_to_dg")
            .formula <- "M_H2O + M_PE <=> M_P-Ethanolamine + M_1,2-DG"
        
        if (reaction == "pe_to_sn1lpe")
            .formula <- "M_H2O + M_PE <=> M_H+ + M_FA + M_ag3pe"
        
        if (reaction == "RHEA:44408") ## pe_to_sn2lpe
            .formula <- "M_H2O + M_PE <=> M_H+ + M_FA + M_2agpe"
        
        if (reaction == "pe_to_nape_sn1")
            .formula <- "M_PE + M_PC <=> M_NAPE + M_2agpc"
        
        if (reaction == "pe_to_nape_sn2")
            .formula <- "M_PE + M_PC <=> M_NAPE + M_ag3pc"
        
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
            .formula <- "M_H2O + M_akac2gpe <=> M_H+ + M_ak2lgpe + M_FA"
        
        if (reaction == "peo_to_napeo_sn1")
            .formula <- "M_akac2gpe + M_PC <=> M_akac2nape + M_2agpc"
        
        if (reaction == "peo_to_napeo_sn2")
            .formula <- "M_akac2gpe + M_PC <=> M_akac2nape + M_ag3pc"
        
        if (reaction == "peo_to_pep")
            .formula <- "M_akac2gpe + M_NADPH + M_H+ + M_O2 <=> M_alkenac2gpe + M_NADP + 2 M_H2O"
        
        if (reaction == "pep_to_lpep")
            .formula <- "M_H2O + M_alkenac2gpe <=> M_H+ + M_alken2gpe + M_FA"
        
        if (reaction == "pep_to_napep_sn1")
            .formula <- "M_alkenac2gpe + M_PC <=> M_alkenac2nape + M_2agpc"
        
        if (reaction == "pep_to_napep_sn2")
            .formula <- "M_alkenac2gpe + M_PC <=> M_alkenac2nape + M_ag3pc"
        
        if (reaction == "pg_to_cl")
            .formula <- "M_cdpdag + M_pg <=> M_H+ + M_CMP + M_clpn"
        
        if (reaction == "pgp_to_pg")
            .formula <- "M_H2O + M_pgp <=> M_Pi + M_pg"
        
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
            .formula <- "M_H2O + M_c17isosphmyln <=> M_P-Choline + M_H+ + M_c17isocrm"
        
        if (reaction == "sphinga_to_dhcer")
            .formula <- "M_AcylCoA + M_c17isosphgn <=> M_CoA + M_c17isodhcrm"
        
        ## tg_to_dg
        if (reaction == "RHEA:33271")
            .formula <- "M_H2O + M_TG = M_H+ + M_FA + M_1,2-DG"
        if (reaction == "RHEA:33272")
            .formula <- "M_H2O + M_TG => M_H+ + M_FA + M_1,2-DG"
        if (reaction == "RHEA:33273")
            .formula <- "M_H2O + M_TG <= M_H+ + M_FA + M_1,2-DG"
        if (reaction == "RHEA:33274")
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
        regex = "[+]")[[1]]
    .formula_substrate <- stringi::stri_replace_all_fixed(str = .formula_substrate, 
        pattern = " ", replacement = "")
    .formula_product <- stringi::stri_split(str = .formula_product, 
        regex = "[+]")[[1]]
    .formula_product <- stringi::stri_replace_all_fixed(str = .formula_product, 
        pattern = " ", replacement = "")

    ## write the substrates and products to the template
    template$reaction_substrate <- .formula_substrate
    template$reaction_product <- .formula_product
    
    ## return the template object
    template
}
