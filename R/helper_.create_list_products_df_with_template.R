#' @name .create_list_products_df_with_template
#' 
#' @title Create the object to return in \code{create_reaction}
#' 
#' @description 
#' Helper function for \code{create_reaction}.
#' 
#' The function \code{.create_list_products_df_with_template} will return a
#' list of length 2 containing the products (in the first entry) and the 
#' updated \code{template} object (in the second entry).
#' 
#' @details
#' Depending on the \code{reaction}, the first list entry will contain one or 
#' multiple entries.
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @param df_reaction \code{data.frame}
#' @param template \code{data.frame}
#' @param reaction \code{character(1)}
#' 
#' @return list of length 2
#' 
#' @examples 
#' FA <- c("FA(14:0(12Me))", "FA(16:0(14Me))", "FA(15:1(9Z)(14Me))",        
#'     "FA(17:0(16Me))", "FA(12:0(11Me))", "FA(13:0(12Me))", "FA(14:0(13Me))",
#'     "FA(15:0(14Me))", "FA(16:0(15Me))", "FA(12:0)", "FA(14:0)")
#' substrates <- list(FA = FA)
#' reaction <- "RHEA:15421"
#' 
#' ## create template
#' template <- LipidNetworkPredictR:::.create_template(template = NA, 
#'     reaction = reaction)
#' 
#' ## create data.frame of substrates
#' df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
#'     substrates = substrates, 
#'     constraints = "", negate = FALSE)
#'     
#' ## add products to data.frame
#' df_reaction <- LipidNetworkPredictR:::.add_products(
#'     substrates = df_substrates, 
#'     reaction = reaction)
#'     
#' ## make new data.frame with reaction template
#' df <- LipidNetworkPredictR:::.create_df_with_template(
#'     df_reaction = df_reaction,
#'     template = template, reaction = reaction)
#' 
#' LipidNetworkPredictR:::.create_list_products_df_with_template(
#'     df_reaction = df_reaction,
#'     template = template,
#'     reaction = reaction)
.create_list_products_df_with_template <- function(df_reaction, template, 
    reaction = "RHEA:15421") {
    
    ## acyldhap_to_alkyldhap
    if (reaction %in% c("RHEA:36171", "RHEA:36172", "RHEA:36173", "RHEA:36174")) 
        l <- list(AlkylDHAP = unique(df_reaction$AlkylDHAP))
    
    ## alkyldhap_to_lpao
    if (reaction %in% c("RHEA:36175", "RHEA:36176", "RHEA:36177", "RHEA:36178"))
        l <- list(LPAO = unique(df_reaction$LPAO))

    if (reaction == "c1p_to_cer")
        l <- list(CER = unique(df_reaction$CER))
    
    ## cdpdg_to_pgp
    if (reaction %in% c("RHEA:12593", "RHEA:12594", "RHEA:12595", "RHEA:12596"))
        l <- list(PGP = unique(df_reaction$PGP))
    
    ## cdpdg_to_pi
    if (reaction %in% c("RHEA:11580", "RHEA:11581", "RHEA:11582", "RHEA:11583"))
        l <- list(PI = unique(df_reaction$PI))
    
    if (reaction == "RHEA:17929") ## cer_to_c1p
        l <- list(C1P = unique(df_reaction$C1P))
    
    if (reaction == "RHEA:12088") ## cer_to_glccer
        l <- list(GLCCER = unique(df_reaction$GLCCER))
    
    if (reaction == "RHEA:18765") ## cer_to_sm
        l <- list(SM = unique(df_reaction$SM))
    
    ## cl_to_lcl
    if (reaction %in% c("RHEA:32935", "RHEA:32936", "RHEA:32937", "RHEA:32938"))
        l <- list(
            LCL = unique(df_reaction$LCL), 
            FA = unique(df_reaction$FA))
       
    
    ## coa_to_acyldhap
    if (reaction %in% c("RHEA:17657", "RHEA:17658", "RHEA:17659", "RHEA:17660"))
        l <- list(AcylDHAP = unique(df_reaction$AcylDHAP))
    
    ## coa_to_fatoh
    if (reaction %in% c("RHEA:52716", "RHEA:52717", "RHEA:52718", "RHEA:52719"))
        l <- list(FATOH = unique(df_reaction$FATOH))
    
    ## coa_to_lpa
    if (reaction %in% 
            c("RHEA:15325", "RHEA:15326", "RHEA:15327", "RHEA:15328"))
        l <- list(LPA = unique(df_reaction$LPA))
    
    ## dg_to_sn1mg
    if (reaction %in% c("RHEA:44712", "RHEA:44713", "RHEA:44714", "RHEA:44715",
            "RHEA:35663", "RHEA:35664", "RHEA:35665", "RHEA:35666"))
        l <- list(sn1MG = unique(df_reaction$sn1MG))
    
    ## dg_to_sn2mg
    if (reaction %in% c("RHEA:33275", "RHEA:33276", "RHEA:33277", "RHEA:33278"))
        l <- list(sn2MG = unique(df_reaction$sn2MG))
    
    ## dg_to_pa
    if (reaction %in% c("RHEA:10272", "RHEA:10273", "RHEA:10274", "RHEA:10275"))
        l <- list(PA = unique(df_reaction$PA))
    
    ## dg_to_pc
    if (reaction %in% c("RHEA:32939", "RHEA:32940", "RHEA:32941", "RHEA:32942"))
        l <- list(PC = unique(df_reaction$PC))
    
    ## dg_to_pe
    if (reaction %in% c("RHEA:32943", "RHEA:32944", "RHEA:32945", "RHEA:32946"))
        l <- list(PE = unique(df_reaction$PE))
    
    ## dg_to_tg
    if (reaction %in% c("RHEA:10868", "RHEA:10869", "RHEA:10870", "RHEA:10871")) 
        l <- list(TG = unique(df_reaction$TG))
    
    ## dgo_to_pco
    if (reaction %in% c("RHEA:36179", "RHEA:36180", "RHEA:36181", "RHEA:36182"))
        l <- list(PCO = unique(df_reaction$PCO))
    
    ## dgo_to_peo
    if (reaction %in% c("RHEA:36187", "RHEA:36188", "RHEA:36189", "RHEA:36190"))
        l <- list(PEO = unique(df_reaction$PEO))
    
    if (reaction == "dhcer_to_cer")
        l <- list(CER = unique(df_reaction$CER))
    
    if (reaction == "RHEA:44620") ## dhcer_to_dhsm
        l <- list(DHSM = unique(df_reaction$DHSM))
    
    if (reaction == "RHEA:19253") ## dhsm_to_dhcer
        l <- list(DHCER = unique(df_reaction$DHCER))
    
    ## fa_to_coa
    if (reaction %in% c("RHEA:15421", "RHEA:15422", "RHEA:15423", "RHEA:15424",
            "RHEA:38883", "RHEA:38884", "RHEA:38885", "RHEA:38886"))
        l <- list(AcylCoA = unique(df_reaction$AcylCoA))
    
    ## lcl_to_cl
    if (reaction %in% c("RHEA:35839", "RHEA:35840", "RHEA:35841", "RHEA:35842"))
        l <- list(CL = unique(df_reaction$CL))
    
    
    if (reaction == "RHEA:45420") ## lnape_to_gpnae
        l <- list(
            GPNAE = unique(df_reaction$GPNAE), 
            FA = unique(df_reaction$FA))
    
    ## lpa_to_pa
    if (reaction %in%
            c("RHEA:19709", "RHEA:19710", "RHEA:19711", "RHEA:19712"))
        l <- list(PA = unique(df_reaction$PA))
    
    ## lpao_to_pao
    if (reaction %in% c("RHEA:36235", "RHEA:36236", "RHEA:36237", "RHEA:36238"))
        l <- list(PAO = unique(df_reaction$PAO))
    
    ## sn1lpc_to_fa
    if (reaction %in% c("RHEA:15177", "RHEA:15178", "RHEA:15179", "RHEA:15180")) 
        l <- list(FA = unique(df_reaction$FA))
    
    ## sn2lpc_to_fa
    if (reaction %in% c("RHEA:44696", "RHEA:44697", "RHEA:44698", "RHEA:44699"))
        l <- list(FA = unique(df_reaction$FA))
    
    ## sn1lpc_to_pc
    if (reaction %in% c("RHEA:12937", "RHEA:12938", "RHEA:12939", "RHEA:12940"))
        l <- list(PC = unique(df_reaction$PC))
    
    ## sn1lpe_to_fa
    if (reaction %in% c("RHEA:32967", "RHEA:32968", "RHEA:32969", "RHEA:32970")) 
        l <- list(FA = unique(df_reaction$FA))
    
    if (reaction == "sn2lpe_to_fa")
        l <- list(FA = unique(df_reaction$FA))
    
    ## sn1lpe_to_pe
    if (reaction %in% c("RHEA:32995", "RHEA:32996", "RHEA:32997", "RHEA:32998")) 
        l <- list(PE = unique(df_reaction$PE))
    
    ## sn1lpi_to_pi
    if (reaction %in% c("RHEA:33195", "RHEA:33196", "RHEA:33197", "RHEA:33198"))
        l <- list(PI = unique(df_reaction$PI))
    
    ## 
    if (reaction == "lpeo_to_peo")
        l <- list(PEO = unique(df_reaction$PEO))
    
    ## lpep_to_pep
    if (reaction == "lpep_to_pep")
        l <- list(PEP = unique(df_reaction$PEP))
    
    ## sn1mg_to_dg
    if (reaction %in% c("RHEA:38463", "RHEA:38464", "RHEA:38465", "RHEA:38466",
            "RHEA:39943", "RHEA:39944", "RHEA:39945", "RHEA:39946"))
        l <- list(DG = unique(df_reaction$DG))
    
    ## sn2mg_to_dg
    if (reaction %in% c("RHEA:32947", "RHEA:32948", "RHEA:32949", "RHEA:32950",
            "RHEA:16741", "RHEA:16742", "RHEA:16743", "RHEA:16744"))
        l <- list(DG = unique(df_reaction$DG))
    
    ## sn1mg_to_fa
    if (reaction %in% c("RHEA:34019", "RHEA:34020", "RHEA:34021", "RHEA:34022"))
        l <- list(FA = unique(df_reaction$FA))
    
    ## sn2mg_to_fa
    if (reaction %in% c("RHEA:32871", "RHEA:32872", "RHEA:32873", "RHEA:32874"))
        l <- list(FA = unique(df_reaction$FA))
    
    if (reaction %in% c("RHEA:33747", "RHEA:33748", "RHEA:33749", "RHEA:33750"))
        l <- list(LPA = unique(df_reaction$LPA))
    
    if (reaction == "sn2mg_to_sn1mg")
        l <- list(sn1MG = unique(df_reaction$sn1MG))
    
    if (reaction == "nae_to_fa")
        l <- list(FA = unique(df_reaction$FA))
    
    if (reaction == "nape_to_lnape")
        l <- list(
            LNAPE = unique(df_reaction$LNAPE), FA = unique(df_reaction$FA))
    
    if (reaction == "nape_to_nae")
        l <- list(
            NAE = unique(df_reaction$NAE), PA = unique(df_reaction$PA))
    
    if (reaction == "nape_to_pnae")
        l <- list(
            PNAE = unique(df_reaction$PNAE), DG = unique(df_reaction$DG))
    
    if (reaction == "napeo_to_nae")
        l <- list(
            NAE = unique(df_reaction$NAE), PAO = unique(df_reaction$PAO))
    
    ## pa_to_cdpdg
    if (reaction %in% c("RHEA:16229", "RHEA:16230", "RHEA:16231", "RHEA:16232"))
        l <- list(CDPDG = unique(df_reaction$CDPDG))
    
    ## pa_to_dg
    if (reaction %in%
            c("RHEA:27429", "RHEA:27430", "RHEA:27431", "RHEA:27432"))
        l <- list(DG = unique(df_reaction$DG))
    
    ## pao_to_dgo
    if (reaction %in% c("RHEA:36239", "RHEA:36240", "RHEA:36241", "RHEA:36242"))
        l <- list(DGO = unique(df_reaction$DGO))
    
    ## pc_to_dg
    if (reaction %in% c("RHEA:10604", "RHEA:10605", "RHEA:10606", "RHEA:10607")) 
        l <- list(DG = unique(df_reaction$DG))
    
    ## pc_to_sn1lpc
    if (reaction %in% c("RHEA:15801", "RHEA:15802", "RHEA:15803", "RHEA:15804")) 
        l <- list(sn1LPC = unique(df_reaction$sn1LPC))
    
    ## pc_to_sn2lpc
    if (reaction %in% c("RHEA:18689", "RHEA:18690","RHEA:18691","RHEA:18692"))
        l <- list(sn2LPC = unique(df_reaction$sn2LPC))
    
    ## pc_to_pa
    if (reaction %in% c("RHEA:14445", "RHEA:14446", "RHEA:14447", "RHEA:14448"))
        l <- list(PA = unique(df_reaction$PA))
    
    ## pc_to_ps
    if (reaction %in% c("RHEA:45088", "RHEA:45089", "RHEA:45090", "RHEA:45091")) 
        l <- list(PS = unique(df_reaction$PS))
    
    if (reaction == "pco_to_lpco")
        l <- list(LPCO = unique(df_reaction$LPCO))
    
    if (reaction == "pe_to_dg")
        l <- list(DG = unique(df_reaction$DG))
    
    ## pe_to_sn1lpe
    if (reaction %in% c("RHEA:44604", "RHEA:44605", "RHEA:44606", "RHEA:44607"))
        l <- list(sn1LPE = unique(df_reaction$sn1LPE))
 
    ## pe_to_sn2lpe
    if (reaction %in% c("RHEA:44408", "RHEA:44409", "RHEA:44410", "RHEA:44411"))
        l <- list(sn2LPE = unique(df_reaction$sn2LPE))
    
    if (reaction == "pe_to_nape_sn1")
        l <- list(
            NAPE = unique(df_reaction$NAPE), 
            sn2LPC = unique(df_reaction$sn2LPC))
    
    if (reaction == "pe_to_nape_sn2")
        l <- list(
            NAPE = unique(df_reaction$NAPE), 
            sn1LPC = unique(df_reaction$sn1LPC))
    
    if (reaction == "pe_to_pa")
        l <- list(PA = unique(df_reaction$PA))
    
        ## pe_to_ps
    if (reaction %in% c("RHEA:27606", "RHEA:27607", "RHEA:27608", "RHEA:27609")) 
        l <- list(PS = unique(df_reaction$PS))
    
    if (reaction == "peo_to_lpeo")
        l <- list(LPEO = unique(df_reaction$LPEO))
    
    if (reaction == "peo_to_napeo_sn1")
        l <- list(
            NAPEO = unique(df_reaction$NAPEO), 
            sn2LPC = unique(df_reaction$sn2LPC))
    
    if (reaction == "peo_to_napeo_sn2")
        l <- list(
            NAPEO = unique(df_reaction$NAPEO),
            sn1LPC = unique(df_reaction$sn1LPC))
    
    ## peo_to_pep
    if (reaction %in% c("RHEA:22956", "RHEA:22957", "RHEA:22958", "RHEA:22959"))
        l <- list(PEP = unique(df_reaction$PEP))
    
    ## pep_to_lpep
    if (reaction %in% c("RHEA:36195", "RHEA:36196", "RHEA:36197", "RHEA:36198"))
        l <- list(LPEP = unique(df_reaction$LPEP))
    
    if (reaction == "pep_to_napep_sn1")
        l <- list(
            NAPEP = unique(df_reaction$NAPEP), 
            sn2LPC = unique(df_reaction$sn2LPC))
    
    if (reaction == "pep_to_napep_sn2")
        l <- list(
            NAPEP = unique(df_reaction$NAPEP), 
            sn1LPC = unique(df_reaction$sn1LPC))
    
    ## pg_to_cl
    if (reaction %in% c("RHEA:32931", "RHEA:32932", "RHEA:32933", "RHEA:32934"))
        l <- list(CL = unique(df_reaction$CL))
    
    ## pgp_to_pg
    if (reaction %in% c("RHEA:33751", "RHEA:33752", "RHEA:33753", "RHEA:33754"))
        l <- list(PG = unique(df_reaction$PG))
    
    ## pi_to_dg
    if (reaction %in% c("RHEA:43484", "RHEA:43485", "RHEA:43486", "RHEA:43487"))
        l <- list(DG = unique(df_reaction$DG))
    
    ## pi_to_sn1lpi
    if (reaction %in% c("RHEA:18001", "RHEA:18002", "RHEA:18003", "RHEA:18004"))
        l <- list(
            sn1LPI = unique(df_reaction$sn1LPI),
            FA = unique(df_reaction$FA)
        )
    
    ## ps_to_pe
    if (reaction %in% c("RHEA:20828", "RHEA:20829", "RHEA:20830", "RHEA:20831"))
        l <- list(PE = unique(df_reaction$PE))
    
    if (reaction == "sm_to_cer")
        l <- list(CER = unique(df_reaction$CER))
    
    if (reaction == "sphinga_to_dhcer")
        l <- list(DHCER = unique(df_reaction$DHCER))
    
    ## tg_to_dg
    if (reaction %in% c("RHEA:33271", "RHEA:33272", "RHEA:33273", "RHEA:33274",
            "RHEA:44864", "RHEA:44865", "RHEA:44866", "RHEA:44867"))
        l <- list(
            sn1Loss_DG = unique(df_reaction$sn1Loss_DG), 
            sn1Loss_FA = unique(df_reaction$sn1Loss_FA),
            sn3Loss_DG = unique(df_reaction$sn3Loss_DG),
            sn3Loss_FA = unique(df_reaction$sn3Loss_FA))
            
    ## return the list together with the template object as a list
    list(l, template)
    
}