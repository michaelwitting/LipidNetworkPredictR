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
    
    if (reaction == "RHEA:36171") ## acdhap_to_alkyldhap
        l <- list(ALKYLDHAP = unique(df_reaction$ALKYLDHAP))
    
    if (reaction == "RHEA:36175") ## alkyldhap_to_lpao
        l <- list(LPAO = unique(df_reaction$LPAO))

    if (reaction == "c1p_to_cer")
        l <- list(CER = unique(df_reaction$CER))
    
    if (reaction == "RHEA:12593") ## cdpdg_to_pgp
        l <- list(PGP = unique(df_reaction$PGP))
    
    if (reaction == "RHEA:11580") ## cdpdg_to_pi
        l <- list(PI = unique(df_reaction$PI))
    
    if (reaction == "RHEA:17929") ## cer_to_c1p
        l <- list(C1P = unique(df_reaction$C1P))
    
    if (reaction == "RHEA:12088") ## cer_to_glccer
        l <- list(GLCCER = unique(df_reaction$GLCCER))
    
    if (reaction == "RHEA:18765") ## cer_to_sm
        l <- list(SM = unique(df_reaction$SM))
    
    if (reaction == "coa_to_acdhap")
        l <- list(DHAP = unique(df_reaction$DHAP))
    
    if (reaction == "RHEA:52716") ## coa_to_fatoh
        l <- list(FATOH = unique(df_reaction$FATOH))
    
    ## coa_to_lpa
    if (reaction %in% 
            c("RHEA:15325", "RHEA:15326", "RHEA:15327", "RHEA:15328"))
        l <- list(LPA = unique(df_reaction$LPA))
    
    ## dg_to_sn1mg
    if (reaction %in% c("RHEA:44712", "RHEA:44713", "RHEA:44714", "RHEA:44715"))
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
    
    if (reaction == "RHEA:36179") ## dgo_to_pco
        l <- list(PCO = unique(df_reaction$PCO))
    
    if (reaction == "RHEA:36187") ## dgo_to_peo
        l <- list(PEO = unique(df_reaction$PEO))
    
    if (reaction == "dhcer_to_cer")
        l <- list(CER = unique(df_reaction$CER))
    
    if (reaction == "RHEA:44620") ## dhcer_to_dhsm
        l <- list(DHSM = unique(df_reaction$DHSM))
    
    if (reaction == "RHEA:19253") ## dhsm_to_dhcer
        l <- list(DHCER = unique(df_reaction$DHCER))
    
    ## fa_to_coa
    if (reaction %in% c("RHEA:15421", "RHEA:15422", "RHEA:15423", "RHEA:15424"))
        l <- list(CoA = unique(df_reaction$CoA))
    
    if (reaction == "RHEA:45420") ## lnape_to_gpnae
        l <- list(
            GPNAE = unique(df_reaction$GPNAE), 
            FA = unique(df_reaction$FA))
    
    ## lpa_to_pa
    if (reaction %in%
            c("RHEA:19709", "RHEA:19710", "RHEA:19711", "RHEA:19712"))
        l <- list(PA = unique(df_reaction$PA))
    
    if (reaction == "RHEA:36235") ## lpao_to_pao
        l <- list(PAO = unique(df_reaction$PAO))
    
    if (reaction == "RHEA:15177") ## sn1lpc_to_fa
        l <- list(FA = unique(df_reaction$FA))
    
    if (reaction == "RHEA:44696") ## sn2lpc_to_fa
        l <- list(FA = unique(df_reaction$FA))
    
    if (reaction == "RHEA:12937") ## sn1lpc_to_pc
        l <- list(PC = unique(df_reaction$PC))
    
    if (reaction == "RHEA:32967") ## sn1lpe_to_fa
        l <- list(FA = unique(df_reaction$FA))
    
    if (reaction == "sn2lpe_to_fa")
        l <- list(FA = unique(df_reaction$FA))
    
    if (reaction == "RHEA:32995") ## sn1lpe_to_pe
        l <- list(PE = unique(df_reaction$PE))
    
    if (reaction == "lpeo_to_peo")
        l <- list(PEO = unique(df_reaction$PEO))
    
    if (reaction == "lpep_to_pep")
        l <- list(PEP = unique(df_reaction$PEP))
    
    ## sn1mg_to_dg
    if (reaction %in% c("RHEA:38463", "RHEA:38464", "RHEA:38465", "RHEA:38466"))
        l <- list(DG = unique(df_reaction$DG))
    
    ## sn2mg_to_dg
    if (reaction %in% c("RHEA:32947", "RHEA:32948", "RHEA:32949", "RHEA:32950"))
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
    
    if (reaction == "pa_to_cdpdg")
        l <- list(CDPDG = unique(df_reaction$CDPDG))
    
    ## pa_to_dg
    if (reaction %in%
            c("RHEA:27429", "RHEA:27430", "RHEA:27431", "RHEA:27432"))
        l <- list(DG = unique(df_reaction$DG))
    
    if (reaction == "pao_to_dgo")
        l <- list(DGO = unique(df_reaction$DGO))
    
    if (reaction == "RHEA:10604") ## pc_to_dg
        l <- list(DG = unique(df_reaction$DG))
    
    if (reaction == "RHEA:15801") ## pc_to_sn1lpc
        l <- list(sn1LPC = unique(df_reaction$sn1LPC))
    
    if (reaction == "RHEA:18689") ## pc_to_sn2lpc
        l <- list(sn2LPC = unique(df_reaction$sn2LPC))
    
    if (reaction == "RHEA:14445") ## pc_to_pa
        l <- list(PA = unique(df_reaction$PA))
    
    ## pc_to_ps
    if (reaction %in% c("RHEA:45088", "RHEA:45089", "RHEA:45090", "RHEA:45091")) 
        l <- list(PS = unique(df_reaction$PS))
    
    if (reaction == "pco_to_lpco")
        l <- list(LPCO = unique(df_reaction$LPCO))
    
    if (reaction == "pe_to_dg")
        l <- list(DG = unique(df_reaction$DG))
    
    if (reaction == "pe_to_sn1lpe")
        l <- list(sn1LPE = unique(df_reaction$sn1LPE))
    
    if (reaction == "RHEA:44408") ## pe_to_sn2lpe
        l <- list(sn2LPE = unique(df_reaction$sn2LPE))
    
    if (reaction == "pe_to_nape_sn1")
        l <- list(
            NAPE = unique(df_reaction$NAPE), LPC = unique(df_reaction$LPC))
    
    if (reaction == "pe_to_nape_sn2")
        l <- list(
            NAPE = unique(df_reaction$NAPE), LPC = unique(df_reaction$LPC))
    
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
            LPC = unique(df_reaction$LPC))
    
    if (reaction == "peo_to_napeo_sn2")
        l <- list(
            NAPEO = unique(df_reaction$NAPEO),
            LPC = unique(df_reaction$LPC))
    
    if (reaction == "peo_to_pep")
        l <- list(PEP = unique(df_reaction$PEP))
    
    if (reaction == "pep_to_lpep")
        l <- list(LPEP = unique(df_reaction$LPEP))
    
    if (reaction == "pep_to_napep_sn1")
        l <- list(
            NAPEP = unique(df_reaction$NAPEP), 
            LPC = unique(df_reaction$LPC))
    
    if (reaction == "pep_to_napep_sn2")
        l <- list(
            NAPEP = unique(df_reaction$NAPEP), 
            LPC = unique(df_reaction$LPC))
    
    if (reaction == "pg_to_cl")
        l <- list(CL = unique(df_reaction$CL))
    
    if (reaction == "pgp_to_pg")
        l <- list(PG = unique(df_reaction$PG))
    
    ## ps_to_pe
    if (reaction %in% c("RHEA:20828", "RHEA:20829", "RHEA:20830", "RHEA:20831"))
        l <- list(PE = unique(df_reaction$PE))
    
    if (reaction == "sm_to_cer")
        l <- list(CER = unique(df_reaction$CER))
    
    if (reaction == "sphinga_to_dhcer")
        l <- list(DHCER = unique(df_reaction$DHCER))
    
    ## tg_to_dg
    if (reaction %in% c("RHEA:33271", "RHEA:33272", "RHEA:33273", "RHEA:33274"))
        l <- list(
            sn1Loss_dg = unique(df_reaction$sn1Loss_dg), 
            sn1Loss_fa = unique(df_reaction$sn1Loss_fa),
            sn3Loss_dg = unique(df_reaction$sn3Loss_dg),
            sn3Loss_fa = unique(df_reaction$sn3Loss_fa))
            
    ## return the list together with the template object as a list
    list(l, template)
    
}