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
#' @param template \code{template}
#' 
#' @return list of length 2
#' 
#' @examples 
#' FA <- c("FA(14:0(12Me))", "FA(16:0(14Me))", "FA(15:1(9Z)(14Me))",        
#'     "FA(17:0(16Me))", "FA(12:0(11Me))", "FA(13:0(12Me))", "FA(14:0(13Me))",
#'     "FA(15:0(14Me))", "FA(16:0(15Me))", "FA(12:0)", "FA(14:0)")
#' substrates <- list(FA = FA)
#' 
#' ## create template
#' template <- wormLipidPredictR:::.create_template(template = NA, 
#'     reaction = "fa_to_coa")
#' 
#' ## create data.frame of substrates
#' df_substrates <- wormLipidPredictR:::.create_substrates_combinations(
#'     substrates = substrates, 
#'     constraints = "", negate = FALSE)
#'     
#' ## add products to data.frame
#' df_reaction <- wormLipidPredictR:::.add_products(
#'     substrates = df_substrates, 
#'     reaction = "fa_to_coa")
#'     
#' ## make new data.frame with reaction template
#' df <- wormLipidPredictR:::.create_df_with_template(
#'     df_reaction = df_reaction,
#'     template = template, reaction = reaction)
#' 
#' wormLipidPredictR:::.create_list_products_df_with_template(
#'     df_reaction = df_reaction,
#'     template = template,
#'     reaction = "fa_to_coa")
.create_list_products_df_with_template <- function(df_reaction, template, 
    reaction = "fa_to_coa") {
    
    if (reaction == "acdhap_to_alkyldhap")
        .l <- list(ALKYLDHAP = unique(df_reaction$ALKYLDHAP))
    
    if (reaction == "alkyldhap_to_lpao")
        .l <- list(LPAO = unique(df_reaction$LPAO))

    if (reaction == "c1p_to_cer")
        .l <- list(CER = unique(df_reaction$CER))
    
    if (reaction == "cdpdg_to_pgp")
        .l <- list(PGP = unique(df_reaction$PGP))
    
    if (reaction == "cdpdg_to_pi")
        .l <- list(PI = unique(df_reaction$PI))
    
    if (reaction == "cer_to_c1p")
        .l <- list(C1P = unique(df_reaction$C1P))
    
    if (reaction == "cer_to_glccer")
        .l <- list(GLCCER = unique(df_reaction$GLCCER))
    
    if (reaction == "cer_to_sm")
        .l <- list(SM = unique(df_reaction$SM))
    
    if (reaction == "coa_to_acdhap")
        .l <- list(DHAP = unique(df_reaction$DHAP))
    
    if (reaction == "coa_to_fatoh")
        .l <- list(FATOH = unique(df_reaction$FATOH))
    
    if (reaction == "coa_to_lpa")
        .l <- list(LPA = unique(df_reaction$LPA))
    
    if (reaction == "dg_to_sn1mg")
        .l <- list(sn1MG = unique(df_reaction$sn1MG))
    
    if (reaction == "dg_to_sn2mg")
        .l <- list(sn2MG = unique(df_reaction$sn2MG))
    
    if (reaction == "dg_to_pa")
        .l <- list(PA = unique(df_reaction$PA))
    
    if (reaction == "dg_to_pc")
        .l <- list(PC = unique(df_reaction$PC))
    
    if (reaction == "dg_to_pe")
        .l <- list(PE = unique(df_reaction$PE))
    
    if (reaction == "dg_to_tg")
        .l <- list(TG = unique(df_reaction$TG))
    
    if (reaction == "dgo_to_pco")
        .l <- list(PCO = unique(df_reaction$PCO))
    
    if (reaction == "dgo_to_peo")
        .l <- list(PEO = unique(df_reaction$PEO))
    
    if (reaction == "dhcer_to_cer")
        .l <- list(CER = unique(df_reaction$CER))
    
    if (reaction == "dhcer_to_dhsm")
        .l <- list(DHSM = unique(df_reaction$DHSM))
    
    if (reaction == "dhsm_to_dhcer")
        .l <- list(DHCER = unique(df_reaction$DHCER))
    
    if (reaction == "fa_to_coa")
        .l <- list(CoA = unique(df_reaction$CoA))
    
    if (reaction == "lnape_to_gpnae")
        .l <- list(
            GPNAE = unique(df_reaction$GPNAE), 
            FA = unique(df_reaction$FA))
    
    if (reaction == "lpa_to_pa")
        .l <- list(PA = unique(df_reaction$PA))
    
    if (reaction == "lpao_to_pao")
        .l <- list(PAO = unique(df_reaction$PAO))
    
    if (reaction == "sn1lpc_to_fa")
        .l <- list(FA = unique(df_reaction$FA))
    
    if (reaction == "sn2lpc_to_fa")
        .l <- list(FA = unique(df_reaction$FA))
    
    if (reaction == "sn1lpc_to_pc")
        .l <- list(PC = unique(df_reaction$PC))
    
    if (reaction == "sn1lpe_to_fa")
        .l <- list(FA = unique(df_reaction$FA))
    
    if (reaction == "sn2lpe_to_fa")
        .l <- list(FA = unique(df_reaction$FA))
    
    if (reaction == "sn1lpe_to_pe")
        .l <- list(PE = unique(df_reaction$PE))
    
    if (reaction == "lpeo_to_peo")
        .l <- list(PEO = unique(df_reaction$PEO))
    
    if (reaction == "lpep_to_pep")
        .l <- list(PEP = unique(df_reaction$PEP))
    
    if (reaction == "sn1mg_to_dg")
        .l <- list(DG = unique(df_reaction$DG))
    
    if (reaction == "sn2mg_to_dg")
        .l <- list(DG = unique(df_reaction$DG))
    
    if (reaction == "sn1mg_to_fa")
        .l <- list(FA = unique(df_reaction$FA))
    
    if (reaction == "sn2mg_to_fa")
        .l <- list(FA = unique(df_reaction$FA))
    
    if (reaction == "sn1mg_to_lpa")
        .l <- list(LPA = unique(df_reaction$LPA))
    
    if (reaction == "sn2mg_to_sn1mg")
        .l <- list(sn1MG = unique(df_reaction$sn1MG))
    
    if (reaction == "nae_to_fa")
        .l <- list(FA = unique(df_reaction$FA))
    
    if (reaction == "nape_to_lnape")
        .l <- list(
            LNAPE = unique(df_reaction$LNAPE), FA = unique(df_reaction$FA))
    
    if (reaction == "nape_to_nae")
        .l <- list(
            NAE = unique(df_reaction$NAE), PA = unique(df_reaction$PA))
    
    if (reaction == "nape_to_pnae")
        .l <- list(
            PNAE = unique(df_reaction$PNAE), DG = unique(df_reaction$DG))
    
    if (reaction == "napeo_to_nae")
        .l <- list(
            NAE = unique(df_reaction$NAE), PAO = unique(df_reaction$PAO))
    
    if (reaction == "pa_to_cdpdg")
        .l <- list(CDPDG = unique(df_reaction$CDPDG))
    
    if (reaction == "pa_to_dg")
        .l <- list(DG = unique(df_reaction$DG))
    
    if (reaction == "pao_to_dgo")
        .l <- list(DGO = unique(df_reaction$DGO))
    
    if (reaction == "pc_to_dg")
        .l <- list(DG = unique(df_reaction$DG))
    
    if (reaction == "pc_to_sn1lpc")
        .l <- list(sn1LPC = unique(df_reaction$sn1LPC))
    
    if (reaction == "pc_to_sn2lpc")
        .l <- list(sn2LPC = unique(df_reaction$sn2LPC))
    
    if (reaction == "pc_to_pa")
        .l <- list(PA = unique(df_reaction$PA))
    
    if (reaction == "pc_to_ps")
        .l <- list(PS = unique(df_reaction$PS))
    
    if (reaction == "pco_to_lpco")
        .l <- list(LPCO = unique(df_reaction$LPCO))
    
    if (reaction == "pe_to_dg")
        .l <- list(DG = unique(df_reaction$DG))
    
    if (reaction == "pe_to_sn1lpe")
        .l <- list(sn1LPE = unique(df_reaction$sn1LPE))
    
    if (reaction == "pe_to_sn2lpe")
        .l <- list(sn2LPE = unique(df_reaction$sn2LPE))
    
    if (reaction == "pe_to_nape_sn1")
        .l <- list(
            NAPE = unique(df_reaction$NAPE), LPC = unique(df_reaction$LPC))
    
    if (reaction == "pe_to_nape_sn2")
        .l <- list(
            NAPE = unique(df_reaction$NAPE), LPC = unique(df_reaction$LPC))
    
    if (reaction == "pe_to_pa")
        .l <- list(PA = unique(df_reaction$PA))
    
    if (reaction == "pe_to_ps")
        .l <- list(PS = unique(df_reaction$PS))
    
    if (reaction == "peo_to_lpeo")
        .l <- list(LPEO = unique(df_reaction$LPEO))
    
    if (reaction == "peo_to_napeo_sn1")
        .l <- list(
            NAPEO = unique(df_reaction$NAPEO), 
            LPC = unique(df_reaction$LPC))
    
    if (reaction == "peo_to_napeo_sn2")
        .l <- list(
            NAPEO = unique(df_reaction$NAPEO),
            LPC = unique(df_reaction$LPC))
    
    if (reaction == "peo_to_pep")
        .l <- list(PEP = unique(df_reaction$PEP))
    
    if (reaction == "pep_to_lpep")
        .l <- list(LPEP = unique(df_reaction$LPEP))
    
    if (reaction == "pep_to_napep_sn1")
        .l <- list(
            NAPEP = unique(df_reaction$NAPEP), 
            LPC = unique(df_reaction$LPC))
    
    if (reaction == "pep_to_napep_sn2")
        .l <- list(
            NAPEP = unique(df_reaction$NAPEP), 
            LPC = unique(df_reaction$LPC))
    
    if (reaction == "pg_to_cl")
        .l <- list(CL = unique(df_reaction$CL))
    
    if (reaction == "pgp_to_pg")
        .l <- list(PG = unique(df_reaction$PG))
    
    if (reaction == "ps_to_pe")
        .l <- list(PE = unique(df_reaction$PE))
    
    if (reaction == "sm_to_cer")
        .l <- list(CER = unique(df_reaction$CER))
    
    if (reaction == "sphinga_to_dhcer")
        .l <- list(DHCER = unique(df_reaction$DHCER))
    
    if (reaction == "tg_to_dg")
        .l <- list(
            sn1Loss_dg = unique(df_reaction$sn1Loss_dg), 
            sn1Loss_fa = unique(df_reaction$sn1Loss_fa),
            sn3Loss_dg = unique(df_reaction$sn3Loss_dg),
            sn3Loss_fa = unique(df_reaction$sn3Loss_fa))
            
    ## return the list together with the template object as a list
    list(.l, template)
    
}