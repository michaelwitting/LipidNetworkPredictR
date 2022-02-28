.create_list_products_df_with_template <- function(df_reactions, df_template, 
    reaction) {
    
    if (reaction == "acdhap_to_alkyldhap")
        .l <- list(list(ALKYLDHAP = unique(df_reactions$ALKYLDHAP)), 
            df_template))
    
    if (reaction == "alkyldhap_to_lpao")
        .l <- list(list(LPAO = unique(df_reactions$LPAO)), df_template)

    if (reaction == "c1p_to_cer")
        .l <- list(list(CER = unique(df_reactions$CER), df_template)
    
    if (reaction == "cdpdg_to_pgp")
        .l <- list(list(PGP = unique(df_reactions$PGP)), df_template)
    
    if (reaction == "cdpdg_to_pi")
        .l <- list(list(PI = unique(df_reactions$PI)), df_template)
    
    if (reaction == "cer_to_c1p")
        .l <- list(list(C1P = unique(df_reactions$C1P)), df_template)
    
    if (reaction == "cer_to_glccer")
        .l <- list(list(GLCCER = unique(df_reactions$GLCCER)), df_template)
    
    if (reaction == "cer_to_sm")
        .l <- list(list(SM = unique(df_reactions$SM)), df_template)
    
    if (reaction == "coa_to_acdhap")
        .l <- list(list(DHAP = unique(df_reactions$DHAP)), df_template)
    
    if (reaction == "coa_to_fatoh")
        .l <- list(list(FATOH = unique(df_reactions$FATOH)), df_template)
    
    if (reaction == "coa_to_lpa")
        .l <- list(list(LPA = unique(df_reactions$LPA)), df_template)
    
    if (reaction == "dg_to_sn1mg")
        .l <- list(list(MG = unique(df_reaction$MG)), df_template)
    
    if (reaction == "dg_to_sn2mg")
        .l <- list(list(MG = unique(df_reactions$MG)), df_template)
    
    if (reaction == "dg_to_pa")
        .l <- list(list(PA = unique(df_reactions$PA)), df_template))
    
    if (reaction == "dg_to_pc")
        .l <- list(list(PC = unique(df_reactions$PC)), df_template)
    
    if (reaction == "dg_to_pe")
        .l <- list(list(PE = unique(df_reactions$PE)), df_template)
    
    if (reaction == "dg_to_tg")
        .l <- list(list(TG = unique(df_reactions$TG)), df_template)
    
    if (reaction == "dgo_to_pco")
        .l <- list(list(PCO = unique(df_reactions$PCO)), df_template)
    
    if (reaction == "dgo_to_peo")
        .l <- list(list(PEO = unique(df_reactions$PEO)), df_template)
    
    if (reaction == "dhcer_to_cer")
        .l <- list(list(CER = unique(df_reactions$CER)), df_template)
    
    if (reaction == "dhcer_to_dhsm")
        .l <- list(list(DHSM = unique(df_reactions$DHSM)), df_template)
    
    if (reaction == "dhsm_to_dhcer")
        .l <- list(list(unique(df_reactions$DHCER)), df_template)
    
    if (reaction == "fa_to_coa")
        .l <- list(list(CoA = unique(df_reactions$CoA)), df_template)
    
    if (reaction == "lnape_to_gpnae")
        .l <- list(
            list(GPNAE = unique(df_reactions$GPNAE), FA = unique(df_reactions$FA)), 
            df_template)
    
    if (reaction == "lpa_to_pa")
        .l <- list(PA = list(unique(df_reactions$PA)), df_template)
    
    if (reaction == "lpao_to_pao")
        .l <- list(PAO = list(unique(df_reactions$PAO)), df_template)
    
    if (reaction == "sn1lpc_to_fa")
        .l <- list(FA = list(unique(df_reactions$FA)), df_template)
    
    if (reaction == "sn21pc_to_fa")
        .l <- list(FA = list(unique(df_reactions$FA)), df_template)
    
    if (reaction == "sn1lpc_to_pc")
        .l <- list(PC = list(unique(df_reactions$PC)), df_template)
    
    if (reaction == "sn1lpe_to_fa")
        .l <- list(FA = list(unique(df_reactions$FA)), df_template)
    
    if (reaction == "sn2lpe_to_fa")
        .l <- list(FA = list(unique(df_reactions$FA)), df_template)
    
    if (reaction == "sn1lpe_to_pe")
        .l <- list(PE = list(unique(df_reactions$PE)), df_template)
    
    if (reaction == "lpeo_to_peo")
        .l <- list(PO = list(unique(df_reactions$PEO)), df_template)
    
    if (reaction == "lpep_to_pep")
        .l <- list(PEP = list(unique(df_reactions$PEP)), df_template)
    
    if (reaction == "sn1mg_to_dg")
        .l <- list(DG = list(unique(df_reactions$DG)), df_template)
    
    if (reaction == "sn2mg_to_dg")
        .l <- list(DG = list(unique(df_reactions$DG)), df_template)
    
    if (reaction == "sn1mg_to_fa")
        .l <- list(FA = list(unique(df_reactions$FA)), df_template)
    
    if (reaction == "sn2mg_to_fa")
        .l <- list(FA = list(unique(df_reactions$FA)), df_template)
    
    if (reaction == "sn1mg_to_lpa")
        .l <- list(LPA = list(unique(df_reactions$LPA)), df_template)
    
    if (reaction == "sn2mg_to_sn1mg")
        .l <- list(sn1MG = list(unique(df_reactions$sn1MG)), df_template)
    
    if (reaction == "nae_to_fa")
        .l <- list(FA = list(unique(df_reactions$FA)), df_template)
    
    if (reaction == "nape_to_lnape")
        .l <- list(
           LNAPE = list(unique(df_reactions$LNAPE), FA = unique(df_reactions$FA)), 
           df_template)
    
    if (reaction == "nape_to_nae")
        .l <- list(
           list(NAE = unique(df_reactions$NAE), PA = unique(df_reactions$PA)), 
           df_template)
    
    if (reaction == "nape_to_pnae")
        .l <- list(
           list(NAE = unique(df_reactions$NAE), DG = unique(df_reactions$DG)), 
           df_template)
    
    if (reaction == "napeo_to_nae")
        .l <- list(
           list(NAE = unique(df_reactions$NAE), PAO = unique(df_reactions$PAO)),
           df_template)
    
    if (reaction == "pa_to_cdpdg")
        .l <- list(list(CDPDG = unique(df_reactions$CDPDG)), df_template)
    
    if (reaction == "pa_to_dg")
        .l <- list(DG = list(unique(df_reactions$DG)), df_template)
    
    if (reaction == "pao_to_dgo")
        .l <- list(list(DGO = unique(df_reactions$DGO)), df_template)
    
    if (reaction == "pc_to_dg")
        .l <- list(list(DG = unique(df_reactions$DG)), df_template)
    
    if (reaction == "pc_to_sn1lpc")
        .l <- list(list(sn1LPC = unique(df_reactions$sn1LPC)), df_template)
    
    if (reaction == "pc_to_sn2lpc")
        .l <- list(list(sn2LPC = unique(df_reactions$sn2LPC)), df_template)
    
    if (reaction == "pc_to_pa")
        .l <- list(list(PA = unique(df_reactions$PA)), df_template)
    
    if (reaction == "pc_to_ps")
        .l <- list(list(PS = unique(df_reactions$PS)), df_template)
    
    if (reaction == "pco_to_lpco")
        .l <- list(list(LPCO = unique(df_reactions$LPCO)), df_template)
    
    if (reaction == "pe_to_dg")
        .l <- list(list(DG = unique(df_reactions$DG)), df_template)
    
    if (reaction == "pe_to_sn1lpe")
        .l <- list(list(sn1LPE = unique(df_reactions$sn1LPE)), df_template)
    
    if (reaction == "pe_to_sn2lpe")
        .l <- list(list(sn2LPE = unique(df_reactions$sn2LPE)), df_template)
    
    if (reaction == "pe_to_nape_sn1")
        .l <- list(
           list(NAPE = unique(df_reactions$NAPE), LPC = unique(NAPEs$LPC)), 
           df_template)
    
    if (reaction == "pe_to_nape_sn2")
        .l <- list(
           list(NAPE = unique(df_reactions$NAPE), LPC = unique(NAPEs$LPC)), 
            df_template)
    
    if (reaction == "pe_to_pa")
        .l <- list(list(PA = unique(df_reactions$PA)), df_template)
    
    if (reaction == "pe_to_ps")
        .l <- list(list(PS = unique(df_reactions$PS)), df_template)
    
    if (reaction == "peo_to_lpeo")
        .l <- list(list(LPEO = unique(df_reactions$LPEO)), df_template)
    
    if (reaction == "peo_to_napeo_sn1")
        .l <- list(
            list(NAPEO = unique(df_reactions$NAPEO), LPC = unique(NAPEOs$LPC)),
            df_template)
    
    if (reaction == "peo_to_napeo_sn2")
        .l <- list(
            list(NAPEO = unique(df_reactions$NAPEO), LPC = unique(NAPEOs$LPC)), 
            df_template)
    
    if (reaction == "peo_to_pep")
        .l <- list(list(PEP = unique(df_reactions$PEP)), df_template)
    
    if (reaction == "pep_to_lpep")
        .l <- list(list(LPEP = unique(df_reactions$LPEP)), df_template)
    
    if (reaction == "pep_to_napep_sn1")
        .l <- list(
            list(NAPEP = unique(df_reactions$NAPEP), LPC = unique(df_reactions$LPC)), 
            df_template)
    
    if (reaction == "pep_to_napep_sn2")
        .l <- list(
            list(NAPEP = unique(df_reactions$NAPEP), LPC = unique(df_reactions$LPC)),
            df_template)
    
    if (reaction == "pg_to_cl")
        .l <- list(
            list(unique(NAPEP = df_reactions$NAPEP), LPC = unique(df_reactions$LPC)), 
            df_template)
    
    if (reaction == "pgp_to_pg")
        .l <- list(list(PG = unique(df_reactions$PG)), df_template)
    
    if (reaction == "ps_to_pe")
        .l <- list(list(PE = unique(df_reactions$PE)), df_template)
    
    if (reaction == "sm_to_cer")
        .l <- list(list(CER = unique(df_reactions$CER)), df_template)
    
    if (reaction == "sphinga_to_dhcer")
        .l <- list(list(DHCER = unique(df_reactions$DHCER)), df_template)
    
    if (reaction == "tg_to_dg")
        .l <- list(
            list(sn1Loss_dg = unique(DGs$sn1Loss_dg, sn1Loss_fa = DGs$sn1Loss_fa)), 
            df_template)
            
    ## return the list 
    .l
    
}