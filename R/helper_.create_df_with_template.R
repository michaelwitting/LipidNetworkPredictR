
#' @name create_df_with_template
#' 
#' @title Create data frame and substitute substrates/products 
.create_df_with_template <- function(df_reactions, template, reaction) {
 
    .check_reaction(reaction)
    
    ## 
    .char <- character(nrow(df_reactions))
    df <- data.frame(Name = .char, ReactionFormula = .char,
        IsReversible = .char, GeneAssociation = .char, Pathway = .char)
    
    ## fill with data
    df[["Name"]] <- template[["reaction_name"]]
    df[["ReactionFormula"]] <- template[["reaction_formula"]]
    df[["IsReversible"]] <- template[["reaction_isReversible"]]
    df[["GeneAssociation"]] <- template[["reaction_geneAssociation"]]
    df[["Pathway"]] <- template[["reaction_pathway"]]
    
    ## get the ReactionFormula
    .formula <- df$ReactionFormula
    
    ## depending on the reaction, adjust fixed and variable substrates and 
    ## products
    if (reaction == "acdhap_to_alkyldhap") {
        .formula <- gsub(x = .formula, pattern = "M_adhap",
            replacement = df_reactions[["ACDHAP"]])
        .formula <- gsub(x = .formula, pattern = "M_alkylR1oh", 
            replacement = df_reactions[["FATOH"]])
        .formula <- gsub(x = .formula, pattern = "M_fatacid", 
            replacement = df_reactions[["FA"]])
        .formula <- gsub(x = .formula, pattern = "M_akdhap", 
            replacement <- df_reactions[["ALKYLDHAP"]])
    }
    
    if (reaction == "alkyldhap_to_lpao") {
        .formula <- gsub(x = .formula, pattern = "M_akdhap", 
            replacement = df_reactions[["ALKYLDHAP"]])
        .formula <- gsub(x = .formula, pattern = "M_alpa_pl", 
            replacement = df_reactions[["LPAO"]])
    }
    
    if (reaction == "c1p_to_cer") {
        .formula <- gsub(x = .formula, pattern = "M_c17isocrm",
            replacement = df_reactions[["CER"]])
        .formula <- gsub(x = .formula, pattern = "M_c17isocrmp",
            replacement = df_reactions[["C1P"]])
    }
    
    if (reaction == "cdpdg_to_pgp") {
        .formula <- gsub(x = .formula, pattern = "M_cdpdag", 
            replacement = df_reactions[["CDPDG"]])
        .formula <- gsub(x = .formula, pattern = "M_pail",
            replacement = df_reactions[["PGP"]])
    }
    
    if (reaction == "cdpdg_to_pi") {
        .formula <- gsub(x = .formula, pattern = "M_cdpdag", 
            replacement = df_reactions[["CDPDG"]])
        .formula <- gsub(x = .formula, pattern = "M_pail", 
            replacement = df_reactions[["PI"]])
    }
    
    if (reaction == "cer_to_c1p") {
        .formula <- gsub(x = .formula, pattern = "M_c17isocrm",
            replacement = df_reactions[["CER"]])
        .formula <- gsub(x = .formula, pattern = "M_c17isocrmp",
            replacement = df_reactions[["C1P"]])
    }
    
    if (reaction == "cer_to_glccer") {
        .formula <- gsub(x = .formula, pattern = "M_c17isocrm", 
            replacement = df_reactions[["CER"]])
        .formula <- gsub(x = .formula, pattern = "M_c17isogluside", 
            replacement = df_reactions[["GLCCER"]])
    }
    
    if (reaction == "cer_to_sm") {
        .formula <- gsub(x = .formula, pattern = "M_c17isocrm", 
            replacement = df_reactions[["CER"]])
        .formula <- gsub(x = .formula, pattern = "M_c17isosphmyln", 
            replacement = df_reactions[["SM"]])
    }
    
    if (reaction == "coa_to_acdhap") {
        .formula <- gsub(x = .formula, pattern = "M_fataccoa",
            replacement = df_reactions[["CoA"]])
        .formula <- gsub(x = .formula, pattern = "M_adhap",
            replacement = df_reactions[["DHAP"]])
    }
    
    if (reaction == "coa_to_fatoh") {
        .formula <- gsub(x = .formula, pattern = "M_fataccoa",
            replacement = df_reactions[["CoA"]])
        .formula <- gsub(x = .formula, pattern = "M_alkylR1oh", 
            replacement = df_reactions[["FATOH"]])
    }
    
    if (reaction == "coa_to_lpa") {
        .formula <- gsub(x = .formula, pattern = "M_fataccoa",
            replacement = df_reactions[["CoA"]])
        .formula <- gsub(x = .formula, pattern = "M_alpa_pl",
            replacement = df_reactions[["LPA"]])
    }
    
    if (reaction == "dg_to_mg") {
        .formula <- gsub(x = .formula, pattern = "M_12dag",
            replacement = df_reactions[["DG"]])
        .formula <- gsub(x = .formula, pattern = "M_1magol",
            replacement = df_reactions[["MG"]])
        .formula <- gsub(x = .formula, pattern = "M_fatacid",
            replacement = df_reactions[["FA"]])
    }
    
    if (reaction == "dg_to_pa") {
        .formula <- gsub(x = .formula, pattern = "M_pa_pl",
            replacement = df_reactions[["PA"]])
        .formula <- gsub(x = .formula, pattern = "M_12dag",
            replacement = df_reactions[["DG"]])
    }
    
    if (reaction == "dg_to_pc") {
        .formula <- gsub(x = .formula, pattern = "M_12dag",
            replacement = df_reactions[["DG"]])
        .formula <- gsub(x = .formula, pattern = "M_pchol",
            replacement = df_reactions[["PC"]])
    }
    
    if (reaction == "dg_to_pe") {
        .formula <- gsub(x = .formula, pattern = "M_12dag",
            replacement = df_reactions[["DG"]])
        .formula <- gsub(x = .formula, pattern = "M_pe",
            replacement = df_reactions[["PE"]])
    }
    
    if (reaction == "dg_to_tg") {
        .formula <- gsub(x = .formula, pattern = "M_fataccoa",
            replacement = df_reactions[["CoA"]])
        .formula <- gsub(x = .formula, pattern = "M_12dag",
            replacement = df_reactions[["DG"]])
        .formula <- gsub(x = .formula, pattern = "M_tag",
            replacement = df_reactions[["TG"]])
    }
    
    if (reaction == "dgo_to_pco") {
        .formula <- gsub(x = .formula, pattern = "M_cdpchol",
            replacement = "CDP-Chol")
        .formula <- gsub(x = .formula, pattern = "M_akac2g",
            replacement = df_reactions[["DGO"]])
        .formula <- gsub(x = .formula, pattern = "M_akac2gchol",
            replacement = df_reactions[["PCO"]])
    }

    if (reaction == "dgo_to_peo") {
        .formula <- gsub(x = .formula, pattern = "M_akac2g",
            replacement = df_reactions[["DGO"]])
        .formula <- gsub(x = .formula, pattern = "M_akac2gpe",
            replacement = df_reactions[["PEO"]])
    }
    
    if (reaction == "dhcer_to_cer") {
        .formula <- gsub(x = .formula, pattern = "M_c17isodhcrm",
            replacement = df_reactions[["DHCER"]])
        .formula <- gsub(x = .formula, pattern = "M_c17isocrm",
            replacement = df_reactions[["CER"]])
    }
    
    if (reaction == "dhcer_to_dhsm") {
        .formula <- gsub(x = .formula, pattern = "M_c17isodhcrm",
            replacement = df_reactions[["DHCER"]])
        .formula <- gsub(x = .formula, pattern = "M_c17isodhsphmyln",
            replacement = df_reactions[["DHSM"]])
    }
    
    if (reaction == "dhsm_to_dhcer") {
        .formula <- gsub(x = .formula, pattern = "M_c17isodhcrm",
            replacement = df_reactions[["DHCER"]])
        .formula <- gsub(x = .formula, pattern = "M_c17isodhsphmyln",
            replacement = df_reactions[["DHSM"]])
    }
    
    if (reaction == "fa_to_coa") {
        .formula <- gsub(x = .formula, pattern = "M_fatacid",
            replacement = df_reactions[["FA"]])
        .formula <- gsub(x = .formula, pattern = "M_fataccoa",
            replacement = df_reactions[["CoA"]])
    }
    
    if (reaction == "lnape_to_gpnae") {
        .formula <- gsub(x = .formula, pattern = "M_lnape",
            replacement = df_reactions[["LNAPE"]])
        .formula <- gsub(x = .formula, pattern = "M_gpnae",
            replacement = df_reactions[["LNAPE"]])
        .formula <- gsub(x = .formula, pattern = "M_fatacid",
            replacement = df_reactions[["FA"]])
    }
    
    if (reaction == "lpa_to_pa") {
        .formula <- gsub(x = .formula, pattern = "M_fataccoa",
            replacement = df_reactions[["CoA"]])
        .formula <- gsub(x = .formula, pattern = "M_alpa_pl",
            replacement = df_reactions[["LPA"]])
        .formula <- gsub(x = .formula, pattern = "M_pa_pl",
            replacement = df_reactions[["PA"]])
    }
    
    if (reaction == "lpao_to_pao") {
        .formula <- gsub(x = .formula, pattern = "M_fataccoa",
            replacement = df_reactions[["CoA"]])
        .formula <- gsub(x = .formula, pattern = "M_alpa_pl",
            replacement = df_reactions[["LPAO"]])
        .formula <- gsub(x = .formula, pattern = "M_pa_pl",
            replacement = df_reactions[["PAO"]])
    }
    
    if (reaction == "sn1lpc_to_fa") {
        .formula <- gsub(x = .formula, pattern = "M_ag3pc",
            replacement = df_reactions[["sn1LPC"]])
        .formula <- gsub(x = .formula, pattern = "M_fatacid",
            replacement = df_reactions[["FA"]])
    }
    
    if (reaction == "sn2lpc_to_fa") {
        .formula <- gsub(x = .formula, pattern = "M_2agpc",
            replacement = df_reactions[["sn2LPC"]])
        .formula <- gsub(x = .formula, pattern = "M_fatacid",
            replacement = df_reactions[["FA"]])
        
    }
    
    if (reaction == "sn1lpc_to_pc") {
        
        .formula <- gsub(x = .formula, pattern = "M_fataccoa",
            replacement = df_reactions[["CoA"]])
        .formula <- gsub(x = .formula, pattern = "M_ag3pc",
            replacement = df_reactions[["sn1LPC"]])
        .formula <- gsub(x = .formula, pattern = "M_pchol",
            replacement = df_reactions[["PC"]])
    }
    
    if (reaction == "sn1lpe_to_fa") {
        .formula <- gsub(x = .formula, pattern = "M_ag3pe",
            replacement = df_reactions[["sn1LPE"]])
        .formula <- gsub(x = .formula, pattern = "M_fatacid",
            replacement = df_reactions[["FA"]])
    }
    
    if (reaction == "sn2lpe_to_fa") {
        .formula <- gsub(x = .formula, pattern = "M_2agpe",
            replacement = df_reactions[["sn2LPE"]])
        .formula <- gsub(x = .formula, pattern = "M_fatacid",
            replacement = df_reactions[["FA"]])
    }
    
    if (reaction == "sn1lpe_to_pe") {
        .formula <- gsub(x = .formula, pattern = "M_fataccoa",
            replacement = df_reactions[["CoA"]])
        .formula <- gsub(x = .formula, pattern = "M_acg3pe",
            replacement = df_reactions[["sn1LPE"]])
        .formula <- gsub(x = .formula, pattern = "M_pe",
            replacement = df_reactions[["PE"]])
    }
    
    if (reaction == "lpeo_to_peo") {
        .formula <- gsub(x = .formula, pattern = "M_fataccoa",
            replacement = df_reactions[["CoA"]])
        .formula <- gsub(x = .formula, pattern = "M_ak2lgpe",
            replacement = df_reactions[["LPEO"]])
        .formula <- gsub(x = .formula, pattern = "M_akac2gpe",
            replacement = df_reactions[["PEO"]])
    }
    
    if (reaction == "lpep_to_pep") {
        .formula <- gsub(x = .formula, pattern = "M_fataccoa",
            replacement = df_reactions[["CoA"]])
        .formula <- gsub(x = .formula, pattern = "M_alken2gpe",
            replacement = df_reactions[["LPEP"]])
        .formula <- gsub(x = .formula, pattern = "M_alkenac2gpe",
            replacement = df_reactions[["PEP"]])
    }
    
    if (reaction == "sn1mg_to_dg") {
        .formula <- gsub(x = .formula, pattern = "M_fataccoa",
            replacement = df_reactions[["CoA"]])
        .formula <- gsub(x = .formula, pattern = "M_1magol",
            replacement = df_reactions[["sn1MG"]])
        .formula <- gsub(x = .formula, pattern = "M_12dag",
            replacement = df_reactions[["DG"]])
    }
    
    if (reaction == "sn2mg_to_dg") {
        
        .formula <- gsub(x = .formula, pattern = "M_fataccoa",
            replacement = df_reactions[["CoA"]])
        .formula <- gsub(x = .formula, pattern = "M_mag",
            replacement = df_reactions[["sn2MG"]])
        .formula <- gsub(x = .formula, pattern = "M_12dag",
            replacement = df_reactions[["DG"]])
    }
    
    if (reaction == "sn1mg_to_fa") {
        .formula <- gsub(x = .formula, pattern = "M_magol",
            replacement = df_reactions[["sn1MG"]])
        .formula <- gsub(x = .formula, pattern = "M_fatacid",
            replacement = df_reactions[["FA"]])   
    }
    
    if (reaction == "sn2mg_to_fa") {
        .formula <- gsub(x = .formula, pattern = "M_mag",
            replacement = df_reactions[["sn2MG"]])
        .formula <- gsub(x = .formula, pattern = "M_fatacid",
            replacement = df_reactions[["FA"]])
    }
    
    if (reaction == "sn1mg_to_lpa") {
        .formula <- gsub(x = .formula, pattern = "M_1magol",
            replacement = df_reactions[["sn1MG"]])
        .formula <- gsub(x = .formula, pattern = "M_alpa_pl",
            replacement = df_reactions[["LPA"]])
    }
    
    if (reaction == "sn2mg_to_sn1mg") {
        .formula <- gsub(x = .formula, pattern = "M_mag",
            replacement = df_reactions[["sn2MG"]])
        .formula <- gsub(x = .formula, pattern = "M_1magol",
            replacement = df_reactions[["sn1MG"]])
    }
    
    if (reaction == "nae_to_fa") {
        .formula <- gsub(x = .formula, pattern = "M_nae",
            replacement = df_reactions[["NAE"]])
        .formula <- gsub(x = .formula, pattern = "M_fatacid",
            replacement = df_reactions[["FA"]])
    }
    
    if (reaction == "nape_to_lnape") {
        .formula <- gsub(x = .formula, pattern = "M_nape",
            replacement = df_reactions[["NAPE"]])
        .formula <- gsub(x = .formula, pattern = "M_lnape",
            replacement = df_reactions[["LNAPE"]])
        .formula <- gsub(x = .formula, pattern = "M_fatacid",
            replacement = df_reactions[["FA"]])
    }
    
    if (reaction == "nape_to_nae") {
        .formula <- gsub(x = .formula, pattern = "M_nape",
            replacement = df_reactions[["NAPE"]])
        .formula <- gsub(x = .formula, pattern = "M_nae",
            replacement = df_reactions[["NAE"]])
        .formula <- gsub(x = .formula, pattern = "M_pa_pl",
            replacement = df_reactions[["PA"]])
    }
    
    if (reaction == "nape_to_pnae") {
        .formula <- gsub(x = .formula, pattern = "M_nape",
            replacement = df_reactions[["NAPE"]])
        .formula <- gsub(x = .formula, pattern = "M_pnae",
            replacement = df_reactions[["PNAE"]])
        .formula <- gsub(x = .formula, pattern = "M_12dag",
            replacement = df_reactions[["DG"]])
    }
    
    if (reaction == "napeo_to_nae") {
        .formula <- gsub(x = .formula, pattern = "M_akac2nape",
            replacement = df_reactions[["NAPEO"]])
        .formula <- gsub(x = .formula, pattern = "M_nae",
            replacement = df_reactions[["NAE"]])
        .formula <- gsub(x = .formula, pattern = "M_akac2gp",
            replacement = df_reactions[["PAO"]])
    }
    
    if (reaction == "pa_to_cdpdg") {
        .formula <- gsub(x = .formula, pattern = "M_pa_pl",
            replacement = df_reactions[["PA"]])
        .formula <- gsub(x = .formula, pattern = "M_cdpdag",
            replacement = df_reactions[["CDPDG"]])
    }
    
    if (reaction == "pa_to_dg") {
        .formula <- gsub(x = .formula, pattern = "M_pa_pl",
            replacement = df_reactions[["PA"]])
        .formula <- gsub(x = .formula, pattern = "M_12dag",
            replacement = df_reactions[["DG"]])
    }
    
    if (reaction == "pao_to_dgo") {
        .formula <- gsub(x = .formula, pattern = "M_akac2gp",
            replacement = df_reactions[["PAO"]])
        .formula <- gsub(x = .formula, pattern = "M_akac2g",
            replacement = df_reactions[["DGO"]])
    }
    
    if (reaction == "pc_to_dg") {
        .formula <- gsub(x = .formula, pattern = "M_pchol", 
            replacement = df_reactions[["PC"]])
        .formula <- gsub(x = .formula, pattern = "M_12dag",
            replacement = df_reactions[["DG"]])
    }
    
    if (reaction == "pc_to_sn1lpc") {
        .formula <- gsub(x = .formula, pattern = "M_pchol",
            replacement = df_reactions[["PC"]])
        .formula <- gsub(x = .formula, pattern = "M_ag3pc",
            replacement = df_reactions[["sn1LPC"]])
        .formula <- gsub(x = .formula, pattern = "M_fatacid",
            replacement = df_reactions[["FA"]])
    }
    
    if (reaction == "pc_to_sn2lpc") {
        .formula <- gsub(x = .formula, pattern = "M_pchol",
            replacement = df_reactions[["PC"]])
        .formula <- gsub(x = .formula, pattern = "M_2agpc", 
            replacement = df_reactions[["sn2LPC"]])
        .formula <- gsub(x = .formula, pattern = "M_fatacid",
            replacement = df_reactions[["FA"]])
    }
    
    if (reaction == "pc_to_pa") {
        .formula <- gsub(x = .formula, pattern = "M_pchol",
            replacement = df_reactions[["PC"]])
        .formula <- gsub(x = .formula, pattern = "M_pa_pl",
            replacement = df_reactions[["PA"]])
    }
    
    if (reaction == "pc_to_ps") {
        .formula <- gsub(x = .formula, pattern = "M_pchol",
            replacement = df_reactions[["PC"]])
        .formula <- gsub(x = .formula, pattern = "M_ps",
            replacement = df_reactions[["PS"]])
    }
    
    if (reaction == "pco_to_lpco"){
        .formula <- gsub(x = .formula, pattern = "M_akac2gchol",
            replacement = df_reactions[["PCO"]])
        .formula <- gsub(x = .formula, pattern = "M_ak2lgchol",
            replacement = df_reactions[["LPCO"]])
        .formula <- gsub(x = .formula, pattern = "M_fatacid",
            replacement = df_reactions[["FA"]])
    } 
    
    if (reaction == "pe_to_dg") {
        .formula <- gsub(x = .formula, pattern = "M_pe",
            replacement = df_reactions[["PE"]])
        .formula <- gsub(x = .formula, pattern = "M_12dag",
            replacement = df_reactions[["DG"]])
    }
    
    if (reaction == "pe_to_sn1lpe") {
        .formula <- gsub(x = .formula, pattern = "M_pe",
            replacement = df_reactions[["PE"]])
        .formula <- gsub(x = .formula, pattern = "M_ag3pe",
            replacement = df_reactions[["sn1LPE"]])
        .formula <- gsub(x = .formula, pattern = "M_fatacid",
            replacement = df_reactions[["FA"]])
    }
    
    if (reaction == "pe_to_sn2lpe") {
        .formula <- gsub(x = .formula, pattern = "M_pe",
            replacement = df_reactions[["PE"]])
        .formula <- gsub(x = .formula, pattern = "M_2agpe",
            replacement = df_reactions[["sn2LPE"]])
        .formula <- gsub(x = .formula, pattern = "M_fatacid",
            replacement = df_reactions[["FA"]]
    }
    
    if (reaction == "pe_to_nape_sn1") {
        .formula <- gsub(x = .formula, pattern = "M_pc",
            replacement = df_reactions[["PC"]])
        .formula <- gsub(x = .formula, pattern = "M_pe",
            replacement = df_reactions[["PE"]])
        .formula <- gsub(x = .formula, pattern = "M_nape",
            replacement = df_reactions[["NAPE"]])
        .formula <- gsub(x = .formula, pattern = "M_2agpc",
            replacement = df_reactions[["LPC"]])
    }
    
    if (reaction == "pe_to_nape_sn2") {
        .formula <- gsub(x = .formula, pattern = "M_pc",
            replacement = df_reactions[["PC"]])
        .formula <- gsub(x = .formula, pattern = "M_pe",
            replacement = df_reactions[["PE"]])
        .formula <- gsub(x = .formula, pattern = "M_nape",
            replacement = df_reactions[["NAPE"]])
        .formula <- gsub(x = .formula, pattern = "M_ag3pc",
            replacement = df_reactions[["LPC"]])   
    }
    
    if (reaction == "pe_to_pa") {
        .formula <- gsub(x = .formula, pattern = "M_pe",
            replacement = df_reactions[["PE"]])
        .formula <- gsub(x = .formula, pattern = "M_pa_pl",
            replacement = df_reactions[["PA"]])
    }
    
    if (reaction == "pe_to_ps") {
        .formula <- gsub(x = .formula, pattern = "M_pchol",
            replacement = df_reactions[["PC"]])
        .formula <- gsub(x = .formula, pattern = "M_ps",
            replacement = df_reactions[["PS"]])
    }
    
    if (reaction == "peo_to_lpeo") {
        # adjust variable substrates and products
        .formula <- gsub(x = .formula, pattern = "M_akac2gpe",
            replacement = df_reactions[["PEO"]])
        .formula <- gsub(x = .formula, pattern = "M_ak2lgpe",
            replacement = df_reactions[["LPEO"]])
        .formula <- gsub(x = .formula, pattern = "M_fatacid",
            replacement = df_reactions[["FA"]])
    }
    
    if (reaction == "peo_to_napeo_sn1") {
        .formula <- gsub(x = .formula, pattern = "M_pc",
            replacement = df_reactions[["PC"]])
        .formula <- gsub(x = .formula, pattern = "M_akac2gpe",
            replacement = df_reactions[["PEO"]])
        .formula <- gsub(x = .formula, pattern = "M_akac2nape",
            replacement = df_reactions[["NAPEO"]])
        .formula <- gsub(x = .formula, pattern = "M_2agpc",
            replacement = df_reactions[["LPC"]])
    }
    
    if (reaction == "peo_to_napeo_sn2") {
        .formula <- gsub(x = .formula, pattern = "M_pc",
            replacement = df_reactions[["PC"]])
        .formula <- gsub(x = .formula, pattern = "M_akac2gpe",
            replacement = df_reactions[["PEO"]])
        .formula <- gsub(x = .formula, pattern = "M_akac2nape",
            replacement = df_reactions[["NAPEO"]])
        .formula <- gsub(x = .formula, pattern = "M_ag3pc", 
            replacement = df_reactions[["LPC"]])
    }
    
    if (reaction == "peo_to_pep") {
        .formula <- gsub(x = .formula, pattern = "M_akac2gpe",
            replacement = df_reactions[["PEO"]])
        .formula <- gsub(x = .formula, pattern = "M_alkenac2gpe_c",
            replacement = df_reactions[["PEP"]])
    }
    
    if (reaction == "pep_to_lpep") {
        .formula <- gsub(x = .formula, pattern = "M_alkenac2gpe",
            replacement = df_reactions[["PEP"]])
        .formula <- gsub(x = .formula, pattern = "M_alken2gpe",
            replacement = df_reactions[["LPEP"]])
        .formula <- gsub(x = .formula, pattern = "M_fatacid",
            replacement = df_reactions[["FA"]])
    }
    
    if (reaction == "pep_to_napep_sn1") {
        .formula <- gsub(x = .formula, pattern = "M_pc",
            replacement = df_reactions[["PC"]])
        .formula <- gsub(x = .formula, pattern = "M_alkenac2gpe",
            replacement = df_reactions[["PEP"]])
        .formula <- gsub(x = .formula, pattern = "M_alkenac2nape",
            replacement = df_reactions[["NAPEP"]])
        .formula <- gsub(x = .formula, pattern = "M_2agpc",
            replacement = df_reactions[["LPC"]])
    }
    
    if (reaction == "pep_to_napep_sn2") {
        .formula <- gsub(x = .formula, pattern = "M_pc",
            replacement = df_reactions[["PC"]])
        .formula <- gsub(x = .formula, pattern = "M_alkenac2gpe",
            replacement = df_reactions[["PEP"]])
        .formula <- gsub(x = .formula, pattern = "M_alkenac2nape",
            replacement = df_reactions[["NAPEP"]])
        .formula <- gsub(x = .formula, pattern = "M_ag3pc",
            replacement = df_reactions[["LPC"]])
    }
    
    if (reaction == "pg_to_cl") {
        .formula <- gsub(x = .formula, pattern = "M_cdpdag",
            replacement = df_reactions[["CDPDG"]])
        .formula <- gsub(x = .formula, pattern = "M_pg",
            replacement = df_reactions[["PG"]])
        .formula <- gsub(x = .formula, pattern = "M_clpn",
            replacement = df_reactions[["CL"]])
    }
    
    if (reaction == "pgp_to_pg") {
        .formula <- gsub(x = .formula, pattern = "M_pgp", 
            replacement = df_reactions[["PGP"]])
        .formula <- gsub(x = .formula, pattern = "M_PG", 
            replacement = df_reactions[["PG"]])
    }
    
    if (reaction == "ps_to_pe") {
        .formula <- gsub(x = .formula, pattern = "M_ps", 
            replacement = df_reactions[["PS"]])
        .formula <- gsub(x = .formula, pattern = "M_pe", 
            replacement = df_reactions[["PE"]])
    }
    
    if (reaction == "sm_to_cer") {
        .formula <- gsub(x = .formula, pattern = "M_c17isocrm", 
            replacement = df_reactions[["CER"]])
        .formula <- gsub(x = .formula, pattern = "M_c17isosphmyln", 
            replacement = df_reactions[["SM"]])
    }
    
    if (reaction == "sphinga_to_dhcer") {
        .formula <- gsub(x = .formula, pattern = "M_fataccoa", 
            replacement = df_reactions[["CoA"]])
        .formula <- gsub(x = .formula, pattern = "M_alpa_pl", 
            replacement = df_reactions[["DHCER"]])
    }
    
    if (reaction == "tg_to_dg") {
        .formula <- gsub(x = .formula, pattern = "M_tag",
            replacement = c(df_reactions[["TG"]], df_reactions[["TG"]]))
        .formula <- gsub(x = .formula, pattern = "M_12dag",
            replacement = c(df_reactions[["sn1Loss_dg"]], df_reactions[["sn3Loss_dg"]]))
        .formula <- gsub(x = .formula, pattern = "M_fatacid", 
            replacement = c(df_reactions[["sn1Loss_fa"]], df_reactions[["sn3Loss_fa"]]) 
    }
    
    ## adjust fixed substrates and products
    .formula <- gsub(x = .formula, pattern = "M_coa ", replacement = "CoA ")
    .formula <- gsub(x = .formula, pattern = "M_coa$", replacement = "CoA")
    .formula <- gsub(x = .formula, pattern = "M_cmp ", replacement = "CMP ")
    .formula <- gsub(x = .formula, pattern = "M_cmp$", replacement = "CMP")
    .formula <- gsub(x = .formula, pattern = "M_ctp ", replacement = "CTP ")
    .formula <- gsub(x = .formula, pattern = "M_ctp$", replacement = "CTP")
    .formula <- gsub(x = .formula, pattern = "M_amp ", replacement = "AMP ")
    .formula <- gsub(x = .formula, pattern = "M_amp$", replacement = "AMP")
    .formula <- gsub(x = .formula, pattern = "M_adp ", replacement = "ADP ")
    .formula <- gsub(x = .formula, pattern = "M_adp$", replacement = "ADP")
    .formula <- gsub(x = .formula, pattern = "M_atp ", replacement = "ATP ")
    .formula <- gsub(x = .formula, pattern = "M_atp$", replacement = "ATP")
    .formula <- gsub(x = .formula, pattern = "M_adp ", replacement = "ADP ")
    .formula <- gsub(x = .formula, pattern = "M_adp$", replacement = "ADP")
    .formula <- gsub(x = .formula, pattern = "M_udp ", replacement = "UDP ")
    .formula <- gsub(x = .formula, pattern = "M_udp$", replacement = "UDP")
    .formula <- gsub(x = .formula, pattern = "M_udpg ", replacement = "UDP-Glucose ")
    .formula <- gsub(x = .formula, pattern = "M_udpg$", replacement = "UDP-Glucose")
    .formula <- gsub(x = .formula, pattern = "M_nadh ", replacement = "NADH ")
    .formula <- gsub(x = .formula, pattern = "M_nadh$", replacement = "NADH")
    .formula <- gsub(x = .formula, pattern = "M_nad ", replacement = "NAD+ ")
    .formula <- gsub(x = .formula, pattern = "M_nad$", replacement = "NAD+")
    .formula <- gsub(x = .formula, pattern = "M_nadph ", replacement = "NADPH ")
    .formula <- gsub(x = .formula, pattern = "M_nadph$", replacement = "NADPH")
    .formula <- gsub(x = .formula, pattern = "M_nadp ", replacement = "NADP+ ")
    .formula <- gsub(x = .formula, pattern = "M_nadp$", replacement = "NADP+")
    .formula <- gsub(x = .formula, pattern = "M_o2 ", replacement = "O2 ")
    .formula <- gsub(x = .formula, pattern = "M_o2$", replacement = "O2")
    .formula <- gsub(x = .formula, pattern = "M_h2o ", replacement = "H2O ")
    .formula <- gsub(x = .formula, pattern = "M_h2o$", replacement = "H2O")
    .formula <- gsub(x = .formula, pattern = "M_co2 ", replacement = "CO2 ")
    .formula <- gsub(x = .formula, pattern = "M_co2$", replacement = "CO2")
    .formula <- gsub(x = .formula, pattern = "M_h ", replacement = "H+ ")
    .formula <- gsub(x = .formula, pattern = "M_h$", replacement = "H+")
    .formula <- gsub(x = .formula, pattern = "M_H ", replacement =  "H+ ")
    .formula <- gsub(x = .formula, pattern = "M_H$", replacement =  "H+")
    .formula <- gsub(x = .formula, pattern = "M_pi ", replacement = "Pi ")
    .formula <- gsub(x = .formula, pattern = "M_pi$", replacement = "Pi")
    .formula <- gsub(x = .formula, pattern = "M_ppi ", replacement = "PPi ")
    .formula <- gsub(x = .formula, pattern = "M_ppi$", replacement = "PPi")
    .formula <- gsub(x = .formula, pattern = "M_glyc ", replacement = "Glycerol ")
    .formula <- gsub(x = .formula, pattern = "M_glyc$", replacement = "Glycerol")
    .formula <- gsub(x = .formula, pattern = "M_glyc3p ", replacement = "Glycerol-3-P ")
    .formula <- gsub(x = .formula, pattern = "M_glyc3p$", replacement = "Glycerol-3-P")
    .formula <- gsub(x = .formula, pattern = "M_inost ", replacement = "Inositol ")
    .formula <- gsub(x = .formula, pattern = "M_inost$", replacement = "Inositol")
    .formula <- gsub(x = .formula, pattern = "M_pchol ", replacement = "PC ")
    .formula <- gsub(x = .formula, pattern = "M_pchol$", replacement = "PC")
    .formula <- gsub(x = .formula, pattern = "M_12dag ", replacement = "DG ")
    .formula <- gsub(x = .formula, pattern = "M_12dag$", replacement = "DG")
    .formula <- gsub(x = .formula, pattern = "M_dhap ", replacement = "Dihydroxyacetone-P ")
    .formula <- gsub(x = .formula, pattern = "M_dhap$", replacement = "Dihydroxyacetone-P")
    .formula <- gsub(x = .formula, pattern = "M_cdpea ", replacement = "CDP-Ethn ")
    .formula <- gsub(x = .formula, pattern = "M_cdpea$", replacement = "CDP-Ethn")
    .formula <- gsub(x = .formula, pattern = "M_cdpchol ", replacement = "CDP-Chol ")
    .formula <- gsub(x = .formula, pattern = "M_cdpchol$", replacement = "CDP-Chol")
    .formula <- gsub(x = .formula, pattern = "M_12dag ", replacement = "DG ")
    .formula <- gsub(x = .formula, pattern = "M_12dag$", replacement = "DG")
    .formula <- gsub(x = .formula, pattern = "M_ser_L ", replacement = "L-Serine ")
    .formula <- gsub(x = .formula, pattern = "M_ser_L$", replacement = "L-Serine")
    .formula <- gsub(x = .formula, pattern = "M_etha ", replacement = "Ethanolamine ")
    .formula <- gsub(x = .formula, pattern = "M_etha$", replacement = "Ethanolamine")
    .formula <- gsub(x = .formula, pattern = "M_c17isosphgn ", replacement = "SPH(d16:0(1OH,3OH)(15Me)) ")
    .formula <- gsub(x = .formula, pattern = "M_c17isosphgn$", replacement = "SPH(d16:0(1OH,3OH)(15Me))")
    .formula <- gsub(x = .formula, pattern = "M_chol ", replacement = "Choline ")
    .formula <- gsub(x = .formula, pattern = "M_chol$", replacement = "Choline")
    .formula <- gsub(x = .formula, pattern = "M_cholp ", replacement = "P-Choline ")
    .formula <- gsub(x = .formula, pattern = "M_cholp$", replacement = "P-Choline")
    .formula <- gsub(x = .formula, pattern = "M_g3pc ", replacement = "Glycerophosphoethanolamine ")
    .formula <- gsub(x = .formula, pattern = "M_g3pc$", replacement = "Glycerophosphoethanolamine")

    
    ## write back to entry ReactionFormula
    df$ReactionFormula <- .formula
    
    ## return the data.frame
    df
    
       
}