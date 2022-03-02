.create_template <- function(template = NA, reaction = "lpa_to_pa") {
    
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
            reaction_isReversible = "", reaction_geneAssociation = "",
            reaction_pathway = "")
        
        ## depending on the reaction, write specific reaction .formula to the
        ## entry reaction_.formula
        if (reaction == "acdhap_to_alkyldhap")
            .formula <- "M_adhap + M_alkylR1oh <=> M_h + M_fatacid + M_akdhap"
        
        if (reaction == "alkyldhap_to_lpao")
            .formula <- "M_h + M_nadph + M_akdhap <=> M_alkgp + M_nadp"
        
        if (reaction == "c1p_to_cer")
            .formula <- "M_h2o + M_c17isocrmp <=> M_pi + M_c17isocrm"
        
        if (reaction == "cdpdg_to_pgp")
            .formula <- "M_glyc3p + M_cdpdag <=> M_h + M_cmp + M_pgp"
        
        if (reaction == "cdpdg_to_pi")
            .formula <- "M_inost + M_cdpdag <=> M_h + M_cmp + M_pail"
        
        if (reaction == "cer_to_c1p")
            .formula <- "M_atp + M_c17isocrm <=> M_h + M_adp + M_c17isocrmp"
        
        if (reaction == "cer_to_glccer")
            .formula <- "M_udpg + M_c17isocrm <=> M_h + M_udp + M_c17isogluside"
        
        if (reaction == "cer_to_sm")
            .formula <- "M_pchol + M_c17isocrm <=> M_12dag + M_c17isosphmyln"
        
        if (reaction == "coa_to_acdhap")
            .formula <- "M_dhap + M_fataccoa <=> M_coa + M_adhap"
        
        if (reaction == "coa_to_fatoh")
            .formula <- "M_fataccoa + 2 M_nadph + 2 M_h <=> M_alkylR1oh + 2 M_nadp + M_coa"
        
        if (reaction == "coa_to_lpa")
            .formula <- "M_glyc3p + M_fataccoa <=> M_coa + M_alpa_pl"
        
        if (reaction == "dg_to_sn1mg")
            .formula <- "M_h2o + M_12dag <=> M_h + M_1magol + M_fatacid"
        
        if (reaction == "dg_to_sn2mg")
            .formula <- "M_h2o + M_12dag <=> M_h + M_mag + M_fatacid"
        
        if (reaction == "dg_to_pa")
            .formula <- "M_atp + M_12dag <=> M_h + M_adp + M_pa_pl"
        
        if (reaction == "dg_to_pc")
            .formula <- "M_cdpchol + M_12dag <=> M_h + M_cmp + M_pchol"
        
        if (reaction == "dg_to_pe")
            .formula <- "M_cdpea + M_12dag <=> M_h + M_cmp + M_pe"
        
        if (reaction == "dg_to_tg")
            .formula <- "M_fataccoa + M_12dag <=> M_coa + M_tag"
        
        if (reaction == "dgo_to_pco")
            .formula <- "M_cdpchol + M_akac2g <=> M_h + M_cmp + M_akac2gchol"
        
        if (reaction == "dgo_to_peo")
            .formula <- "M_cdpea + M_akac2g <=> M_h + M_cmp + M_akac2gpe"
        
        if (reaction == "dhcer_to_cer")
            .formula <- "M_h + M_nadh + M_o2 + M_c17isodhcrm <=> 2.0 M_h2o + M_nad + M_c17isocrm"
        
        if (reaction == "dhcer_to_dhsm")
            .formula <- "M_pchol + M_c17isodhcrm <=> M_12dag + M_c17isodhsphmyln"
        
        if (reaction == "dhsm_to_dhcer")
            .formula <- "M_h2o + M_c17isodhsphmyln <=> M_cholp + M_h + M_c17isodhcrm"
        
        if (reaction == "fa_to_coa")
            .formula <- "M_atp + M_coa + M_fatacid <=> M_ppi + M_amp + M_fataccoa"
        
        if (reaction == "lnape_to_gpnae")
            .formula <- "M_lnape + M_h2o <=> M_gpnae + M_fatacid"
        
        if (reaction == "lpa_to_pa")
            .formula <- "M_fataccoa + M_alpa_pl <=> M_coa + M_pa_pl"
        
        if (reaction == "lpao_to_pao")
            .formula <- "M_alkgp + M_fataccoa <=> M_akac2gp + M_coa"
        
        if (reaction == "sn1lpc_to_fa")
            .formula <- "M_h2o + M_ag3pc <=> M_g3pc + M_h + M_fatacid"
        
        if (reaction == "sn21pc_to_fa")
            .formula <- "M_h2o +  M_2agpc <=> M_g3pc + M_h + M_fatacid"
        
        if (reaction == "sn1lpc_to_pc")
            .formula <- "M_ag3pc + M_fataccoa <=> M_pchol + M_coa"
        
        if (reaction == "sn1lpe_to_fa")
            .formula <- "M_h2o + M_ag3pe <=> M_h + M_fatacid + M_g3pe"
        
        if (reaction == "sn2lpe_to_fa")
            .formula <- "M_h2o +  M_2agpe <=> M_h + M_fatacid + M_g3pe"
        
        if (reaction == "sn1lpe_to_pe")
            .formula <- "M_fataccoa + M_acg3pe <=> M_coa + M_pe"
        
        if (reaction == "lpeo_to_peo")
            .formula <- "M_fataccoa + M_ak2lgpe <=> M_coa + M_akac2gpe"
        
        if (reaction == "lpep_to_pep")
            .formula <- "M_fataccoa + M_alken2gpe <=> M_coa + M_alkenac2gpe"
        
        if (reaction == "sn1mg_to_dg")
            .formula <- "M_1magol + M_fataccoa <=> M_coa + M_12dag"
        
        if (reaction == "sn2mg_to_dg")
            .formula <- "M_fataccoa + M_mag <=> M_coa + M_12dag"
        
        if (reaction == "sn1mg_to_fa")
            .formula <- "M_h2o + M_magol <=> M_glyc + M_h + M_fatacid"
        
        if (reaction == "sn2mg_to_fa")
            .formula <-  "M_h2o + M_mag <=> M_glyc + M_h + M_fatacid"
            
        if (reaction == "sn1mg_to_lpa")
            .formula <- "M_atp + M_1magol <=> M_h + M_adp + M_alpa_pl"
        
        if (reaction == "sn2mg_to_sn1mg")
            .formula <- "M_mag <=> M_1magol"
        
        if (reaction == "nae_to_fa")
            .formula <- "M_h2o + M_nae <=> M_etha + M_h + M_fatacid"
        
        if (reaction == "nape_to_lnape")
            .formula <- "M_nape + M_h2o <=> M_lnape + M_fatacid"
        
        if (reaction == "nape_to_nae")
            .formula <- "M_nape + M_h2o <=> M_nae + M_pa_pl"
        
        if (reaction == "nape_to_pnae")
            .formula <- "M_nape + M_h2o <=> M_pnae + M_12dag"
        
        if (reaction == "napeo_to_nae")
            .formula <- "M_akac2nape + M_h2o <=> M_nae + M_akac2gp"
        
        if (reaction == "pa_to_cdpdg")
            .formula <- "M_ctp + M_pa_pl <=> M_ppi + M_cdpdag"
        
        if (reaction == "pa_to_dg")
            .formula <- "M_h2o + M_pa_pl <=> M_pi + M_12dag"
        
        if (reaction == "pao_to_dgo")
            .formula <- "M_h2o + M_akac2gp <=> M_pi + M_akac2g"
        
        if (reaction == "pc_to_dg")
            .formula <- "M_h2o + M_pchol <=> M_cholp + M_12dag"
        
        if (reaction == "pc_to_sn1lpc")
            .formula <- "M_h2o +  M_pchol <=> M_ag3pc + M_fatacid"
        
        if (reaction == "pc_to_sn2lpc")
            .formula <- "M_h2o + M_pchol <=> M_2agpc + M_fatacid"
        
        if (reaction == "pc_to_pa")
            .formula <- "M_h2o + M_pchol <=> M_chol + M_pa_pl"
        
        if (reaction == "pc_to_ps")
            .formula <- "M_ser_L + M_pchol <=> M_chol + M_ps"
        
        if (reaction == "pco_to_lpco")
            .formula <- "M_h2o + M_akac2gchol <=> M_h + M_ak2lgchol + M_fatacid"
        
        if (reaction == "pe_to_dg")
            .formula <- "M_h2o + M_pe <=> M_ethamp + M_12dag"
        
        if (reaction == "pe_to_sn1lpe")
            .formula <- "M_h2o + M_pe <=> M_h + M_fatacid + M_ag3pe"
        
        if (reaction == "pe_to_sn2lpe")
            .formula <- "M_h2o + M_pe <=> M_h + M_fatacid + M_2agpe"
        
        if (reaction == "pe_to_nape_sn1")
            .formula <- "M_pe + M_pc <=> M_nape + M_2agpc"
        
        if (reaction == "pe_to_nape_sn2")
            .formula <- "M_pe + M_pc <=> M_nape + M_ag3pc"
        
        if (reaction == "pe_to_pa")
            .formula <- "M_h2o + M_pe <=> M_etha + M_h + M_pa_pl"
        
        if (reaction == "pe_to_ps")
            .formula <- "M_ser_L + M_pe <=> M_etha + M_ps"
        
        if (reaction == "peo_to_lpeo")
            .formula <- "M_h2o + M_akac2gpe <=> M_h + M_ak2lgpe + M_fatacid"
        
        if (reaction == "peo_to_napeo_sn1")
            .formula <- "M_akac2gpe + M_pc <=> M_akac2nape + M_2agpc"
        
        if (reaction == "peo_to_napeo_sn2")
            .formula <- "M_akac2gpe + M_pc <=> M_akac2nape + M_ag3pc"
        
        if (reaction == "peo_to_pep")
            .formula <- "M_akac2gpe + M_nadph + M_h + M_o2 <=> M_alkenac2gpe + M_nadp + 2 M_h2o"
        
        if (reaction == "pep_to_lpep")
            .formula <- "M_h2o + M_alkenac2gpe <=> M_h + M_alken2gpe + M_fatacid"
        
        if (reaction == "pep_to_napep_sn1")
            .formula <- "M_alkenac2gpe + M_pc <=> M_alkenac2nape + M_2agpc"
        
        if (reaction == "pep_to_napep_sn2")
            .formula <- "M_alkenac2gpe + M_pc <=> M_alkenac2nape + M_ag3pc"
        
        if (reaction == "pg_to_cl")
            .formula <- "M_cdpdag + M_pg <=> M_h + M_cmp + M_clpn"
        
        if (reaction == "pgp_to_pg")
            .formula <- "M_h2o + M_pgp <=> M_pi + M_pg"
        
        if (reaction == "ps_to_pe")
            .formula <- "M_h + M_ps <=> M_co2 + M_pe"
        
        if (reaction == "sm_to_cer")
            .formula <- "M_h2o + M_c17isosphmyln <=> M_cholp + M_h + M_c17isocrm"
        
        if (reaction == "sphinga_to_dhcer")
            .formula <- "M_fataccoa + M_c17isosphgn <=> M_coa + M_c17isodhcrm"
        
        if (reaction == "tg_to_dg")
            .formula <- "M_h2o + M_tag <=> M_h + M_fatacid + M_12dag"
            
        template$reaction_formula <- .formula
    
    }
    
    ## return the template object
    template
}
