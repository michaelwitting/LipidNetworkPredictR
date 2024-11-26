## function .create_list_reactions_with_template
test_that(".create_list_reactions_with_template works", {
    
    ## acyldhap_to_alkyldhap
    acyldhap <- "DHAP(18:0)"
    fao <- "FAO(16:0)"
    substrates <- list(AcylDHAP = acyldhap, FAO = fao)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36171"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    
    expect_equal(l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(l[[1]]$FAO, "FAO(16:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AlkylDHAP, "DHAP(O-16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylDHAP + M_FAO = M_AlkylDHAP + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylDHAP", "M_FAO"))
    expect_equal(l[[2]]$reaction_product, c("M_AlkylDHAP", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57534 + CHEBI:17135 = CHEBI:73315 + CHEBI:57560 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57534", "CHEBI:17135"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:73315", "CHEBI:57560", "CHEBI:15378"))
    
    reaction <- "RHEA:36172"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(l[[1]]$FAO, "FAO(16:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AlkylDHAP, "DHAP(O-16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylDHAP + M_FAO => M_AlkylDHAP + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylDHAP", "M_FAO"))
    expect_equal(l[[2]]$reaction_product, c("M_AlkylDHAP", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57534 + CHEBI:17135 => CHEBI:73315 + CHEBI:57560 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57534", "CHEBI:17135"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:73315", "CHEBI:57560", "CHEBI:15378"))
    
    reaction <- "RHEA:36173"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(l[[1]]$FAO, "FAO(16:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AlkylDHAP, "DHAP(O-16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylDHAP + M_FAO <= M_AlkylDHAP + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylDHAP", "M_FAO"))
    expect_equal(l[[2]]$reaction_product, c("M_AlkylDHAP", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57534 + CHEBI:17135 <= CHEBI:73315 + CHEBI:57560 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57534", "CHEBI:17135"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:73315", "CHEBI:57560", "CHEBI:15378"))
    
    reaction <- "RHEA:36174"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(l[[1]]$FAO, "FAO(16:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AlkylDHAP, "DHAP(O-16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylDHAP + M_FAO <=> M_AlkylDHAP + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylDHAP", "M_FAO"))
    expect_equal(l[[2]]$reaction_product, c("M_AlkylDHAP", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57534 + CHEBI:17135 <=> CHEBI:73315 + CHEBI:57560 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57534", "CHEBI:17135"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:73315", "CHEBI:57560", "CHEBI:15378"))
    
    ## alkyldhap_to_lpao
    alkyldhap <- "DHAP(O-18:0)"
    substrates = list(AlkylDHAP = alkyldhap)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36175"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AlkylDHAP, "DHAP(O-18:0)")
    expect_equal(l[[1]]$sn1LPAO, "PA(O-18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AlkylDHAP + M_H+ + M_NADPH = M_LPA-O + M_NADP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AlkylDHAP", "M_H+", "M_NADPH"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA-O", "M_NADP"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:73315 + CHEBI:15378 + CHEBI:57783 = CHEBI:58014 + CHEBI:58349")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:73315", "CHEBI:15378", "CHEBI:57783"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58014", "CHEBI:58349"))
    
    reaction <- "RHEA:36176"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AlkylDHAP, "DHAP(O-18:0)")
    expect_equal(l[[1]]$sn1LPAO, "PA(O-18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AlkylDHAP + M_H+ + M_NADPH => M_LPA-O + M_NADP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AlkylDHAP", "M_H+", "M_NADPH"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA-O", "M_NADP"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:73315 + CHEBI:15378 + CHEBI:57783 => CHEBI:58014 + CHEBI:58349")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:73315", "CHEBI:15378", "CHEBI:57783"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58014", "CHEBI:58349"))
    
    reaction <- "RHEA:36177"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AlkylDHAP, "DHAP(O-18:0)")
    expect_equal(l[[1]]$sn1LPAO, "PA(O-18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AlkylDHAP + M_H+ + M_NADPH <= M_LPA-O + M_NADP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AlkylDHAP", "M_H+", "M_NADPH"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA-O", "M_NADP"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:73315 + CHEBI:15378 + CHEBI:57783 <= CHEBI:58014 + CHEBI:58349")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:73315", "CHEBI:15378", "CHEBI:57783"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58014", "CHEBI:58349"))
    
    reaction <- "RHEA:36178"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AlkylDHAP, "DHAP(O-18:0)")
    expect_equal(l[[1]]$sn1LPAO, "PA(O-18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AlkylDHAP + M_H+ + M_NADPH <=> M_LPA-O + M_NADP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AlkylDHAP", "M_H+", "M_NADPH"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA-O", "M_NADP"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:73315 + CHEBI:15378 + CHEBI:57783 <=> CHEBI:58014 + CHEBI:58349")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:73315", "CHEBI:15378", "CHEBI:57783"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58014", "CHEBI:58349"))
    
    ## cerp_to_cer
    cerp <- "CerP(16:0(3OH,4OH,15Me)/18:0)"
    substrates <- list(CerP = cerp)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:33743"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CerP, "CerP(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_CerP = M_Pi + M_Cer")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_CerP"))
    expect_equal(l[[2]]$reaction_product, c("M_Pi", "M_Cer"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:57674 = CHEBI:43474 + CHEBI:52639")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:57674"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:43474", "CHEBI:52639"))
    
    reaction <- "RHEA:33744"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CerP, "CerP(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_CerP => M_Pi + M_Cer")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_CerP"))
    expect_equal(l[[2]]$reaction_product, c("M_Pi", "M_Cer"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:57674 => CHEBI:43474 + CHEBI:52639")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:57674"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:43474", "CHEBI:52639"))
    
    reaction <- "RHEA:33745"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CerP, "CerP(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_CerP <= M_Pi + M_Cer")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_CerP"))
    expect_equal(l[[2]]$reaction_product, c("M_Pi", "M_Cer"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:57674 <= CHEBI:43474 + CHEBI:52639")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:57674"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:43474", "CHEBI:52639"))
    
    reaction <- "RHEA:33746"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CerP, "CerP(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_CerP <=> M_Pi + M_Cer")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_CerP"))
    expect_equal(l[[2]]$reaction_product, c("M_Pi", "M_Cer"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:57674 <=> CHEBI:43474 + CHEBI:52639")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:57674"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:43474", "CHEBI:52639"))
    
    ## cdpdg_to_pgp
    cdpdg <- "CDP-DG(12:0(11Me)/14:0)"
    substrates <- list(CDPDG = cdpdg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:12593"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(l[[1]]$PGP, "PGP(12:0(11Me)/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-DG + M_Glycerol-3-P = M_PGP + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-DG", "M_Glycerol-3-P"))
    expect_equal(l[[2]]$reaction_product, c("M_PGP", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58332 + CHEBI:57597 = CHEBI:60110 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58332", "CHEBI:57597"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:60110", "CHEBI:60377", "CHEBI:15378"))
    
    reaction <- "RHEA:12594"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(l[[1]]$PGP, "PGP(12:0(11Me)/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-DG + M_Glycerol-3-P => M_PGP + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-DG", "M_Glycerol-3-P"))
    expect_equal(l[[2]]$reaction_product, c("M_PGP", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58332 + CHEBI:57597 => CHEBI:60110 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58332", "CHEBI:57597"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:60110", "CHEBI:60377", "CHEBI:15378"))
    
    reaction <- "RHEA:12595"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(l[[1]]$PGP, "PGP(12:0(11Me)/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-DG + M_Glycerol-3-P <= M_PGP + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-DG", "M_Glycerol-3-P"))
    expect_equal(l[[2]]$reaction_product, c("M_PGP", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58332 + CHEBI:57597 <= CHEBI:60110 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58332", "CHEBI:57597"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:60110", "CHEBI:60377", "CHEBI:15378"))
    
    reaction <- "RHEA:12596"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(l[[1]]$PGP, "PGP(12:0(11Me)/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-DG + M_Glycerol-3-P <=> M_PGP + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-DG", "M_Glycerol-3-P"))
    expect_equal(l[[2]]$reaction_product, c("M_PGP", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58332 + CHEBI:57597 <=> CHEBI:60110 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58332", "CHEBI:57597"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:60110", "CHEBI:60377", "CHEBI:15378"))
    
    ## cdpdg_to_pi
    cdgdg <- "CDP-DG(12:0(11Me)/14:0)"
    substrates <- list(CDPDG = cdpdg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:11580"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(l[[1]]$PI, "PI(12:0(11Me)/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-DG + M_myo-Inositol = M_PI + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-DG", "M_myo-Inositol"))
    expect_equal(l[[2]]$reaction_product, c("M_PI", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58332 + CHEBI:17268 = CHEBI:57880 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58332", "CHEBI:17268"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57880", "CHEBI:60377", "CHEBI:15378"))
    
    reaction <- "RHEA:11581"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(l[[1]]$PI, "PI(12:0(11Me)/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-DG + M_myo-Inositol => M_PI + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-DG", "M_myo-Inositol"))
    expect_equal(l[[2]]$reaction_product, c("M_PI", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58332 + CHEBI:17268 => CHEBI:57880 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58332", "CHEBI:17268"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57880", "CHEBI:60377", "CHEBI:15378"))
    
    reaction <- "RHEA:11582"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(l[[1]]$PI, "PI(12:0(11Me)/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-DG + M_myo-Inositol <= M_PI + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-DG", "M_myo-Inositol"))
    expect_equal(l[[2]]$reaction_product, c("M_PI", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58332 + CHEBI:17268 <= CHEBI:57880 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58332", "CHEBI:17268"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57880", "CHEBI:60377", "CHEBI:15378"))
    
    reaction <- "RHEA:11583"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(l[[1]]$PI, "PI(12:0(11Me)/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-DG + M_myo-Inositol <=> M_PI + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-DG", "M_myo-Inositol"))
    expect_equal(l[[2]]$reaction_product, c("M_PI", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58332 + CHEBI:17268 <=> CHEBI:57880 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58332", "CHEBI:17268"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57880", "CHEBI:60377", "CHEBI:15378"))
    
    ## cer_to_cerp
    cer <- "Cer(16:0(3OH,4OH,15Me)/18:0)"
    substrates <- list(Cer = cer)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:17929"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[1]]$CerP, "CerP(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_Cer + M_ATP = M_ADP + M_CerP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_Cer", "M_ATP"))
    expect_equal(l[[2]]$reaction_product, c("M_ADP", "M_CerP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:52639 + CHEBI:30616 = CHEBI:456216 + CHEBI:57674 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:52639", "CHEBI:30616"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:456216", "CHEBI:57674", "CHEBI:15378"))
    
    reaction <- "RHEA:17930"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[1]]$CerP, "CerP(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_Cer + M_ATP => M_ADP + M_CerP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_Cer", "M_ATP"))
    expect_equal(l[[2]]$reaction_product, c("M_ADP", "M_CerP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:52639 + CHEBI:30616 => CHEBI:456216 + CHEBI:57674 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:52639", "CHEBI:30616"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:456216", "CHEBI:57674", "CHEBI:15378"))
    
    reaction <- "RHEA:17931"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[1]]$CerP, "CerP(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_Cer + M_ATP <= M_ADP + M_CerP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_Cer", "M_ATP"))
    expect_equal(l[[2]]$reaction_product, c("M_ADP", "M_CerP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:52639 + CHEBI:30616 <= CHEBI:456216 + CHEBI:57674 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:52639", "CHEBI:30616"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:456216", "CHEBI:57674", "CHEBI:15378"))
    
    reaction <- "RHEA:17932"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[1]]$CerP, "CerP(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_Cer + M_ATP <=> M_ADP + M_CerP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_Cer", "M_ATP"))
    expect_equal(l[[2]]$reaction_product, c("M_ADP", "M_CerP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:52639 + CHEBI:30616 <=> CHEBI:456216 + CHEBI:57674 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:52639", "CHEBI:30616"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:456216", "CHEBI:57674", "CHEBI:15378"))
    
    ## cer_to_glccer
    cer <- "Cer(16:0(3OH,4OH,15Me)/18:0)"
    substrates <- list(Cer = cer)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:12088"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[1]]$GlcCer, "GlcCer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_Cer + M_UDP-Glucose = M_GlcCer + M_H+ + M_UDP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_Cer", "M_UDP-Glucose"))
    expect_equal(l[[2]]$reaction_product, c("M_GlcCer",  "M_H+", "M_UDP"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:52639 + CHEBI:58885 = CHEBI:22801 + CHEBI:15378 + CHEBI:58223")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:52639", "CHEBI:58885"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:22801", "CHEBI:15378", "CHEBI:58223"))
    
    reaction <- "RHEA:12089"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[1]]$GlcCer, "GlcCer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_Cer + M_UDP-Glucose => M_GlcCer + M_H+ + M_UDP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_Cer", "M_UDP-Glucose"))
    expect_equal(l[[2]]$reaction_product, c("M_GlcCer",  "M_H+", "M_UDP"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:52639 + CHEBI:58885 => CHEBI:22801 + CHEBI:15378 + CHEBI:58223")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:52639", "CHEBI:58885"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:22801", "CHEBI:15378", "CHEBI:58223"))
    
    reaction <- "RHEA:12090"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[1]]$GlcCer, "GlcCer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_Cer + M_UDP-Glucose <= M_GlcCer + M_H+ + M_UDP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_Cer", "M_UDP-Glucose"))
    expect_equal(l[[2]]$reaction_product, c("M_GlcCer",  "M_H+", "M_UDP"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:52639 + CHEBI:58885 <= CHEBI:22801 + CHEBI:15378 + CHEBI:58223")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:52639", "CHEBI:58885"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:22801", "CHEBI:15378", "CHEBI:58223"))
    
    reaction <- "RHEA:12091"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[1]]$GlcCer, "GlcCer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_Cer + M_UDP-Glucose <=> M_GlcCer + M_H+ + M_UDP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_Cer", "M_UDP-Glucose"))
    expect_equal(l[[2]]$reaction_product, c("M_GlcCer",  "M_H+", "M_UDP"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:52639 + CHEBI:58885 <=> CHEBI:22801 + CHEBI:15378 + CHEBI:58223")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:52639", "CHEBI:58885"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:22801", "CHEBI:15378", "CHEBI:58223"))
    
    ## cer_to_sm
    cer <- "Cer(16:0(3OH,4OH,15Me)/18:0)"
    substrates <- list(Cer = cer)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <-  "RHEA:18765"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[1]]$SM, "SM(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_Cer = M_1,2-DG + M_SM")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_Cer"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_SM"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:52639 = CHEBI:17815 + CHEBI:17636")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:52639"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:17636"))
    
    reaction <-  "RHEA:18766"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[1]]$SM, "SM(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_Cer => M_1,2-DG + M_SM")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_Cer"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_SM"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:52639 => CHEBI:17815 + CHEBI:17636")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:52639"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:17636"))
    
    reaction <-  "RHEA:18767"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[1]]$SM, "SM(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_Cer <= M_1,2-DG + M_SM")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_Cer"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_SM"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:52639 <= CHEBI:17815 + CHEBI:17636")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:52639"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:17636"))
    
    reaction <-  "RHEA:18768"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[1]]$SM, "SM(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_Cer <=> M_1,2-DG + M_SM")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_Cer"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_SM"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:52639 <=> CHEBI:17815 + CHEBI:17636")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:52639"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:17636"))
    
    ## cl_to_lcl
    cl <- "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])"
    substrates <- list(CL = cl)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:32935"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[1]]$LCL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[1]]$FA, "FA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CL + M_H2O = M_1,2,4-LCL + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CL", "M_H2O" ))
    expect_equal(l[[2]]$reaction_product, c("M_1,2,4-LCL", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:62237 + CHEBI:15377 = CHEBI:64743 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:62237", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64743", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:32936"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[1]]$LCL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[1]]$FA, "FA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CL + M_H2O => M_1,2,4-LCL + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CL", "M_H2O" ))
    expect_equal(l[[2]]$reaction_product, c("M_1,2,4-LCL", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:62237 + CHEBI:15377 => CHEBI:64743 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:62237", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64743", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:32937"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[1]]$LCL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[1]]$FA, "FA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CL + M_H2O <= M_1,2,4-LCL + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CL", "M_H2O" ))
    expect_equal(l[[2]]$reaction_product, c("M_1,2,4-LCL", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:62237 + CHEBI:15377 <= CHEBI:64743 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:62237", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64743", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:32938"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[1]]$LCL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[1]]$FA, "FA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CL + M_H2O <=> M_1,2,4-LCL + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CL", "M_H2O" ))
    expect_equal(l[[2]]$reaction_product, c("M_1,2,4-LCL", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:62237 + CHEBI:15377 <=> CHEBI:64743 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:62237", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64743", "CHEBI:28868", "CHEBI:15378"))
    
    ## coa_to_acyldhap
    acylcoa <- "CoA(18:0)"
    substrates <- list(AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:17657"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_Dihydroxyacetone-P = M_AcylDHAP + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_Dihydroxyacetone-P"))
    expect_equal(l[[2]]$reaction_product, c("M_AcylDHAP", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57642 = CHEBI:57534 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58342", "CHEBI:57642"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57534", "CHEBI:57287"))
    
    reaction <- "RHEA:17658"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_Dihydroxyacetone-P => M_AcylDHAP + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_Dihydroxyacetone-P"))
    expect_equal(l[[2]]$reaction_product, c("M_AcylDHAP", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57642 => CHEBI:57534 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58342", "CHEBI:57642"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57534", "CHEBI:57287"))
    
    reaction <- "RHEA:17659"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_Dihydroxyacetone-P <= M_AcylDHAP + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_Dihydroxyacetone-P"))
    expect_equal(l[[2]]$reaction_product, c("M_AcylDHAP", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57642 <= CHEBI:57534 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58342", "CHEBI:57642"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57534", "CHEBI:57287"))
    
    reaction <- "RHEA:17660"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_Dihydroxyacetone-P <=> M_AcylDHAP + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_Dihydroxyacetone-P"))
    expect_equal(l[[2]]$reaction_product, c("M_AcylDHAP", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57642 <=> CHEBI:57534 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58342", "CHEBI:57642"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57534", "CHEBI:57287"))
    
    ## coa_to_FAO
    acylcoa <- "CoA(18:0)"
    substrates <- list(AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:52716"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$FAO, "FAO(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + 2 M_H+ + 2 M_NADPH = M_FAO + M_CoA + 2 M_NADP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "2 M_H+", "2 M_NADPH"))
    expect_equal(l[[2]]$reaction_product, c("M_FAO", "M_CoA", "2 M_NADP"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:83139 + 2 CHEBI:15378 + 2 CHEBI:57783 = CHEBI:77396 + CHEBI:57287 + 2 CHEBI:58349")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:83139", "2 CHEBI:15378", "2 CHEBI:57783"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77396", "CHEBI:57287", "2 CHEBI:58349"))
    
    reaction <- "RHEA:52717"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$FAO, "FAO(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + 2 M_H+ + 2 M_NADPH => M_FAO + M_CoA + 2 M_NADP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "2 M_H+", "2 M_NADPH"))
    expect_equal(l[[2]]$reaction_product, c("M_FAO", "M_CoA", "2 M_NADP"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:83139 + 2 CHEBI:15378 + 2 CHEBI:57783 => CHEBI:77396 + CHEBI:57287 + 2 CHEBI:58349")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:83139", "2 CHEBI:15378", "2 CHEBI:57783"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77396", "CHEBI:57287", "2 CHEBI:58349"))
    
    reaction <- "RHEA:52718"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$FAO, "FAO(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + 2 M_H+ + 2 M_NADPH <= M_FAO + M_CoA + 2 M_NADP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "2 M_H+", "2 M_NADPH"))
    expect_equal(l[[2]]$reaction_product, c("M_FAO", "M_CoA", "2 M_NADP"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:83139 + 2 CHEBI:15378 + 2 CHEBI:57783 <= CHEBI:77396 + CHEBI:57287 + 2 CHEBI:58349")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:83139", "2 CHEBI:15378", "2 CHEBI:57783"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77396", "CHEBI:57287", "2 CHEBI:58349"))
    
    reaction <- "RHEA:52719"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$FAO, "FAO(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + 2 M_H+ + 2 M_NADPH <=> M_FAO + M_CoA + 2 M_NADP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "2 M_H+", "2 M_NADPH"))
    expect_equal(l[[2]]$reaction_product, c("M_FAO", "M_CoA", "2 M_NADP"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:83139 + 2 CHEBI:15378 + 2 CHEBI:57783 <=> CHEBI:77396 + CHEBI:57287 + 2 CHEBI:58349")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:83139", "2 CHEBI:15378", "2 CHEBI:57783"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77396", "CHEBI:57287", "2 CHEBI:58349"))
    
    ## coa_to_lpa
    acylcoa <- "CoA(18:0)"
    substrates <- list(AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:15325"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_Glycerol-3-P = M_LPA + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_Glycerol-3-P"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57597 = CHEBI:57970 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58342", "CHEBI:57597"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57970", "CHEBI:57287"))
    
    reaction <- "RHEA:15326"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_Glycerol-3-P => M_LPA + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_Glycerol-3-P"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57597 => CHEBI:57970 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58342", "CHEBI:57597"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57970", "CHEBI:57287"))
    
    reaction <- "RHEA:15327"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_Glycerol-3-P <= M_LPA + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_Glycerol-3-P"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57597 <= CHEBI:57970 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58342", "CHEBI:57597"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57970", "CHEBI:57287"))
    
    reaction <- "RHEA:15328"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_Glycerol-3-P <=> M_LPA + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_Glycerol-3-P"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57597 <=> CHEBI:57970 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58342", "CHEBI:57597"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57970", "CHEBI:57287"))
    
    ## dg_to_sn1mg
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:44712"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_H2O = M_1-MG + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-MG", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:49172 + CHEBI:15377 = CHEBI:35759 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:49172", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:35759", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:44713"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_H2O => M_1-MG + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-MG", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:49172 + CHEBI:15377 => CHEBI:35759 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:49172", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:35759", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:44714"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_H2O <= M_1-MG + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-MG", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:49172 + CHEBI:15377 <= CHEBI:35759 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:49172", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:35759", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:44715"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_H2O <=> M_1-MG + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-MG", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:49172 + CHEBI:15377 <=> CHEBI:35759 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:49172", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:35759", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:35663"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_H2O = M_1-MG + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-MG", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:15377 = CHEBI:64683 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17815", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64683", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:35664"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_H2O => M_1-MG + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-MG", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:15377 => CHEBI:64683 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17815", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64683", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:35665"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), 
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_H2O <= M_1-MG + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-MG", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:15377 <= CHEBI:64683 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17815", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64683", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:35666"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_H2O <=> M_1-MG + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-MG", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:15377 <=> CHEBI:64683 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17815", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64683", "CHEBI:28868", "CHEBI:15378"))
    
    ## dg_to_sn2mg
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:33275"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn2MG, "MG(0:0/16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_H2O = M_2-MG + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_2-MG", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:17815 = CHEBI:17389 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:17815"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17389", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:33276"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn2MG, "MG(0:0/16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_H2O => M_2-MG + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_2-MG", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:17815 => CHEBI:17389 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:17815"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17389", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:33277"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn2MG, "MG(0:0/16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_H2O <= M_2-MG + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_2-MG", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:17815 <= CHEBI:17389 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:17815"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17389", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:33278"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn2MG, "MG(0:0/16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_H2O <=> M_2-MG + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_2-MG", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:17815 <=> CHEBI:17389 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:17815"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17389", "CHEBI:28868", "CHEBI:15378"))
    
    ## dg_to_pa
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:10272"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PA, "PA(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_ATP = M_PA + M_ADP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_ATP"))
    expect_equal(l[[2]]$reaction_product, c("M_PA", "M_ADP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:30616 = CHEBI:58608 + CHEBI:456216 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17815", "CHEBI:30616"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58608", "CHEBI:456216", "CHEBI:15378"))
    
    reaction <- "RHEA:10273"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PA, "PA(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_ATP => M_PA + M_ADP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_ATP"))
    expect_equal(l[[2]]$reaction_product, c("M_PA", "M_ADP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:30616 => CHEBI:58608 + CHEBI:456216 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17815", "CHEBI:30616"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58608", "CHEBI:456216", "CHEBI:15378"))
    
    reaction <- "RHEA:10274"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PA, "PA(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_ATP <= M_PA + M_ADP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_ATP"))
    expect_equal(l[[2]]$reaction_product, c("M_PA", "M_ADP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:30616 <= CHEBI:58608 + CHEBI:456216 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17815", "CHEBI:30616"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58608", "CHEBI:456216", "CHEBI:15378"))
    
    reaction <- "RHEA:10275"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PA, "PA(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_ATP <=> M_PA + M_ADP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_ATP"))
    expect_equal(l[[2]]$reaction_product, c("M_PA", "M_ADP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:30616 <=> CHEBI:58608 + CHEBI:456216 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17815", "CHEBI:30616"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58608", "CHEBI:456216", "CHEBI:15378"))
    
    ## dg_to_pc
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:32939"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PC, "PC(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_CDP-Choline = M_PC + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_CDP-Choline"))
    expect_equal(l[[2]]$reaction_product, c("M_PC", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58779 = CHEBI:57643 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17815", "CHEBI:58779"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57643", "CHEBI:60377", "CHEBI:15378"))
    
    reaction <- "RHEA:32940"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PC, "PC(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_CDP-Choline => M_PC + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_CDP-Choline"))
    expect_equal(l[[2]]$reaction_product, c("M_PC", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58779 => CHEBI:57643 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17815", "CHEBI:58779"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57643", "CHEBI:60377", "CHEBI:15378"))
    
    reaction <- "RHEA:32941"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PC, "PC(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_CDP-Choline <= M_PC + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_CDP-Choline"))
    expect_equal(l[[2]]$reaction_product, c("M_PC", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58779 <= CHEBI:57643 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17815", "CHEBI:58779"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57643", "CHEBI:60377", "CHEBI:15378"))
    
    reaction <- "RHEA:32942"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PC, "PC(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_CDP-Choline <=> M_PC + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_CDP-Choline"))
    expect_equal(l[[2]]$reaction_product, c("M_PC", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58779 <=> CHEBI:57643 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17815", "CHEBI:58779"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57643", "CHEBI:60377", "CHEBI:15378"))
    
    ## dg_to_pe
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:32943"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PE, "PE(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_CDP-Ethanolamine = M_PE + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_CDP-Ethanolamine"))
    expect_equal(l[[2]]$reaction_product, c("M_PE", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:57876 = CHEBI:64612 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17815", "CHEBI:57876"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64612", "CHEBI:60377", "CHEBI:15378"))
    
    reaction <- "RHEA:32944"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PE, "PE(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_CDP-Ethanolamine => M_PE + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_CDP-Ethanolamine"))
    expect_equal(l[[2]]$reaction_product, c("M_PE", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:57876 => CHEBI:64612 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17815", "CHEBI:57876"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64612", "CHEBI:60377", "CHEBI:15378"))
    
    reaction <- "RHEA:32945"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PE, "PE(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_CDP-Ethanolamine <= M_PE + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_CDP-Ethanolamine"))
    expect_equal(l[[2]]$reaction_product, c("M_PE", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:57876 <= CHEBI:64612 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17815", "CHEBI:57876"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64612", "CHEBI:60377", "CHEBI:15378"))
    
    reaction <- "RHEA:32946"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PE, "PE(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_CDP-Ethanolamine <=> M_PE + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_CDP-Ethanolamine"))
    expect_equal(l[[2]]$reaction_product, c("M_PE", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:57876 <=> CHEBI:64612 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17815", "CHEBI:57876"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64612", "CHEBI:60377", "CHEBI:15378"))
    
    ## dg_to_tg
    dg <- "DG(18:0/16:0/0:0)"
    acylcoa <- "CoA(14:0)"
    substrates <- list(DG = dg, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:10868"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_AcylCoA = M_TG + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_TG", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58342 = CHEBI:64615 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17815", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64615", "CHEBI:57287"))
    
    reaction <- "RHEA:10869"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_AcylCoA => M_TG + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_TG", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58342 => CHEBI:64615 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17815", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64615", "CHEBI:57287"))
    
    reaction <- "RHEA:10870"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_AcylCoA <= M_TG + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_TG", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58342 <= CHEBI:64615 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17815", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64615", "CHEBI:57287"))
    
    reaction <- "RHEA:10871"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1,2-DG + M_AcylCoA <=> M_TG + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2-DG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_TG", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58342 <=> CHEBI:64615 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17815", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64615", "CHEBI:57287"))
    
    ## dgo_to_pco
    dgo <- "DG(O-18:0/16:0/0:0)"
    substrates <- list(DGO = dgo)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36179"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(l[[1]]$PCO, "PC(O-18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_DG-O + M_CDP-Choline = M_PC-O + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_DG-O", "M_CDP-Choline"))
    expect_equal(l[[2]]$reaction_product, c("M_PC-O", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:52595 + CHEBI:58779 = CHEBI:36702 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:52595", "CHEBI:58779"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:36702", "CHEBI:60377", "CHEBI:15378"))
    
    reaction <- "RHEA:36180"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(l[[1]]$PCO, "PC(O-18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_DG-O + M_CDP-Choline => M_PC-O + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_DG-O", "M_CDP-Choline"))
    expect_equal(l[[2]]$reaction_product, c("M_PC-O", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:52595 + CHEBI:58779 => CHEBI:36702 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:52595", "CHEBI:58779"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:36702", "CHEBI:60377", "CHEBI:15378"))
    
    reaction <- "RHEA:36181"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(l[[1]]$PCO, "PC(O-18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_DG-O + M_CDP-Choline <= M_PC-O + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_DG-O", "M_CDP-Choline"))
    expect_equal(l[[2]]$reaction_product, c("M_PC-O", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:52595 + CHEBI:58779 <= CHEBI:36702 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:52595", "CHEBI:58779"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:36702", "CHEBI:60377", "CHEBI:15378"))
    
    reaction <- "RHEA:36182"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(l[[1]]$PCO, "PC(O-18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_DG-O + M_CDP-Choline <=> M_PC-O + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_DG-O", "M_CDP-Choline"))
    expect_equal(l[[2]]$reaction_product, c("M_PC-O", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:52595 + CHEBI:58779 <=> CHEBI:36702 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:52595", "CHEBI:58779"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:36702", "CHEBI:60377", "CHEBI:15378"))
    
    ## dgo_to_peo
    dgo <- "DG(O-18:0/16:0/0:0)"
    substrates <- list(DGO = dgo)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36187"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
       reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(l[[1]]$PEO, "PE(O-18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_DG-O + M_CDP-Ethanolamine = M_PE-O + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_DG-O", "M_CDP-Ethanolamine"))
    expect_equal(l[[2]]$reaction_product, c("M_PE-O", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:52595 + CHEBI:57876 = CHEBI:60520 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:52595", "CHEBI:57876"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:60520", "CHEBI:60377", "CHEBI:15378"))
    
    reaction <- "RHEA:36188"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(l[[1]]$PEO, "PE(O-18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_DG-O + M_CDP-Ethanolamine => M_PE-O + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_DG-O", "M_CDP-Ethanolamine"))
    expect_equal(l[[2]]$reaction_product, c("M_PE-O", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:52595 + CHEBI:57876 => CHEBI:60520 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:52595", "CHEBI:57876"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:60520", "CHEBI:60377", "CHEBI:15378"))
    
    reaction <- "RHEA:36189"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(l[[1]]$PEO, "PE(O-18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_DG-O + M_CDP-Ethanolamine <= M_PE-O + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_DG-O", "M_CDP-Ethanolamine"))
    expect_equal(l[[2]]$reaction_product, c("M_PE-O", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:52595 + CHEBI:57876 <= CHEBI:60520 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:52595", "CHEBI:57876"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:60520", "CHEBI:60377", "CHEBI:15378"))
    
    reaction <- "RHEA:36190"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(l[[1]]$PEO, "PE(O-18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_DG-O + M_CDP-Ethanolamine <=> M_PE-O + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_DG-O", "M_CDP-Ethanolamine"))
    expect_equal(l[[2]]$reaction_product, c("M_PE-O", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:52595 + CHEBI:57876 <=> CHEBI:60520 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:52595", "CHEBI:57876"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:60520", "CHEBI:60377", "CHEBI:15378"))
    
    ## dhcer_to_cer
    dhcer <- "Cer(d16:0(3OH,4OH)(15Me)/12:0)" ###################################################
    substrates <- list(DhCer = dhcer)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:46544"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DhCer, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(l[[1]]$Cer, "Cer(d16:1(4E)(3OH,4OH)(15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_DhCer + 2 M_Fe2+-cytochrome_b5 + 2 M_H+ + M_O2 = 2 M_Fe3+-cytochrome_b5 + M_Cer + 2 M_H2O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_DhCer", "2 M_Fe2+-cytochrome_b5", "2 M_H+", "M_O2"))
    expect_equal(l[[2]]$reaction_product, c("2 M_Fe3+-cytochrome_b5", "M_Cer", "2 M_H2O"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:31488 + 2 CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 = 2 CHEBI:57540 + CHEBI:52639 + 2 CHEBI:15377")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:31488", "2 CHEBI:29033", "2 CHEBI:15378", "CHEBI:15379"))
    expect_equal(l[[2]]$reaction_product_chebi, c("2 CHEBI:57540", "CHEBI:52639", "2 CHEBI:15377"))
    
    reaction <- "RHEA:46545"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DhCer, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(l[[1]]$Cer, "Cer(d16:1(4E)(3OH,4OH)(15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_DhCer + 2 M_Fe2+-cytochrome_b5 + 2 M_H+ + M_O2 => 2 M_Fe3+-cytochrome_b5 + M_Cer + 2 M_H2O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_DhCer", "2 M_Fe2+-cytochrome_b5", "2 M_H+", "M_O2"))
    expect_equal(l[[2]]$reaction_product, c("2 M_Fe3+-cytochrome_b5", "M_Cer", "2 M_H2O"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:31488 + 2 CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 => 2 CHEBI:57540 + CHEBI:52639 + 2 CHEBI:15377")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:31488", "2 CHEBI:29033", "2 CHEBI:15378", "CHEBI:15379"))
    expect_equal(l[[2]]$reaction_product_chebi, c("2 CHEBI:57540", "CHEBI:52639", "2 CHEBI:15377"))
    
    reaction <- "RHEA:46546"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DhCer, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(l[[1]]$Cer, "Cer(d16:1(4E)(3OH,4OH)(15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_DhCer + 2 M_Fe2+-cytochrome_b5 + 2 M_H+ + M_O2 <= 2 M_Fe3+-cytochrome_b5 + M_Cer + 2 M_H2O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_DhCer", "2 M_Fe2+-cytochrome_b5", "2 M_H+", "M_O2"))
    expect_equal(l[[2]]$reaction_product, c("2 M_Fe3+-cytochrome_b5", "M_Cer", "2 M_H2O"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:31488 + 2 CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 <= 2 CHEBI:57540 + CHEBI:52639 + 2 CHEBI:15377")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:31488", "2 CHEBI:29033", "2 CHEBI:15378", "CHEBI:15379"))
    expect_equal(l[[2]]$reaction_product_chebi, c("2 CHEBI:57540", "CHEBI:52639", "2 CHEBI:15377"))
    
    reaction <- "RHEA:46547"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DhCer, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(l[[1]]$Cer, "Cer(d16:1(4E)(3OH,4OH)(15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_DhCer + 2 M_Fe2+-cytochrome_b5 + 2 M_H+ + M_O2 <=> 2 M_Fe3+-cytochrome_b5 + M_Cer + 2 M_H2O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_DhCer", "2 M_Fe2+-cytochrome_b5", "2 M_H+", "M_O2"))
    expect_equal(l[[2]]$reaction_product, c("2 M_Fe3+-cytochrome_b5", "M_Cer", "2 M_H2O"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:31488 + 2 CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 <=> 2 CHEBI:57540 + CHEBI:52639 + 2 CHEBI:15377")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:31488", "2 CHEBI:29033", "2 CHEBI:15378", "CHEBI:15379"))
    expect_equal(l[[2]]$reaction_product_chebi, c("2 CHEBI:57540", "CHEBI:52639", "2 CHEBI:15377"))
    
    ## dhcer_to_dhsm
    dhcer <- "Cer(16:1(3OH,4OH,15Me)/12:0)"
    substrates <- list(DhCer = dhcer)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:44620"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DhCer, "Cer(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[1]]$DhSM, "SM(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_DhCer = M_1,2-DG + M_DhSM")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_DhCer"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_DhSM"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:31488 = CHEBI:17815 + CHEBI:67090")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:31488"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:67090"))
    
    reaction <- "RHEA:44621"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DhCer, "Cer(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[1]]$DhSM, "SM(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_DhCer => M_1,2-DG + M_DhSM")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_DhCer"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_DhSM"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:31488 => CHEBI:17815 + CHEBI:67090")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:31488"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:67090"))
    
    reaction <- "RHEA:44622"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DhCer, "Cer(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[1]]$DhSM, "SM(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_DhCer <= M_1,2-DG + M_DhSM")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_DhCer"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_DhSM"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:31488 <= CHEBI:17815 + CHEBI:67090")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:31488"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:67090"))
    
    reaction <- "RHEA:44623"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DhCer, "Cer(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[1]]$DhSM, "SM(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_DhCer <=> M_1,2-DG + M_DhSM")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_DhCer"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_DhSM"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:31488 <=> CHEBI:17815 + CHEBI:67090")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:31488"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:67090"))
    
    ## dhsm_to_dhcer
    dhsm <- "SM(16:1(3OH,4OH,15Me)/12:0)"
    substrates <- list(DhSM = dhsm)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:45300"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DhSM, "SM(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[1]]$DhCer, "Cer(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_DhSM + M_H2O = M_DhCer + M_H+ + M_Phosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_DhSM", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_DhCer", "M_H+", "M_Phosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64583 + CHEBI:15377 = CHEBI:83273 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64583", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:83273", "CHEBI:15378", "CHEBI:295975"))
    
    reaction <- "RHEA:45301"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DhSM, "SM(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[1]]$DhCer, "Cer(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_DhSM + M_H2O => M_DhCer + M_H+ + M_Phosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_DhSM", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_DhCer", "M_H+", "M_Phosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64583 + CHEBI:15377 => CHEBI:83273 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64583", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:83273", "CHEBI:15378", "CHEBI:295975"))
    
    reaction <- "RHEA:45302"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DhSM, "SM(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[1]]$DhCer, "Cer(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_DhSM + M_H2O <= M_DhCer + M_H+ + M_Phosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_DhSM", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_DhCer", "M_H+", "M_Phosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64583 + CHEBI:15377 <= CHEBI:83273 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64583", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:83273", "CHEBI:15378", "CHEBI:295975"))
    
    reaction <- "RHEA:45303"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DhSM, "SM(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[1]]$DhCer, "Cer(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_DhSM + M_H2O <=> M_DhCer + M_H+ + M_Phosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_DhSM", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_DhCer", "M_H+", "M_Phosphocholine"))    
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64583 + CHEBI:15377 <=> CHEBI:83273 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64583", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:83273", "CHEBI:15378", "CHEBI:295975"))
    
    ## fa_to_coa
    fa <- "FA(18:0)"
    substrates = list(FA = fa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:15421"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_FA + M_ATP + M_CoA = M_AcylCoA + M_AMP + M_PPi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_FA", "M_ATP", "M_CoA"))
    expect_equal(l[[2]]$reaction_product, c("M_AcylCoA", "M_AMP", "M_PPi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57560 + CHEBI:30616 + CHEBI:57287 = CHEBI:83139 + CHEBI:456215 + CHEBI:33019")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57560", "CHEBI:30616", "CHEBI:57287"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:83139", "CHEBI:456215", "CHEBI:33019"))
    
    reaction <- "RHEA:15422"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_FA + M_ATP + M_CoA => M_AcylCoA + M_AMP + M_PPi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_FA", "M_ATP", "M_CoA"))
    expect_equal(l[[2]]$reaction_product, c("M_AcylCoA", "M_AMP", "M_PPi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57560 + CHEBI:30616 + CHEBI:57287 => CHEBI:83139 + CHEBI:456215 + CHEBI:33019")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57560", "CHEBI:30616", "CHEBI:57287"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:83139", "CHEBI:456215", "CHEBI:33019"))
    
    reaction <- "RHEA:15423"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_FA + M_ATP + M_CoA <= M_AcylCoA + M_AMP + M_PPi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_FA", "M_ATP", "M_CoA"))
    expect_equal(l[[2]]$reaction_product, c("M_AcylCoA", "M_AMP", "M_PPi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57560 + CHEBI:30616 + CHEBI:57287 <= CHEBI:83139 + CHEBI:456215 + CHEBI:33019")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57560", "CHEBI:30616", "CHEBI:57287"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:83139", "CHEBI:456215", "CHEBI:33019"))
    
    reaction <- "RHEA:15424"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_FA + M_ATP + M_CoA <=> M_AcylCoA + M_AMP + M_PPi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_FA", "M_ATP", "M_CoA"))
    expect_equal(l[[2]]$reaction_product, c("M_AcylCoA", "M_AMP", "M_PPi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57560 + CHEBI:30616 + CHEBI:57287 <=> CHEBI:83139 + CHEBI:456215 + CHEBI:33019")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57560", "CHEBI:30616", "CHEBI:57287"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:83139", "CHEBI:456215", "CHEBI:33019"))
    
    reaction <- "RHEA:38883"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_FA + M_ATP + M_CoA = M_AcylCoA + M_AMP + M_PPi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_FA", "M_ATP", "M_CoA"))
    expect_equal(l[[2]]$reaction_product, c("M_AcylCoA", "M_AMP", "M_PPi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:28868 + CHEBI:30616 + CHEBI:57287 = CHEBI:77636 + CHEBI:456215 + CHEBI:33019")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:28868", "CHEBI:30616", "CHEBI:57287"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77636", "CHEBI:456215", "CHEBI:33019"))
    
    reaction <- "RHEA:38884"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_FA + M_ATP + M_CoA => M_AcylCoA + M_AMP + M_PPi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_FA", "M_ATP", "M_CoA"))
    expect_equal(l[[2]]$reaction_product, c("M_AcylCoA", "M_AMP", "M_PPi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:28868 + CHEBI:30616 + CHEBI:57287 => CHEBI:77636 + CHEBI:456215 + CHEBI:33019")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:28868", "CHEBI:30616", "CHEBI:57287"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77636", "CHEBI:456215", "CHEBI:33019"))
    
    reaction <- "RHEA:38885"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_FA + M_ATP + M_CoA <= M_AcylCoA + M_AMP + M_PPi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_FA", "M_ATP", "M_CoA"))
    expect_equal(l[[2]]$reaction_product, c("M_AcylCoA", "M_AMP", "M_PPi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:28868 + CHEBI:30616 + CHEBI:57287 <= CHEBI:77636 + CHEBI:456215 + CHEBI:33019")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:28868", "CHEBI:30616", "CHEBI:57287"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77636", "CHEBI:456215", "CHEBI:33019"))
    
    reaction <- "RHEA:38886"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_FA + M_ATP + M_CoA <=> M_AcylCoA + M_AMP + M_PPi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_FA", "M_ATP", "M_CoA"))
    expect_equal(l[[2]]$reaction_product, c("M_AcylCoA", "M_AMP", "M_PPi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:28868 + CHEBI:30616 + CHEBI:57287 <=> CHEBI:77636 + CHEBI:456215 + CHEBI:33019")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:28868", "CHEBI:30616", "CHEBI:57287"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77636", "CHEBI:456215", "CHEBI:33019"))
    
    ## lcl_to_cl
    lcl <- "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])"
    acylcoa <- "CoA(18:4(6Z,9Z,12Z,15Z))"
    substrates <- list(LCL = lcl, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:35839"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$LCL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_1,2,4-LCL + M_AcylCoA = M_CL + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2,4-LCL", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CL", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64743 + CHEBI:58342 = CHEBI:62237 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64743", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:62237", "CHEBI:57287"))
    
    reaction <- "RHEA:35840"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$LCL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_1,2,4-LCL + M_AcylCoA => M_CL + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2,4-LCL", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CL", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64743 + CHEBI:58342 => CHEBI:62237 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64743", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:62237", "CHEBI:57287"))
    
    reaction <- "RHEA:35841"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$LCL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_1,2,4-LCL + M_AcylCoA <= M_CL + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2,4-LCL", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CL", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64743 + CHEBI:58342 <= CHEBI:62237 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64743", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:62237", "CHEBI:57287"))
    
    reaction <- "RHEA:35842"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$LCL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_1,2,4-LCL + M_AcylCoA <=> M_CL + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1,2,4-LCL", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CL", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64743 + CHEBI:58342 <=> CHEBI:62237 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64743", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:62237", "CHEBI:57287"))
    
    ## lnape_to_gpnae
    lnape <- "NAPE(14:0/0:0/0:0)"
    substrates <- list(LNAPE = lnape)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:45420"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$LNAPE, "NAPE(14:0/0:0/0:0)")
    expect_equal(l[[1]]$GPNAE, "GPNAE(0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_H2O + M_LNAPE = M_FA + M_H+ + M_GPNAE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_LNAPE"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_H+", "M_GPNAE"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:85216 = CHEBI:28868 + CHEBI:15378 + CHEBI:85225")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:85216"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:15378", "CHEBI:85225"))
    
    reaction <- "RHEA:45421"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$LNAPE, "NAPE(14:0/0:0/0:0)")
    expect_equal(l[[1]]$GPNAE, "GPNAE(0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_H2O + M_LNAPE => M_FA + M_H+ + M_GPNAE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_LNAPE"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_H+", "M_GPNAE"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:85216 => CHEBI:28868 + CHEBI:15378 + CHEBI:85225")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:85216"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:15378", "CHEBI:85225"))
    
    reaction <- "RHEA:45422"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$LNAPE, "NAPE(14:0/0:0/0:0)")
    expect_equal(l[[1]]$GPNAE, "GPNAE(0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_H2O + M_LNAPE <= M_FA + M_H+ + M_GPNAE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_LNAPE"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_H+", "M_GPNAE"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:85216 <= CHEBI:28868 + CHEBI:15378 + CHEBI:85225")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:85216"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:15378", "CHEBI:85225"))
    
    reaction <- "RHEA:45423"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$LNAPE, "NAPE(14:0/0:0/0:0)")
    expect_equal(l[[1]]$GPNAE, "GPNAE(0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_H2O + M_LNAPE <=> M_FA + M_H+ + M_GPNAE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_LNAPE"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_H+", "M_GPNAE"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:85216 <=> CHEBI:28868 + CHEBI:15378 + CHEBI:85225")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:85216"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:15378", "CHEBI:85225"))
    
    ## lpa_to_pa
    lpa <- "PA(18:0/0:0)"
    acylcoa <- "CoA(14:0)"
    substrates <- list(sn1LPA = lpa, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:19709"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(l[[1]]$PA, "PA(18:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_LPA + M_AcylCoA = M_PA + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_LPA", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PA", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57970 + CHEBI:58342 = CHEBI:58608 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57970", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58608", "CHEBI:57287"))
    
    reaction <- "RHEA:19710"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(l[[1]]$PA, "PA(18:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_LPA + M_AcylCoA => M_PA + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_LPA", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PA", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57970 + CHEBI:58342 => CHEBI:58608 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57970", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58608", "CHEBI:57287"))
    
    reaction <- "RHEA:19711"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(l[[1]]$PA, "PA(18:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_LPA + M_AcylCoA <= M_PA + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_LPA", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PA", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57970 + CHEBI:58342 <= CHEBI:58608 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57970", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58608", "CHEBI:57287"))
    
    reaction <- "RHEA:19712"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(l[[1]]$PA, "PA(18:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_LPA + M_AcylCoA <=> M_PA + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_LPA", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PA", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57970 + CHEBI:58342 <=> CHEBI:58608 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57970", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58608", "CHEBI:57287"))
    
    ## lpao_to_pao
    lpao <- "PA(O-18:0/0:0)"
    acylcoa <- "CoA(14:0)"
    substrates <- list(sn1LPAO = lpao, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36235"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPAO, "PA(O-18:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(l[[1]]$PAO, "PA(O-18:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_LPA-O + M_AcylCoA = M_PA-O + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_LPA-O", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PA-O", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58014 + CHEBI:58342 = CHEBI:73332 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58014", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:73332", "CHEBI:57287"))
    
    reaction <- "RHEA:36236"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPAO, "PA(O-18:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(l[[1]]$PAO, "PA(O-18:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_LPA-O + M_AcylCoA => M_PA-O + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_LPA-O", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PA-O", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58014 + CHEBI:58342 => CHEBI:73332 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58014", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:73332", "CHEBI:57287"))
    
    reaction <- "RHEA:36237"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPAO, "PA(O-18:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(l[[1]]$PAO, "PA(O-18:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_LPA-O + M_AcylCoA <= M_PA-O + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_LPA-O", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PA-O", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58014 + CHEBI:58342 <= CHEBI:73332 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58014", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:73332", "CHEBI:57287"))
    
    reaction <- "RHEA:36238"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPAO, "PA(O-18:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(l[[1]]$PAO, "PA(O-18:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_LPA-O + M_AcylCoA <=> M_PA-O + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_LPA-O", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PA-O", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58014 + CHEBI:58342 <=> CHEBI:73332 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58014", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:73332", "CHEBI:57287"))
    
    ## sn1lpc_to_fa
    sn1lpc <- "PC(14:0/0:0)"
    substrates <- list(sn1LPC = sn1lpc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:15177"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPC, "PC(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPC + M_H2O = M_FA + M_H+ + M_Glycerophosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_H+", "M_Glycerophosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58168 + CHEBI:15377 = CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58168", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:15378", "CHEBI:16870"))
    
    reaction <- "RHEA:15178"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPC, "PC(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPC + M_H2O => M_FA + M_H+ + M_Glycerophosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_H+", "M_Glycerophosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58168 + CHEBI:15377 => CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58168", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:15378", "CHEBI:16870"))
    
    reaction <- "RHEA:15179"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPC, "PC(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPC + M_H2O <= M_FA + M_H+ + M_Glycerophosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_H+", "M_Glycerophosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58168 + CHEBI:15377 <= CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58168", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:15378", "CHEBI:16870"))
    
    reaction <- "RHEA:15180"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPC, "PC(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPC + M_H2O <=> M_FA + M_H+ + M_Glycerophosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_H+", "M_Glycerophosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58168 + CHEBI:15377 <=> CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58168", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:15378", "CHEBI:16870"))
    
    ## sn2lpc_to_fa
    sn2lpc <- "PC(0:0/14:0)"
    substrates <- list(sn2LPC = sn2lpc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:44696"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/14:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_2-LPC + M_H2O = M_FA + M_H+ + M_Glycerophosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_2-LPC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_H+", "M_Glycerophosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57875 + CHEBI:15377 = CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57875", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:15378", "CHEBI:16870"))
    
    reaction <- "RHEA:44697"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/14:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_2-LPC + M_H2O => M_FA + M_H+ + M_Glycerophosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_2-LPC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_H+", "M_Glycerophosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57875 + CHEBI:15377 => CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57875", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:15378", "CHEBI:16870"))
    
    reaction <- "RHEA:44698"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/14:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_2-LPC + M_H2O <= M_FA + M_H+ + M_Glycerophosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_2-LPC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_H+", "M_Glycerophosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57875 + CHEBI:15377 <= CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57875", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:15378", "CHEBI:16870"))
    
    reaction <- "RHEA:44699"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/14:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_2-LPC + M_H2O <=> M_FA + M_H+ + M_Glycerophosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_2-LPC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_H+", "M_Glycerophosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57875 + CHEBI:15377 <=> CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57875", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:15378", "CHEBI:16870"))
    
    ## sn1lpc_to_pc
    sn1lpc <- "PC(14:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPC = sn1lpc, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:12937"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPC, "PC(14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PC, "PC(14:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPC + M_AcylCoA = M_PC + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPC", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PC", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58168 + CHEBI:58342 = CHEBI:57643 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58168", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57643", "CHEBI:57287"))
    
    reaction <- "RHEA:12938"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPC, "PC(14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PC, "PC(14:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPC + M_AcylCoA => M_PC + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPC", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PC", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58168 + CHEBI:58342 => CHEBI:57643 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58168", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57643", "CHEBI:57287"))
    
    reaction <- "RHEA:12939"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPC, "PC(14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PC, "PC(14:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPC + M_AcylCoA <= M_PC + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPC", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PC", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58168 + CHEBI:58342 <= CHEBI:57643 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58168", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57643", "CHEBI:57287"))
    
    reaction <- "RHEA:12940"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPC, "PC(14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PC, "PC(14:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPC + M_AcylCoA <=> M_PC + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPC", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PC", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58168 + CHEBI:58342 <=> CHEBI:57643 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58168", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57643", "CHEBI:57287"))
    
    ## sn1lpe_to_fa
    pe <- "PE(14:0/0:0)"
    substrates <- list(sn1LPE = pe)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:32967"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPE + M_H2O = M_FA + M_H+ + M_Glycerophosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPE", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_H+", "M_Glycerophosphoethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64381 + CHEBI:15377 = CHEBI:28868 + CHEBI:15378 + CHEBI:143890")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64381", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:15378", "CHEBI:143890"))
    
    reaction <- "RHEA:32968"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPE + M_H2O => M_FA + M_H+ + M_Glycerophosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPE", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_H+", "M_Glycerophosphoethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64381 + CHEBI:15377 => CHEBI:28868 + CHEBI:15378 + CHEBI:143890")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64381", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:15378", "CHEBI:143890"))
    
    reaction <- "RHEA:32969"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPE + M_H2O <= M_FA + M_H+ + M_Glycerophosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPE", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_H+", "M_Glycerophosphoethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64381 + CHEBI:15377 <= CHEBI:28868 + CHEBI:15378 + CHEBI:143890")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64381", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:15378", "CHEBI:143890"))
    
    reaction <- "RHEA:32970"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPE + M_H2O <=> M_FA + M_H+ + M_Glycerophosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPE", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_H+", "M_Glycerophosphoethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64381 + CHEBI:15377 <=> CHEBI:28868 + CHEBI:15378 + CHEBI:143890")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64381", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:15378", "CHEBI:143890"))
    
    ## sn2lpe_to_fa
    pe <- "PE(0:0/14:0)"
    substrates <- list(sn2LPE = pe)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "sn2lpe_to_fa"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2LPE, "PE(0:0/14:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_2-LPE + M_H2O <=> M_FA + M_H+ + M_Glycerophosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_2-LPE", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_H+", "M_Glycerophosphoethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:65213 + CHEBI:15377 <=> CHEBI:28868 + CHEBI:15378 + CHEBI:143890")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:65213", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:15378", "CHEBI:143890"))
    
    ## sn1lpe_to_pe
    pe <- "PE(14:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPE = pe, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:32995"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PE, "PE(14:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPE + M_AcylCoA = M_PE + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPE", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PE", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64381 + CHEBI:58342 = CHEBI:64612 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64381", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64612", "CHEBI:57287"))
    
    reaction <- "RHEA:32996"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PE, "PE(14:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPE + M_AcylCoA => M_PE + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPE", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PE", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64381 + CHEBI:58342 => CHEBI:64612 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64381", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64612", "CHEBI:57287"))
    
    reaction <- "RHEA:32997"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PE, "PE(14:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPE + M_AcylCoA <= M_PE + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPE", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PE", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64381 + CHEBI:58342 <= CHEBI:64612 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64381", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64612", "CHEBI:57287"))
    
    reaction <- "RHEA:32998"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PE, "PE(14:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPE + M_AcylCoA <=> M_PE + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPE", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PE", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64381 + CHEBI:58342 <=> CHEBI:64612 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64381", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64612", "CHEBI:57287"))
    
    ## sn1lpi_to_pi
    sn1lpi <- "PI(16:0/0:0)"
    acylcoa <- "CoA(18:1(9Z))"
    substrates = list(sn1LPI = sn1lpi, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:33195"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:1(9Z))")
    expect_equal(l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPI + M_AcylCoA = M_PI + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPI", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PI", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64771 + CHEBI:58342 = CHEBI:57880 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64771", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57880", "CHEBI:57287"))
    
    reaction <- "RHEA:33196"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:1(9Z))")
    expect_equal(l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPI + M_AcylCoA => M_PI + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPI", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PI", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64771 + CHEBI:58342 => CHEBI:57880 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64771", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57880", "CHEBI:57287"))
    
    reaction <- "RHEA:33197"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:1(9Z))")
    expect_equal(l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPI + M_AcylCoA <= M_PI + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPI", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PI", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64771 + CHEBI:58342 <= CHEBI:57880 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64771", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57880", "CHEBI:57287"))
    
    reaction <- "RHEA:33198"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:1(9Z))")
    expect_equal(l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPI + M_AcylCoA <=> M_PI + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPI", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PI", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64771 + CHEBI:58342 <=> CHEBI:57880 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64771", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57880", "CHEBI:57287"))
    
    ## lpeo_to_peo
    lpeo <- "PE(O-16:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPEO = lpeo, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "lpeo_to_peo"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEO, "PE(O-16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PEO, "PE(O-16:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_LPE-O <=> M_CoA + M_PE-O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_LPE-O"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_PE-O"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58342 +  <=> CHEBI:57287 + CHEBI:75028") #############################
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58342", ""))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57287", "CHEBI:75028"))
    
    ## lpep_to_pep
    lpep <- "PE(P-16:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPEP = lpep, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:16245"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PEP, "PE(P-16:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_LPE-P + M_AcylCoA = M_PE-P + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_LPE-P", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PE-P", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:58342 = CHEBI:77290 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77288", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77290", "CHEBI:57287"))
    
    reaction <- "RHEA:16246"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PEP, "PE(P-16:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_LPE-P + M_AcylCoA => M_PE-P + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_LPE-P", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PE-P", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:58342 => CHEBI:77290 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77288", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77290", "CHEBI:57287"))
    
    reaction <- "RHEA:16247"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PEP, "PE(P-16:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_LPE-P + M_AcylCoA <= M_PE-P + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_LPE-P", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PE-P", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:58342 <= CHEBI:77290 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77288", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77290", "CHEBI:57287"))
    
    reaction <- "RHEA:16248"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PEP, "PE(P-16:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_LPE-P + M_AcylCoA <=> M_PE-P + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_LPE-P", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PE-P", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:58342 <=> CHEBI:77290 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77288", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77290", "CHEBI:57287"))
    
    ## sn1mg_to_dg
    sn1mg <- "MG(14:0/0:0/0:0)"
    acylcoa <- "CoA(16:0)"
    substrates <- list(sn1MG = sn1mg, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:38463"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_AcylCoA = M_1,2-DG + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64683 + CHEBI:58342 = CHEBI:17815 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64683", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:57287"))
    
    reaction <- "RHEA:38464"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_AcylCoA => M_1,2-DG + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64683 + CHEBI:58342 => CHEBI:17815 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64683", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:57287"))
    
    reaction <- "RHEA:38465"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_AcylCoA <= M_1,2-DG + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64683 + CHEBI:58342 <= CHEBI:17815 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64683", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:57287"))
    
    reaction <- "RHEA:38466"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_AcylCoA <=> M_1,2-DG + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64683 + CHEBI:58342 <=> CHEBI:17815 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64683", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:57287"))
    
    reaction <- "RHEA:39943"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_AcylCoA = M_1,2-DG + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:35759 + CHEBI:58342 = CHEBI:49172 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:35759", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:49172", "CHEBI:57287"))
    
    reaction <- "RHEA:39944"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_AcylCoA => M_1,2-DG + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:35759 + CHEBI:58342 => CHEBI:49172 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:35759", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:49172", "CHEBI:57287"))
    
    reaction <- "RHEA:39945"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_AcylCoA <= M_1,2-DG + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:35759 + CHEBI:58342 <= CHEBI:49172 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:35759", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:49172", "CHEBI:57287"))
    
    reaction <- "RHEA:39946"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_AcylCoA <=> M_1,2-DG + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:35759 + CHEBI:58342 <=> CHEBI:49172 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:35759", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:49172", "CHEBI:57287"))
    
    ## sn2mg_to_dg
    sn2mg <- "MG(0:0/14:0/0:0)"
    acylcoa <- "CoA(16:0)"
    substrates <- list(sn2MG = sn2mg, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:32947"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_2-MG + M_AcylCoA = M_1,2-DG + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_2-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 = CHEBI:17815 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17389", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:57287"))
    
    reaction <- "RHEA:32948"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_2-MG + M_AcylCoA => M_1,2-DG + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_2-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 => CHEBI:17815 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17389", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:57287"))
    
    reaction <- "RHEA:32949"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_2-MG + M_AcylCoA <= M_1,2-DG + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_2-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 <= CHEBI:17815 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17389", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:57287"))
    
    reaction <- "RHEA:32950"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_2-MG + M_AcylCoA <=> M_1,2-DG + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_2-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 <=> CHEBI:17815 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17389", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:57287"))
    
    reaction <- "RHEA:16741"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_2-MG + M_AcylCoA = M_1,2-DG + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_2-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 = CHEBI:49172 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17389", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:49172", "CHEBI:57287"))
    
    reaction <- "RHEA:16742"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_2-MG + M_AcylCoA => M_1,2-DG + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_2-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 => CHEBI:49172 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17389", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:49172", "CHEBI:57287"))
    
    reaction <- "RHEA:16743"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_2-MG + M_AcylCoA <= M_1,2-DG + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_2-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 <= CHEBI:49172 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17389", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:49172", "CHEBI:57287"))
    
    reaction <- "RHEA:16744"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_2-MG + M_AcylCoA <=> M_1,2-DG + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_2-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 <=> CHEBI:49172 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:17389", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:49172", "CHEBI:57287"))
    
    ## sn1mg_to_fa
    sn1mg <- "MG(14:0/0:0/0:0)"
    substrates <- list(sn1MG = sn1mg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:34019"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_H2O = M_FA + M_Glycerol + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_Glycerol", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:35759 + CHEBI:15377 = CHEBI:28868 + CHEBI:17754 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:35759", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:17754", "CHEBI:15378"))
    
    reaction <- "RHEA:34020"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_H2O => M_FA + M_Glycerol + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_Glycerol", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:35759 + CHEBI:15377 => CHEBI:28868 + CHEBI:17754 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:35759", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:17754", "CHEBI:15378"))
    
    reaction <- "RHEA:34021"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_H2O <= M_FA + M_Glycerol + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_Glycerol", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:35759 + CHEBI:15377 <= CHEBI:28868 + CHEBI:17754 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:35759", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:17754", "CHEBI:15378"))
    
    reaction <- "RHEA:34022"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_H2O <=> M_FA + M_Glycerol + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_Glycerol", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:35759 + CHEBI:15377 <=> CHEBI:28868 + CHEBI:17754 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:35759", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:17754", "CHEBI:15378"))
    
    ## sn2mg_to_fa
    sn2mg <- "MG(0:0/14:0/0:0)"
    substrates <- list(sn2MG = sn2mg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:32871"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_2-MG = M_FA + M_Glycerol + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_2-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_Glycerol", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:17389 = CHEBI:29067 + CHEBI:17754 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:17389"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:29067", "CHEBI:17754", "CHEBI:15378"))
    
    reaction <- "RHEA:32872"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_2-MG => M_FA + M_Glycerol + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_2-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_Glycerol", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:17389 => CHEBI:29067 + CHEBI:17754 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:17389"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:29067", "CHEBI:17754", "CHEBI:15378"))
    
    reaction <- "RHEA:32873"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_2-MG <= M_FA + M_Glycerol + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_2-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_Glycerol", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:17389 <= CHEBI:29067 + CHEBI:17754 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:17389"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:29067", "CHEBI:17754", "CHEBI:15378"))
    
    reaction <- "RHEA:32874"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_2-MG <=> M_FA + M_Glycerol + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_2-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_Glycerol", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:17389 <=> CHEBI:29067 + CHEBI:17754 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:17389"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:29067", "CHEBI:17754", "CHEBI:15378"))
    
    ## sn1mg_to_lpa
    sn1mg <- "MG(14:0/0:0/0:0)"
    substrates <- list(sn1MG = sn1mg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:33747"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$sn1LPA, "PA(14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_ATP = M_LPA + M_ADP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_ATP"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA", "M_ADP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64683 + CHEBI:30616 = CHEBI:57970 + CHEBI:456216 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64683", "CHEBI:30616"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57970", "CHEBI:456216", "CHEBI:15378"))
    
    reaction <- "RHEA:33748"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$sn1LPA, "PA(14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_ATP => M_LPA + M_ADP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_ATP"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA", "M_ADP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64683 + CHEBI:30616 => CHEBI:57970 + CHEBI:456216 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64683", "CHEBI:30616"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57970", "CHEBI:456216", "CHEBI:15378"))
    
    reaction <- "RHEA:33749"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$sn1LPA, "PA(14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_ATP <= M_LPA + M_ADP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_ATP"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA", "M_ADP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64683 + CHEBI:30616 <= CHEBI:57970 + CHEBI:456216 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64683", "CHEBI:30616"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57970", "CHEBI:456216", "CHEBI:15378"))
    
    reaction <- "RHEA:33750"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$sn1LPA, "PA(14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_ATP <=> M_LPA + M_ADP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_ATP"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA", "M_ADP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64683 + CHEBI:30616 <=> CHEBI:57970 + CHEBI:456216 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64683", "CHEBI:30616"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57970", "CHEBI:456216", "CHEBI:15378"))
    
    ## sn2mg_to_sn1mg
    sn2mg <- "MG(0:0/14:0/0:0)"
    substrates = list(sn2MG = sn2mg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    reaction <- "sn2mg_to_sn1mg"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_2-MG <=> M_1-MG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_2-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_1-MG"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:35759 <=> CHEBI:17389")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:35759"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17389"))
    
    ## nae_to_fa
    nae <- "NAE(18:0)"
    substrates <- list(NAE = nae)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:17505"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$NAE, "NAE(18:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_NAE = M_FA + M_Ethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_NAE"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_Ethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:15897 = CHEBI:57560 + CHEBI:57603")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:15897"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57560", "CHEBI:57603"))
    
    reaction <- "RHEA:17506"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$NAE, "NAE(18:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_NAE => M_FA + M_Ethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_NAE"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_Ethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:15897 => CHEBI:57560 + CHEBI:57603")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:15897"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57560", "CHEBI:57603"))
    
    reaction <- "RHEA:17507"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$NAE, "NAE(18:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_NAE <= M_FA + M_Ethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_NAE"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_Ethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:15897 <= CHEBI:57560 + CHEBI:57603")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:15897"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57560", "CHEBI:57603"))
    
    reaction <- "RHEA:17508"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$NAE, "NAE(18:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_NAE <=> M_FA + M_Ethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_NAE"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_Ethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:15897 <=> CHEBI:57560 + CHEBI:57603")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:15897"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57560", "CHEBI:57603"))
    
    reaction <- "RHEA:39995"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$NAE, "NAE(18:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_NAE = M_FA + M_Ethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_NAE"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_Ethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:52640 = CHEBI:28868 + CHEBI:57603")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:52640"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:57603"))
    
    reaction <- "RHEA:39996"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$NAE, "NAE(18:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_NAE => M_FA + M_Ethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_NAE"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_Ethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:52640 => CHEBI:28868 + CHEBI:57603")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:52640"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:57603"))
    
    reaction <- "RHEA:39997"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$NAE, "NAE(18:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_NAE <= M_FA + M_Ethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_NAE"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_Ethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:52640 <= CHEBI:28868 + CHEBI:57603")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:52640"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:57603"))
    
    reaction <- "RHEA:39998"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$NAE, "NAE(18:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_NAE <=> M_FA + M_Ethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_NAE"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_Ethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:52640 <=> CHEBI:28868 + CHEBI:57603")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:52640"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:57603"))
    
    ## nape_to_lnape
    nape <- "NAPE(14:0/16:0/18:0)"
    substrates = list(NAPE = nape)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:45460"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction, template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(l[[1]]$LNAPE, "NAPE(14:0/0:0/18:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_NAPE = M_FA + M_H+ + M_LNAPE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_NAPE"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_H+", "M_LNAPE"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 = CHEBI:28868 + CHEBI:15378 + CHEBI:85216")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:62537"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:15378", "CHEBI:85216"))
    
    reaction <- "RHEA:45461"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction, template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(l[[1]]$LNAPE, "NAPE(14:0/0:0/18:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_NAPE => M_FA + M_H+ + M_LNAPE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_NAPE"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_H+", "M_LNAPE"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 => CHEBI:28868 + CHEBI:15378 + CHEBI:85216")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:62537"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:15378", "CHEBI:85216"))
    
    reaction <- "RHEA:45462"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction, template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(l[[1]]$LNAPE, "NAPE(14:0/0:0/18:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_NAPE <= M_FA + M_H+ + M_LNAPE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_NAPE"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_H+", "M_LNAPE"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 <= CHEBI:28868 + CHEBI:15378 + CHEBI:85216")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:62537"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:15378", "CHEBI:85216"))
    
    reaction <- "RHEA:45463"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction, template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(l[[1]]$LNAPE, "NAPE(14:0/0:0/18:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_NAPE <=> M_FA + M_H+ + M_LNAPE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_NAPE"))
    expect_equal(l[[2]]$reaction_product, c("M_FA", "M_H+", "M_LNAPE"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 <=> CHEBI:28868 + CHEBI:15378 + CHEBI:85216")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:62537"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:28868", "CHEBI:15378", "CHEBI:85216"))
    
    ## nape_to_nae
    nape <- "NAPE(14:0/16:0/18:0)"
    substrates = list(NAPE = nape)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:33159"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction, template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(l[[1]]$NAE, "NAE(18:0)")
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_NAPE = M_PA + M_NAE + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_NAPE"))
    expect_equal(l[[2]]$reaction_product, c("M_PA", "M_NAE", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 = CHEBI:58608 + CHEBI:52640 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:62537"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58608", "CHEBI:52640", "CHEBI:15378"))
    
    reaction <- "RHEA:33160"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction, template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(l[[1]]$NAE, "NAE(18:0)")
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_NAPE => M_PA + M_NAE + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_NAPE"))
    expect_equal(l[[2]]$reaction_product, c("M_PA", "M_NAE", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 => CHEBI:58608 + CHEBI:52640 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:62537"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58608", "CHEBI:52640", "CHEBI:15378"))
    
    reaction <- "RHEA:33161"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction, template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(l[[1]]$NAE, "NAE(18:0)")
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_NAPE <= M_PA + M_NAE + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_NAPE"))
    expect_equal(l[[2]]$reaction_product, c("M_PA", "M_NAE", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 <= CHEBI:58608 + CHEBI:52640 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:62537"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58608", "CHEBI:52640", "CHEBI:15378"))
    
    reaction <- "RHEA:33162"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction, template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(l[[1]]$NAE, "NAE(18:0)")
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_NAPE <=> M_PA + M_NAE + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_NAPE"))
    expect_equal(l[[2]]$reaction_product, c("M_PA", "M_NAE", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 <=> CHEBI:58608 + CHEBI:52640 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:62537"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58608", "CHEBI:52640", "CHEBI:15378"))
    
    ## nape_to_pnae
    nape <- "NAPE(14:0/16:0/18:0)"
    substrates <- list(NAPE = nape)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "nape_to_pnae"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(l[[1]]$PNAE, "PNAE(18:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_NAPE + M_H2O <=> M_PNAE + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_NAPE", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_PNAE", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:62537 + CHEBI:15377 <=> CHEBI:145538 + CHEBI:17815")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:62537", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:145538", "CHEBI:17815"))
    
    ## napeo_to_nae
    napeo <- "NAPE(O-18:0/16:0/14:0)"
    substrates = list(NAPEO = napeo)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "napeo_to_nae"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$NAPEO, "NAPE(O-18:0/16:0/14:0)")
    expect_equal(l[[1]]$NAE, "NAE(14:0)")
    expect_equal(l[[1]]$PAO, "PA(O-18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_NAPEO + M_H2O <=> M_NAE + M_PA-O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_NAPEO", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_NAE", "M_PA-O"))
    expect_equal(l[[2]]$reaction_formula_chebi, " + CHEBI:15377 <=> CHEBI:52640 + CHEBI:73332") ################
    expect_equal(l[[2]]$reaction_substrate_chebi, c("", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:52640", "CHEBI:73332"))
    
    ## pa_to_cdpdg
    pa <- "PA(14:0/16:0)"
    substrates <- list(PA = pa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:16229"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[1]]$CDPDG, "CDP-DG(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PA + M_CTP + M_H+ = M_CDP-DG + M_PPi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PA", "M_CTP", "M_H+"))
    expect_equal(l[[2]]$reaction_product, c("M_CDP-DG", "M_PPi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58608 + CHEBI:37563 + CHEBI:15378 = CHEBI:58332 + CHEBI:33019")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58608", "CHEBI:37563", "CHEBI:15378"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58332", "CHEBI:33019"))
    
    reaction <- "RHEA:16230"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[1]]$CDPDG, "CDP-DG(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PA + M_CTP + M_H+ => M_CDP-DG + M_PPi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PA", "M_CTP", "M_H+"))
    expect_equal(l[[2]]$reaction_product, c("M_CDP-DG", "M_PPi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58608 + CHEBI:37563 + CHEBI:15378 => CHEBI:58332 + CHEBI:33019")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58608", "CHEBI:37563", "CHEBI:15378"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58332", "CHEBI:33019"))
    
    reaction <- "RHEA:16231"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[1]]$CDPDG, "CDP-DG(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PA + M_CTP + M_H+ <= M_CDP-DG + M_PPi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PA", "M_CTP", "M_H+"))
    expect_equal(l[[2]]$reaction_product, c("M_CDP-DG", "M_PPi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58608 + CHEBI:37563 + CHEBI:15378 <= CHEBI:58332 + CHEBI:33019")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58608", "CHEBI:37563", "CHEBI:15378"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58332", "CHEBI:33019"))
    
    reaction <- "RHEA:16232"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[1]]$CDPDG, "CDP-DG(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PA + M_CTP + M_H+ <=> M_CDP-DG + M_PPi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PA", "M_CTP", "M_H+"))
    expect_equal(l[[2]]$reaction_product, c("M_CDP-DG", "M_PPi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58608 + CHEBI:37563 + CHEBI:15378 <=> CHEBI:58332 + CHEBI:33019")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58608", "CHEBI:37563", "CHEBI:15378"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58332", "CHEBI:33019"))
    
    ## pa_to_dg
    pa <- "PA(14:0/16:0)"
    substrates = list(PA = pa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:27429"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PA + M_H2O = M_1,2-DG + M_Pi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PA", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_Pi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58608 + CHEBI:15377 = CHEBI:17815 + CHEBI:43474")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58608", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:43474"))
    
    reaction <- "RHEA:27430"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PA + M_H2O => M_1,2-DG + M_Pi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PA", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_Pi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58608 + CHEBI:15377 => CHEBI:17815 + CHEBI:43474")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58608", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:43474"))
    
    reaction <- "RHEA:27431"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PA + M_H2O <= M_1,2-DG + M_Pi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PA", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_Pi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58608 + CHEBI:15377 <= CHEBI:17815 + CHEBI:43474")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58608", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:43474"))
    
    reaction <- "RHEA:27432"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PA + M_H2O <=> M_1,2-DG + M_Pi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PA", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_Pi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:58608 + CHEBI:15377 <=> CHEBI:17815 + CHEBI:43474")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:58608", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:43474"))
    
    ## pao_to_dgo
    pao <- "PA(O-14:0/16:0)"
    substrates = list(PAO = pao)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36239"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PAO, "PA(O-14:0/16:0)")
    expect_equal(l[[1]]$DGO, "DG(O-14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PA-O + M_H2O = M_DG-O + M_Pi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PA-O", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_DG-O", "M_Pi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:73332 + CHEBI:15377 = CHEBI:52595 + CHEBI:43474")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:73332", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:52595", "CHEBI:43474"))
    
    reaction <- "RHEA:36240"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PAO, "PA(O-14:0/16:0)")
    expect_equal(l[[1]]$DGO, "DG(O-14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PA-O + M_H2O => M_DG-O + M_Pi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PA-O", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_DG-O", "M_Pi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:73332 + CHEBI:15377 => CHEBI:52595 + CHEBI:43474")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:73332", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:52595", "CHEBI:43474"))
    
    reaction <- "RHEA:36241"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PAO, "PA(O-14:0/16:0)")
    expect_equal(l[[1]]$DGO, "DG(O-14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PA-O + M_H2O <= M_DG-O + M_Pi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PA-O", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_DG-O", "M_Pi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:73332 + CHEBI:15377 <= CHEBI:52595 + CHEBI:43474")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:73332", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:52595", "CHEBI:43474"))
    
    reaction <- "RHEA:36242"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PAO, "PA(O-14:0/16:0)")
    expect_equal(l[[1]]$DGO, "DG(O-14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PA-O + M_H2O <=> M_DG-O + M_Pi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PA-O", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_DG-O", "M_Pi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:73332 + CHEBI:15377 <=> CHEBI:52595 + CHEBI:43474")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:73332", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:52595", "CHEBI:43474"))
    
    ## pc_to_dg
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:10604"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$DG, "DG(20:0/18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_H2O = M_1,2-DG + M_H+ + M_Phosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_H+", "M_Phosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 = CHEBI:17815 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:15378", "CHEBI:295975"))
    
    reaction <- "RHEA:10605"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$DG, "DG(20:0/18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_H2O => M_1,2-DG + M_H+ + M_Phosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_H+", "M_Phosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 => CHEBI:17815 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:15378", "CHEBI:295975"))
    
    reaction <- "RHEA:10606"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$DG, "DG(20:0/18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_H2O <= M_1,2-DG + M_H+ + M_Phosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_H+", "M_Phosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <= CHEBI:17815 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:15378", "CHEBI:295975"))
    
    reaction <- "RHEA:10607"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$DG, "DG(20:0/18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_H2O <=> M_1,2-DG + M_H+ + M_Phosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_H+", "M_Phosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <=> CHEBI:17815 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:15378", "CHEBI:295975"))
    
    ## pc_to_sn1lpc
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:15801"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn1LPC, "PC(20:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_H2O = M_1-LPC + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPC", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 = CHEBI:58168 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58168", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:15802"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn1LPC, "PC(20:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_H2O => M_1-LPC + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPC", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 => CHEBI:58168 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58168", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:15803"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn1LPC, "PC(20:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_H2O <= M_1-LPC + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPC", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <= CHEBI:58168 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58168", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:15804"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn1LPC, "PC(20:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_H2O <=> M_1-LPC + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPC", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <=> CHEBI:58168 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58168", "CHEBI:28868", "CHEBI:15378"))
    
    ## pc_to_sn2lpc
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:18689"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/18:0)")
    expect_equal(l[[1]]$FA, "FA(20:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_H2O = M_2-LPC + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_2-LPC", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 = CHEBI:57875 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57875", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:18690"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/18:0)")
    expect_equal(l[[1]]$FA, "FA(20:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_H2O => M_2-LPC + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_2-LPC", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 => CHEBI:57875 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57875", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:18691"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/18:0)")
    expect_equal(l[[1]]$FA, "FA(20:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_H2O <= M_2-LPC + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_2-LPC", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <= CHEBI:57875 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57875", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:18692"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/18:0)")
    expect_equal(l[[1]]$FA, "FA(20:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_H2O <=> M_2-LPC + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_2-LPC", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <=> CHEBI:57875 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57875", "CHEBI:28868", "CHEBI:15378"))
    
    ## pc_to_pa
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:14445"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$PA, "PA(20:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_H2O = M_PA + M_Choline + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_PA", "M_Choline", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 = CHEBI:58608 + CHEBI:15354 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58608", "CHEBI:15354", "CHEBI:15378"))
    
    reaction <- "RHEA:14446"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$PA, "PA(20:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_H2O => M_PA + M_Choline + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_PA", "M_Choline", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 => CHEBI:58608 + CHEBI:15354 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58608", "CHEBI:15354", "CHEBI:15378"))
    
    reaction <- "RHEA:14447"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$PA, "PA(20:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_H2O <= M_PA + M_Choline + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_PA", "M_Choline", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <= CHEBI:58608 + CHEBI:15354 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58608", "CHEBI:15354", "CHEBI:15378"))
    
    reaction <- "RHEA:14448"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$PA, "PA(20:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_H2O <=> M_PA + M_Choline + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_PA", "M_Choline", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <=> CHEBI:58608 + CHEBI:15354 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58608", "CHEBI:15354", "CHEBI:15378"))
    
    ## pc_to_ps
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:45088"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$PS, "PS(20:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_L-Serine = M_PS + M_Choline")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_L-Serine"))
    expect_equal(l[[2]]$reaction_product, c("M_PS", "M_Choline"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:33384 = CHEBI:57262 + CHEBI:15354")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:33384"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57262", "CHEBI:15354"))
    
    reaction <- "RHEA:45089"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$PS, "PS(20:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_L-Serine => M_PS + M_Choline")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_L-Serine"))
    expect_equal(l[[2]]$reaction_product, c("M_PS", "M_Choline"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:33384 => CHEBI:57262 + CHEBI:15354")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:33384"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57262", "CHEBI:15354"))
    
    reaction <- "RHEA:45090"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$PS, "PS(20:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_L-Serine <= M_PS + M_Choline")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_L-Serine"))
    expect_equal(l[[2]]$reaction_product, c("M_PS", "M_Choline"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:33384 <= CHEBI:57262 + CHEBI:15354")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:33384"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57262", "CHEBI:15354"))
    
    reaction <- "RHEA:45091"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$PS, "PS(20:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_L-Serine <=> M_PS + M_Choline")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_L-Serine"))
    expect_equal(l[[2]]$reaction_product, c("M_PS", "M_Choline"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:33384 <=> CHEBI:57262 + CHEBI:15354")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:33384"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57262", "CHEBI:15354"))
    
    ## pco_to_lpco
    pco <- "PC(O-16:0/14:0)"
    substrates <- list(PCO = pco)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36231"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PCO, "PC(O-16:0/14:0)")
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC-O + M_H2O = M_LPC-O + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC-O", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_LPC-O", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:36702 + CHEBI:15377 = CHEBI:30909 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:36702", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:30909", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:36232"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PCO, "PC(O-16:0/14:0)")
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC-O + M_H2O => M_LPC-O + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC-O", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_LPC-O", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:36702 + CHEBI:15377 => CHEBI:30909 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:36702", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:30909", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:36233"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PCO, "PC(O-16:0/14:0)")
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC-O + M_H2O <= M_LPC-O + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC-O", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_LPC-O", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:36702 + CHEBI:15377 <= CHEBI:30909 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:36702", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:30909", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:36234"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PCO, "PC(O-16:0/14:0)")
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC-O + M_H2O <=> M_LPC-O + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC-O", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_LPC-O", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:36702 + CHEBI:15377 <=> CHEBI:30909 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:36702", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:30909", "CHEBI:28868", "CHEBI:15378"))
    
    ## lpco_to_lpao
    sn1lpco <- "PC(O-16:0/0:0)"
    substrates <- list(sn1LPCO = sn1lpco)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:39927"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$sn1LPAO, "PA(O-16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPC-O + M_H2O = M_1-LPA-O + M_Choline + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPC-O", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPA-O", "M_Choline", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 = CHEBI:58014 + CHEBI:15354 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:30909", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58014", "CHEBI:15354", "CHEBI:15378"))
    
    reaction <- "RHEA:39928"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$sn1LPAO, "PA(O-16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPC-O + M_H2O => M_1-LPA-O + M_Choline + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPC-O", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPA-O", "M_Choline", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 => CHEBI:58014 + CHEBI:15354 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:30909", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58014", "CHEBI:15354", "CHEBI:15378"))
    
    reaction <- "RHEA:39929"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$sn1LPAO, "PA(O-16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPC-O + M_H2O <= M_1-LPA-O + M_Choline + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPC-O", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPA-O", "M_Choline", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 <= CHEBI:58014 + CHEBI:15354 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:30909", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58014", "CHEBI:15354", "CHEBI:15378"))
    
    reaction <- "RHEA:39930"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$sn1LPAO, "PA(O-16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPC-O + M_H2O <=> M_1-LPA-O + M_Choline + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPC-O", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPA-O", "M_Choline", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 <=> CHEBI:58014 + CHEBI:15354 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:30909", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58014", "CHEBI:15354", "CHEBI:15378"))
    
    ## lpco_to_mgo
    sn1lpco <- "PC(O-16:0/0:0)"
    substrates <- list(sn1LPCO = sn1lpco)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36083"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$sn1MGO, "MG(O-16:0/0:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPC-O + M_H2O = M_1-MG-O + M_H+ + M_Phosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPC-O", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-MG-O", "M_H+", "M_Phosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 = CHEBI:15850 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:30909", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:15850", "CHEBI:15378", "CHEBI:295975"))
    
    reaction <- "RHEA:36084"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$sn1MGO, "MG(O-16:0/0:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPC-O + M_H2O => M_1-MG-O + M_H+ + M_Phosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPC-O", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-MG-O", "M_H+", "M_Phosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 => CHEBI:15850 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:30909", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:15850", "CHEBI:15378", "CHEBI:295975"))
    
    reaction <- "RHEA:36085"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$sn1MGO, "MG(O-16:0/0:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPC-O + M_H2O <= M_1-MG-O + M_H+ + M_Phosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPC-O", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-MG-O", "M_H+", "M_Phosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 <= CHEBI:15850 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:30909", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:15850", "CHEBI:15378", "CHEBI:295975"))
    
    reaction <- "RHEA:36086"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$sn1MGO, "MG(O-16:0/0:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPC-O + M_H2O <=> M_1-MG-O + M_H+ + M_Phosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPC-O", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-MG-O", "M_H+", "M_Phosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 <=> CHEBI:15850 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:30909", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:15850", "CHEBI:15378", "CHEBI:295975"))
    
    ## lpco_to_pco
    sn1lpco <- "PC(O-16:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPCO = sn1lpco, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:23992"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PCO, "PC(O-16:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPC-O + M_AcylCoA = M_PC-O + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPC-O", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PC-O", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:58342 = CHEBI:36702 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:30909", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:36702", "CHEBI:57287"))
    
    reaction <- "RHEA:23993"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PCO, "PC(O-16:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPC-O + M_AcylCoA => M_PC-O + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPC-O", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PC-O", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:58342 => CHEBI:36702 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:30909", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:36702", "CHEBI:57287"))
    
    reaction <- "RHEA:23994"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PCO, "PC(O-16:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPC-O + M_AcylCoA <= M_PC-O + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPC-O", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PC-O", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:58342 <= CHEBI:36702 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:30909", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:36702", "CHEBI:57287"))
    
    reaction <- "RHEA:23995"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PCO, "PC(O-16:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPC-O + M_AcylCoA <=> M_PC-O + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPC-O", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_PC-O", "M_CoA"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:58342 <=> CHEBI:36702 + CHEBI:57287")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:30909", "CHEBI:58342"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:36702", "CHEBI:57287"))
    
    ## pe_to_dg
    pe <- "PE(14:0/16:0)"
    substrates <- list(PE = pe)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:78951"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PE = M_P-Ethanolamine + M_1,2-DG + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_P-Ethanolamine", "M_1,2-DG", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:64612 = CHEBI:58190 + CHEBI:17815 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:64612"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58190", "CHEBI:17815", "CHEBI:15378"))
    
    reaction <- "RHEA:78952"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PE => M_P-Ethanolamine + M_1,2-DG + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_P-Ethanolamine", "M_1,2-DG", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:64612 => CHEBI:58190 + CHEBI:17815 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:64612"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58190", "CHEBI:17815", "CHEBI:15378"))
    
    reaction <- "RHEA:78953"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PE <= M_P-Ethanolamine + M_1,2-DG + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_P-Ethanolamine", "M_1,2-DG", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:64612 <= CHEBI:58190 + CHEBI:17815 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:64612"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58190", "CHEBI:17815", "CHEBI:15378"))
    
    reaction <- "RHEA:78954"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PE <=> M_P-Ethanolamine + M_1,2-DG + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_P-Ethanolamine", "M_1,2-DG", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:64612 <=> CHEBI:58190 + CHEBI:17815 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:64612"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58190", "CHEBI:17815", "CHEBI:15378"))
    
    ## pe_to_sn1lpe
    pe <- "PE(14:0/16:0)"
    substrates = list(PE = pe)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:44604"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE + M_H2O = M_1-LPE + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPE", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 = CHEBI:64381 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64612", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64381", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:44605"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE + M_H2O => M_1-LPE + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPE", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 => CHEBI:64381 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64612", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64381", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:44606"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE + M_H2O <= M_1-LPE + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPE", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 <= CHEBI:64381 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64612", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64381", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:44607"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE + M_H2O <=> M_1-LPE + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPE", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 <=> CHEBI:64381 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64612", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64381", "CHEBI:28868", "CHEBI:15378"))
    
    ## pe_to_sn2lpe
    pe <- "PE(14:0/16:0)"
    substrates = list(PE = pe)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:44408"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$sn2LPE, "PE(0:0/14:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE + M_H2O = M_2-LPE + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_2-LPE", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 = CHEBI:65213 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64612", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:65213", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:44409"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$sn2LPE, "PE(0:0/14:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE + M_H2O => M_2-LPE + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_2-LPE", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 => CHEBI:65213 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64612", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:65213", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:44410"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$sn2LPE, "PE(0:0/14:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE + M_H2O <= M_2-LPE + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_2-LPE", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 <= CHEBI:65213 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64612", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:65213", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:44411"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$sn2LPE, "PE(0:0/14:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE + M_H2O <=> M_2-LPE + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_2-LPE", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 <=> CHEBI:65213 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64612", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:65213", "CHEBI:28868", "CHEBI:15378"))
    
    ## pe_to_nape_sn1
    pe <- "PE(14:0/16:0)"
    pc <- "PC(18:0/20:0)"
    substrates <- list(PE = pe, PC = pc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:45188"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$PC, "PC(18:0/20:0)")
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/20:0)")
    expect_equal(l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_PE = M_2-LPC + M_H+ + M_NAPE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_2-LPC", "M_H+", "M_NAPE"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 = CHEBI:57875 + CHEBI:15378 + CHEBI:62537")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:64612"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57875", "CHEBI:15378", "CHEBI:62537"))
    
    reaction <- "RHEA:45189"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$PC, "PC(18:0/20:0)")
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/20:0)")
    expect_equal(l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_PE => M_2-LPC + M_H+ + M_NAPE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_2-LPC", "M_H+", "M_NAPE"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 => CHEBI:57875 + CHEBI:15378 + CHEBI:62537")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:64612"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57875", "CHEBI:15378", "CHEBI:62537"))
    
    reaction <- "RHEA:45190"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$PC, "PC(18:0/20:0)")
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/20:0)")
    expect_equal(l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_PE <= M_2-LPC + M_H+ + M_NAPE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_2-LPC", "M_H+", "M_NAPE"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 <= CHEBI:57875 + CHEBI:15378 + CHEBI:62537")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:64612"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57875", "CHEBI:15378", "CHEBI:62537"))
    
    reaction <- "RHEA:45191"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$PC, "PC(18:0/20:0)")
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/20:0)")
    expect_equal(l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_PE <=> M_2-LPC + M_H+ + M_NAPE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_2-LPC", "M_H+", "M_NAPE"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 <=> CHEBI:57875 + CHEBI:15378 + CHEBI:62537")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:64612"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57875", "CHEBI:15378", "CHEBI:62537"))
    
    ## pe_to_nape_sn2
    pe <- "PE(14:0/16:0)"
    pc <- "PC(18:0/20:0)"
    substrates <- list(PE = pe, PC = pc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:45192"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$PC, "PC(18:0/20:0)")
    expect_equal(l[[1]]$sn1LPC, "PC(18:0/0:0)")
    expect_equal(l[[1]]$NAPE, "NAPE(14:0/16:0/20:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_PE = M_1-LPC + M_H+ + M_NAPE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPC", "M_H+", "M_NAPE"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 = CHEBI:58168 + CHEBI:15378 + CHEBI:62537")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:64612"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58168", "CHEBI:15378", "CHEBI:62537"))
    
    reaction <- "RHEA:45193"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$PC, "PC(18:0/20:0)")
    expect_equal(l[[1]]$sn1LPC, "PC(18:0/0:0)")
    expect_equal(l[[1]]$NAPE, "NAPE(14:0/16:0/20:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_PE => M_1-LPC + M_H+ + M_NAPE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPC", "M_H+", "M_NAPE"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 => CHEBI:58168 + CHEBI:15378 + CHEBI:62537")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:64612"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58168", "CHEBI:15378", "CHEBI:62537"))
    
    reaction <- "RHEA:45194"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$PC, "PC(18:0/20:0)")
    expect_equal(l[[1]]$sn1LPC, "PC(18:0/0:0)")
    expect_equal(l[[1]]$NAPE, "NAPE(14:0/16:0/20:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_PE <= M_1-LPC + M_H+ + M_NAPE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPC", "M_H+", "M_NAPE"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 <= CHEBI:58168 + CHEBI:15378 + CHEBI:62537")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:64612"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58168", "CHEBI:15378", "CHEBI:62537"))
    
    reaction <- "RHEA:45195"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$PC, "PC(18:0/20:0)")
    expect_equal(l[[1]]$sn1LPC, "PC(18:0/0:0)")
    expect_equal(l[[1]]$NAPE, "NAPE(14:0/16:0/20:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PC + M_PE <=> M_1-LPC + M_H+ + M_NAPE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PC", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPC", "M_H+", "M_NAPE"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 <=> CHEBI:58168 + CHEBI:15378 + CHEBI:62537")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57643", "CHEBI:64612"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58168", "CHEBI:15378", "CHEBI:62537"))
    
    ## pe_to_pa
    pe <- "PE(14:0/16:0)"
    substrates <- list(PE = pe)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "pe_to_pa"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE + M_H2O <=> M_PA + M_Ethanolamine + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_PA", "M_Ethanolamine", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 <=> CHEBI:58608 + CHEBI:57603 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64612", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58608", "CHEBI:57603", "CHEBI:15378"))
    
    ## pe_to_ps
    pe <- "PE(14:0/16:0)"
    substrates <- list(PE = pe)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:27606" 
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$PS, "PS(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE + M_L-Serine = M_PS + M_Ethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE", "M_L-Serine"))
    expect_equal(l[[2]]$reaction_product, c("M_PS", "M_Ethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:33384 = CHEBI:57262 + CHEBI:57603")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64612", "CHEBI:33384"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57262", "CHEBI:57603"))
    
    reaction <- "RHEA:27607"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$PS, "PS(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE + M_L-Serine => M_PS + M_Ethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE", "M_L-Serine"))
    expect_equal(l[[2]]$reaction_product, c("M_PS", "M_Ethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:33384 => CHEBI:57262 + CHEBI:57603")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64612", "CHEBI:33384"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57262", "CHEBI:57603"))
    
    reaction <- "RHEA:27608"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$PS, "PS(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE + M_L-Serine <= M_PS + M_Ethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE", "M_L-Serine"))
    expect_equal(l[[2]]$reaction_product, c("M_PS", "M_Ethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:33384 <= CHEBI:57262 + CHEBI:57603")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64612", "CHEBI:33384"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57262", "CHEBI:57603"))
    
    reaction <- "RHEA:27609"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$PS, "PS(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE + M_L-Serine <=> M_PS + M_Ethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE", "M_L-Serine"))
    expect_equal(l[[2]]$reaction_product, c("M_PS", "M_Ethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:33384 <=> CHEBI:57262 + CHEBI:57603")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64612", "CHEBI:33384"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57262", "CHEBI:57603"))
    
    ## peo_to_lpeo
    peo <- "PE(O-16:0/14:0)"
    substrates <- list(PEO = peo)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "peo_to_lpeo"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEO, "PE(O-16:0/14:0)")
    expect_equal(l[[1]]$sn1LPEO, "PE(O-16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE-O + M_H2O <=> M_LPE-O + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE-O", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_LPE-O", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:75028 + CHEBI:15377 <=>  + CHEBI:28868 + CHEBI:15378") ##################
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:75028", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("", "CHEBI:28868", "CHEBI:15378"))
    
    ## peo_to_napeo_sn1
    peo <- "PE(O-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    substrates <- list(PEO = peo, PC = pc) ################
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "peo_to_napeo_sn1"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEO, "PE(O-16:0/14:0)")
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/18:0)")
    expect_equal(l[[1]]$NAPEO, "NAPE(O-16:0/14:0/20:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE-O + M_PC <=> M_NAPEO + M_2-LPC")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE-O", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_NAPEO", "M_2-LPC"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:75028 + CHEBI:57643 <=>  + CHEBI:58168") ##################
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:75028", "CHEBI:57643"))
    expect_equal(l[[2]]$reaction_product_chebi, c("", "CHEBI:58168"))
    
    ## peo_to_napeo_sn2
    peo <- "PE(O-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    substrates <- list(PEO = peo, PC = pc) ################
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "peo_to_napeo_sn2"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEO, "PE(O-16:0/14:0)")
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn1LPC, "PC(20:0/0:0)")
    expect_equal(l[[1]]$NAPEO, "NAPE(O-16:0/14:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE-O + M_PC <=> M_NAPEO + M_1-LPC")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE-O", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_NAPEO", "M_1-LPC"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:75028 + CHEBI:57643 <=>  + CHEBI:58168") ##################
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:75028", "CHEBI:57643"))
    expect_equal(l[[2]]$reaction_product_chebi, c("", "CHEBI:58168"))
    
    ## peo_to_pep
    peo <- "PE(O-16:0/14:0)"
    substrates <- list(PEO = peo)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:22956"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEO, "PE(O-16:0/14:0)")
    expect_equal(l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE-O + M_Fe2+-cytochrome_b5 + 2 M_H+ + M_O2 = M_PE-P + M_Fe3+-cytochrome_b5 + 2 M_H2O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE-O", "M_Fe2+-cytochrome_b5", "2 M_H+", "M_O2"))
    expect_equal(l[[2]]$reaction_product, c("M_PE-P", "M_Fe3+-cytochrome_b5", "2 M_H2O"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:75028 + CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 = CHEBI:77290 + CHEBI:29034 + 2 CHEBI:15377")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:75028", "CHEBI:29033", "2 CHEBI:15378", "CHEBI:15379"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77290", "CHEBI:29034", "2 CHEBI:15377"))
    
    reaction <- "RHEA:22957"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEO, "PE(O-16:0/14:0)")
    expect_equal(l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE-O + M_Fe2+-cytochrome_b5 + 2 M_H+ + M_O2 => M_PE-P + M_Fe3+-cytochrome_b5 + 2 M_H2O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE-O", "M_Fe2+-cytochrome_b5", "2 M_H+", "M_O2"))
    expect_equal(l[[2]]$reaction_product, c("M_PE-P", "M_Fe3+-cytochrome_b5", "2 M_H2O"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:75028 + CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 => CHEBI:77290 + CHEBI:29034 + 2 CHEBI:15377")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:75028", "CHEBI:29033", "2 CHEBI:15378", "CHEBI:15379"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77290", "CHEBI:29034", "2 CHEBI:15377"))
    
    reaction <- "RHEA:22958"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEO, "PE(O-16:0/14:0)")
    expect_equal(l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE-O + M_Fe2+-cytochrome_b5 + 2 M_H+ + M_O2 <= M_PE-P + M_Fe3+-cytochrome_b5 + 2 M_H2O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE-O", "M_Fe2+-cytochrome_b5", "2 M_H+", "M_O2"))
    expect_equal(l[[2]]$reaction_product, c("M_PE-P", "M_Fe3+-cytochrome_b5", "2 M_H2O"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:75028 + CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 <= CHEBI:77290 + CHEBI:29034 + 2 CHEBI:15377")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:75028", "CHEBI:29033", "2 CHEBI:15378", "CHEBI:15379"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77290", "CHEBI:29034", "2 CHEBI:15377"))
    
    reaction <- "RHEA:22959"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEO, "PE(O-16:0/14:0)")
    expect_equal(l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE-O + M_Fe2+-cytochrome_b5 + 2 M_H+ + M_O2 <=> M_PE-P + M_Fe3+-cytochrome_b5 + 2 M_H2O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE-O", "M_Fe2+-cytochrome_b5", "2 M_H+", "M_O2"))
    expect_equal(l[[2]]$reaction_product, c("M_PE-P", "M_Fe3+-cytochrome_b5", "2 M_H2O"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:75028 + CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 <=> CHEBI:77290 + CHEBI:29034 + 2 CHEBI:15377")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:75028", "CHEBI:29033", "2 CHEBI:15378", "CHEBI:15379"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77290", "CHEBI:29034", "2 CHEBI:15377"))
    
    ## pep_to_lpep
    pep <- "PE(P-16:0/14:0)"
    substrates <- list(PEP = pep)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36195"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE-P + M_H2O = M_LPE-P + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE-P", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_LPE-P", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77290 + CHEBI:15377 = CHEBI:77288 + CHEBI:29067 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77290", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77288", "CHEBI:29067", "CHEBI:15378"))
    
    reaction <- "RHEA:36196"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE-P + M_H2O => M_LPE-P + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE-P", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_LPE-P", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77290 + CHEBI:15377 => CHEBI:77288 + CHEBI:29067 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77290", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77288", "CHEBI:29067", "CHEBI:15378"))
    
    reaction <- "RHEA:36197"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE-P + M_H2O <= M_LPE-P + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE-P", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_LPE-P", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77290 + CHEBI:15377 <= CHEBI:77288 + CHEBI:29067 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77290", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77288", "CHEBI:29067", "CHEBI:15378"))
    
    reaction <- "RHEA:36198"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE-P + M_H2O <=> M_LPE-P + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE-P", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_LPE-P", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77290 + CHEBI:15377 <=> CHEBI:77288 + CHEBI:29067 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77290", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77288", "CHEBI:29067", "CHEBI:15378"))
    
    ## lpep_to_fal
    sn1lpep <- "PE(P-16:0/0:0)"
    substrates = list(sn1LPEP = sn1lpep)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:16905"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$FAL, "FAL(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,"M_1-LPE-P + M_H2O = M_FAL + M_Glycerophosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPE-P", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_FAL", "M_Glycerophosphoethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 = CHEBI:73359 + CHEBI:143890")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77288", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:73359", "CHEBI:143890"))
    
    reaction <- "RHEA:16906"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$FAL, "FAL(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,"M_1-LPE-P + M_H2O => M_FAL + M_Glycerophosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPE-P", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_FAL", "M_Glycerophosphoethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 => CHEBI:73359 + CHEBI:143890")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77288", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:73359", "CHEBI:143890"))
    
    reaction <- "RHEA:16907"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$FAL, "FAL(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,"M_1-LPE-P + M_H2O <= M_FAL + M_Glycerophosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPE-P", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_FAL", "M_Glycerophosphoethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 <= CHEBI:73359 + CHEBI:143890")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77288", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:73359", "CHEBI:143890"))
    
    reaction <- "RHEA:16908"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$FAL, "FAL(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,"M_1-LPE-P + M_H2O <=> M_FAL + M_Glycerophosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPE-P", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_FAL", "M_Glycerophosphoethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 <=> CHEBI:73359 + CHEBI:143890")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77288", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:73359", "CHEBI:143890"))
    
    ## lpep_to_lpap
    sn1lpep <- "PE(P-16:0/0:0)"
    substrates <- list(sn1LPEP = sn1lpep)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36203"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$sn1LPAP, "PA(P-16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPE-P + M_H2O = M_LPA-P + M_Ethanolamine + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPE-P", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA-P", "M_Ethanolamine", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 = CHEBI:77283 + CHEBI:57603 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77288", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77283", "CHEBI:57603", "CHEBI:15378"))
    
    reaction <- "RHEA:36204"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$sn1LPAP, "PA(P-16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPE-P + M_H2O => M_LPA-P + M_Ethanolamine + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPE-P", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA-P", "M_Ethanolamine", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 => CHEBI:77283 + CHEBI:57603 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77288", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77283", "CHEBI:57603", "CHEBI:15378"))
    
    reaction <- "RHEA:36205"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$sn1LPAP, "PA(P-16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPE-P + M_H2O <= M_LPA-P + M_Ethanolamine + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPE-P", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA-P", "M_Ethanolamine", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 <= CHEBI:77283 + CHEBI:57603 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77288", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77283", "CHEBI:57603", "CHEBI:15378"))
    
    reaction <- "RHEA:36206"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$sn1LPAP, "PA(P-16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPE-P + M_H2O <=> M_LPA-P + M_Ethanolamine + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPE-P", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA-P", "M_Ethanolamine", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 <=> CHEBI:77283 + CHEBI:57603 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77288", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77283", "CHEBI:57603", "CHEBI:15378"))
    
    ## lpep_to_mgp
    sn1lpep <- "PE(P-16:0/0:0)"
    substrates <- list(sn1LPEP = sn1lpep)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36199"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$sn1MGP, "MG(P-16:0/0:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPE-P + M_H2O = M_1-MG-P + M_H+ + M_Phosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPE-P", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-MG-P", "M_H+", "M_Phosphoethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 = CHEBI:77297 + CHEBI:15378 + CHEBI:58190")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77288", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77297", "CHEBI:15378", "CHEBI:58190"))
    
    reaction <- "RHEA:36200"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$sn1MGP, "MG(P-16:0/0:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPE-P + M_H2O => M_1-MG-P + M_H+ + M_Phosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPE-P", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-MG-P", "M_H+", "M_Phosphoethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 => CHEBI:77297 + CHEBI:15378 + CHEBI:58190")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77288", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77297", "CHEBI:15378", "CHEBI:58190"))
    
    reaction <- "RHEA:36201"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$sn1MGP, "MG(P-16:0/0:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPE-P + M_H2O <= M_1-MG-P + M_H+ + M_Phosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPE-P", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-MG-P", "M_H+", "M_Phosphoethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 <= CHEBI:77297 + CHEBI:15378 + CHEBI:58190")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77288", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77297", "CHEBI:15378", "CHEBI:58190"))
    
    reaction <- "RHEA:36202"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$sn1MGP, "MG(P-16:0/0:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-LPE-P + M_H2O <=> M_1-MG-P + M_H+ + M_Phosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-LPE-P", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-MG-P", "M_H+", "M_Phosphoethanolamine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 <=> CHEBI:77297 + CHEBI:15378 + CHEBI:58190")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77288", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:77297", "CHEBI:15378", "CHEBI:58190"))
    
    ## pep_to_napep_sn1
    pep <- "PE(P-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    substrates <- list(PEP = pep, PC = pc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:63596"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/18:0)")
    expect_equal(l[[1]]$NAPEP, "NAPE(P-16:0/14:0/20:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_PE-P + M_PC = M_2-LPC + M_H+ + M_NAPEP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE-P", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_2-LPC", "M_H+", "M_NAPEP"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77290 + CHEBI:57643 = CHEBI:57875 + CHEBI:15378 + CHEBI:140451")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77290", "CHEBI:57643"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57875", "CHEBI:15378", "CHEBI:140451"))
    
    reaction <- "RHEA:63597"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/18:0)")
    expect_equal(l[[1]]$NAPEP, "NAPE(P-16:0/14:0/20:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_PE-P + M_PC => M_2-LPC + M_H+ + M_NAPEP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE-P", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_2-LPC", "M_H+", "M_NAPEP"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77290 + CHEBI:57643 => CHEBI:57875 + CHEBI:15378 + CHEBI:140451")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77290", "CHEBI:57643"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57875", "CHEBI:15378", "CHEBI:140451"))
    
    reaction <- "RHEA:63598"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/18:0)")
    expect_equal(l[[1]]$NAPEP, "NAPE(P-16:0/14:0/20:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_PE-P + M_PC <= M_2-LPC + M_H+ + M_NAPEP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE-P", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_2-LPC", "M_H+", "M_NAPEP"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77290 + CHEBI:57643 <= CHEBI:57875 + CHEBI:15378 + CHEBI:140451")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77290", "CHEBI:57643"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57875", "CHEBI:15378", "CHEBI:140451"))
    
    reaction <- "RHEA:63599"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(), reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/18:0)")
    expect_equal(l[[1]]$NAPEP, "NAPE(P-16:0/14:0/20:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_PE-P + M_PC <=> M_2-LPC + M_H+ + M_NAPEP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE-P", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_2-LPC", "M_H+", "M_NAPEP"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77290 + CHEBI:57643 <=> CHEBI:57875 + CHEBI:15378 + CHEBI:140451")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77290", "CHEBI:57643"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:57875", "CHEBI:15378", "CHEBI:140451"))
    
    ## pep_to_napep_sn2
    pep <- "PE(P-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    substrates <- list(PEP = pep, PC = pc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "pep_to_napep_sn2"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn1LPC, "PC(20:0/0:0)")
    expect_equal(l[[1]]$NAPEP, "NAPE(P-16:0/14:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE-P + M_PC <=> M_NAPEP + M_1-LPC")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE-P", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_NAPEP", "M_1-LPC"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77290 + CHEBI:57643 <=> CHEBI:140451 + CHEBI:58168")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77290", "CHEBI:57643"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:140451", "CHEBI:58168"))
    
    ## pg_to_cl
    pg <- "PG(18:4(6Z,9Z,12Z,15Z)/14:0)"
    cdpdg <- "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)"
    substrates = list(PG = pg, CDPDG = cdpdg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:32931"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PG, "PG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(l[[1]]$CDPDG, "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(l[[1]]$PGs2, "18:4(6Z,9Z,12Z,15Z)/14:0") ##
    expect_equal(l[[1]]$CDPDGs2, "18:4(6Z,9Z,12Z,15Z)/14:0") ##
    expect_equal(l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PG + M_CDP-DG = M_CL + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PG", "M_CDP-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_CL", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64716 + CHEBI:58332 = CHEBI:62237 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64716", "CHEBI:58332"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:62237", "CHEBI:60377", "CHEBI:15378"))
    
    reaction <- "RHEA:32932"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PG, "PG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(l[[1]]$CDPDG, "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(l[[1]]$PGs2, "18:4(6Z,9Z,12Z,15Z)/14:0") ##
    expect_equal(l[[1]]$CDPDGs2, "18:4(6Z,9Z,12Z,15Z)/14:0") ##
    expect_equal(l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PG + M_CDP-DG => M_CL + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PG", "M_CDP-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_CL", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64716 + CHEBI:58332 => CHEBI:62237 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64716", "CHEBI:58332"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:62237", "CHEBI:60377", "CHEBI:15378"))
    
    reaction <- "RHEA:32933"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PG, "PG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(l[[1]]$CDPDG, "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(l[[1]]$PGs2, "18:4(6Z,9Z,12Z,15Z)/14:0") ##
    expect_equal(l[[1]]$CDPDGs2, "18:4(6Z,9Z,12Z,15Z)/14:0") ##
    expect_equal(l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PG + M_CDP-DG <= M_CL + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PG", "M_CDP-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_CL", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64716 + CHEBI:58332 <= CHEBI:62237 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64716", "CHEBI:58332"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:62237", "CHEBI:60377", "CHEBI:15378"))
    
    reaction <- "RHEA:32934"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PG, "PG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(l[[1]]$CDPDG, "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(l[[1]]$PGs2, "18:4(6Z,9Z,12Z,15Z)/14:0") ##
    expect_equal(l[[1]]$CDPDGs2, "18:4(6Z,9Z,12Z,15Z)/14:0") ##
    expect_equal(l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PG + M_CDP-DG <=> M_CL + M_CMP + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PG", "M_CDP-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_CL", "M_CMP", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64716 + CHEBI:58332 <=> CHEBI:62237 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64716", "CHEBI:58332"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:62237", "CHEBI:60377", "CHEBI:15378"))
    
    ## pgp_to_pg
    pgp <- "PGP(16:0/14:0)"
    substrates <- list(PGP = pgp)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:33751"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PGP, "PGP(16:0/14:0)")
    expect_equal(l[[1]]$PG, "PG(16:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PGP + M_H2O = M_PG + M_Pi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PGP", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_PG", "M_Pi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:60110 + CHEBI:15377 = CHEBI:64716 + CHEBI:43474")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:60110", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64716", "CHEBI:43474"))
    
    reaction <- "RHEA:33752"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PGP, "PGP(16:0/14:0)")
    expect_equal(l[[1]]$PG, "PG(16:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PGP + M_H2O => M_PG + M_Pi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PGP", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_PG", "M_Pi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:60110 + CHEBI:15377 => CHEBI:64716 + CHEBI:43474")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:60110", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64716", "CHEBI:43474"))
    
    reaction <- "RHEA:33753"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PGP, "PGP(16:0/14:0)")
    expect_equal(l[[1]]$PG, "PG(16:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PGP + M_H2O <= M_PG + M_Pi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PGP", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_PG", "M_Pi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:60110 + CHEBI:15377 <= CHEBI:64716 + CHEBI:43474")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:60110", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64716", "CHEBI:43474"))
    
    reaction <- "RHEA:33754"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PGP, "PGP(16:0/14:0)")
    expect_equal(l[[1]]$PG, "PG(16:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PGP + M_H2O <=> M_PG + M_Pi")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PGP", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_PG", "M_Pi"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:60110 + CHEBI:15377 <=> CHEBI:64716 + CHEBI:43474")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:60110", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64716", "CHEBI:43474"))
    
    ## pi_to_dg
    pi <- "PI(16:0/18:1(9Z))"
    substrates <- list(PI = pi)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:43484"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(l[[1]]$DG, "DG(16:0/18:1(9Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PI + M_H2O = M_myo-Inositol-1-P + M_1,2-DG + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PI", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_myo-Inositol-1-P", "M_1,2-DG", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 = CHEBI:58433 + CHEBI:17815 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57880", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58433", "CHEBI:17815", "CHEBI:15378"))
    
    reaction <- "RHEA:43485"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(l[[1]]$DG, "DG(16:0/18:1(9Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PI + M_H2O => M_myo-Inositol-1-P + M_1,2-DG + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PI", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_myo-Inositol-1-P", "M_1,2-DG", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 => CHEBI:58433 + CHEBI:17815 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57880", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58433", "CHEBI:17815", "CHEBI:15378"))
    
    reaction <- "RHEA:43486"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(l[[1]]$DG, "DG(16:0/18:1(9Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PI + M_H2O <= M_myo-Inositol-1-P + M_1,2-DG + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PI", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_myo-Inositol-1-P", "M_1,2-DG", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 <= CHEBI:58433 + CHEBI:17815 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57880", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58433", "CHEBI:17815", "CHEBI:15378"))
    
    reaction <- "RHEA:43487"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(l[[1]]$DG, "DG(16:0/18:1(9Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PI + M_H2O <=> M_myo-Inositol-1-P + M_1,2-DG + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PI", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_myo-Inositol-1-P", "M_1,2-DG", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 <=> CHEBI:58433 + CHEBI:17815 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57880", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:58433", "CHEBI:17815", "CHEBI:15378"))
    
    ## pi_to_sn1lpi
    pi <- "PI(16:0/18:1(9Z))"
    substrates <- list(PI = pi)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:18001"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:1(9Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PI + M_H2O = M_1-LPI + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PI", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPI", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 = CHEBI:64771 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57880", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64771", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:18002"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:1(9Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PI + M_H2O => M_1-LPI + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PI", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPI", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 => CHEBI:64771 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57880", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64771", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:18003"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:1(9Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PI + M_H2O <= M_1-LPI + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PI", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPI", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 <= CHEBI:64771 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57880", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64771", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:18004"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:1(9Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PI + M_H2O <=> M_1-LPI + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PI", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPI", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 <=> CHEBI:64771 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57880", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64771", "CHEBI:28868", "CHEBI:15378"))
    
    ## ps_to_pe
    ps <- "PS(14:0/14:0)"
    substrates = list(PS = ps)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:20828"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PS, "PS(14:0/14:0)")
    expect_equal(l[[1]]$PE, "PE(14:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PS + M_H+ = M_PE + M_CO2")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PS", "M_H+"))
    expect_equal(l[[2]]$reaction_product, c("M_PE", "M_CO2"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57262 + CHEBI:15378 = CHEBI:64612 + CHEBI:16526")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57262", "CHEBI:15378"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64612", "CHEBI:16526"))
    
    reaction <- "RHEA:20829"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PS, "PS(14:0/14:0)")
    expect_equal(l[[1]]$PE, "PE(14:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PS + M_H+ => M_PE + M_CO2")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PS", "M_H+"))
    expect_equal(l[[2]]$reaction_product, c("M_PE", "M_CO2"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57262 + CHEBI:15378 => CHEBI:64612 + CHEBI:16526")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57262", "CHEBI:15378"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64612", "CHEBI:16526"))
    
    reaction <- "RHEA:20830"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PS, "PS(14:0/14:0)")
    expect_equal(l[[1]]$PE, "PE(14:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PS + M_H+ <= M_PE + M_CO2")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PS", "M_H+"))
    expect_equal(l[[2]]$reaction_product, c("M_PE", "M_CO2"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57262 + CHEBI:15378 <= CHEBI:64612 + CHEBI:16526")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57262", "CHEBI:15378"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64612", "CHEBI:16526"))
    
    reaction <- "RHEA:20831"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PS, "PS(14:0/14:0)")
    expect_equal(l[[1]]$PE, "PE(14:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PS + M_H+ <=> M_PE + M_CO2")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PS", "M_H+"))
    expect_equal(l[[2]]$reaction_product, c("M_PE", "M_CO2"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:57262 + CHEBI:15378 <=> CHEBI:64612 + CHEBI:16526")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:57262", "CHEBI:15378"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:64612", "CHEBI:16526"))
    
    ## sm_to_cer
    sm <- "SM(16:0(3OH,4OH,15Me)/12:0)"
    substrates <- list(SM = sm)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:45644"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$SM, "SM(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_SM = M_H+ + M_Cer + M_Phosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_SM"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_Cer", "M_Phosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:78646 = CHEBI:15378 + CHEBI:72959 + CHEBI:295975")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:78646"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:15378", "CHEBI:72959", "CHEBI:295975"))
    
    reaction <- "RHEA:45645"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$SM, "SM(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_SM => M_H+ + M_Cer + M_Phosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_SM"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_Cer", "M_Phosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:78646 => CHEBI:15378 + CHEBI:72959 + CHEBI:295975")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:78646"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:15378", "CHEBI:72959", "CHEBI:295975"))
    
    reaction <- "RHEA:45646"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$SM, "SM(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_SM <= M_H+ + M_Cer + M_Phosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_SM"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_Cer", "M_Phosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:78646 <= CHEBI:15378 + CHEBI:72959 + CHEBI:295975")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:78646"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:15378", "CHEBI:72959", "CHEBI:295975"))
    
    reaction <- "RHEA:45647"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$SM, "SM(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_SM <=> M_H+ + M_Cer + M_Phosphocholine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_SM"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_Cer", "M_Phosphocholine"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:78646 <=> CHEBI:15378 + CHEBI:72959 + CHEBI:295975")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:15377", "CHEBI:78646"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:15378", "CHEBI:72959", "CHEBI:295975"))
    
    ## sphinga_to_dhcer
    acylcoa <- "CoA(12:0)"
    substrates <- list(AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:53424"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(12:0)")
    expect_equal(l[[1]]$DhCer, "Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_Sphinganine = M_DhCer + M_CoA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_Sphinganine"))
    expect_equal(l[[2]]$reaction_product, c("M_DhCer", "M_CoA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77636 + CHEBI:84410 = CHEBI:83273 + CHEBI:57287 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77636", "CHEBI:84410"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:83273", "CHEBI:57287", "CHEBI:15378"))
    
    reaction <- "RHEA:53425"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(12:0)")
    expect_equal(l[[1]]$DhCer, "Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_Sphinganine => M_DhCer + M_CoA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_Sphinganine"))
    expect_equal(l[[2]]$reaction_product, c("M_DhCer", "M_CoA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77636 + CHEBI:84410 => CHEBI:83273 + CHEBI:57287 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77636", "CHEBI:84410"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:83273", "CHEBI:57287", "CHEBI:15378"))
    
    reaction <- "RHEA:53426"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(12:0)")
    expect_equal(l[[1]]$DhCer, "Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_Sphinganine <= M_DhCer + M_CoA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_Sphinganine"))
    expect_equal(l[[2]]$reaction_product, c("M_DhCer", "M_CoA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77636 + CHEBI:84410 <= CHEBI:83273 + CHEBI:57287 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77636", "CHEBI:84410"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:83273", "CHEBI:57287", "CHEBI:15378"))
    
    reaction <- "RHEA:53427"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(12:0)")
    expect_equal(l[[1]]$DhCer, "Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_Sphinganine <=> M_DhCer + M_CoA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_Sphinganine"))
    expect_equal(l[[2]]$reaction_product, c("M_DhCer", "M_CoA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:77636 + CHEBI:84410 <=> CHEBI:83273 + CHEBI:57287 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:77636", "CHEBI:84410"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:83273", "CHEBI:57287", "CHEBI:15378"))
    
    ## tg_to_dg
    tg <- "TG(18:0/16:0/14:0)"
    substrates = list(TG = tg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:33271"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_TG + M_H2O = M_1,2-DG + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_TG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64615 + CHEBI:15377 = CHEBI:17815 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64615", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:33272"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction) 
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_TG + M_H2O => M_1,2-DG + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_TG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64615 + CHEBI:15377 => CHEBI:17815 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64615", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:33273"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_TG + M_H2O <= M_1,2-DG + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_TG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64615 + CHEBI:15377 <= CHEBI:17815 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64615", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:33274"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_TG + M_H2O <=> M_1,2-DG + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_TG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64615 + CHEBI:15377 <=> CHEBI:17815 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64615", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:17815", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:44864"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_TG + M_H2O = M_1,2-DG + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_TG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64615 + CHEBI:15377 = CHEBI:49172 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64615", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:49172", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:44865"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_TG + M_H2O => M_1,2-DG + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_TG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64615 + CHEBI:15377 => CHEBI:49172 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64615", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:49172", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:44866"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_TG + M_H2O <= M_1,2-DG + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_TG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64615 + CHEBI:15377 <= CHEBI:49172 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64615", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:49172", "CHEBI:28868", "CHEBI:15378"))
    
    reaction <- "RHEA:44867"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- .create_template(template = list(),
        reaction = reaction)
    df <- .create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_TG + M_H2O <=> M_1,2-DG + M_FA + M_H+")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_TG", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_1,2-DG", "M_FA", "M_H+"))
    expect_equal(l[[2]]$reaction_formula_chebi, "CHEBI:64615 + CHEBI:15377 <=> CHEBI:49172 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l[[2]]$reaction_substrate_chebi, c("CHEBI:64615", "CHEBI:15377"))
    expect_equal(l[[2]]$reaction_product_chebi, c("CHEBI:49172", "CHEBI:28868", "CHEBI:15378"))
})

