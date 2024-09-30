## function .create_list_reactions_with_template
test_that(".create_list_reactions_with_template works", {
    
    ## acyldhap_to_alkyldhap
    acyldhap <- "DHAP(18:0)"
    fao <- "FAO(16:0)"
    substrates <- list(AcylDHAP = acyldhap, FAO = fao)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36171"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    
    expect_equal(l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(l[[1]]$FAO, "FAO(16:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AlkylDHAP, "DHAP(O-16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylDHAP + M_FAO = M_H+ + M_FA + M_AlkylDHAP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylDHAP", "M_FAO"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_AlkylDHAP"))
    
    reaction <- "RHEA:36172"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(l[[1]]$FAO, "FAO(16:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AlkylDHAP, "DHAP(O-16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylDHAP + M_FAO => M_H+ + M_FA + M_AlkylDHAP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylDHAP", "M_FAO"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_AlkylDHAP"))
    
    reaction <- "RHEA:36173"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(l[[1]]$FAO, "FAO(16:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AlkylDHAP, "DHAP(O-16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylDHAP + M_FAO <= M_H+ + M_FA + M_AlkylDHAP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylDHAP", "M_FAO"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_AlkylDHAP"))
    
    reaction <- "RHEA:36174"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(l[[1]]$FAO, "FAO(16:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AlkylDHAP, "DHAP(O-16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylDHAP + M_FAO <=> M_H+ + M_FA + M_AlkylDHAP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylDHAP", "M_FAO"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_AlkylDHAP"))
    
    ## alkyldhap_to_lpao
    alkyldhap <- "DHAP(O-18:0)"
    substrates = list(AlkylDHAP = alkyldhap)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36175"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AlkylDHAP, "DHAP(O-18:0)")
    expect_equal(l[[1]]$sn1LPAO, "PA(O-18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H+ + M_NADPH + M_AlkylDHAP = M_LPA-O + M_NADP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H+", "M_NADPH", "M_AlkylDHAP"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA-O", "M_NADP"))
    
    reaction <- "RHEA:36176"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AlkylDHAP, "DHAP(O-18:0)")
    expect_equal(l[[1]]$sn1LPAO, "PA(O-18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H+ + M_NADPH + M_AlkylDHAP => M_LPA-O + M_NADP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H+", "M_NADPH", "M_AlkylDHAP"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA-O", "M_NADP"))
    
    reaction <- "RHEA:36177"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AlkylDHAP, "DHAP(O-18:0)")
    expect_equal(l[[1]]$sn1LPAO, "PA(O-18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H+ + M_NADPH + M_AlkylDHAP <= M_LPA-O + M_NADP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H+", "M_NADPH", "M_AlkylDHAP"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA-O", "M_NADP"))
    
    reaction <- "RHEA:36178"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AlkylDHAP, "DHAP(O-18:0)")
    expect_equal(l[[1]]$sn1LPAO, "PA(O-18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H+ + M_NADPH + M_AlkylDHAP <=> M_LPA-O + M_NADP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H+", "M_NADPH", "M_AlkylDHAP"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA-O", "M_NADP"))
    
    ## cerp_to_cer
    cerp <- "CerP(16:0(3OH,4OH,15Me)/18:0)"
    substrates <- list(CerP = cerp)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    reaction <- "cerp_to_cer"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    ## cdpdg_to_pgp
    cdpdg <- "CDP-DG(12:0(11Me)/14:0)"
    substrates <- list(CDPDG = cdpdg)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:12593"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(l[[1]]$PGP, "PGP(12:0(11Me)/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_Glycerol-3-P + M_CDP-DG = M_H+ + M_CMP + M_PGP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_Glycerol-3-P", "M_CDP-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PGP"))
    
    reaction <- "RHEA:12594"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(l[[1]]$PGP, "PGP(12:0(11Me)/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_Glycerol-3-P + M_CDP-DG => M_H+ + M_CMP + M_PGP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_Glycerol-3-P", "M_CDP-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PGP"))
    
    reaction <- "RHEA:12595"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(l[[1]]$PGP, "PGP(12:0(11Me)/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_Glycerol-3-P + M_CDP-DG <= M_H+ + M_CMP + M_PGP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_Glycerol-3-P", "M_CDP-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PGP"))
    
    reaction <- "RHEA:12596"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(l[[1]]$PGP, "PGP(12:0(11Me)/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_Glycerol-3-P + M_CDP-DG <=> M_H+ + M_CMP + M_PGP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_Glycerol-3-P", "M_CDP-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PGP"))
    
    ## cdpdg_to_pi
    cdgdg <- "CDP-DG(12:0(11Me)/14:0)"
    substrates <- list(CDPDG = cdpdg)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:11580"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(l[[1]]$PI, "PI(12:0(11Me)/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_myo-Inositol + M_CDP-DG = M_H+ + M_CMP + M_PI")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_myo-Inositol", "M_CDP-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PI"))
    
    reaction <- "RHEA:11581"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(l[[1]]$PI, "PI(12:0(11Me)/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_myo-Inositol + M_CDP-DG => M_H+ + M_CMP + M_PI")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_myo-Inositol", "M_CDP-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PI"))
    
    reaction <- "RHEA:11582"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(l[[1]]$PI, "PI(12:0(11Me)/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_myo-Inositol + M_CDP-DG <= M_H+ + M_CMP + M_PI")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_myo-Inositol", "M_CDP-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PI"))
    
    reaction <- "RHEA:11583"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(l[[1]]$PI, "PI(12:0(11Me)/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_myo-Inositol + M_CDP-DG <=> M_H+ + M_CMP + M_PI")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_myo-Inositol", "M_CDP-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PI"))
    
    ## cer_to_cerp
    cer <- "Cer(16:0(3OH,4OH,15Me)/18:0)"
    substrates <- list(Cer = cer)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    reaction <- "RHEA:17929"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[1]]$CerP, "CerP(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_ATP + M_Cer <=> M_H+ + M_ADP + M_CerP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_ATP", "M_Cer"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_ADP", "M_CerP"))
    
    ## cer_to_glccer
    cer <- "Cer(16:0(3OH,4OH,15Me)/18:0)"
    substrates <- list(Cer = cer)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    reaction <- "RHEA:12088"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[1]]$GlcCer, "GlcCer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_UDP-Glucose + M_Cer <=> M_H+ + M_UDP + M_GlcCer")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_UDP-Glucose", "M_Cer"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_UDP", "M_GlcCer"))
    
    ## cer_to_sm
    cer <- "Cer(16:0(3OH,4OH,15Me)/18:0)"
    substrates <- list(Cer = cer)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    reaction <-  "RHEA:18765"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    ## cl_to_lcl
    cl <- "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])"
    substrates <- list(CL = cl)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:32935"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[1]]$LCL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[1]]$FA, "FA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_CL = M_1,2,4-LCL + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_CL" ))
    expect_equal(l[[2]]$reaction_product, c("M_1,2,4-LCL", "M_H+", "M_FA"))
    
    reaction <- "RHEA:32936"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[1]]$LCL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[1]]$FA, "FA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_CL => M_1,2,4-LCL + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_CL" ))
    expect_equal(l[[2]]$reaction_product, c("M_1,2,4-LCL", "M_H+", "M_FA"))
    
    reaction <- "RHEA:32937"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[1]]$LCL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[1]]$FA, "FA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_CL <= M_1,2,4-LCL + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_CL" ))
    expect_equal(l[[2]]$reaction_product, c("M_1,2,4-LCL", "M_H+", "M_FA"))
    
    reaction <- "RHEA:32938"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[1]]$LCL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(l[[1]]$FA, "FA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_CL <=> M_1,2,4-LCL + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_CL" ))
    expect_equal(l[[2]]$reaction_product, c("M_1,2,4-LCL", "M_H+", "M_FA"))
    
    ## coa_to_acyldhap
    acylcoa <- "CoA(18:0)"
    substrates <- list(AcylCoA = acylcoa)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:17657"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_Dihydroxyacetone-P + M_AcylCoA = M_CoA + M_AcylDHAP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_Dihydroxyacetone-P", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_AcylDHAP"))
    
    reaction <- "RHEA:17658"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_Dihydroxyacetone-P + M_AcylCoA => M_CoA + M_AcylDHAP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_Dihydroxyacetone-P", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_AcylDHAP"))
    
    reaction <- "RHEA:17659"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_Dihydroxyacetone-P + M_AcylCoA <= M_CoA + M_AcylDHAP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_Dihydroxyacetone-P", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_AcylDHAP"))
    
    reaction <- "RHEA:17660"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_Dihydroxyacetone-P + M_AcylCoA <=> M_CoA + M_AcylDHAP")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_Dihydroxyacetone-P", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_AcylDHAP"))
    
    ## coa_to_FAO
    acylcoa <- "CoA(18:0)"
    substrates <- list(AcylCoA = acylcoa)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:52716"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$FAO, "FAO(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + 2 M_NADPH + 2 M_H+ = M_FAO + 2 M_NADP + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "2 M_NADPH", "2 M_H+"))
    expect_equal(l[[2]]$reaction_product, c("M_FAO", "2 M_NADP", "M_CoA" ))
    
    reaction <- "RHEA:52717"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$FAO, "FAO(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + 2 M_NADPH + 2 M_H+ => M_FAO + 2 M_NADP + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "2 M_NADPH", "2 M_H+"))
    expect_equal(l[[2]]$reaction_product, c("M_FAO", "2 M_NADP", "M_CoA" ))
    
    reaction <- "RHEA:52718"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$FAO, "FAO(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + 2 M_NADPH + 2 M_H+ <= M_FAO + 2 M_NADP + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "2 M_NADPH", "2 M_H+"))
    expect_equal(l[[2]]$reaction_product, c("M_FAO", "2 M_NADP", "M_CoA" )) 
    
    reaction <- "RHEA:52719"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$FAO, "FAO(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + 2 M_NADPH + 2 M_H+ <=> M_FAO + 2 M_NADP + M_CoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "2 M_NADPH", "2 M_H+"))
    expect_equal(l[[2]]$reaction_product, c("M_FAO", "2 M_NADP", "M_CoA" ))
    
    ## coa_to_lpa
    acylcoa <- "CoA(18:0)"
    substrates <- list(AcylCoA = acylcoa)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:15325"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_Glycerol-3-P + M_AcylCoA = M_CoA + M_LPA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_Glycerol-3-P", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_LPA"))
    
    reaction <- "RHEA:15326"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_Glycerol-3-P + M_AcylCoA => M_CoA + M_LPA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_Glycerol-3-P", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_LPA"))
    
    reaction <- "RHEA:15327"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_Glycerol-3-P + M_AcylCoA <= M_CoA + M_LPA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_Glycerol-3-P", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_LPA"))
    
    reaction <- "RHEA:15328"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_Glycerol-3-P + M_AcylCoA <=> M_CoA + M_LPA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_Glycerol-3-P", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_LPA"))
    
    ## dg_to_sn1mg
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:44712"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1,2-DG = M_H+ + M_1-MG + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_1-MG", "M_FA"))
    
    reaction <- "RHEA:44713"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1,2-DG => M_H+ + M_1-MG + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_1-MG", "M_FA")) 
    
    reaction <- "RHEA:44714"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1,2-DG <= M_H+ + M_1-MG + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_1-MG", "M_FA"))
    
    reaction <- "RHEA:44715"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1,2-DG <=> M_H+ + M_1-MG + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_1-MG", "M_FA"))
    
    reaction <- "RHEA:35663"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1,2-DG = M_H+ + M_1-MG + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_1-MG", "M_FA"))
    
    reaction <- "RHEA:35664"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1,2-DG => M_H+ + M_1-MG + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_1-MG", "M_FA"))
    
    reaction <- "RHEA:35665"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1,2-DG <= M_H+ + M_1-MG + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_1-MG", "M_FA"))
    
    reaction <- "RHEA:35666"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1,2-DG <=> M_H+ + M_1-MG + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_1-MG", "M_FA"))
    
    ## dg_to_sn2mg
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:33275"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn2MG, "MG(0:0/16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1,2-DG = M_H+ + M_2-MG + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_2-MG", "M_FA"))
    
    reaction <- "RHEA:33276"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn2MG, "MG(0:0/16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1,2-DG => M_H+ + M_2-MG + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_2-MG", "M_FA"))
    
    reaction <- "RHEA:33277"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn2MG, "MG(0:0/16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1,2-DG <= M_H+ + M_2-MG + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_2-MG", "M_FA"))
    
    reaction <- "RHEA:33278"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$sn2MG, "MG(0:0/16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1,2-DG <=> M_H+ + M_2-MG + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_2-MG", "M_FA"))
    
    ## dg_to_pa
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:10272"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PA, "PA(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_ATP + M_1,2-DG = M_H+ + M_ADP + M_PA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_ATP", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_ADP", "M_PA"))
    
    reaction <- "RHEA:10273"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PA, "PA(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_ATP + M_1,2-DG => M_H+ + M_ADP + M_PA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_ATP", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_ADP", "M_PA"))
    
    reaction <- "RHEA:10274"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PA, "PA(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_ATP + M_1,2-DG <= M_H+ + M_ADP + M_PA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_ATP", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_ADP", "M_PA"))
    
    reaction <- "RHEA:10275"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PA, "PA(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_ATP + M_1,2-DG <=> M_H+ + M_ADP + M_PA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_ATP", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_ADP", "M_PA"))
    
    ## dg_to_pc
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:32939"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PC, "PC(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-Choline + M_1,2-DG = M_H+ + M_CMP + M_PC")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-Choline", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PC"))
    
    reaction <- "RHEA:32940"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PC, "PC(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-Choline + M_1,2-DG => M_H+ + M_CMP + M_PC")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-Choline", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PC"))
    
    reaction <- "RHEA:32941"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PC, "PC(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-Choline + M_1,2-DG <= M_H+ + M_CMP + M_PC")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-Choline", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PC"))
    
    reaction <- "RHEA:32942"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PC, "PC(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-Choline + M_1,2-DG <=> M_H+ + M_CMP + M_PC")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-Choline", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PC"))
    
    ## dg_to_pe
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:32943"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PE, "PE(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-Ethanolamine + M_1,2-DG = M_H+ + M_CMP + M_PE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-Ethanolamine", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PE"))
    
    reaction <- "RHEA:32944"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PE, "PE(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-Ethanolamine + M_1,2-DG => M_H+ + M_CMP + M_PE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-Ethanolamine", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PE"))
    
    reaction <- "RHEA:32945"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PE, "PE(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-Ethanolamine + M_1,2-DG <= M_H+ + M_CMP + M_PE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-Ethanolamine", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PE"))
    
    reaction <- "RHEA:32946"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$PE, "PE(18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-Ethanolamine + M_1,2-DG <=> M_H+ + M_CMP + M_PE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-Ethanolamine", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PE"))
    
    ## dg_to_tg
    dg <- "DG(18:0/16:0/0:0)"
    acylcoa <- "CoA(14:0)"
    substrates <- list(DG = dg, AcylCoA = acylcoa)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:10868"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_1,2-DG = M_CoA + M_TG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_TG"))
    
    reaction <- "RHEA:10869"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_1,2-DG => M_CoA + M_TG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_TG"))
    
    reaction <- "RHEA:10870"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_1,2-DG <= M_CoA + M_TG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_TG"))
    
    reaction <- "RHEA:10871"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_1,2-DG <=> M_CoA + M_TG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_1,2-DG"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_TG"))
    
    ## dgo_to_pco
    dgo <- "DG(O-18:0/16:0/0:0)"
    substrates <- list(DGO = dgo)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36179"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(l[[1]]$PCO, "PC(O-18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-Choline + M_DG-O = M_H+ + M_CMP + M_PC-O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-Choline", "M_DG-O"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PC-O"))
    
    reaction <- "RHEA:36180"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(l[[1]]$PCO, "PC(O-18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-Choline + M_DG-O => M_H+ + M_CMP + M_PC-O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-Choline", "M_DG-O"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PC-O"))
    
    reaction <- "RHEA:36181"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(l[[1]]$PCO, "PC(O-18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-Choline + M_DG-O <= M_H+ + M_CMP + M_PC-O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-Choline", "M_DG-O"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PC-O"))
    
    reaction <- "RHEA:36182"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(l[[1]]$PCO, "PC(O-18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-Choline + M_DG-O <=> M_H+ + M_CMP + M_PC-O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-Choline", "M_DG-O"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PC-O"))
    
    ## dgo_to_peo
    dgo <- "DG(O-18:0/16:0/0:0)"
    substrates <- list(DGO = dgo)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36187"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
       reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(l[[1]]$PEO, "PE(O-18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-Ethanolamine + M_DG-O = M_H+ + M_CMP + M_PE-O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-Ethanolamine", "M_DG-O"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PE-O"))
    
    reaction <- "RHEA:36188"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(l[[1]]$PEO, "PE(O-18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-Ethanolamine + M_DG-O => M_H+ + M_CMP + M_PE-O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-Ethanolamine", "M_DG-O"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PE-O"))
    
    reaction <- "RHEA:36189"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(l[[1]]$PEO, "PE(O-18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-Ethanolamine + M_DG-O <= M_H+ + M_CMP + M_PE-O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-Ethanolamine", "M_DG-O"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PE-O"))
    
    reaction <- "RHEA:36190"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(l[[1]]$PEO, "PE(O-18:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CDP-Ethanolamine + M_DG-O <=> M_H+ + M_CMP + M_PE-O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-Ethanolamine", "M_DG-O"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_PE-O"))
    
    ## dhcer_to_cer
    dhcer <- "Cer(d16:0(3OH,4OH)(15Me)/12:0)" ###################################################
    substrates <- list(DhCer = dhcer)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "dhcer_to_cer"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DhCer, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(l[[1]]$Cer, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H+ + M_NADH + M_O2 + M_DhCer <=> 2 M_H2O + M_NAD + M_Cer")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H+", "M_NADH", "M_O2", "M_DhCer"))
    expect_equal(l[[2]]$reaction_product, c("2 M_H2O", "M_NAD", "M_Cer"))
    
    ## dhcer_to_dhsm
    dhcer <- "Cer(16:1(3OH,4OH,15Me)/12:0)"
    substrates <- list(DhCer = dhcer)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:44620"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    ## dhsm_to_dhcer
    dhsm <- "SM(16:1(3OH,4OH,15Me)/12:0)" ###################################################
    substrates <- list(DhSM = dhsm)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:45300"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DhSM, "SM(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[1]]$DhCer, "Cer(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_DhSM = M_Phosphocholine + M_H+ + M_DhCer")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_DhSM"))
    expect_equal(l[[2]]$reaction_product, c("M_Phosphocholine", "M_H+", "M_DhCer"))
    
    reaction <- "RHEA:45301"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DhSM, "SM(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[1]]$DhCer, "Cer(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_DhSM => M_Phosphocholine + M_H+ + M_DhCer")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_DhSM"))
    expect_equal(l[[2]]$reaction_product, c("M_Phosphocholine", "M_H+", "M_DhCer"))
    
    reaction <- "RHEA:45302"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DhSM, "SM(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[1]]$DhCer, "Cer(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_DhSM <= M_Phosphocholine + M_H+ + M_DhCer")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_DhSM"))
    expect_equal(l[[2]]$reaction_product, c("M_Phosphocholine", "M_H+", "M_DhCer"))
    
    reaction <- "RHEA:45303"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$DhSM, "SM(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[1]]$DhCer, "Cer(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_DhSM <=> M_Phosphocholine + M_H+ + M_DhCer")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_DhSM"))
    expect_equal(l[[2]]$reaction_product, c("M_Phosphocholine", "M_H+", "M_DhCer"))    
    
    ## fa_to_coa
    fa <- "FA(18:0)"
    substrates = list(FA = fa)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:15421"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_ATP", "M_CoA", "M_FA"))
    expect_equal(l[[2]]$reaction_product, c("M_PPi", "M_AMP", "M_AcylCoA")) 
    
    reaction <- "RHEA:15422"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_ATP + M_CoA + M_FA => M_PPi + M_AMP + M_AcylCoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_ATP", "M_CoA", "M_FA"))
    expect_equal(l[[2]]$reaction_product, c("M_PPi", "M_AMP", "M_AcylCoA"))
    
    reaction <- "RHEA:15423"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_ATP + M_CoA + M_FA <= M_PPi + M_AMP + M_AcylCoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_ATP", "M_CoA", "M_FA"))
    expect_equal(l[[2]]$reaction_product, c("M_PPi", "M_AMP", "M_AcylCoA"))
    
    reaction <- "RHEA:15424"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_ATP + M_CoA + M_FA <=> M_PPi + M_AMP + M_AcylCoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_ATP", "M_CoA", "M_FA"))
    expect_equal(l[[2]]$reaction_product, c("M_PPi", "M_AMP", "M_AcylCoA"))
    
    reaction <- "RHEA:38883"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_ATP", "M_CoA", "M_FA"))
    expect_equal(l[[2]]$reaction_product, c("M_PPi", "M_AMP", "M_AcylCoA"))
    
    reaction <- "RHEA:38884"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_ATP + M_CoA + M_FA => M_PPi + M_AMP + M_AcylCoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_ATP", "M_CoA", "M_FA"))
    expect_equal(l[[2]]$reaction_product, c("M_PPi", "M_AMP", "M_AcylCoA"))
    
    reaction <- "RHEA:38885"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_ATP + M_CoA + M_FA <= M_PPi + M_AMP + M_AcylCoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_ATP", "M_CoA", "M_FA"))
    expect_equal(l[[2]]$reaction_product, c("M_PPi", "M_AMP", "M_AcylCoA"))
    
    reaction <- "RHEA:38886"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_ATP + M_CoA + M_FA <=> M_PPi + M_AMP + M_AcylCoA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_ATP", "M_CoA", "M_FA"))
    expect_equal(l[[2]]$reaction_product, c("M_PPi", "M_AMP", "M_AcylCoA"))
    
    ## lcl_to_cl
    lcl <- "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])"
    acylcoa <- "CoA(18:4(6Z,9Z,12Z,15Z))"
    substrates <- list(LCL = lcl, AcylCoA = acylcoa)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:35839"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    reaction <- "RHEA:35840"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    reaction <- "RHEA:35841"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    reaction <- "RHEA:35842"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    ## lnape_to_gpnae
    lnape <- "NAPE(14:0/0:0/0:0)"
    substrates <- list(LNAPE = lnape)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:45420"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$LNAPE, "NAPE(14:0/0:0/0:0)")
    expect_equal(l[[1]]$GPNAE, "GPNAE(0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_LNAPE + M_H2O <=> M_GPNAE + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_LNAPE", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_GPNAE", "M_FA"))
    
    ## lpa_to_pa
    lpa <- "PA(18:0/0:0)"
    acylcoa <- "CoA(14:0)"
    substrates <- list(sn1LPA = lpa, AcylCoA = acylcoa)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:19709"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(l[[1]]$PA, "PA(18:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_LPA + M_AcylCoA = M_CoA + M_PA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_LPA", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_PA"))
    
    reaction <- "RHEA:19710"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(l[[1]]$PA, "PA(18:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_LPA + M_AcylCoA => M_CoA + M_PA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_LPA", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_PA"))
    
    reaction <- "RHEA:19711"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(l[[1]]$PA, "PA(18:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_LPA + M_AcylCoA <= M_CoA + M_PA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_LPA", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_PA"))
    
    reaction <- "RHEA:19712"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(l[[1]]$PA, "PA(18:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,  "M_LPA + M_AcylCoA <=> M_CoA + M_PA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_LPA", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_PA"))
    
    ## lpao_to_pao
    lpao <- "PA(O-18:0/0:0)"
    acylcoa <- "CoA(14:0)"
    substrates <- list(sn1LPAO = lpao, AcylCoA = acylcoa)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36235"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    reaction <- "RHEA:36236"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    reaction <- "RHEA:36237"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    reaction <- "RHEA:36238"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    ## sn1lpc_to_fa
    sn1lpc <- "PC(14:0/0:0)"
    substrates <- list(sn1LPC = sn1lpc)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:15177"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPC, "PC(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPC = M_Glycerophosphocholine + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPC"))
    expect_equal(l[[2]]$reaction_product, c("M_Glycerophosphocholine", "M_H+", "M_FA"))
    
    reaction <- "RHEA:15178"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPC, "PC(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPC => M_Glycerophosphocholine + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPC"))
    expect_equal(l[[2]]$reaction_product, c("M_Glycerophosphocholine", "M_H+", "M_FA"))
    
    reaction <- "RHEA:15179"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPC, "PC(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPC <= M_Glycerophosphocholine + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPC"))
    expect_equal(l[[2]]$reaction_product, c("M_Glycerophosphocholine", "M_H+", "M_FA"))
    
    reaction <- "RHEA:15180"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPC, "PC(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPC <=> M_Glycerophosphocholine + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPC"))
    expect_equal(l[[2]]$reaction_product, c("M_Glycerophosphocholine", "M_H+", "M_FA"))
    
    ## sn2lpc_to_fa
    sn2lpc <- "PC(0:0/14:0)"
    substrates <- list(sn2LPC = sn2lpc)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:44696"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/14:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_2-LPC = M_Glycerophosphocholine + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_2-LPC"))
    expect_equal(l[[2]]$reaction_product, c("M_Glycerophosphocholine", "M_H+", "M_FA"))
    
    reaction <- "RHEA:44697"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/14:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_2-LPC => M_Glycerophosphocholine + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_2-LPC"))
    expect_equal(l[[2]]$reaction_product, c("M_Glycerophosphocholine", "M_H+", "M_FA"))
    
    reaction <- "RHEA:44698"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/14:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_2-LPC <= M_Glycerophosphocholine + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_2-LPC"))
    expect_equal(l[[2]]$reaction_product, c("M_Glycerophosphocholine", "M_H+", "M_FA"))
    
    reaction <- "RHEA:44699"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/14:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_2-LPC <=> M_Glycerophosphocholine + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_2-LPC"))
    expect_equal(l[[2]]$reaction_product, c("M_Glycerophosphocholine", "M_H+", "M_FA"))
    
    ## sn1lpc_to_pc
    sn1lpc <- "PC(14:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPC = sn1lpc, AcylCoA = acylcoa)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:12937"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    reaction <- "RHEA:12938"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    reaction <- "RHEA:12939"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    reaction <- "RHEA:12940"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    ## sn1lpe_to_fa
    pe <- "PE(14:0/0:0)"
    substrates <- list(sn1LPE = pe)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:32967"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPE = M_H+ + M_FA + M_Glycerophosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPE"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_Glycerophosphoethanolamine"))
    
    reaction <- "RHEA:32968"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPE => M_H+ + M_FA + M_Glycerophosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPE"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_Glycerophosphoethanolamine"))
    
    reaction <- "RHEA:32969"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPE <= M_H+ + M_FA + M_Glycerophosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPE"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_Glycerophosphoethanolamine"))
    
    reaction <- "RHEA:32970"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPE <=> M_H+ + M_FA + M_Glycerophosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPE"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_Glycerophosphoethanolamine"))
    
    ## sn2lpe_to_fa
    pe <- "PE(0:0/14:0)"
    substrates <- list(sn2LPE = pe)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "sn2lpe_to_fa"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2LPE, "PE(0:0/14:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_2-LPE <=> M_H+ + M_FA + M_Glycerophosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_2-LPE"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_Glycerophosphoethanolamine"))
    
    ## sn1lpe_to_pe
    pe <- "PE(14:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPE = pe, AcylCoA = acylcoa)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:32995"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PE, "PE(14:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_1-LPE = M_CoA + M_PE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_1-LPE"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_PE"))
    
    reaction <- "RHEA:32996"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PE, "PE(14:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_1-LPE => M_CoA + M_PE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_1-LPE"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_PE"))
    
    reaction <- "RHEA:32997"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PE, "PE(14:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_1-LPE <= M_CoA + M_PE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_1-LPE"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_PE"))
    
    reaction <- "RHEA:32998"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PE, "PE(14:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_1-LPE <=> M_CoA + M_PE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_1-LPE"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_PE"))
    
    ## sn1lpi_to_pi
    sn1lpi <- "PI(16:0/0:0)"
    acylcoa <- "CoA(18:1(9Z))"
    substrates = list(sn1LPI = sn1lpi, AcylCoA = acylcoa)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:33195"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    reaction <- "RHEA:33196"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    reaction <- "RHEA:33197"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    reaction <- "RHEA:33198"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    ## lpeo_to_peo
    lpeo <- "PE(O-16:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPEO = lpeo, AcylCoA = acylcoa)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "lpeo_to_peo"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEO, "PE(O-16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PEO, "PE(O-16:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_ak2lgpe <=> M_CoA + M_PE-O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_ak2lgpe"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_PE-O"))
    
    ## lpep_to_pep
    lpep <- "PE(P-16:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPEP = lpep, AcylCoA = acylcoa)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:16245"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PEP, "PE(P-16:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_LPE-P = M_CoA + M_PE-P")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_LPE-P"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_PE-P"))
    
    reaction <- "RHEA:16246"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PEP, "PE(P-16:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_LPE-P => M_CoA + M_PE-P")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_LPE-P"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_PE-P"))
    
    reaction <- "RHEA:16247"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PEP, "PE(P-16:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_LPE-P <= M_CoA + M_PE-P")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_LPE-P"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_PE-P"))
    
    reaction <- "RHEA:16248"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(l[[1]]$PEP, "PE(P-16:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_LPE-P <=> M_CoA + M_PE-P")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_LPE-P"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_PE-P"))
    
    ## sn1mg_to_dg
    sn1mg <- "MG(14:0/0:0/0:0)"
    acylcoa <- "CoA(16:0)"
    substrates <- list(sn1MG = sn1mg, AcylCoA = acylcoa)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:38463"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_AcylCoA = M_CoA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:38464"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_AcylCoA => M_CoA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:38465"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_AcylCoA <= M_CoA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:38466"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_AcylCoA <=> M_CoA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:39943"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_AcylCoA = M_CoA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:39944"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_AcylCoA => M_CoA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:39945"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_AcylCoA <= M_CoA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:39946"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_1-MG + M_AcylCoA <=> M_CoA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    ## sn2mg_to_dg
    sn2mg <- "MG(0:0/14:0/0:0)"
    acylcoa <- "CoA(16:0)"
    substrates <- list(sn2MG = sn2mg, AcylCoA = acylcoa)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:32947"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_1-MG = M_CoA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_1-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:32948"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_1-MG => M_CoA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_1-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:32949"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_1-MG <= M_CoA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_1-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:32950"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_1-MG <=> M_CoA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_1-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:16741"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_1-MG = M_CoA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_1-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:16742"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_1-MG => M_CoA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_1-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:16743"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_1-MG <= M_CoA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_1-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:16744"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_1-MG <=> M_CoA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_1-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    ## sn1mg_to_fa
    sn1mg <- "MG(14:0/0:0/0:0)"
    substrates <- list(sn1MG = sn1mg)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:34019"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-MG = M_Glycerol + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_Glycerol", "M_H+", "M_FA"))
    
    reaction <- "RHEA:34020"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-MG => M_Glycerol + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_Glycerol", "M_H+", "M_FA"))
    
    reaction <- "RHEA:34021"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-MG <= M_Glycerol + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_Glycerol", "M_H+", "M_FA"))
    
    reaction <- "RHEA:34022"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-MG <=> M_Glycerol + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_Glycerol", "M_H+", "M_FA"))
    
    ## sn2mg_to_fa
    sn2mg <- "MG(0:0/14:0/0:0)"
    substrates <- list(sn2MG = sn2mg)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:32871"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_2-MG = M_Glycerol + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_2-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_Glycerol", "M_H+", "M_FA"))
    
    reaction <- "RHEA:32872"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_2-MG => M_Glycerol + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_2-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_Glycerol", "M_H+", "M_FA"))
    
    reaction <- "RHEA:32873"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_2-MG <= M_Glycerol + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_2-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_Glycerol", "M_H+", "M_FA"))
    
    reaction <- "RHEA:32874"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_2-MG <=> M_Glycerol + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_2-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_Glycerol", "M_H+", "M_FA"))
    
    ## sn1mg_to_lpa
    sn1mg <- "MG(14:0/0:0/0:0)"
    substrates <- list(sn1MG = sn1mg)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:33747"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$sn1LPA, "PA(14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_ATP + M_1-MG = M_H+ + M_ADP + M_LPA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_ATP", "M_1-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_ADP", "M_LPA"))
    
    reaction <- "RHEA:33748"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$sn1LPA, "PA(14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_ATP + M_1-MG => M_H+ + M_ADP + M_LPA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_ATP", "M_1-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_ADP", "M_LPA"))
    
    reaction <- "RHEA:33749"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$sn1LPA, "PA(14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_ATP + M_1-MG <= M_H+ + M_ADP + M_LPA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_ATP", "M_1-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_ADP", "M_LPA"))
    
    reaction <- "RHEA:33750"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(l[[1]]$sn1LPA, "PA(14:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_ATP + M_1-MG <=> M_H+ + M_ADP + M_LPA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_ATP", "M_1-MG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_ADP", "M_LPA"))
    
    ## sn2mg_to_sn1mg
    sn2mg <- "MG(0:0/14:0/0:0)"
    substrates = list(sn2MG = sn2mg)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    reaction <- "sn2mg_to_sn1mg"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    ## nae_to_fa
    nae <- "NAE(18:0)"
    substrates <- list(NAE = nae)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "nae_to_fa"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$NAE, "NAE(18:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_NAE <=> M_Ethanolamine + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_NAE"))
    expect_equal(l[[2]]$reaction_product, c("M_Ethanolamine", "M_H+", "M_FA"))
    
    ## nape_to_lnape
    nape <- "NAPE(14:0/16:0/18:0)"
    substrates = list(NAPE = nape)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    reaction <- "nape_to_lnape"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(l[[1]]$LNAPE, "NAPE(14:0/0:0/18:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_NAPE + M_H2O <=> M_LNAPE + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_NAPE", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_LNAPE", "M_FA"))
    
    ## nape_to_nae
    nape <- "NAPE(14:0/16:0/18:0)"
    substrates = list(NAPE = nape)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    reaction <- "nape_to_nae"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(l[[1]]$NAE, "NAE(18:0)")
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_NAPE + M_H2O <=> M_NAE + M_PA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_NAPE", "M_H2O"))
    expect_equal(l[[2]]$reaction_product, c("M_NAE", "M_PA"))
    
    ## nape_to_pnae
    nape <- "NAPE(14:0/16:0/18:0)"
    substrates <- list(NAPE = nape)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    reaction <- "nape_to_pnae"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    ## napeo_to_nae
    napeo <- "NAPE(O-18:0/16:0/14:0)"
    substrates = list(NAPEO = napeo)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    reaction <- "napeo_to_nae"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    ## pa_to_cdpdg
    pa <- "PA(14:0/16:0)"
    substrates <- list(PA = pa)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:16229"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[1]]$CDPDG, "CDP-DG(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CTP + M_PA = M_PPi + M_CDP-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CTP", "M_PA"))
    expect_equal(l[[2]]$reaction_product, c("M_PPi", "M_CDP-DG"))
    
    reaction <- "RHEA:16230"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[1]]$CDPDG, "CDP-DG(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CTP + M_PA => M_PPi + M_CDP-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CTP", "M_PA"))
    expect_equal(l[[2]]$reaction_product, c("M_PPi", "M_CDP-DG"))
    
    reaction <- "RHEA:16231"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[1]]$CDPDG, "CDP-DG(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CTP + M_PA <= M_PPi + M_CDP-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CTP", "M_PA"))
    expect_equal(l[[2]]$reaction_product, c("M_PPi", "M_CDP-DG"))
    
    reaction <- "RHEA:16232"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[1]]$CDPDG, "CDP-DG(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_CTP + M_PA <=> M_PPi + M_CDP-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CTP", "M_PA"))
    expect_equal(l[[2]]$reaction_product, c("M_PPi", "M_CDP-DG"))
    
    ## pa_to_dg
    pa <- "PA(14:0/16:0)"
    substrates = list(PA = pa)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:27429"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PA = M_Pi + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PA"))
    expect_equal(l[[2]]$reaction_product, c("M_Pi", "M_1,2-DG"))
    
    reaction <- "RHEA:27430"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PA => M_Pi + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PA"))
    expect_equal(l[[2]]$reaction_product, c("M_Pi", "M_1,2-DG"))
    
    reaction <- "RHEA:27431"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PA <= M_Pi + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PA"))
    expect_equal(l[[2]]$reaction_product, c("M_Pi", "M_1,2-DG"))
    
    reaction <- "RHEA:27432"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PA <=> M_Pi + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PA"))
    expect_equal(l[[2]]$reaction_product, c("M_Pi", "M_1,2-DG"))
    
    ## pao_to_dgo
    pao <- "PA(O-14:0/16:0)"
    substrates = list(PAO = pao)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36239"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PAO, "PA(O-14:0/16:0)")
    expect_equal(l[[1]]$DGO, "DG(O-14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PA-O = M_Pi + M_DG-O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PA-O"))
    expect_equal(l[[2]]$reaction_product, c("M_Pi", "M_DG-O"))
    
    reaction <- "RHEA:36240"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PAO, "PA(O-14:0/16:0)")
    expect_equal(l[[1]]$DGO, "DG(O-14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PA-O => M_Pi + M_DG-O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PA-O"))
    expect_equal(l[[2]]$reaction_product, c("M_Pi", "M_DG-O"))
    
    
    reaction <- "RHEA:36241"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PAO, "PA(O-14:0/16:0)")
    expect_equal(l[[1]]$DGO, "DG(O-14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PA-O <= M_Pi + M_DG-O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PA-O"))
    expect_equal(l[[2]]$reaction_product, c("M_Pi", "M_DG-O"))
    
    reaction <- "RHEA:36242"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PAO, "PA(O-14:0/16:0)")
    expect_equal(l[[1]]$DGO, "DG(O-14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PA-O <=> M_Pi + M_DG-O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PA-O"))
    expect_equal(l[[2]]$reaction_product, c("M_Pi", "M_DG-O"))
    
    ## pc_to_dg
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:10604"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$DG, "DG(20:0/18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PC = M_Phosphocholine + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_Phosphocholine", "M_1,2-DG"))
    
    reaction <- "RHEA:10605"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$DG, "DG(20:0/18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PC => M_Phosphocholine + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_Phosphocholine", "M_1,2-DG"))
    
    reaction <- "RHEA:10606"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$DG, "DG(20:0/18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PC <= M_Phosphocholine + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_Phosphocholine", "M_1,2-DG"))
    
    reaction <- "RHEA:10607"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$DG, "DG(20:0/18:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PC <=> M_Phosphocholine + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_Phosphocholine", "M_1,2-DG"))
    
    ## pc_to_sn1lpc
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:15801"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn1LPC, "PC(20:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PC = M_1-LPC + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPC", "M_FA"))
    
    reaction <- "RHEA:15802"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn1LPC, "PC(20:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PC => M_1-LPC + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPC", "M_FA"))
    
    reaction <- "RHEA:15803"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn1LPC, "PC(20:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PC <= M_1-LPC + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPC", "M_FA"))
    
    reaction <- "RHEA:15804"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn1LPC, "PC(20:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PC <=> M_1-LPC + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPC", "M_FA"))
    
    ## pc_to_sn2lpc
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:18689"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/18:0)")
    expect_equal(l[[1]]$FA, "FA(20:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PC = M_2-LPC + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_2-LPC", "M_FA"))
    
    reaction <- "RHEA:18690"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/18:0)")
    expect_equal(l[[1]]$FA, "FA(20:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PC => M_2-LPC + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_2-LPC", "M_FA"))
    
    reaction <- "RHEA:18691"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/18:0)")
    expect_equal(l[[1]]$FA, "FA(20:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PC <= M_2-LPC + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_2-LPC", "M_FA"))
    
    reaction <- "RHEA:18692"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/18:0)")
    expect_equal(l[[1]]$FA, "FA(20:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PC <=> M_2-LPC + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_2-LPC", "M_FA"))
    
    ## pc_to_pa
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:14445"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$PA, "PA(20:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PC = M_Choline + M_PA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_Choline", "M_PA"))
    
    reaction <- "RHEA:14446"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$PA, "PA(20:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PC => M_Choline + M_PA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_Choline", "M_PA"))
    
    reaction <- "RHEA:14447"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$PA, "PA(20:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PC <= M_Choline + M_PA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_Choline", "M_PA"))
    
    reaction <- "RHEA:14448"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$PA, "PA(20:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PC <=> M_Choline + M_PA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_Choline", "M_PA"))
    
    ## pc_to_ps
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:45088"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$PS, "PS(20:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_L-Serine + M_PC = M_Choline + M_PS")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_L-Serine", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_Choline", "M_PS"))
    
    reaction <- "RHEA:45089"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$PS, "PS(20:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_L-Serine + M_PC => M_Choline + M_PS")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_L-Serine", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_Choline", "M_PS"))
    
    reaction <- "RHEA:45090"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$PS, "PS(20:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_L-Serine + M_PC <= M_Choline + M_PS")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_L-Serine", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_Choline", "M_PS"))
    
    reaction <- "RHEA:45091"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$PS, "PS(20:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_L-Serine + M_PC <=> M_Choline + M_PS")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_L-Serine", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_Choline", "M_PS"))
    
    ## pco_to_lpco
    pco <- "PC(O-16:0/14:0)"
    substrates <- list(PCO = pco)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36231"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PCO, "PC(O-16:0/14:0)")
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PC-O = M_H+ + M_LPC-O + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PC-O"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_LPC-O", "M_FA"))
    
    reaction <- "RHEA:36232"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PCO, "PC(O-16:0/14:0)")
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PC-O => M_H+ + M_LPC-O + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PC-O"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_LPC-O", "M_FA"))
    
    reaction <- "RHEA:36233"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PCO, "PC(O-16:0/14:0)")
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PC-O <= M_H+ + M_LPC-O + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PC-O"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_LPC-O", "M_FA"))
    
    reaction <- "RHEA:36234"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PCO, "PC(O-16:0/14:0)")
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PC-O <=> M_H+ + M_LPC-O + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PC-O"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_LPC-O", "M_FA"))
    
    ## lpco_to_lpao
    sn1lpco <- "PC(O-16:0/0:0)"
    substrates <- list(sn1LPCO = sn1lpco)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:39927"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$sn1LPAO, "PA(O-16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPC-O = M_1-LPA-O + M_H+ + M_Choline")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPC-O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPA-O", "M_H+", "M_Choline"))
    
    reaction <- "RHEA:39928"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$sn1LPAO, "PA(O-16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPC-O => M_1-LPA-O + M_H+ + M_Choline")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPC-O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPA-O", "M_H+", "M_Choline"))
    
    reaction <- "RHEA:39929"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$sn1LPAO, "PA(O-16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPC-O <= M_1-LPA-O + M_H+ + M_Choline")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPC-O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPA-O", "M_H+", "M_Choline"))
    
    reaction <- "RHEA:39930"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$sn1LPAO, "PA(O-16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPC-O <=> M_1-LPA-O + M_H+ + M_Choline")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPC-O"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPA-O", "M_H+", "M_Choline"))
    
    ## lpco_to_mgo
    sn1lpco <- "PC(O-16:0/0:0)"
    substrates <- list(sn1LPCO = sn1lpco)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36083"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$sn1MGO, "MG(O-16:0/0:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPC-O = M_H+ + M_Phosphocholine + M_1-MG-O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPC-O"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_Phosphocholine", "M_1-MG-O"))
    
    reaction <- "RHEA:36084"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$sn1MGO, "MG(O-16:0/0:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPC-O => M_H+ + M_Phosphocholine + M_1-MG-O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPC-O"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_Phosphocholine", "M_1-MG-O"))
    
    reaction <- "RHEA:36085"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$sn1MGO, "MG(O-16:0/0:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPC-O <= M_H+ + M_Phosphocholine + M_1-MG-O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPC-O"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_Phosphocholine", "M_1-MG-O"))
    
    reaction <- "RHEA:36086"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(l[[1]]$sn1MGO, "MG(O-16:0/0:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPC-O <=> M_H+ + M_Phosphocholine + M_1-MG-O")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPC-O"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_Phosphocholine", "M_1-MG-O"))
    
    ## lpco_to_pco
    sn1lpco <- "PC(O-16:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPCO = sn1lpco, AcylCoA = acylcoa)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:23992"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    reaction <- "RHEA:23993"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    reaction <- "RHEA:23994"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    reaction <- "RHEA:23995"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    ## pe_to_dg
    pe <- "PE(14:0/16:0)"
    substrates <- list(PE = pe)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "pe_to_dg"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PE <=> M_P-Ethanolamine + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_P-Ethanolamine", "M_1,2-DG"))
    
    ## pe_to_sn1lpe
    pe <- "PE(14:0/16:0)"
    substrates = list(PE = pe)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:44604"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PE = M_H+ + M_FA + M_1-LPE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_1-LPE"))
    
    reaction <- "RHEA:44605"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PE => M_H+ + M_FA + M_1-LPE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_1-LPE"))
    
    reaction <- "RHEA:44606"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PE <= M_H+ + M_FA + M_1-LPE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_1-LPE"))
    
    reaction <- "RHEA:44607"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PE <=> M_H+ + M_FA + M_1-LPE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_1-LPE"))
    
    ## pe_to_sn2lpe
    pe <- "PE(14:0/16:0)"
    substrates = list(PE = pe)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:44408"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$sn2LPE, "PE(0:0/14:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PE = M_H+ + M_FA + M_2-LPE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_2-LPE"))
    
    reaction <- "RHEA:44409"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$sn2LPE, "PE(0:0/14:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PE => M_H+ + M_FA + M_2-LPE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_2-LPE"))
    
    reaction <- "RHEA:44410"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$sn2LPE, "PE(0:0/14:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PE <= M_H+ + M_FA + M_2-LPE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_2-LPE"))
    
    reaction <- "RHEA:44411"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$sn2LPE, "PE(0:0/14:0)")
    expect_equal(l[[1]]$FA, "FA(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PE <=> M_H+ + M_FA + M_2-LPE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_2-LPE"))
    
    ## pe_to_nape_sn1
    pe <- "PE(14:0/16:0)"
    pc <- "PC(18:0/20:0)"
    substrates <- list(PE = pe, PC = pc)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "pe_to_nape_sn1"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$PC, "PC(18:0/20:0)")
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/20:0)")
    expect_equal(l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE + M_PC <=> M_NAPE + M_2-LPC")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_NAPE", "M_2-LPC"))
    
    ## pe_to_nape_sn2
    pe <- "PE(14:0/16:0)"
    pc <- "PC(18:0/20:0)"
    substrates <- list(PE = pe, PC = pc)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "pe_to_nape_sn2"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$PC, "PC(18:0/20:0)")
    expect_equal(l[[1]]$sn1LPC, "PC(18:0/0:0)")
    expect_equal(l[[1]]$NAPE, "NAPE(14:0/16:0/20:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE + M_PC <=> M_NAPE + M_1-LPC")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_NAPE", "M_1-LPC"))
    
    ## pe_to_pa
    pe <- "PE(14:0/16:0)"
    substrates <- list(PE = pe)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "pe_to_pa"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PE <=> M_Ethanolamine + M_H+ + M_PA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_Ethanolamine", "M_H+", "M_PA"))
    
    ## pe_to_ps
    pe <- "PE(14:0/16:0)"
    substrates <- list(PE = pe)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:27606" 
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$PS, "PS(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_L-Serine + M_PE = M_Ethanolamine + M_PS")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_L-Serine", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_Ethanolamine", "M_PS"))
    
    reaction <- "RHEA:27607"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$PS, "PS(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_L-Serine + M_PE => M_Ethanolamine + M_PS")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_L-Serine", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_Ethanolamine", "M_PS"))
    
    reaction <- "RHEA:27608"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$PS, "PS(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_L-Serine + M_PE <= M_Ethanolamine + M_PS")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_L-Serine", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_Ethanolamine", "M_PS"))
    
    reaction <- "RHEA:27609"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(l[[1]]$PS, "PS(14:0/16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_L-Serine + M_PE <=> M_Ethanolamine + M_PS")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_L-Serine", "M_PE"))
    expect_equal(l[[2]]$reaction_product, c("M_Ethanolamine", "M_PS"))
    
    ## peo_to_lpeo
    peo <- "PE(O-16:0/14:0)"
    substrates <- list(PEO = peo)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "peo_to_lpeo"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEO, "PE(O-16:0/14:0)")
    expect_equal(l[[1]]$sn1LPEO, "PE(O-16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PE-O <=> M_H+ + M_LPE-O + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PE-O"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_LPE-O", "M_FA"))
    
    ## peo_to_napeo_sn1
    peo <- "PE(O-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    substrates <- list(PEO = peo, PC = pc) ################
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "peo_to_napeo_sn1"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    ## peo_to_napeo_sn2
    peo <- "PE(O-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    substrates <- list(PEO = peo, PC = pc) ################
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "peo_to_napeo_sn2"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    ## peo_to_pep
    peo <- "PE(O-16:0/14:0)"
    substrates <- list(PEO = peo)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:22956"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    reaction <- "RHEA:22957"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    reaction <- "RHEA:22958"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    reaction <- "RHEA:22959"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    ## pep_to_lpep
    pep <- "PE(P-16:0/14:0)"
    substrates <- list(PEP = pep)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36195"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PE-P = M_H+ + M_LPE-P + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PE-P"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_LPE-P", "M_FA"))
    
    reaction <- "RHEA:36196"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PE-P => M_H+ + M_LPE-P + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PE-P"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_LPE-P", "M_FA"))
    
    reaction <- "RHEA:36197"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PE-P <= M_H+ + M_LPE-P + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PE-P"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_LPE-P", "M_FA"))
    
    reaction <- "RHEA:36198"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PE-P <=> M_H+ + M_LPE-P + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PE-P"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_LPE-P", "M_FA"))
    
    ## lpep_to_fal
    sn1lpep <- "PE(P-16:0/0:0)"
    substrates = list(sn1LPEP = sn1lpep)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:16905"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$FAL, "FAL(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,"M_H2O + M_1-LPE-P = M_FAL + M_Glycerophosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(l[[2]]$reaction_product, c("M_FAL", "M_Glycerophosphoethanolamine"))
    
    reaction <- "RHEA:16906"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$FAL, "FAL(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,"M_H2O + M_1-LPE-P => M_FAL + M_Glycerophosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(l[[2]]$reaction_product, c("M_FAL", "M_Glycerophosphoethanolamine"))
    
    reaction <- "RHEA:16907"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$FAL, "FAL(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,"M_H2O + M_1-LPE-P <= M_FAL + M_Glycerophosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(l[[2]]$reaction_product, c("M_FAL", "M_Glycerophosphoethanolamine"))
    
    reaction <- "RHEA:16908"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$FAL, "FAL(16:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula,"M_H2O + M_1-LPE-P <=> M_FAL + M_Glycerophosphoethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(l[[2]]$reaction_product, c("M_FAL", "M_Glycerophosphoethanolamine"))
    
    ## lpep_to_lpap
    sn1lpep <- "PE(P-16:0/0:0)"
    substrates <- list(sn1LPEP = sn1lpep)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36203"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$sn1LPAP, "PA(P-16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPE-P = M_LPA-P + M_H+ + M_Ethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA-P", "M_H+", "M_Ethanolamine"))
    
    reaction <- "RHEA:36204"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$sn1LPAP, "PA(P-16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPE-P => M_LPA-P + M_H+ + M_Ethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA-P", "M_H+", "M_Ethanolamine"))
    
    reaction <- "RHEA:36205"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$sn1LPAP, "PA(P-16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPE-P <= M_LPA-P + M_H+ + M_Ethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA-P", "M_H+", "M_Ethanolamine"))
    
    reaction <- "RHEA:36206"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$sn1LPAP, "PA(P-16:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPE-P <=> M_LPA-P + M_H+ + M_Ethanolamine")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(l[[2]]$reaction_product, c("M_LPA-P", "M_H+", "M_Ethanolamine"))
    
    ## lpep_to_mgp
    sn1lpep <- "PE(P-16:0/0:0)"
    substrates <- list(sn1LPEP = sn1lpep)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:36199"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$sn1MGP, "MG(P-16:0/0:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPE-P = M_Phosphoethanolamine + M_H+ + M_1-MG-P")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(l[[2]]$reaction_product, c("M_Phosphoethanolamine", "M_H+", "M_1-MG-P"))
    
    reaction <- "RHEA:36200"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$sn1MGP, "MG(P-16:0/0:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPE-P => M_Phosphoethanolamine + M_H+ + M_1-MG-P")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(l[[2]]$reaction_product, c("M_Phosphoethanolamine", "M_H+", "M_1-MG-P"))
    
    reaction <- "RHEA:36201"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$sn1MGP, "MG(P-16:0/0:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPE-P <= M_Phosphoethanolamine + M_H+ + M_1-MG-P")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(l[[2]]$reaction_product, c("M_Phosphoethanolamine", "M_H+", "M_1-MG-P"))
    
    reaction <- "RHEA:36202"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(l[[1]]$sn1MGP, "MG(P-16:0/0:0/0:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_1-LPE-P <=> M_Phosphoethanolamine + M_H+ + M_1-MG-P")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(l[[2]]$reaction_product, c("M_Phosphoethanolamine", "M_H+", "M_1-MG-P"))
    
    ## pep_to_napep_sn1
    pep <- "PE(P-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    substrates <- list(PEP = pep, PC = pc)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "pep_to_napep_sn1"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(l[[1]]$sn2LPC, "PC(0:0/18:0)")
    expect_equal(l[[1]]$NAPEP, "NAPE(P-16:0/14:0/20:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_PE-P + M_PC <=> M_NAPEP + M_2-LPC")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_PE-P", "M_PC"))
    expect_equal(l[[2]]$reaction_product, c("M_NAPEP", "M_2-LPC"))
    
    ## pep_to_napep_sn2
    pep <- "PE(P-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    substrates <- list(PEP = pep, PC = pc)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "pep_to_napep_sn2"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    
    ## pg_to_cl
    pg <- "PG(18:4(6Z,9Z,12Z,15Z)/14:0)"
    cdpdg <- "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)"
    substrates = list(PG = pg, CDPDG = cdpdg)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:32931"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    expect_equal(l[[2]]$reaction_formula, "M_CDP-DG + M_PG = M_H+ + M_CMP + M_CL")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-DG", "M_PG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_CL"))
    
    reaction <- "RHEA:32932"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    expect_equal(l[[2]]$reaction_formula, "M_CDP-DG + M_PG => M_H+ + M_CMP + M_CL")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-DG", "M_PG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_CL"))
    
    reaction <- "RHEA:32933"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    expect_equal(l[[2]]$reaction_formula, "M_CDP-DG + M_PG <= M_H+ + M_CMP + M_CL")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-DG", "M_PG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_CL"))
    
    reaction <- "RHEA:32934"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    expect_equal(l[[2]]$reaction_formula, "M_CDP-DG + M_PG <=> M_H+ + M_CMP + M_CL")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_CDP-DG", "M_PG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_CMP", "M_CL"))
    
    ## pgp_to_pg
    pgp <- "PGP(16:0/14:0)"
    substrates <- list(PGP = pgp)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:33751"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PGP, "PGP(16:0/14:0)")
    expect_equal(l[[1]]$PG, "PG(16:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PGP = M_Pi + M_PG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PGP"))
    expect_equal(l[[2]]$reaction_product, c("M_Pi", "M_PG"))
    
    reaction <- "RHEA:33752"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PGP, "PGP(16:0/14:0)")
    expect_equal(l[[1]]$PG, "PG(16:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PGP => M_Pi + M_PG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PGP"))
    expect_equal(l[[2]]$reaction_product, c("M_Pi", "M_PG"))
    
    reaction <- "RHEA:33753"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PGP, "PGP(16:0/14:0)")
    expect_equal(l[[1]]$PG, "PG(16:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PGP <= M_Pi + M_PG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PGP"))
    expect_equal(l[[2]]$reaction_product, c("M_Pi", "M_PG"))
    
    reaction <- "RHEA:33754"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PGP, "PGP(16:0/14:0)")
    expect_equal(l[[1]]$PG, "PG(16:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PGP <=> M_Pi + M_PG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PGP"))
    expect_equal(l[[2]]$reaction_product, c("M_Pi", "M_PG"))
    
    ## pi_to_dg
    pi <- "PI(16:0/18:1(9Z))"
    substrates <- list(PI = pi)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:43484"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(l[[1]]$DG, "DG(16:0/18:1(9Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PI = M_myo-Inositol-1-P + M_H+ + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PI"))
    expect_equal(l[[2]]$reaction_product, c("M_myo-Inositol-1-P", "M_H+", "M_1,2-DG"))
    
    reaction <- "RHEA:43485"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(l[[1]]$DG, "DG(16:0/18:1(9Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PI => M_myo-Inositol-1-P + M_H+ + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PI"))
    expect_equal(l[[2]]$reaction_product, c("M_myo-Inositol-1-P", "M_H+", "M_1,2-DG"))
    
    reaction <- "RHEA:43486"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(l[[1]]$DG, "DG(16:0/18:1(9Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PI <= M_myo-Inositol-1-P + M_H+ + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PI"))
    expect_equal(l[[2]]$reaction_product, c("M_myo-Inositol-1-P", "M_H+", "M_1,2-DG"))
    
    reaction <- "RHEA:43487"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(l[[1]]$DG, "DG(16:0/18:1(9Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PI <=> M_myo-Inositol-1-P + M_H+ + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PI"))
    expect_equal(l[[2]]$reaction_product, c("M_myo-Inositol-1-P", "M_H+", "M_1,2-DG"))
    
    ## pi_to_sn1lpi
    pi <- "PI(16:0/18:1(9Z))"
    substrates <- list(PI = pi)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:18001"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:1(9Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PI = M_1-LPI + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PI"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPI", "M_H+", "M_FA"))
    
    reaction <- "RHEA:18002"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:1(9Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PI => M_1-LPI + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PI"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPI", "M_H+", "M_FA"))
    
    reaction <- "RHEA:18003"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:1(9Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PI <= M_1-LPI + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PI"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPI", "M_H+", "M_FA"))
    
    reaction <- "RHEA:18004"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(l[[1]]$FA, "FA(18:1(9Z))")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_PI <=> M_1-LPI + M_H+ + M_FA")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_PI"))
    expect_equal(l[[2]]$reaction_product, c("M_1-LPI", "M_H+", "M_FA"))
    
    ## ps_to_pe
    ps <- "PS(14:0/14:0)"
    substrates = list(PS = ps)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:20828"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PS, "PS(14:0/14:0)")
    expect_equal(l[[1]]$PE, "PE(14:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H+ + M_PS = M_CO2 + M_PE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H+", "M_PS"))
    expect_equal(l[[2]]$reaction_product, c("M_CO2", "M_PE"))
    
    reaction <- "RHEA:20829"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PS, "PS(14:0/14:0)")
    expect_equal(l[[1]]$PE, "PE(14:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H+ + M_PS => M_CO2 + M_PE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H+", "M_PS"))
    expect_equal(l[[2]]$reaction_product, c("M_CO2", "M_PE"))
    
    reaction <- "RHEA:20830"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PS, "PS(14:0/14:0)")
    expect_equal(l[[1]]$PE, "PE(14:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H+ + M_PS <= M_CO2 + M_PE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H+", "M_PS"))
    expect_equal(l[[2]]$reaction_product, c("M_CO2", "M_PE"))
    
    reaction <- "RHEA:20831"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$PS, "PS(14:0/14:0)")
    expect_equal(l[[1]]$PE, "PE(14:0/14:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H+ + M_PS <=> M_CO2 + M_PE")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H+", "M_PS"))
    expect_equal(l[[2]]$reaction_product, c("M_CO2", "M_PE"))
    
    ## sm_to_cer
    sm <- "SM(16:0(3OH,4OH,15Me)/12:0)"
    substrates <- list(SM = sm)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "sm_to_cer"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$SM, "SM(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_SM <=> M_Phosphocholine + M_H+ + M_Cer")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_SM"))
    expect_equal(l[[2]]$reaction_product, c("M_Phosphocholine", "M_H+", "M_Cer"))
    
    ## sphinga_to_dhcer
    acylcoa <- "CoA(12:0)" 
    substrates <- list(AcylCoA = acylcoa) #######################################
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:53424"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(12:0)")
    expect_equal(l[[1]]$DhCer, "Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_Sphinganine = M_CoA + M_DhCer")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_Sphinganine"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_DhCer"))
    
    reaction <- "RHEA:53425"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(12:0)")
    expect_equal(l[[1]]$DhCer, "Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_Sphinganine => M_CoA + M_DhCer")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_Sphinganine"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_DhCer"))
    
    reaction <- "RHEA:53426"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(12:0)")
    expect_equal(l[[1]]$DhCer, "Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_Sphinganine <= M_CoA + M_DhCer")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_Sphinganine"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_DhCer"))
    
    reaction <- "RHEA:53427"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    l <- .create_list_reactants_with_template(df_reaction = df_reaction,
        template = template)
    expect_equal(l[[1]]$AcylCoA, "CoA(12:0)")
    expect_equal(l[[1]]$DhCer, "Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l[[2]]$reaction_name, "")
    expect_equal(l[[2]]$reaction_formula, "M_AcylCoA + M_Sphinganine <=> M_CoA + M_DhCer")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_AcylCoA", "M_Sphinganine"))
    expect_equal(l[[2]]$reaction_product, c("M_CoA", "M_DhCer"))
    
    ## tg_to_dg
    tg <- "TG(18:0/16:0/14:0)"
    substrates = list(TG = tg)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    
    reaction <- "RHEA:33271"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_TG = M_H+ + M_FA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_TG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_1,2-DG"))
    
    reaction <- "RHEA:33272"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction) 
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_TG => M_H+ + M_FA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_TG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_1,2-DG"))
    
    reaction <- "RHEA:33273"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_TG <= M_H+ + M_FA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_TG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_1,2-DG"))
    
    reaction <- "RHEA:33274"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_TG <=> M_H+ + M_FA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_TG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_1,2-DG"))
    
    reaction <- "RHEA:44864"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_TG = M_H+ + M_FA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_TG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_1,2-DG"))
    
    reaction <- "RHEA:44865"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_TG => M_H+ + M_FA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_TG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_1,2-DG"))
    
    reaction <- "RHEA:44866"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_TG <= M_H+ + M_FA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_TG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_1,2-DG"))
    
    reaction <- "RHEA:44867"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
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
    expect_equal(l[[2]]$reaction_formula, "M_H2O + M_TG <=> M_H+ + M_FA + M_1,2-DG")
    expect_equal(l[[2]]$reaction_isReversible, "")
    expect_equal(l[[2]]$reaction_geneAssociation, "")
    expect_equal(l[[2]]$reaction_pathway, "")
    expect_equal(l[[2]]$reaction_substrate, c("M_H2O", "M_TG"))
    expect_equal(l[[2]]$reaction_product, c("M_H+", "M_FA", "M_1,2-DG"))
})


