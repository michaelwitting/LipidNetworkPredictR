## function .create_reaction
test_that(".create_df_with_template works", {
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, "DHAP(18:0) + FAO(16:0) = H+ + FA(18:0) + DHAP(O-16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "DHAP(18:0) + FAO(16:0)")
    expect_equal(df$reaction_product, "H+ + FA(18:0) + DHAP(O-16:0)")
    
    reaction <- "RHEA:36172"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, "DHAP(18:0) + FAO(16:0) => H+ + FA(18:0) + DHAP(O-16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "DHAP(18:0) + FAO(16:0)")
    expect_equal(df$reaction_product, "H+ + FA(18:0) + DHAP(O-16:0)")
    
    reaction <- "RHEA:36173"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, "DHAP(18:0) + FAO(16:0) <= H+ + FA(18:0) + DHAP(O-16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "DHAP(18:0) + FAO(16:0)")
    expect_equal(df$reaction_product, "H+ + FA(18:0) + DHAP(O-16:0)")
    
    reaction <- "RHEA:36174"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, "DHAP(18:0) + FAO(16:0) <=> H+ + FA(18:0) + DHAP(O-16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "DHAP(18:0) + FAO(16:0)")
    expect_equal(df$reaction_product, "H+ + FA(18:0) + DHAP(O-16:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, "H+ + NADPH + DHAP(O-18:0) = PA(O-18:0/0:0) + NADP+")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H+ + NADPH + DHAP(O-18:0)")
    expect_equal(df$reaction_product, "PA(O-18:0/0:0) + NADP+")
    
    reaction <- "RHEA:36176"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, "H+ + NADPH + DHAP(O-18:0) => PA(O-18:0/0:0) + NADP+")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H+ + NADPH + DHAP(O-18:0)")
    expect_equal(df$reaction_product, "PA(O-18:0/0:0) + NADP+")
    
    reaction <- "RHEA:36177"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, "H+ + NADPH + DHAP(O-18:0) <= PA(O-18:0/0:0) + NADP+")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H+ + NADPH + DHAP(O-18:0)")
    expect_equal(df$reaction_product, "PA(O-18:0/0:0) + NADP+")
    
    reaction <- "RHEA:36178"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, "H+ + NADPH + DHAP(O-18:0) <=> PA(O-18:0/0:0) + NADP+")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H+ + NADPH + DHAP(O-18:0)")
    expect_equal(df$reaction_product, "PA(O-18:0/0:0) + NADP+")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, "H2O + CerP(16:0(3OH,4OH,15Me)/18:0) <=> Pi + Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + CerP(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(df$reaction_product, "Pi + Cer(16:0(3OH,4OH,15Me)/18:0)")

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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, "Glycerol-3-P + CDP-DG(12:0(11Me)/14:0) = H+ + CMP + PGP(12:0(11Me)/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "Glycerol-3-P + CDP-DG(12:0(11Me)/14:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PGP(12:0(11Me)/14:0)")
    
    reaction <- "RHEA:12594"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, "Glycerol-3-P + CDP-DG(12:0(11Me)/14:0) => H+ + CMP + PGP(12:0(11Me)/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "Glycerol-3-P + CDP-DG(12:0(11Me)/14:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PGP(12:0(11Me)/14:0)")
    
    reaction <- "RHEA:12595"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, "Glycerol-3-P + CDP-DG(12:0(11Me)/14:0) <= H+ + CMP + PGP(12:0(11Me)/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "Glycerol-3-P + CDP-DG(12:0(11Me)/14:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PGP(12:0(11Me)/14:0)")
    
    reaction <- "RHEA:12596"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, "Glycerol-3-P + CDP-DG(12:0(11Me)/14:0) <=> H+ + CMP + PGP(12:0(11Me)/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "Glycerol-3-P + CDP-DG(12:0(11Me)/14:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PGP(12:0(11Me)/14:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, "myo-Inositol + CDP-DG(12:0(11Me)/14:0) = H+ + CMP + PI(12:0(11Me)/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "myo-Inositol + CDP-DG(12:0(11Me)/14:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PI(12:0(11Me)/14:0)")
    
    reaction <- "RHEA:11581"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, "myo-Inositol + CDP-DG(12:0(11Me)/14:0) => H+ + CMP + PI(12:0(11Me)/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "myo-Inositol + CDP-DG(12:0(11Me)/14:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PI(12:0(11Me)/14:0)")
    
    reaction <- "RHEA:11582"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, "myo-Inositol + CDP-DG(12:0(11Me)/14:0) <= H+ + CMP + PI(12:0(11Me)/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "myo-Inositol + CDP-DG(12:0(11Me)/14:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PI(12:0(11Me)/14:0)")
    
    reaction <- "RHEA:11583"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, "myo-Inositol + CDP-DG(12:0(11Me)/14:0) <=> H+ + CMP + PI(12:0(11Me)/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "myo-Inositol + CDP-DG(12:0(11Me)/14:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PI(12:0(11Me)/14:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "ATP + Cer(16:0(3OH,4OH,15Me)/18:0) <=> H+ + ADP + CerP(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "ATP + Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(df$reaction_product, "H+ + ADP + CerP(16:0(3OH,4OH,15Me)/18:0)")                  
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "UDP-Glucose + Cer(16:0(3OH,4OH,15Me)/18:0) <=> H+ + UDP + GlcCer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "UDP-Glucose + Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(df$reaction_product, "H+ + UDP + GlcCer(16:0(3OH,4OH,15Me)/18:0)")     
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PC + Cer(16:0(3OH,4OH,15Me)/18:0) <=> DG + SM(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PC + Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(df$reaction_product, "DG + SM(16:0(3OH,4OH,15Me)/18:0)") 
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) = CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + H+ + FA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(df$reaction_product, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + H+ + FA(18:4(6Z,9Z,12Z,15Z))") 
    
    reaction <- "RHEA:32936"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) => CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + H+ + FA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(df$reaction_product, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + H+ + FA(18:4(6Z,9Z,12Z,15Z))") 
    
    reaction <- "RHEA:32937"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) <= CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + H+ + FA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(df$reaction_product, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + H+ + FA(18:4(6Z,9Z,12Z,15Z))") 
    
    reaction <- "RHEA:32938"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) <=> CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + H+ + FA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(df$reaction_product, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + H+ + FA(18:4(6Z,9Z,12Z,15Z))") 
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "Dihydroxyacetone-P + CoA(18:0) = CoA + DHAP(18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "Dihydroxyacetone-P + CoA(18:0)")
    expect_equal(df$reaction_product, "CoA + DHAP(18:0)") 
    
    reaction <- "RHEA:17658"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "Dihydroxyacetone-P + CoA(18:0) => CoA + DHAP(18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "Dihydroxyacetone-P + CoA(18:0)")
    expect_equal(df$reaction_product, "CoA + DHAP(18:0)") 
    
    reaction <- "RHEA:17659"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "Dihydroxyacetone-P + CoA(18:0) <= CoA + DHAP(18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "Dihydroxyacetone-P + CoA(18:0)")
    expect_equal(df$reaction_product, "CoA + DHAP(18:0)") 
    
    reaction <- "RHEA:17660"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "Dihydroxyacetone-P + CoA(18:0) <=> CoA + DHAP(18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "Dihydroxyacetone-P + CoA(18:0)")
    expect_equal(df$reaction_product, "CoA + DHAP(18:0)") 
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(18:0) + 2 NADPH + 2 H+ = FAO(18:0) + 2 NADP+ + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(18:0) + 2 NADPH + 2 H+")
    expect_equal(df$reaction_product, "FAO(18:0) + 2 NADP+ + CoA") 
    
    reaction <- "RHEA:52717"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(18:0) + 2 NADPH + 2 H+ => FAO(18:0) + 2 NADP+ + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(18:0) + 2 NADPH + 2 H+")
    expect_equal(df$reaction_product, "FAO(18:0) + 2 NADP+ + CoA") 
    
    reaction <- "RHEA:52718"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(18:0) + 2 NADPH + 2 H+ <= FAO(18:0) + 2 NADP+ + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(18:0) + 2 NADPH + 2 H+")
    expect_equal(df$reaction_product, "FAO(18:0) + 2 NADP+ + CoA") 
    
    reaction <- "RHEA:52719"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(18:0) + 2 NADPH + 2 H+ <=> FAO(18:0) + 2 NADP+ + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(18:0) + 2 NADPH + 2 H+")
    expect_equal(df$reaction_product, "FAO(18:0) + 2 NADP+ + CoA") 
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "Glycerol-3-P + CoA(18:0) = CoA + PA(18:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "Glycerol-3-P + CoA(18:0)")
    expect_equal(df$reaction_product, "CoA + PA(18:0/0:0)") 
    
    reaction <- "RHEA:15326"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "Glycerol-3-P + CoA(18:0) => CoA + PA(18:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "Glycerol-3-P + CoA(18:0)")
    expect_equal(df$reaction_product, "CoA + PA(18:0/0:0)") 
    
    reaction <- "RHEA:15327"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "Glycerol-3-P + CoA(18:0) <= CoA + PA(18:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "Glycerol-3-P + CoA(18:0)")
    expect_equal(df$reaction_product, "CoA + PA(18:0/0:0)") 
    
    reaction <- "RHEA:15328"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "Glycerol-3-P + CoA(18:0) <=> CoA + PA(18:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "Glycerol-3-P + CoA(18:0)")
    expect_equal(df$reaction_product, "CoA + PA(18:0/0:0)") 
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) = H+ + MG(18:0/0:0/0:0) + FA(16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + MG(18:0/0:0/0:0) + FA(16:0)") 
    
    reaction <- "RHEA:44713"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) => H+ + MG(18:0/0:0/0:0) + FA(16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + MG(18:0/0:0/0:0) + FA(16:0)") 
    
    reaction <- "RHEA:44714"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) <= H+ + MG(18:0/0:0/0:0) + FA(16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + MG(18:0/0:0/0:0) + FA(16:0)") 
    
    reaction <- "RHEA:44715"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) <=> H+ + MG(18:0/0:0/0:0) + FA(16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + MG(18:0/0:0/0:0) + FA(16:0)") 
    
    reaction <- "RHEA:35663"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) = H+ + MG(18:0/0:0/0:0) + FA(16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + MG(18:0/0:0/0:0) + FA(16:0)") 
    
    reaction <- "RHEA:35664"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) => H+ + MG(18:0/0:0/0:0) + FA(16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + MG(18:0/0:0/0:0) + FA(16:0)") 
    
    reaction <- "RHEA:35665"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA, 
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) <= H+ + MG(18:0/0:0/0:0) + FA(16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + MG(18:0/0:0/0:0) + FA(16:0)") 
    
    reaction <- "RHEA:35666"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) <=> H+ + MG(18:0/0:0/0:0) + FA(16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + MG(18:0/0:0/0:0) + FA(16:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) = H+ + MG(0:0/16:0/0:0) + FA(18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + MG(0:0/16:0/0:0) + FA(18:0)")
    
    reaction <- "RHEA:33276"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) => H+ + MG(0:0/16:0/0:0) + FA(18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + MG(0:0/16:0/0:0) + FA(18:0)")
    
    reaction <- "RHEA:33277"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) <= H+ + MG(0:0/16:0/0:0) + FA(18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + MG(0:0/16:0/0:0) + FA(18:0)")
    
    reaction <- "RHEA:33278"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) <=> H+ + MG(0:0/16:0/0:0) + FA(18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + MG(0:0/16:0/0:0) + FA(18:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "ATP + DG(18:0/16:0/0:0) = H+ + ADP + PA(18:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "ATP + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + ADP + PA(18:0/16:0)")
    
    reaction <- "RHEA:10273"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "ATP + DG(18:0/16:0/0:0) => H+ + ADP + PA(18:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "ATP + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + ADP + PA(18:0/16:0)")
    
    reaction <- "RHEA:10274"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "ATP + DG(18:0/16:0/0:0) <= H+ + ADP + PA(18:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "ATP + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + ADP + PA(18:0/16:0)")
    
    reaction <- "RHEA:10275"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "ATP + DG(18:0/16:0/0:0) <=> H+ + ADP + PA(18:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "ATP + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + ADP + PA(18:0/16:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CDP-Choline + DG(18:0/16:0/0:0) = H+ + CMP + PC(18:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CDP-Choline + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PC(18:0/16:0)")
    
    reaction <- "RHEA:32940"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CDP-Choline + DG(18:0/16:0/0:0) => H+ + CMP + PC(18:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CDP-Choline + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PC(18:0/16:0)")
    
    reaction <- "RHEA:32941"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CDP-Choline + DG(18:0/16:0/0:0) <= H+ + CMP + PC(18:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CDP-Choline + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PC(18:0/16:0)")
    
    reaction <- "RHEA:32942"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CDP-Choline + DG(18:0/16:0/0:0) <=> H+ + CMP + PC(18:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CDP-Choline + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PC(18:0/16:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CDP-Ethanolamine + DG(18:0/16:0/0:0) = H+ + CMP + PE(18:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CDP-Ethanolamine + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PE(18:0/16:0)")
    
    reaction <- "RHEA:32944"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CDP-Ethanolamine + DG(18:0/16:0/0:0) => H+ + CMP + PE(18:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CDP-Ethanolamine + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PE(18:0/16:0)")
    
    
    reaction <- "RHEA:32945"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CDP-Ethanolamine + DG(18:0/16:0/0:0) <= H+ + CMP + PE(18:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CDP-Ethanolamine + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PE(18:0/16:0)")
    
    reaction <- "RHEA:32946"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CDP-Ethanolamine + DG(18:0/16:0/0:0) <=> H+ + CMP + PE(18:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CDP-Ethanolamine + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PE(18:0/16:0)")

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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(14:0) + DG(18:0/16:0/0:0) = CoA + TG(18:0/16:0/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(14:0) + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "CoA + TG(18:0/16:0/14:0)")
    
    reaction <- "RHEA:10869"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(14:0) + DG(18:0/16:0/0:0) => CoA + TG(18:0/16:0/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(14:0) + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "CoA + TG(18:0/16:0/14:0)")
    
    reaction <- "RHEA:10870"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(14:0) + DG(18:0/16:0/0:0) <= CoA + TG(18:0/16:0/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(14:0) + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "CoA + TG(18:0/16:0/14:0)")
    
    reaction <- "RHEA:10871"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(14:0) + DG(18:0/16:0/0:0) <=> CoA + TG(18:0/16:0/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(14:0) + DG(18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "CoA + TG(18:0/16:0/14:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CDP-Choline + DG(O-18:0/16:0/0:0) = H+ + CMP + PC(O-18:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CDP-Choline + DG(O-18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PC(O-18:0/16:0)")
    
    reaction <- "RHEA:36180"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CDP-Choline + DG(O-18:0/16:0/0:0) => H+ + CMP + PC(O-18:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CDP-Choline + DG(O-18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PC(O-18:0/16:0)")
    
    reaction <- "RHEA:36181"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CDP-Choline + DG(O-18:0/16:0/0:0) <= H+ + CMP + PC(O-18:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CDP-Choline + DG(O-18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PC(O-18:0/16:0)")
    
    reaction <- "RHEA:36182"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CDP-Choline + DG(O-18:0/16:0/0:0) <=> H+ + CMP + PC(O-18:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CDP-Choline + DG(O-18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PC(O-18:0/16:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CDP-Ethanolamine + DG(O-18:0/16:0/0:0) = H+ + CMP + PE(O-18:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CDP-Ethanolamine + DG(O-18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PE(O-18:0/16:0)")
    
    reaction <- "RHEA:36188"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CDP-Ethanolamine + DG(O-18:0/16:0/0:0) => H+ + CMP + PE(O-18:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CDP-Ethanolamine + DG(O-18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PE(O-18:0/16:0)")
    
    reaction <- "RHEA:36189"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CDP-Ethanolamine + DG(O-18:0/16:0/0:0) <= H+ + CMP + PE(O-18:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CDP-Ethanolamine + DG(O-18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PE(O-18:0/16:0)")
    
    reaction <- "RHEA:36190"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CDP-Ethanolamine + DG(O-18:0/16:0/0:0) <=> H+ + CMP + PE(O-18:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CDP-Ethanolamine + DG(O-18:0/16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + CMP + PE(O-18:0/16:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H+ + NADH + O2 + Cer(d16:0(3OH,4OH)(15Me)/12:0) <=> 2 H2O + NAD+ + Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H+ + NADH + O2 + Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(df$reaction_product, "2 H2O + NAD+ + Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    
    ## dhcer_to_dhsm
    dhcer <- "Cer(16:1(3OH,4OH,15Me)/12:0)" ###################################################
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PC + Cer(16:1(3OH,4OH,15Me)/12:0) <=> DG + SM(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PC + Cer(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(df$reaction_product, "DG + SM(16:1(3OH,4OH,15Me)/12:0)")            
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + SM(16:1(3OH,4OH,15Me)/12:0) = Phosphocholine + H+ + Cer(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + SM(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(df$reaction_product, "Phosphocholine + H+ + Cer(16:1(3OH,4OH,15Me)/12:0)")            
    
    reaction <- "RHEA:45301"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + SM(16:1(3OH,4OH,15Me)/12:0) => Phosphocholine + H+ + Cer(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + SM(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(df$reaction_product, "Phosphocholine + H+ + Cer(16:1(3OH,4OH,15Me)/12:0)")  
    
    reaction <- "RHEA:45302"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + SM(16:1(3OH,4OH,15Me)/12:0) <= Phosphocholine + H+ + Cer(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + SM(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(df$reaction_product, "Phosphocholine + H+ + Cer(16:1(3OH,4OH,15Me)/12:0)")  
    
    reaction <- "RHEA:45303"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + SM(16:1(3OH,4OH,15Me)/12:0) <=> Phosphocholine + H+ + Cer(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + SM(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(df$reaction_product, "Phosphocholine + H+ + Cer(16:1(3OH,4OH,15Me)/12:0)")  
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "ATP + CoA + FA(18:0) = PPi + AMP + CoA(18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "ATP + CoA + FA(18:0)")
    expect_equal(df$reaction_product, "PPi + AMP + CoA(18:0)") 
    
    reaction <- "RHEA:15422"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "ATP + CoA + FA(18:0) => PPi + AMP + CoA(18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "ATP + CoA + FA(18:0)")
    expect_equal(df$reaction_product, "PPi + AMP + CoA(18:0)") 
    
    reaction <- "RHEA:15423"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "ATP + CoA + FA(18:0) <= PPi + AMP + CoA(18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "ATP + CoA + FA(18:0)")
    expect_equal(df$reaction_product, "PPi + AMP + CoA(18:0)") 
    
    reaction <- "RHEA:15424"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "ATP + CoA + FA(18:0) <=> PPi + AMP + CoA(18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "ATP + CoA + FA(18:0)")
    expect_equal(df$reaction_product, "PPi + AMP + CoA(18:0)") 
    
    reaction <- "RHEA:38883"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "ATP + CoA + FA(18:0) = PPi + AMP + CoA(18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "ATP + CoA + FA(18:0)")
    expect_equal(df$reaction_product, "PPi + AMP + CoA(18:0)") 
    
    reaction <- "RHEA:38884"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "ATP + CoA + FA(18:0) => PPi + AMP + CoA(18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "ATP + CoA + FA(18:0)")
    expect_equal(df$reaction_product, "PPi + AMP + CoA(18:0)") 
    
    reaction <- "RHEA:38885"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "ATP + CoA + FA(18:0) <= PPi + AMP + CoA(18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "ATP + CoA + FA(18:0)")
    expect_equal(df$reaction_product, "PPi + AMP + CoA(18:0)") 
    
    reaction <- "RHEA:38886"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "ATP + CoA + FA(18:0) <=> PPi + AMP + CoA(18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "ATP + CoA + FA(18:0)")
    expect_equal(df$reaction_product, "PPi + AMP + CoA(18:0)") 
    
    ## lcl_to_cl
    lcl <- "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])" ###################
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z)) = CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(df$reaction_product, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA")
    
    reaction <- "RHEA:35840"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z)) => CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(df$reaction_product, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA")
    
    reaction <- "RHEA:35841"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z)) <= CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(df$reaction_product, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA")
    
    reaction <- "RHEA:35842"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z)) <=> CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(df$reaction_product, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "NAPE(14:0/0:0/0:0) + H2O <=> GPNAE(0:0) + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "NAPE(14:0/0:0/0:0) + H2O")
    expect_equal(df$reaction_product, "GPNAE(0:0) + FA(14:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PA(18:0/0:0) + CoA(14:0) = CoA + PA(18:0/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PA(18:0/0:0) + CoA(14:0)")
    expect_equal(df$reaction_product, "CoA + PA(18:0/14:0)")
    
    reaction <- "RHEA:19710"
    df_reaction <- .add_products(substrates = df_substrates, 
                        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
                                                        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
                 "PA(18:0/0:0) + CoA(14:0) => CoA + PA(18:0/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PA(18:0/0:0) + CoA(14:0)")
    expect_equal(df$reaction_product, "CoA + PA(18:0/14:0)")
    
    reaction <- "RHEA:19711"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PA(18:0/0:0) + CoA(14:0) <= CoA + PA(18:0/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PA(18:0/0:0) + CoA(14:0)")
    expect_equal(df$reaction_product, "CoA + PA(18:0/14:0)")
    
    reaction <- "RHEA:19712"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PA(18:0/0:0) + CoA(14:0) <=> CoA + PA(18:0/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PA(18:0/0:0) + CoA(14:0)")
    expect_equal(df$reaction_product, "CoA + PA(18:0/14:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PA(O-18:0/0:0) + CoA(14:0) = PA(O-18:0/14:0) + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PA(O-18:0/0:0) + CoA(14:0)")
    expect_equal(df$reaction_product, "PA(O-18:0/14:0) + CoA")
    
    reaction <- "RHEA:36236"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PA(O-18:0/0:0) + CoA(14:0) => PA(O-18:0/14:0) + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PA(O-18:0/0:0) + CoA(14:0)")
    expect_equal(df$reaction_product, "PA(O-18:0/14:0) + CoA")
    
    reaction <- "RHEA:36237"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PA(O-18:0/0:0) + CoA(14:0) <= PA(O-18:0/14:0) + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PA(O-18:0/0:0) + CoA(14:0)")
    expect_equal(df$reaction_product, "PA(O-18:0/14:0) + CoA")
    
    reaction <- "RHEA:36238"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PA(O-18:0/0:0) + CoA(14:0) <=> PA(O-18:0/14:0) + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PA(O-18:0/0:0) + CoA(14:0)")
    expect_equal(df$reaction_product, "PA(O-18:0/14:0) + CoA")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(14:0/0:0) = Glycerophosphocholine + H+ + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(14:0/0:0)")
    expect_equal(df$reaction_product, "Glycerophosphocholine + H+ + FA(14:0)")
    
    reaction <- "RHEA:15178"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(14:0/0:0) => Glycerophosphocholine + H+ + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(14:0/0:0)")
    expect_equal(df$reaction_product, "Glycerophosphocholine + H+ + FA(14:0)")
    
    reaction <- "RHEA:15179"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(14:0/0:0) <= Glycerophosphocholine + H+ + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(14:0/0:0)")
    expect_equal(df$reaction_product, "Glycerophosphocholine + H+ + FA(14:0)")
    
    reaction <- "RHEA:15180"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(14:0/0:0) <=> Glycerophosphocholine + H+ + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(14:0/0:0)")
    expect_equal(df$reaction_product, "Glycerophosphocholine + H+ + FA(14:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(0:0/14:0) = Glycerophosphocholine + H+ + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(0:0/14:0)")
    expect_equal(df$reaction_product, "Glycerophosphocholine + H+ + FA(14:0)")
    
    reaction <- "RHEA:44697"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(0:0/14:0) => Glycerophosphocholine + H+ + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(0:0/14:0)")
    expect_equal(df$reaction_product, "Glycerophosphocholine + H+ + FA(14:0)")
    
    reaction <- "RHEA:44698"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(0:0/14:0) <= Glycerophosphocholine + H+ + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(0:0/14:0)")
    expect_equal(df$reaction_product, "Glycerophosphocholine + H+ + FA(14:0)")
    
    reaction <- "RHEA:44699"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(0:0/14:0) <=> Glycerophosphocholine + H+ + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(0:0/14:0)")
    expect_equal(df$reaction_product, "Glycerophosphocholine + H+ + FA(14:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PC(14:0/0:0) + CoA(18:0) = PC(14:0/18:0) + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PC(14:0/0:0) + CoA(18:0)")
    expect_equal(df$reaction_product, "PC(14:0/18:0) + CoA")
    
    reaction <- "RHEA:12938"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PC(14:0/0:0) + CoA(18:0) => PC(14:0/18:0) + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PC(14:0/0:0) + CoA(18:0)")
    expect_equal(df$reaction_product, "PC(14:0/18:0) + CoA")
    
    reaction <- "RHEA:12939"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PC(14:0/0:0) + CoA(18:0) <= PC(14:0/18:0) + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PC(14:0/0:0) + CoA(18:0)")
    expect_equal(df$reaction_product, "PC(14:0/18:0) + CoA")
    
    reaction <- "RHEA:12940"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PC(14:0/0:0) + CoA(18:0) <=> PC(14:0/18:0) + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PC(14:0/0:0) + CoA(18:0)")
    expect_equal(df$reaction_product, "PC(14:0/18:0) + CoA")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(14:0/0:0) = H+ + FA(14:0) + Glycerophosphoethanolamine")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(14:0/0:0)")
    expect_equal(df$reaction_product, "H+ + FA(14:0) + Glycerophosphoethanolamine")
    
    reaction <- "RHEA:32968"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(14:0/0:0) => H+ + FA(14:0) + Glycerophosphoethanolamine")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(14:0/0:0)")
    expect_equal(df$reaction_product, "H+ + FA(14:0) + Glycerophosphoethanolamine")
    
    reaction <- "RHEA:32969"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(14:0/0:0) <= H+ + FA(14:0) + Glycerophosphoethanolamine")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(14:0/0:0)")
    expect_equal(df$reaction_product, "H+ + FA(14:0) + Glycerophosphoethanolamine")
    
    reaction <- "RHEA:32970"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(14:0/0:0) <=> H+ + FA(14:0) + Glycerophosphoethanolamine")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(14:0/0:0)")
    expect_equal(df$reaction_product, "H+ + FA(14:0) + Glycerophosphoethanolamine")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(0:0/14:0) <=> H+ + FA(14:0) + Glycerophosphoethanolamine")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(0:0/14:0)")
    expect_equal(df$reaction_product, "H+ + FA(14:0) + Glycerophosphoethanolamine")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(18:0) + PE(14:0/0:0) = CoA + PE(14:0/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(18:0) + PE(14:0/0:0)")
    expect_equal(df$reaction_product, "CoA + PE(14:0/18:0)")
    
    reaction <- "RHEA:32996"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(18:0) + PE(14:0/0:0) => CoA + PE(14:0/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(18:0) + PE(14:0/0:0)")
    expect_equal(df$reaction_product, "CoA + PE(14:0/18:0)")
    
    reaction <- "RHEA:32997"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(18:0) + PE(14:0/0:0) <= CoA + PE(14:0/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(18:0) + PE(14:0/0:0)")
    expect_equal(df$reaction_product, "CoA + PE(14:0/18:0)")
    
    reaction <- "RHEA:32998"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(18:0) + PE(14:0/0:0) <=> CoA + PE(14:0/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(18:0) + PE(14:0/0:0)")
    expect_equal(df$reaction_product, "CoA + PE(14:0/18:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PI(16:0/0:0) + CoA(18:1(9Z)) = PI(16:0/18:1(9Z)) + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PI(16:0/0:0) + CoA(18:1(9Z))")
    expect_equal(df$reaction_product, "PI(16:0/18:1(9Z)) + CoA")
    
    reaction <- "RHEA:33196"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PI(16:0/0:0) + CoA(18:1(9Z)) => PI(16:0/18:1(9Z)) + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PI(16:0/0:0) + CoA(18:1(9Z))")
    expect_equal(df$reaction_product, "PI(16:0/18:1(9Z)) + CoA")
    
    reaction <- "RHEA:33197"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PI(16:0/0:0) + CoA(18:1(9Z)) <= PI(16:0/18:1(9Z)) + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PI(16:0/0:0) + CoA(18:1(9Z))")
    expect_equal(df$reaction_product, "PI(16:0/18:1(9Z)) + CoA")
    
    reaction <- "RHEA:33198"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PI(16:0/0:0) + CoA(18:1(9Z)) <=> PI(16:0/18:1(9Z)) + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PI(16:0/0:0) + CoA(18:1(9Z))")
    expect_equal(df$reaction_product, "PI(16:0/18:1(9Z)) + CoA")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(18:0) + PE(O-16:0/0:0) <=> CoA + PE(O-16:0/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(18:0) + PE(O-16:0/0:0)")
    expect_equal(df$reaction_product, "CoA + PE(O-16:0/18:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(18:0) + PE(P-16:0/0:0) = CoA + PE(P-16:0/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(18:0) + PE(P-16:0/0:0)")
    expect_equal(df$reaction_product, "CoA + PE(P-16:0/18:0)")
    
    reaction <- "RHEA:16246"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(18:0) + PE(P-16:0/0:0) => CoA + PE(P-16:0/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(18:0) + PE(P-16:0/0:0)")
    expect_equal(df$reaction_product, "CoA + PE(P-16:0/18:0)")
    
    reaction <- "RHEA:16247"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(18:0) + PE(P-16:0/0:0) <= CoA + PE(P-16:0/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(18:0) + PE(P-16:0/0:0)")
    expect_equal(df$reaction_product, "CoA + PE(P-16:0/18:0)")
    
    reaction <- "RHEA:16248"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(18:0) + PE(P-16:0/0:0) <=> CoA + PE(P-16:0/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(18:0) + PE(P-16:0/0:0)")
    expect_equal(df$reaction_product, "CoA + PE(P-16:0/18:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) = CoA + DG(14:0/16:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "MG(14:0/0:0/0:0) + CoA(16:0)")
    expect_equal(df$reaction_product, "CoA + DG(14:0/16:0/0:0)")
    
    reaction <- "RHEA:38464"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) => CoA + DG(14:0/16:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "MG(14:0/0:0/0:0) + CoA(16:0)")
    expect_equal(df$reaction_product, "CoA + DG(14:0/16:0/0:0)")
    
    reaction <- "RHEA:38465"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) <= CoA + DG(14:0/16:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "MG(14:0/0:0/0:0) + CoA(16:0)")
    expect_equal(df$reaction_product, "CoA + DG(14:0/16:0/0:0)")
    
    reaction <- "RHEA:38466"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) <=> CoA + DG(14:0/16:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "MG(14:0/0:0/0:0) + CoA(16:0)")
    expect_equal(df$reaction_product, "CoA + DG(14:0/16:0/0:0)")
    
    reaction <- "RHEA:39943"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) = CoA + DG(14:0/16:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "MG(14:0/0:0/0:0) + CoA(16:0)")
    expect_equal(df$reaction_product, "CoA + DG(14:0/16:0/0:0)")
    
    reaction <- "RHEA:39944"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) => CoA + DG(14:0/16:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "MG(14:0/0:0/0:0) + CoA(16:0)")
    expect_equal(df$reaction_product, "CoA + DG(14:0/16:0/0:0)")
    
    reaction <- "RHEA:39945"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) <= CoA + DG(14:0/16:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "MG(14:0/0:0/0:0) + CoA(16:0)")
    expect_equal(df$reaction_product, "CoA + DG(14:0/16:0/0:0)")
    
    reaction <- "RHEA:39946"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) <=> CoA + DG(14:0/16:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "MG(14:0/0:0/0:0) + CoA(16:0)")
    expect_equal(df$reaction_product, "CoA + DG(14:0/16:0/0:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(16:0) + MG(0:0/14:0/0:0) = CoA + DG(16:0/14:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(16:0) + MG(0:0/14:0/0:0)")
    expect_equal(df$reaction_product, "CoA + DG(16:0/14:0/0:0)")
    
    reaction <- "RHEA:32948"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(16:0) + MG(0:0/14:0/0:0) => CoA + DG(16:0/14:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(16:0) + MG(0:0/14:0/0:0)")
    expect_equal(df$reaction_product, "CoA + DG(16:0/14:0/0:0)")
    
    reaction <- "RHEA:32949"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(16:0) + MG(0:0/14:0/0:0) <= CoA + DG(16:0/14:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(16:0) + MG(0:0/14:0/0:0)")
    expect_equal(df$reaction_product, "CoA + DG(16:0/14:0/0:0)")
    
    reaction <- "RHEA:32950"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(16:0) + MG(0:0/14:0/0:0) <=> CoA + DG(16:0/14:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(16:0) + MG(0:0/14:0/0:0)")
    expect_equal(df$reaction_product, "CoA + DG(16:0/14:0/0:0)")
    
    reaction <- "RHEA:16741"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(16:0) + MG(0:0/14:0/0:0) = CoA + DG(16:0/14:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(16:0) + MG(0:0/14:0/0:0)")
    expect_equal(df$reaction_product, "CoA + DG(16:0/14:0/0:0)")
    
    reaction <- "RHEA:16742"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(16:0) + MG(0:0/14:0/0:0) => CoA + DG(16:0/14:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(16:0) + MG(0:0/14:0/0:0)")
    expect_equal(df$reaction_product, "CoA + DG(16:0/14:0/0:0)")
    
    reaction <- "RHEA:16743"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(16:0) + MG(0:0/14:0/0:0) <= CoA + DG(16:0/14:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(16:0) + MG(0:0/14:0/0:0)")
    expect_equal(df$reaction_product, "CoA + DG(16:0/14:0/0:0)")
    
    reaction <- "RHEA:16744"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(16:0) + MG(0:0/14:0/0:0) <=> CoA + DG(16:0/14:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(16:0) + MG(0:0/14:0/0:0)")
    expect_equal(df$reaction_product, "CoA + DG(16:0/14:0/0:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + MG(14:0/0:0/0:0) = Glycerol + H+ + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + MG(14:0/0:0/0:0)")
    expect_equal(df$reaction_product, "Glycerol + H+ + FA(14:0)")
    
    reaction <- "RHEA:34020"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + MG(14:0/0:0/0:0) => Glycerol + H+ + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + MG(14:0/0:0/0:0)")
    expect_equal(df$reaction_product, "Glycerol + H+ + FA(14:0)")
    
    reaction <- "RHEA:34021"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + MG(14:0/0:0/0:0) <= Glycerol + H+ + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + MG(14:0/0:0/0:0)")
    expect_equal(df$reaction_product, "Glycerol + H+ + FA(14:0)")
    
    reaction <- "RHEA:34022"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + MG(14:0/0:0/0:0) <=> Glycerol + H+ + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + MG(14:0/0:0/0:0)")
    expect_equal(df$reaction_product, "Glycerol + H+ + FA(14:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + MG(0:0/14:0/0:0) = Glycerol + H+ + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + MG(0:0/14:0/0:0)")
    expect_equal(df$reaction_product, "Glycerol + H+ + FA(14:0)")
    
    reaction <- "RHEA:32872"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + MG(0:0/14:0/0:0) => Glycerol + H+ + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + MG(0:0/14:0/0:0)")
    expect_equal(df$reaction_product, "Glycerol + H+ + FA(14:0)")
    
    reaction <- "RHEA:32873"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + MG(0:0/14:0/0:0) <= Glycerol + H+ + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + MG(0:0/14:0/0:0)")
    expect_equal(df$reaction_product, "Glycerol + H+ + FA(14:0)")
    
    reaction <- "RHEA:32874"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + MG(0:0/14:0/0:0) <=> Glycerol + H+ + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + MG(0:0/14:0/0:0)")
    expect_equal(df$reaction_product, "Glycerol + H+ + FA(14:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "ATP + MG(14:0/0:0/0:0) = H+ + ADP + PA(14:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "ATP + MG(14:0/0:0/0:0)")
    expect_equal(df$reaction_product, "H+ + ADP + PA(14:0/0:0)")
    
    reaction <- "RHEA:33748"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "ATP + MG(14:0/0:0/0:0) => H+ + ADP + PA(14:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "ATP + MG(14:0/0:0/0:0)")
    expect_equal(df$reaction_product, "H+ + ADP + PA(14:0/0:0)")
    
    reaction <- "RHEA:33749"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "ATP + MG(14:0/0:0/0:0) <= H+ + ADP + PA(14:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "ATP + MG(14:0/0:0/0:0)")
    expect_equal(df$reaction_product, "H+ + ADP + PA(14:0/0:0)")
    
    reaction <- "RHEA:33750"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "ATP + MG(14:0/0:0/0:0) <=> H+ + ADP + PA(14:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "ATP + MG(14:0/0:0/0:0)")
    expect_equal(df$reaction_product, "H+ + ADP + PA(14:0/0:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "MG(0:0/14:0/0:0) <=> MG(14:0/0:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "MG(0:0/14:0/0:0)")
    expect_equal(df$reaction_product, "MG(14:0/0:0/0:0)")

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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + NAE(18:0) <=> Ethanolamine + H+ + FA(18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + NAE(18:0)")
    expect_equal(df$reaction_product, "Ethanolamine + H+ + FA(18:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "NAPE(14:0/16:0/18:0) + H2O <=> NAPE(14:0/0:0/18:0) + FA(16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "NAPE(14:0/16:0/18:0) + H2O")
    expect_equal(df$reaction_product, "NAPE(14:0/0:0/18:0) + FA(16:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "NAPE(14:0/16:0/18:0) + H2O <=> NAE(18:0) + PA(14:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "NAPE(14:0/16:0/18:0) + H2O")
    expect_equal(df$reaction_product, "NAE(18:0) + PA(14:0/16:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "NAPE(14:0/16:0/18:0) + H2O <=> PNAE(18:0) + DG(14:0/16:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "NAPE(14:0/16:0/18:0) + H2O")
    expect_equal(df$reaction_product, "PNAE(18:0) + DG(14:0/16:0/0:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "NAPE(O-18:0/16:0/14:0) + H2O <=> NAE(14:0) + PA(O-18:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "NAPE(O-18:0/16:0/14:0) + H2O")
    expect_equal(df$reaction_product, "NAE(14:0) + PA(O-18:0/16:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CTP + PA(14:0/16:0) = PPi + CDP-DG(14:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CTP + PA(14:0/16:0)")
    expect_equal(df$reaction_product, "PPi + CDP-DG(14:0/16:0)")
    
    reaction <- "RHEA:16230"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CTP + PA(14:0/16:0) => PPi + CDP-DG(14:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CTP + PA(14:0/16:0)")
    expect_equal(df$reaction_product, "PPi + CDP-DG(14:0/16:0)")
    
    reaction <- "RHEA:16231"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CTP + PA(14:0/16:0) <= PPi + CDP-DG(14:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CTP + PA(14:0/16:0)")
    expect_equal(df$reaction_product, "PPi + CDP-DG(14:0/16:0)")
    
    reaction <- "RHEA:16232"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CTP + PA(14:0/16:0) <=> PPi + CDP-DG(14:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CTP + PA(14:0/16:0)")
    expect_equal(df$reaction_product, "PPi + CDP-DG(14:0/16:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PA(14:0/16:0) = Pi + DG(14:0/16:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PA(14:0/16:0)")
    expect_equal(df$reaction_product, "Pi + DG(14:0/16:0/0:0)")
    
    reaction <- "RHEA:27430"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PA(14:0/16:0) => Pi + DG(14:0/16:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PA(14:0/16:0)")
    expect_equal(df$reaction_product, "Pi + DG(14:0/16:0/0:0)")
    
    reaction <- "RHEA:27431"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PA(14:0/16:0) <= Pi + DG(14:0/16:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PA(14:0/16:0)")
    expect_equal(df$reaction_product, "Pi + DG(14:0/16:0/0:0)")
    
    reaction <- "RHEA:27432"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PA(14:0/16:0) <=> Pi + DG(14:0/16:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PA(14:0/16:0)")
    expect_equal(df$reaction_product, "Pi + DG(14:0/16:0/0:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PA(O-14:0/16:0) = Pi + DG(O-14:0/16:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PA(O-14:0/16:0)")
    expect_equal(df$reaction_product, "Pi + DG(O-14:0/16:0/0:0)")
    
    reaction <- "RHEA:36240"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PA(O-14:0/16:0) => Pi + DG(O-14:0/16:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PA(O-14:0/16:0)")
    expect_equal(df$reaction_product, "Pi + DG(O-14:0/16:0/0:0)")
    
    
    reaction <- "RHEA:36241"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PA(O-14:0/16:0) <= Pi + DG(O-14:0/16:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PA(O-14:0/16:0)")
    expect_equal(df$reaction_product, "Pi + DG(O-14:0/16:0/0:0)")
    
    reaction <- "RHEA:36242"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PA(O-14:0/16:0) <=> Pi + DG(O-14:0/16:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PA(O-14:0/16:0)")
    expect_equal(df$reaction_product, "Pi + DG(O-14:0/16:0/0:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(20:0/18:0) = Phosphocholine + DG(20:0/18:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "Phosphocholine + DG(20:0/18:0/0:0)")
    
    reaction <- "RHEA:10605"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(20:0/18:0) => Phosphocholine + DG(20:0/18:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "Phosphocholine + DG(20:0/18:0/0:0)")
    
    reaction <- "RHEA:10606"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(20:0/18:0) <= Phosphocholine + DG(20:0/18:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "Phosphocholine + DG(20:0/18:0/0:0)")
    
    reaction <- "RHEA:10607"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(20:0/18:0) <=> Phosphocholine + DG(20:0/18:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "Phosphocholine + DG(20:0/18:0/0:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(20:0/18:0) = PC(20:0/0:0) + FA(18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "PC(20:0/0:0) + FA(18:0)")
    
    reaction <- "RHEA:15802"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(20:0/18:0) => PC(20:0/0:0) + FA(18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "PC(20:0/0:0) + FA(18:0)")
    
    reaction <- "RHEA:15803"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(20:0/18:0) <= PC(20:0/0:0) + FA(18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "PC(20:0/0:0) + FA(18:0)")
    
    reaction <- "RHEA:15804"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(20:0/18:0) <=> PC(20:0/0:0) + FA(18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "PC(20:0/0:0) + FA(18:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(20:0/18:0) = PC(0:0/18:0) + FA(20:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "PC(0:0/18:0) + FA(20:0)")
    
    reaction <- "RHEA:18690"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(20:0/18:0) => PC(0:0/18:0) + FA(20:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "PC(0:0/18:0) + FA(20:0)")
    
    reaction <- "RHEA:18691"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(20:0/18:0) <= PC(0:0/18:0) + FA(20:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "PC(0:0/18:0) + FA(20:0)")
    
    reaction <- "RHEA:18692"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(20:0/18:0) <=> PC(0:0/18:0) + FA(20:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "PC(0:0/18:0) + FA(20:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(20:0/18:0) = Choline + PA(20:0/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "Choline + PA(20:0/18:0)")
    
    reaction <- "RHEA:14446"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(20:0/18:0) => Choline + PA(20:0/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "Choline + PA(20:0/18:0)")
    
    reaction <- "RHEA:14447"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(20:0/18:0) <= Choline + PA(20:0/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "Choline + PA(20:0/18:0)")
    
    reaction <- "RHEA:14448"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(20:0/18:0) <=> Choline + PA(20:0/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "Choline + PA(20:0/18:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "L-Serine + PC(20:0/18:0) = Choline + PS(20:0/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "L-Serine + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "Choline + PS(20:0/18:0)")
    
    reaction <- "RHEA:45089"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "L-Serine + PC(20:0/18:0) => Choline + PS(20:0/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "L-Serine + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "Choline + PS(20:0/18:0)")
    
    reaction <- "RHEA:45090"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "L-Serine + PC(20:0/18:0) <= Choline + PS(20:0/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "L-Serine + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "Choline + PS(20:0/18:0)")
    
    reaction <- "RHEA:45091"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "L-Serine + PC(20:0/18:0) <=> Choline + PS(20:0/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "L-Serine + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "Choline + PS(20:0/18:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(O-16:0/14:0) = H+ + PC(O-16:0/0:0) + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(O-16:0/14:0)")
    expect_equal(df$reaction_product, "H+ + PC(O-16:0/0:0) + FA(14:0)")
    
    reaction <- "RHEA:36232"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(O-16:0/14:0) => H+ + PC(O-16:0/0:0) + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(O-16:0/14:0)")
    expect_equal(df$reaction_product, "H+ + PC(O-16:0/0:0) + FA(14:0)")
    
    reaction <- "RHEA:36233"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(O-16:0/14:0) <= H+ + PC(O-16:0/0:0) + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(O-16:0/14:0)")
    expect_equal(df$reaction_product, "H+ + PC(O-16:0/0:0) + FA(14:0)")
    
    reaction <- "RHEA:36234"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(O-16:0/14:0) <=> H+ + PC(O-16:0/0:0) + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(O-16:0/14:0)")
    expect_equal(df$reaction_product, "H+ + PC(O-16:0/0:0) + FA(14:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(O-16:0/0:0) = PA(O-16:0/0:0) + H+ + Choline")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(O-16:0/0:0)")
    expect_equal(df$reaction_product, "PA(O-16:0/0:0) + H+ + Choline")
    
    reaction <- "RHEA:39928"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(O-16:0/0:0) => PA(O-16:0/0:0) + H+ + Choline")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(O-16:0/0:0)")
    expect_equal(df$reaction_product, "PA(O-16:0/0:0) + H+ + Choline")
    
    reaction <- "RHEA:39929"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(O-16:0/0:0) <= PA(O-16:0/0:0) + H+ + Choline")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(O-16:0/0:0)")
    expect_equal(df$reaction_product, "PA(O-16:0/0:0) + H+ + Choline")
    
    reaction <- "RHEA:39930"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(O-16:0/0:0) <=> PA(O-16:0/0:0) + H+ + Choline")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(O-16:0/0:0)")
    expect_equal(df$reaction_product, "PA(O-16:0/0:0) + H+ + Choline")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(O-16:0/0:0) = H+ + Phosphocholine + MG(O-16:0/0:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(O-16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + Phosphocholine + MG(O-16:0/0:0/0:0)")
    
    reaction <- "RHEA:36084"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(O-16:0/0:0) => H+ + Phosphocholine + MG(O-16:0/0:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(O-16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + Phosphocholine + MG(O-16:0/0:0/0:0)")
    
    reaction <- "RHEA:36085"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(O-16:0/0:0) <= H+ + Phosphocholine + MG(O-16:0/0:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(O-16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + Phosphocholine + MG(O-16:0/0:0/0:0)")
    
    reaction <- "RHEA:36086"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PC(O-16:0/0:0) <=> H+ + Phosphocholine + MG(O-16:0/0:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PC(O-16:0/0:0)")
    expect_equal(df$reaction_product, "H+ + Phosphocholine + MG(O-16:0/0:0/0:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PC(O-16:0/0:0) + CoA(18:0) = PC(O-16:0/18:0) + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PC(O-16:0/0:0) + CoA(18:0)")
    expect_equal(df$reaction_product, "PC(O-16:0/18:0) + CoA")
    
    reaction <- "RHEA:23993"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PC(O-16:0/0:0) + CoA(18:0) => PC(O-16:0/18:0) + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PC(O-16:0/0:0) + CoA(18:0)")
    expect_equal(df$reaction_product, "PC(O-16:0/18:0) + CoA")
    
    reaction <- "RHEA:23994"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PC(O-16:0/0:0) + CoA(18:0) <= PC(O-16:0/18:0) + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PC(O-16:0/0:0) + CoA(18:0)")
    expect_equal(df$reaction_product, "PC(O-16:0/18:0) + CoA")
    
    reaction <- "RHEA:23995"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PC(O-16:0/0:0) + CoA(18:0) <=> PC(O-16:0/18:0) + CoA")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PC(O-16:0/0:0) + CoA(18:0)")
    expect_equal(df$reaction_product, "PC(O-16:0/18:0) + CoA")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(14:0/16:0) <=> P-Ethanolamine + DG(14:0/16:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(14:0/16:0)")
    expect_equal(df$reaction_product, "P-Ethanolamine + DG(14:0/16:0/0:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(14:0/16:0) = H+ + FA(16:0) + PE(14:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(14:0/16:0)")
    expect_equal(df$reaction_product, "H+ + FA(16:0) + PE(14:0/0:0)")
    
    reaction <- "RHEA:44605"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(14:0/16:0) => H+ + FA(16:0) + PE(14:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(14:0/16:0)")
    expect_equal(df$reaction_product, "H+ + FA(16:0) + PE(14:0/0:0)")
    
    reaction <- "RHEA:44606"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(14:0/16:0) <= H+ + FA(16:0) + PE(14:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(14:0/16:0)")
    expect_equal(df$reaction_product, "H+ + FA(16:0) + PE(14:0/0:0)")
    
    reaction <- "RHEA:44607"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(14:0/16:0) <=> H+ + FA(16:0) + PE(14:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(14:0/16:0)")
    expect_equal(df$reaction_product, "H+ + FA(16:0) + PE(14:0/0:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(14:0/16:0) = H+ + FA(16:0) + PE(0:0/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(14:0/16:0)")
    expect_equal(df$reaction_product, "H+ + FA(16:0) + PE(0:0/14:0)")
    
    reaction <- "RHEA:44409"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(14:0/16:0) => H+ + FA(16:0) + PE(0:0/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(14:0/16:0)")
    expect_equal(df$reaction_product, "H+ + FA(16:0) + PE(0:0/14:0)")
    
    reaction <- "RHEA:44410"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(14:0/16:0) <= H+ + FA(16:0) + PE(0:0/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(14:0/16:0)")
    expect_equal(df$reaction_product, "H+ + FA(16:0) + PE(0:0/14:0)")
    
    reaction <- "RHEA:44411"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(14:0/16:0) <=> H+ + FA(16:0) + PE(0:0/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(14:0/16:0)")
    expect_equal(df$reaction_product, "H+ + FA(16:0) + PE(0:0/14:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PE(14:0/16:0) + PC(18:0/20:0) <=> NAPE(14:0/16:0/18:0) + PC(0:0/20:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PE(14:0/16:0) + PC(18:0/20:0)")
    expect_equal(df$reaction_product, "NAPE(14:0/16:0/18:0) + PC(0:0/20:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PE(14:0/16:0) + PC(18:0/20:0) <=> NAPE(14:0/16:0/20:0) + PC(18:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PE(14:0/16:0) + PC(18:0/20:0)")
    expect_equal(df$reaction_product, "NAPE(14:0/16:0/20:0) + PC(18:0/0:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(14:0/16:0) <=> Ethanolamine + H+ + PA(14:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(14:0/16:0)")
    expect_equal(df$reaction_product, "Ethanolamine + H+ + PA(14:0/16:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "L-Serine + PE(14:0/16:0) = Ethanolamine + PS(14:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "L-Serine + PE(14:0/16:0)")
    expect_equal(df$reaction_product, "Ethanolamine + PS(14:0/16:0)")
    
    reaction <- "RHEA:27607"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "L-Serine + PE(14:0/16:0) => Ethanolamine + PS(14:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "L-Serine + PE(14:0/16:0)")
    expect_equal(df$reaction_product, "Ethanolamine + PS(14:0/16:0)")
    
    reaction <- "RHEA:27608"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "L-Serine + PE(14:0/16:0) <= Ethanolamine + PS(14:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "L-Serine + PE(14:0/16:0)")
    expect_equal(df$reaction_product, "Ethanolamine + PS(14:0/16:0)")
    
    reaction <- "RHEA:27609"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "L-Serine + PE(14:0/16:0) <=> Ethanolamine + PS(14:0/16:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "L-Serine + PE(14:0/16:0)")
    expect_equal(df$reaction_product, "Ethanolamine + PS(14:0/16:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(O-16:0/14:0) <=> H+ + PE(O-16:0/0:0) + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(O-16:0/14:0)")
    expect_equal(df$reaction_product, "H+ + PE(O-16:0/0:0) + FA(14:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PE(O-16:0/14:0) + PC(20:0/18:0) <=> NAPE(O-16:0/14:0/20:0) + PC(0:0/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PE(O-16:0/14:0) + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "NAPE(O-16:0/14:0/20:0) + PC(0:0/18:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PE(O-16:0/14:0) + PC(20:0/18:0) <=> NAPE(O-16:0/14:0/18:0) + PC(20:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PE(O-16:0/14:0) + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "NAPE(O-16:0/14:0/18:0) + PC(20:0/0:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2 = PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2")
    expect_equal(df$reaction_product, "PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    
    reaction <- "RHEA:22957"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2 => PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2")
    expect_equal(df$reaction_product, "PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    
    reaction <- "RHEA:22958"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2 <= PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2")
    expect_equal(df$reaction_product, "PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    
    reaction <- "RHEA:22959"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2 <=> PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2")
    expect_equal(df$reaction_product, "PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(P-16:0/14:0) = H+ + PE(P-16:0/0:0) + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(P-16:0/14:0)")
    expect_equal(df$reaction_product, "H+ + PE(P-16:0/0:0) + FA(14:0)")
    
    reaction <- "RHEA:36196"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(P-16:0/14:0) => H+ + PE(P-16:0/0:0) + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(P-16:0/14:0)")
    expect_equal(df$reaction_product, "H+ + PE(P-16:0/0:0) + FA(14:0)")
    
    reaction <- "RHEA:36197"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(P-16:0/14:0) <= H+ + PE(P-16:0/0:0) + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(P-16:0/14:0)")
    expect_equal(df$reaction_product, "H+ + PE(P-16:0/0:0) + FA(14:0)")
    
    reaction <- "RHEA:36198"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(P-16:0/14:0) <=> H+ + PE(P-16:0/0:0) + FA(14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(P-16:0/14:0)")
    expect_equal(df$reaction_product, "H+ + PE(P-16:0/0:0) + FA(14:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(P-16:0/0:0) = FAL(16:0) + Glycerophosphoethanolamine")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(P-16:0/0:0)")
    expect_equal(df$reaction_product, "FAL(16:0) + Glycerophosphoethanolamine")
    
    reaction <- "RHEA:16906"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(P-16:0/0:0) => FAL(16:0) + Glycerophosphoethanolamine")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(P-16:0/0:0)")
    expect_equal(df$reaction_product, "FAL(16:0) + Glycerophosphoethanolamine")
    
    reaction <- "RHEA:16907"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(P-16:0/0:0) <= FAL(16:0) + Glycerophosphoethanolamine")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(P-16:0/0:0)")
    expect_equal(df$reaction_product, "FAL(16:0) + Glycerophosphoethanolamine")
    
    reaction <- "RHEA:16908"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(P-16:0/0:0) <=> FAL(16:0) + Glycerophosphoethanolamine")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(P-16:0/0:0)")
    expect_equal(df$reaction_product, "FAL(16:0) + Glycerophosphoethanolamine")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(P-16:0/0:0) = PA(P-16:0/0:0) + H+ + Ethanolamine")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(P-16:0/0:0)")
    expect_equal(df$reaction_product, "PA(P-16:0/0:0) + H+ + Ethanolamine")
    
    reaction <- "RHEA:36204"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(P-16:0/0:0) => PA(P-16:0/0:0) + H+ + Ethanolamine")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(P-16:0/0:0)")
    expect_equal(df$reaction_product, "PA(P-16:0/0:0) + H+ + Ethanolamine")
    
    reaction <- "RHEA:36205"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(P-16:0/0:0) <= PA(P-16:0/0:0) + H+ + Ethanolamine")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(P-16:0/0:0)")
    expect_equal(df$reaction_product, "PA(P-16:0/0:0) + H+ + Ethanolamine")
    
    reaction <- "RHEA:36206"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(P-16:0/0:0) <=> PA(P-16:0/0:0) + H+ + Ethanolamine")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(P-16:0/0:0)")
    expect_equal(df$reaction_product, "PA(P-16:0/0:0) + H+ + Ethanolamine")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(P-16:0/0:0) = Phosphoethanolamine + H+ + MG(P-16:0/0:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(P-16:0/0:0)")
    expect_equal(df$reaction_product, "Phosphoethanolamine + H+ + MG(P-16:0/0:0/0:0)")
    
    reaction <- "RHEA:36200"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(P-16:0/0:0) => Phosphoethanolamine + H+ + MG(P-16:0/0:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(P-16:0/0:0)")
    expect_equal(df$reaction_product, "Phosphoethanolamine + H+ + MG(P-16:0/0:0/0:0)")
    
    reaction <- "RHEA:36201"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(P-16:0/0:0) <= Phosphoethanolamine + H+ + MG(P-16:0/0:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(P-16:0/0:0)")
    expect_equal(df$reaction_product, "Phosphoethanolamine + H+ + MG(P-16:0/0:0/0:0)")
    
    
    reaction <- "RHEA:36202"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PE(P-16:0/0:0) <=> Phosphoethanolamine + H+ + MG(P-16:0/0:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PE(P-16:0/0:0)")
    expect_equal(df$reaction_product, "Phosphoethanolamine + H+ + MG(P-16:0/0:0/0:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PE(P-16:0/14:0) + PC(20:0/18:0) <=> NAPE(P-16:0/14:0/20:0) + PC(0:0/18:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PE(P-16:0/14:0) + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "NAPE(P-16:0/14:0/20:0) + PC(0:0/18:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "PE(P-16:0/14:0) + PC(20:0/18:0) <=> NAPE(P-16:0/14:0/18:0) + PC(20:0/0:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "PE(P-16:0/14:0) + PC(20:0/18:0)")
    expect_equal(df$reaction_product, "NAPE(P-16:0/14:0/18:0) + PC(20:0/0:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0) + PG(18:4(6Z,9Z,12Z,15Z)/14:0) = H+ + CMP + CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0) + PG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(df$reaction_product, "H+ + CMP + CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    
    reaction <- "RHEA:32932"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0) + PG(18:4(6Z,9Z,12Z,15Z)/14:0) => H+ + CMP + CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0) + PG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(df$reaction_product, "H+ + CMP + CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    
    reaction <- "RHEA:32933"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0) + PG(18:4(6Z,9Z,12Z,15Z)/14:0) <= H+ + CMP + CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0) + PG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(df$reaction_product, "H+ + CMP + CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    
    reaction <- "RHEA:32934"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0) + PG(18:4(6Z,9Z,12Z,15Z)/14:0) <=> H+ + CMP + CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0) + PG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(df$reaction_product, "H+ + CMP + CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    
    ## pgp_to_pg
    pgp <- "PGP(16:0/14:0)"
    substrates <- list(PGP = pgp)
    df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
        substrates = substrates)
    .names <- c("PGP", "PG")
    .values <- c("PGP(16:0/14:0)", "PG(16:0/14:0)")
    
    reaction <- "RHEA:33751"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PGP(16:0/14:0) = Pi + PG(16:0/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PGP(16:0/14:0)")
    expect_equal(df$reaction_product, "Pi + PG(16:0/14:0)")
    
    reaction <- "RHEA:33752"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PGP(16:0/14:0) => Pi + PG(16:0/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PGP(16:0/14:0)")
    expect_equal(df$reaction_product, "Pi + PG(16:0/14:0)")
    
    reaction <- "RHEA:33753"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PGP(16:0/14:0) <= Pi + PG(16:0/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PGP(16:0/14:0)")
    expect_equal(df$reaction_product, "Pi + PG(16:0/14:0)")
    
    reaction <- "RHEA:33754"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PGP(16:0/14:0) <=> Pi + PG(16:0/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PGP(16:0/14:0)")
    expect_equal(df$reaction_product, "Pi + PG(16:0/14:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PI(16:0/18:1(9Z)) = myo-Inositol-1-P + H+ + DG(16:0/18:1(9Z))")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PI(16:0/18:1(9Z))")
    expect_equal(df$reaction_product, "myo-Inositol-1-P + H+ + DG(16:0/18:1(9Z))")
    
    reaction <- "RHEA:43485"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PI(16:0/18:1(9Z)) => myo-Inositol-1-P + H+ + DG(16:0/18:1(9Z))")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PI(16:0/18:1(9Z))")
    expect_equal(df$reaction_product, "myo-Inositol-1-P + H+ + DG(16:0/18:1(9Z))")
    
    
    reaction <- "RHEA:43486"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PI(16:0/18:1(9Z)) <= myo-Inositol-1-P + H+ + DG(16:0/18:1(9Z))")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PI(16:0/18:1(9Z))")
    expect_equal(df$reaction_product, "myo-Inositol-1-P + H+ + DG(16:0/18:1(9Z))")
    
    reaction <- "RHEA:43487"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PI(16:0/18:1(9Z)) <=> myo-Inositol-1-P + H+ + DG(16:0/18:1(9Z))")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PI(16:0/18:1(9Z))")
    expect_equal(df$reaction_product, "myo-Inositol-1-P + H+ + DG(16:0/18:1(9Z))")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PI(16:0/18:1(9Z)) = PI(16:0/0:0) + H+ + FA(18:1(9Z))")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PI(16:0/18:1(9Z))")
    expect_equal(df$reaction_product, "PI(16:0/0:0) + H+ + FA(18:1(9Z))")
    
    reaction <- "RHEA:18002"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PI(16:0/18:1(9Z)) => PI(16:0/0:0) + H+ + FA(18:1(9Z))")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PI(16:0/18:1(9Z))")
    expect_equal(df$reaction_product, "PI(16:0/0:0) + H+ + FA(18:1(9Z))")
    
    reaction <- "RHEA:18003"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PI(16:0/18:1(9Z)) <= PI(16:0/0:0) + H+ + FA(18:1(9Z))")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PI(16:0/18:1(9Z))")
    expect_equal(df$reaction_product, "PI(16:0/0:0) + H+ + FA(18:1(9Z))")
    
    reaction <- "RHEA:18004"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + PI(16:0/18:1(9Z)) <=> PI(16:0/0:0) + H+ + FA(18:1(9Z))")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + PI(16:0/18:1(9Z))")
    expect_equal(df$reaction_product, "PI(16:0/0:0) + H+ + FA(18:1(9Z))")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H+ + PS(14:0/14:0) = CO2 + PE(14:0/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H+ + PS(14:0/14:0)")
    expect_equal(df$reaction_product, "CO2 + PE(14:0/14:0)")
    
    reaction <- "RHEA:20829"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H+ + PS(14:0/14:0) => CO2 + PE(14:0/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H+ + PS(14:0/14:0)")
    expect_equal(df$reaction_product, "CO2 + PE(14:0/14:0)")
    
    reaction <- "RHEA:20830"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H+ + PS(14:0/14:0) <= CO2 + PE(14:0/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H+ + PS(14:0/14:0)")
    expect_equal(df$reaction_product, "CO2 + PE(14:0/14:0)")
    
    reaction <- "RHEA:20831"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H+ + PS(14:0/14:0) <=> CO2 + PE(14:0/14:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H+ + PS(14:0/14:0)")
    expect_equal(df$reaction_product, "CO2 + PE(14:0/14:0)")

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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "H2O + SM(16:0(3OH,4OH,15Me)/12:0) <=> Phosphocholine + H+ + Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "H2O + SM(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(df$reaction_product, "Phosphocholine + H+ + Cer(16:0(3OH,4OH,15Me)/12:0)")
    
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
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me)) = CoA + Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me))")
    expect_equal(df$reaction_product, "CoA + Cer(16:0(3OH,4OH,15Me)/12:0)")
    
    reaction <- "RHEA:53425"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me)) => CoA + Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me))")
    expect_equal(df$reaction_product, "CoA + Cer(16:0(3OH,4OH,15Me)/12:0)")
    
    reaction <- "RHEA:53426"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me)) <= CoA + Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me))")
    expect_equal(df$reaction_product, "CoA + Cer(16:0(3OH,4OH,15Me)/12:0)")
    
    reaction <- "RHEA:53427"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, "")
    expect_equal(df$reaction_formula, 
        "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me)) <=> CoA + Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(df$reaction_isReversible, "")
    expect_equal(df$reaction_geneAssociation, "")
    expect_equal(df$reaction_pathway, "")
    expect_equal(df$reaction_substrate, "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me))")
    expect_equal(df$reaction_product, "CoA + Cer(16:0(3OH,4OH,15Me)/12:0)")
    
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
    expect_equal(df$reaction_name, c("", ""))
    expect_equal(df$reaction_formula, 
        c("H2O + TG(18:0/16:0/14:0) = H+ + FA(18:0) + DG(14:0/16:0/0:0)", 
            "H2O + TG(18:0/16:0/14:0) = H+ + FA(14:0) + DG(18:0/16:0/0:0)"))
    expect_equal(df$reaction_isReversible, c("", ""))
    expect_equal(df$reaction_geneAssociation, c("", ""))
    expect_equal(df$reaction_pathway, c("", ""))
    expect_equal(df$reaction_substrate, 
        c("H2O + TG(18:0/16:0/14:0)", "H2O + TG(18:0/16:0/14:0)"))
    expect_equal(df$reaction_product, 
        c("H+ + FA(18:0) + DG(14:0/16:0/0:0)", "H+ + FA(14:0) + DG(18:0/16:0/0:0)"))
    
    reaction <- "RHEA:33272"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction) 
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, c("", ""))
    expect_equal(df$reaction_formula, 
        c("H2O + TG(18:0/16:0/14:0) => H+ + FA(18:0) + DG(14:0/16:0/0:0)", 
            "H2O + TG(18:0/16:0/14:0) => H+ + FA(14:0) + DG(18:0/16:0/0:0)"))
    expect_equal(df$reaction_isReversible, c("", ""))
    expect_equal(df$reaction_geneAssociation, c("", ""))
    expect_equal(df$reaction_pathway, c("", ""))
    expect_equal(df$reaction_substrate, 
        c("H2O + TG(18:0/16:0/14:0)", "H2O + TG(18:0/16:0/14:0)"))
    expect_equal(df$reaction_product, 
        c("H+ + FA(18:0) + DG(14:0/16:0/0:0)", "H+ + FA(14:0) + DG(18:0/16:0/0:0)"))
    
    reaction <- "RHEA:33273"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, c("", ""))
    expect_equal(df$reaction_formula, 
        c("H2O + TG(18:0/16:0/14:0) <= H+ + FA(18:0) + DG(14:0/16:0/0:0)", 
            "H2O + TG(18:0/16:0/14:0) <= H+ + FA(14:0) + DG(18:0/16:0/0:0)"))
    expect_equal(df$reaction_isReversible, c("", ""))
    expect_equal(df$reaction_geneAssociation, c("", ""))
    expect_equal(df$reaction_pathway, c("", ""))
    expect_equal(df$reaction_substrate, 
        c("H2O + TG(18:0/16:0/14:0)", "H2O + TG(18:0/16:0/14:0)"))
    expect_equal(df$reaction_product, 
        c("H+ + FA(18:0) + DG(14:0/16:0/0:0)", "H+ + FA(14:0) + DG(18:0/16:0/0:0)"))
    
    reaction <- "RHEA:33274"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, c("", ""))
    expect_equal(df$reaction_formula, 
        c("H2O + TG(18:0/16:0/14:0) <=> H+ + FA(18:0) + DG(14:0/16:0/0:0)", 
            "H2O + TG(18:0/16:0/14:0) <=> H+ + FA(14:0) + DG(18:0/16:0/0:0)"))
    expect_equal(df$reaction_isReversible, c("", ""))
    expect_equal(df$reaction_geneAssociation, c("", ""))
    expect_equal(df$reaction_pathway, c("", ""))
    expect_equal(df$reaction_substrate, 
        c("H2O + TG(18:0/16:0/14:0)", "H2O + TG(18:0/16:0/14:0)"))
    expect_equal(df$reaction_product, 
        c("H+ + FA(18:0) + DG(14:0/16:0/0:0)", "H+ + FA(14:0) + DG(18:0/16:0/0:0)"))
     
    reaction <- "RHEA:44864"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, c("", ""))
    expect_equal(df$reaction_formula, 
        c("H2O + TG(18:0/16:0/14:0) = H+ + FA(18:0) + DG(14:0/16:0/0:0)", 
            "H2O + TG(18:0/16:0/14:0) = H+ + FA(14:0) + DG(18:0/16:0/0:0)"))
    expect_equal(df$reaction_isReversible, c("", ""))
    expect_equal(df$reaction_geneAssociation, c("", ""))
    expect_equal(df$reaction_pathway, c("", ""))
    expect_equal(df$reaction_substrate, 
        c("H2O + TG(18:0/16:0/14:0)", "H2O + TG(18:0/16:0/14:0)"))
    expect_equal(df$reaction_product, 
        c("H+ + FA(18:0) + DG(14:0/16:0/0:0)", "H+ + FA(14:0) + DG(18:0/16:0/0:0)"))
    
    reaction <- "RHEA:44865"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, c("", ""))
    expect_equal(df$reaction_formula, 
        c("H2O + TG(18:0/16:0/14:0) => H+ + FA(18:0) + DG(14:0/16:0/0:0)", 
            "H2O + TG(18:0/16:0/14:0) => H+ + FA(14:0) + DG(18:0/16:0/0:0)"))
    expect_equal(df$reaction_isReversible, c("", ""))
    expect_equal(df$reaction_geneAssociation, c("", ""))
    expect_equal(df$reaction_pathway, c("", ""))
    expect_equal(df$reaction_substrate, 
        c("H2O + TG(18:0/16:0/14:0)", "H2O + TG(18:0/16:0/14:0)"))
    expect_equal(df$reaction_product, 
        c("H+ + FA(18:0) + DG(14:0/16:0/0:0)", "H+ + FA(14:0) + DG(18:0/16:0/0:0)"))
    
    reaction <- "RHEA:44866"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, c("", ""))
    expect_equal(df$reaction_formula, 
        c("H2O + TG(18:0/16:0/14:0) <= H+ + FA(18:0) + DG(14:0/16:0/0:0)", 
            "H2O + TG(18:0/16:0/14:0) <= H+ + FA(14:0) + DG(18:0/16:0/0:0)"))
    expect_equal(df$reaction_isReversible, c("", ""))
    expect_equal(df$reaction_geneAssociation, c("", ""))
    expect_equal(df$reaction_pathway, c("", ""))
    expect_equal(df$reaction_substrate, 
        c("H2O + TG(18:0/16:0/14:0)", "H2O + TG(18:0/16:0/14:0)"))
    expect_equal(df$reaction_product, 
        c("H+ + FA(18:0) + DG(14:0/16:0/0:0)", "H+ + FA(14:0) + DG(18:0/16:0/0:0)"))
    
    reaction <- "RHEA:44867"
    df_reaction <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    template <- LipidNetworkPredictR:::.create_template(template = NA,
        reaction = reaction)
    df <- LipidNetworkPredictR:::.create_df_with_template(
        df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(df$reaction_name, c("", ""))
    expect_equal(df$reaction_formula, 
        c("H2O + TG(18:0/16:0/14:0) <=> H+ + FA(18:0) + DG(14:0/16:0/0:0)", 
            "H2O + TG(18:0/16:0/14:0) <=> H+ + FA(14:0) + DG(18:0/16:0/0:0)"))
    expect_equal(df$reaction_isReversible, c("", ""))
    expect_equal(df$reaction_geneAssociation, c("", ""))
    expect_equal(df$reaction_pathway, c("", ""))
    expect_equal(df$reaction_substrate, 
        c("H2O + TG(18:0/16:0/14:0)", "H2O + TG(18:0/16:0/14:0)"))
    expect_equal(df$reaction_product, 
        c("H+ + FA(18:0) + DG(14:0/16:0/0:0)", "H+ + FA(14:0) + DG(18:0/16:0/0:0)"))
     
})
