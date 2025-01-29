## function .create_reaction
test_that(".create_list_with_template works", {
    
    ## dg_to_pa
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    
    reaction <- "RHEA:10272"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + ATP = PA(18:0/16:0) + ADP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + ATP")
    expect_equal(l$reaction_product, "PA(18:0/16:0) + ADP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:30616 = CHEBI:58608 + CHEBI:456216 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:30616")
    expect_equal(l$reaction_product_chebi, "CHEBI:58608 + CHEBI:456216 + CHEBI:15378")
    
    reaction <- "RHEA:10273"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + ATP => PA(18:0/16:0) + ADP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + ATP")
    expect_equal(l$reaction_product, "PA(18:0/16:0) + ADP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:30616 => CHEBI:58608 + CHEBI:456216 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:30616")
    expect_equal(l$reaction_product_chebi, "CHEBI:58608 + CHEBI:456216 + CHEBI:15378")
    
    reaction <- "RHEA:10274"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + ATP <= PA(18:0/16:0) + ADP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + ATP")
    expect_equal(l$reaction_product, "PA(18:0/16:0) + ADP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:30616 <= CHEBI:58608 + CHEBI:456216 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:30616")
    expect_equal(l$reaction_product_chebi, "CHEBI:58608 + CHEBI:456216 + CHEBI:15378")
    
    reaction <- "RHEA:10275"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + ATP <=> PA(18:0/16:0) + ADP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + ATP")
    expect_equal(l$reaction_product, "PA(18:0/16:0) + ADP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:30616 <=> CHEBI:58608 + CHEBI:456216 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:30616")
    expect_equal(l$reaction_product_chebi, "CHEBI:58608 + CHEBI:456216 + CHEBI:15378")
    
    
    ## pc_to_dg
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    
    reaction <- "RHEA:10604"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(20:0/18:0) + H2O = DG(20:0/18:0/0:0) + H+ + Phosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(20:0/18:0) + H2O")
    expect_equal(l$reaction_product, "DG(20:0/18:0/0:0) + H+ + Phosphocholine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 = CHEBI:17815 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:15378 + CHEBI:295975")
    
    reaction <- "RHEA:10605"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(20:0/18:0) + H2O => DG(20:0/18:0/0:0) + H+ + Phosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(20:0/18:0) + H2O")
    expect_equal(l$reaction_product, "DG(20:0/18:0/0:0) + H+ + Phosphocholine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 => CHEBI:17815 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:15378 + CHEBI:295975")
    
    reaction <- "RHEA:10606"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(20:0/18:0) + H2O <= DG(20:0/18:0/0:0) + H+ + Phosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(20:0/18:0) + H2O")
    expect_equal(l$reaction_product, "DG(20:0/18:0/0:0) + H+ + Phosphocholine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <= CHEBI:17815 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:15378 + CHEBI:295975")
    
    reaction <- "RHEA:10607"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(20:0/18:0) + H2O <=> DG(20:0/18:0/0:0) + H+ + Phosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(20:0/18:0) + H2O")
    expect_equal(l$reaction_product, "DG(20:0/18:0/0:0) + H+ + Phosphocholine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <=> CHEBI:17815 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:15378 + CHEBI:295975")
    
    ## dg_to_tg
    dg <- "DG(18:0/16:0/0:0)"
    acylcoa <- "CoA(14:0)"
    substrates <- list(DG = dg, AcylCoA = acylcoa)
    
    reaction <- "RHEA:10868"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + CoA(14:0) = TG(18:0/16:0/14:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + CoA(14:0)")
    expect_equal(l$reaction_product, "TG(18:0/16:0/14:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58342 = CHEBI:64615 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:64615 + CHEBI:57287")
    
    reaction <- "RHEA:10869"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + CoA(14:0) => TG(18:0/16:0/14:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + CoA(14:0)")
    expect_equal(l$reaction_product, "TG(18:0/16:0/14:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58342 => CHEBI:64615 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:64615 + CHEBI:57287")
    
    reaction <- "RHEA:10870"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + CoA(14:0) <= TG(18:0/16:0/14:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + CoA(14:0)")
    expect_equal(l$reaction_product, "TG(18:0/16:0/14:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58342 <= CHEBI:64615 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:64615 + CHEBI:57287")
    
    reaction <- "RHEA:10871"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + CoA(14:0) <=> TG(18:0/16:0/14:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + CoA(14:0)")
    expect_equal(l$reaction_product, "TG(18:0/16:0/14:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58342 <=> CHEBI:64615 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:64615 + CHEBI:57287")
    
    ## cdpdg_to_pi
    cdpdg <- "CDP-DG(12:0(11Me)/14:0)"
    substrates <- list(CDPDG = cdpdg)
    
    reaction <- "RHEA:11580"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "CDP-DG(12:0(11Me)/14:0) + myo-Inositol = PI(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CDP-DG(12:0(11Me)/14:0) + myo-Inositol")
    expect_equal(l$reaction_product, "PI(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58332 + CHEBI:17268 = CHEBI:57880 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58332 + CHEBI:17268")
    expect_equal(l$reaction_product_chebi, "CHEBI:57880 + CHEBI:60377 + CHEBI:15378")
    
    reaction <- "RHEA:11581"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "CDP-DG(12:0(11Me)/14:0) + myo-Inositol => PI(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CDP-DG(12:0(11Me)/14:0) + myo-Inositol")
    expect_equal(l$reaction_product, "PI(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58332 + CHEBI:17268 => CHEBI:57880 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58332 + CHEBI:17268")
    expect_equal(l$reaction_product_chebi, "CHEBI:57880 + CHEBI:60377 + CHEBI:15378")
    
    reaction <- "RHEA:11582"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "CDP-DG(12:0(11Me)/14:0) + myo-Inositol <= PI(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CDP-DG(12:0(11Me)/14:0) + myo-Inositol")
    expect_equal(l$reaction_product, "PI(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58332 + CHEBI:17268 <= CHEBI:57880 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58332 + CHEBI:17268")
    expect_equal(l$reaction_product_chebi, "CHEBI:57880 + CHEBI:60377 + CHEBI:15378")
    
    reaction <- "RHEA:11583"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "CDP-DG(12:0(11Me)/14:0) + myo-Inositol <=> PI(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CDP-DG(12:0(11Me)/14:0) + myo-Inositol")
    expect_equal(l$reaction_product, "PI(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58332 + CHEBI:17268 <=> CHEBI:57880 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58332 + CHEBI:17268")
    expect_equal(l$reaction_product_chebi, "CHEBI:57880 + CHEBI:60377 + CHEBI:15378")
    
    ## cer_to_glccer
    cer <- "Cer(d16:1(3OH,4OH)(15Me)/18:0)"
    substrates <- list(Cer = cer)
    
    reaction <- "RHEA:12088"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "Cer(d16:1(3OH,4OH)(15Me)/18:0) + UDP-Glucose = GlcCer(d16:1(3OH,4OH)(15Me)/18:0) + H+ + UDP")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "Cer(d16:1(3OH,4OH)(15Me)/18:0) + UDP-Glucose")
    expect_equal(l$reaction_product, "GlcCer(d16:1(3OH,4OH)(15Me)/18:0) + H+ + UDP")
    expect_equal(l$reaction_formula_chebi, "CHEBI:52639 + CHEBI:58885 = CHEBI:22801 + CHEBI:15378 + CHEBI:58223")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:52639 + CHEBI:58885")
    expect_equal(l$reaction_product_chebi, "CHEBI:22801 + CHEBI:15378 + CHEBI:58223")
    
    reaction <- "RHEA:12089"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "Cer(d16:1(3OH,4OH)(15Me)/18:0) + UDP-Glucose => GlcCer(d16:1(3OH,4OH)(15Me)/18:0) + H+ + UDP")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "Cer(d16:1(3OH,4OH)(15Me)/18:0) + UDP-Glucose")
    expect_equal(l$reaction_product, "GlcCer(d16:1(3OH,4OH)(15Me)/18:0) + H+ + UDP")
    expect_equal(l$reaction_formula_chebi, "CHEBI:52639 + CHEBI:58885 => CHEBI:22801 + CHEBI:15378 + CHEBI:58223")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:52639 + CHEBI:58885")
    expect_equal(l$reaction_product_chebi, "CHEBI:22801 + CHEBI:15378 + CHEBI:58223")
    
    reaction <- "RHEA:12090"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "Cer(d16:1(3OH,4OH)(15Me)/18:0) + UDP-Glucose <= GlcCer(d16:1(3OH,4OH)(15Me)/18:0) + H+ + UDP")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "Cer(d16:1(3OH,4OH)(15Me)/18:0) + UDP-Glucose")
    expect_equal(l$reaction_product, "GlcCer(d16:1(3OH,4OH)(15Me)/18:0) + H+ + UDP")
    expect_equal(l$reaction_formula_chebi, "CHEBI:52639 + CHEBI:58885 <= CHEBI:22801 + CHEBI:15378 + CHEBI:58223")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:52639 + CHEBI:58885")
    expect_equal(l$reaction_product_chebi, "CHEBI:22801 + CHEBI:15378 + CHEBI:58223")
    
    reaction <- "RHEA:12091"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "Cer(d16:1(3OH,4OH)(15Me)/18:0) + UDP-Glucose <=> GlcCer(d16:1(3OH,4OH)(15Me)/18:0) + H+ + UDP")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "Cer(d16:1(3OH,4OH)(15Me)/18:0) + UDP-Glucose")
    expect_equal(l$reaction_product, "GlcCer(d16:1(3OH,4OH)(15Me)/18:0) + H+ + UDP")
    expect_equal(l$reaction_formula_chebi, "CHEBI:52639 + CHEBI:58885 <=> CHEBI:22801 + CHEBI:15378 + CHEBI:58223")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:52639 + CHEBI:58885")
    expect_equal(l$reaction_product_chebi, "CHEBI:22801 + CHEBI:15378 + CHEBI:58223")
    
    ## cdpdg_to_pgp
    cdpdg <- "CDP-DG(12:0(11Me)/14:0)"
    substrates <- list(CDPDG = cdpdg)
    
    reaction <- "RHEA:12593"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "CDP-DG(12:0(11Me)/14:0) + Glycerol-3-P = PGP(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CDP-DG(12:0(11Me)/14:0) + Glycerol-3-P")
    expect_equal(l$reaction_product, "PGP(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58332 + CHEBI:57597 = CHEBI:60110 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58332 + CHEBI:57597")
    expect_equal(l$reaction_product_chebi, "CHEBI:60110 + CHEBI:60377 + CHEBI:15378")
    
    reaction <- "RHEA:12594"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "CDP-DG(12:0(11Me)/14:0) + Glycerol-3-P => PGP(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CDP-DG(12:0(11Me)/14:0) + Glycerol-3-P")
    expect_equal(l$reaction_product, "PGP(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58332 + CHEBI:57597 => CHEBI:60110 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58332 + CHEBI:57597")
    expect_equal(l$reaction_product_chebi, "CHEBI:60110 + CHEBI:60377 + CHEBI:15378")
    
    reaction <- "RHEA:12595"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "CDP-DG(12:0(11Me)/14:0) + Glycerol-3-P <= PGP(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CDP-DG(12:0(11Me)/14:0) + Glycerol-3-P")
    expect_equal(l$reaction_product, "PGP(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58332 + CHEBI:57597 <= CHEBI:60110 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58332 + CHEBI:57597")
    expect_equal(l$reaction_product_chebi, "CHEBI:60110 + CHEBI:60377 + CHEBI:15378")
    
    reaction <- "RHEA:12596"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "CDP-DG(12:0(11Me)/14:0) + Glycerol-3-P <=> PGP(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CDP-DG(12:0(11Me)/14:0) + Glycerol-3-P")
    expect_equal(l$reaction_product, "PGP(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58332 + CHEBI:57597 <=> CHEBI:60110 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58332 + CHEBI:57597")
    expect_equal(l$reaction_product_chebi, "CHEBI:60110 + CHEBI:60377 + CHEBI:15378")
    
    ## sn1lpc_to_pc
    sn1lpc <- "PC(14:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPC = sn1lpc, AcylCoA = acylcoa)
    
    reaction <- "RHEA:12937"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(14:0/0:0) + CoA(18:0) = PC(14:0/18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(14:0/0:0) + CoA(18:0)")
    expect_equal(l$reaction_product, "PC(14:0/18:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58168 + CHEBI:58342 = CHEBI:57643 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58168 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:57643 + CHEBI:57287")
    
    reaction <- "RHEA:12938"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(14:0/0:0) + CoA(18:0) => PC(14:0/18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(14:0/0:0) + CoA(18:0)")
    expect_equal(l$reaction_product, "PC(14:0/18:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58168 + CHEBI:58342 => CHEBI:57643 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58168 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:57643 + CHEBI:57287")
    
    reaction <- "RHEA:12939"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(14:0/0:0) + CoA(18:0) <= PC(14:0/18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(14:0/0:0) + CoA(18:0)")
    expect_equal(l$reaction_product, "PC(14:0/18:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58168 + CHEBI:58342 <= CHEBI:57643 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58168 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:57643 + CHEBI:57287")
    
    reaction <- "RHEA:12940"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(14:0/0:0) + CoA(18:0) <=> PC(14:0/18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(14:0/0:0) + CoA(18:0)")
    expect_equal(l$reaction_product, "PC(14:0/18:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58168 + CHEBI:58342 <=> CHEBI:57643 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58168 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:57643 + CHEBI:57287")
    
    ## pc_to_pa
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    
    reaction <- "RHEA:14445"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(20:0/18:0) + H2O = PA(20:0/18:0) + Choline + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(20:0/18:0) + H2O")
    expect_equal(l$reaction_product, "PA(20:0/18:0) + Choline + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 = CHEBI:58608 + CHEBI:15354 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:58608 + CHEBI:15354 + CHEBI:15378")
    
    reaction <- "RHEA:14446"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(20:0/18:0) + H2O => PA(20:0/18:0) + Choline + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(20:0/18:0) + H2O")
    expect_equal(l$reaction_product, "PA(20:0/18:0) + Choline + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 => CHEBI:58608 + CHEBI:15354 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:58608 + CHEBI:15354 + CHEBI:15378")
    
    reaction <- "RHEA:14447"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(20:0/18:0) + H2O <= PA(20:0/18:0) + Choline + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(20:0/18:0) + H2O")
    expect_equal(l$reaction_product, "PA(20:0/18:0) + Choline + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <= CHEBI:58608 + CHEBI:15354 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:58608 + CHEBI:15354 + CHEBI:15378")
    
    reaction <- "RHEA:14448"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(20:0/18:0) + H2O <=> PA(20:0/18:0) + Choline + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(20:0/18:0) + H2O")
    expect_equal(l$reaction_product, "PA(20:0/18:0) + Choline + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <=> CHEBI:58608 + CHEBI:15354 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:58608 + CHEBI:15354 + CHEBI:15378")
    
    ## sn1lpc_to_fa
    sn1lpc <- "PC(14:0/0:0)"
    substrates <- list(sn1LPC = sn1lpc)
    
    reaction <- "RHEA:15177"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(14:0/0:0) + H2O = FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(14:0/0:0) + H2O")
    expect_equal(l$reaction_product, "FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58168 + CHEBI:15377 = CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58168 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    
    reaction <- "RHEA:15178"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(14:0/0:0) + H2O => FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(14:0/0:0) + H2O")
    expect_equal(l$reaction_product, "FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58168 + CHEBI:15377 => CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58168 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    
    reaction <- "RHEA:15179"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(14:0/0:0) + H2O <= FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(14:0/0:0) + H2O")
    expect_equal(l$reaction_product, "FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58168 + CHEBI:15377 <= CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58168 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    
    reaction <- "RHEA:15180"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(14:0/0:0) + H2O <=> FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(14:0/0:0) + H2O")
    expect_equal(l$reaction_product, "FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58168 + CHEBI:15377 <=> CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58168 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    
    ## coa_to_lpa
    acylcoa <- "CoA(18:0)"
    substrates <- list(AcylCoA = acylcoa)
    
    reaction <- "RHEA:15325"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CoA(18:0) + Glycerol-3-P = PA(18:0/0:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CoA(18:0) + Glycerol-3-P")
    expect_equal(l$reaction_product, "PA(18:0/0:0) + CoA") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57597 = CHEBI:57970 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58342 + CHEBI:57597")
    expect_equal(l$reaction_product_chebi, "CHEBI:57970 + CHEBI:57287")
    
    reaction <- "RHEA:15326"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CoA(18:0) + Glycerol-3-P => PA(18:0/0:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CoA(18:0) + Glycerol-3-P")
    expect_equal(l$reaction_product, "PA(18:0/0:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57597 => CHEBI:57970 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58342 + CHEBI:57597")
    expect_equal(l$reaction_product_chebi, "CHEBI:57970 + CHEBI:57287")
    
    reaction <- "RHEA:15327"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CoA(18:0) + Glycerol-3-P <= PA(18:0/0:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CoA(18:0) + Glycerol-3-P")
    expect_equal(l$reaction_product, "PA(18:0/0:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57597 <= CHEBI:57970 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58342 + CHEBI:57597")
    expect_equal(l$reaction_product_chebi, "CHEBI:57970 + CHEBI:57287")
    
    reaction <- "RHEA:15328"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CoA(18:0) + Glycerol-3-P <=> PA(18:0/0:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CoA(18:0) + Glycerol-3-P")
    expect_equal(l$reaction_product, "PA(18:0/0:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57597 <=> CHEBI:57970 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58342 + CHEBI:57597")
    expect_equal(l$reaction_product_chebi, "CHEBI:57970 + CHEBI:57287")
    
    ## fa_to_coa
    fa <- "FA(18:0)"
    substrates = list(FA = fa)
    
    reaction <- "RHEA:15421"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "FA(18:0) + ATP + CoA = CoA(18:0) + AMP + PPi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "FA(18:0) + ATP + CoA")
    expect_equal(l$reaction_product, "CoA(18:0) + AMP + PPi") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:57560 + CHEBI:30616 + CHEBI:57287 = CHEBI:83139 + CHEBI:456215 + CHEBI:33019")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57560 + CHEBI:30616 + CHEBI:57287")
    expect_equal(l$reaction_product_chebi, "CHEBI:83139 + CHEBI:456215 + CHEBI:33019")
    
    reaction <- "RHEA:15422"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "FA(18:0) + ATP + CoA => CoA(18:0) + AMP + PPi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "FA(18:0) + ATP + CoA")
    expect_equal(l$reaction_product, "CoA(18:0) + AMP + PPi") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:57560 + CHEBI:30616 + CHEBI:57287 => CHEBI:83139 + CHEBI:456215 + CHEBI:33019")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57560 + CHEBI:30616 + CHEBI:57287")
    expect_equal(l$reaction_product_chebi, "CHEBI:83139 + CHEBI:456215 + CHEBI:33019")
    
    reaction <- "RHEA:15423"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "FA(18:0) + ATP + CoA <= CoA(18:0) + AMP + PPi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "FA(18:0) + ATP + CoA")
    expect_equal(l$reaction_product, "CoA(18:0) + AMP + PPi")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57560 + CHEBI:30616 + CHEBI:57287 <= CHEBI:83139 + CHEBI:456215 + CHEBI:33019")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57560 + CHEBI:30616 + CHEBI:57287")
    expect_equal(l$reaction_product_chebi, "CHEBI:83139 + CHEBI:456215 + CHEBI:33019")
    
    reaction <- "RHEA:15424"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "FA(18:0) + ATP + CoA <=> CoA(18:0) + AMP + PPi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "FA(18:0) + ATP + CoA")
    expect_equal(l$reaction_product, "CoA(18:0) + AMP + PPi") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:57560 + CHEBI:30616 + CHEBI:57287 <=> CHEBI:83139 + CHEBI:456215 + CHEBI:33019")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57560 + CHEBI:30616 + CHEBI:57287")
    expect_equal(l$reaction_product_chebi, "CHEBI:83139 + CHEBI:456215 + CHEBI:33019")
    
    reaction <- "RHEA:38883"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "FA(18:0) + ATP + CoA = CoA(18:0) + AMP + PPi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "FA(18:0) + ATP + CoA")
    expect_equal(l$reaction_product, "CoA(18:0) + AMP + PPi")
    expect_equal(l$reaction_formula_chebi, "CHEBI:28868 + CHEBI:30616 + CHEBI:57287 = CHEBI:77636 + CHEBI:456215 + CHEBI:33019")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:28868 + CHEBI:30616 + CHEBI:57287")
    expect_equal(l$reaction_product_chebi, "CHEBI:77636 + CHEBI:456215 + CHEBI:33019")
    
    reaction <- "RHEA:38884"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "FA(18:0) + ATP + CoA => CoA(18:0) + AMP + PPi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "FA(18:0) + ATP + CoA")
    expect_equal(l$reaction_product, "CoA(18:0) + AMP + PPi") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:28868 + CHEBI:30616 + CHEBI:57287 => CHEBI:77636 + CHEBI:456215 + CHEBI:33019")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:28868 + CHEBI:30616 + CHEBI:57287")
    expect_equal(l$reaction_product_chebi, "CHEBI:77636 + CHEBI:456215 + CHEBI:33019")
    
    reaction <- "RHEA:38885"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "FA(18:0) + ATP + CoA <= CoA(18:0) + AMP + PPi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "FA(18:0) + ATP + CoA")
    expect_equal(l$reaction_product, "CoA(18:0) + AMP + PPi")
    expect_equal(l$reaction_formula_chebi, "CHEBI:28868 + CHEBI:30616 + CHEBI:57287 <= CHEBI:77636 + CHEBI:456215 + CHEBI:33019")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:28868 + CHEBI:30616 + CHEBI:57287")
    expect_equal(l$reaction_product_chebi, "CHEBI:77636 + CHEBI:456215 + CHEBI:33019")
    
    reaction <- "RHEA:38886"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "FA(18:0) + ATP + CoA <=> CoA(18:0) + AMP + PPi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "FA(18:0) + ATP + CoA")
    expect_equal(l$reaction_product, "CoA(18:0) + AMP + PPi")
    expect_equal(l$reaction_formula_chebi, "CHEBI:28868 + CHEBI:30616 + CHEBI:57287 <=> CHEBI:77636 + CHEBI:456215 + CHEBI:33019")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:28868 + CHEBI:30616 + CHEBI:57287")
    expect_equal(l$reaction_product_chebi, "CHEBI:77636 + CHEBI:456215 + CHEBI:33019")
    
    ## pc_to_sn1lpc
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    
    reaction <- "RHEA:15801"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(20:0/18:0) + H2O = PC(20:0/0:0) + FA(18:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(20:0/18:0) + H2O")
    expect_equal(l$reaction_product, "PC(20:0/0:0) + FA(18:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 = CHEBI:58168 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:58168 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:15802"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(20:0/18:0) + H2O => PC(20:0/0:0) + FA(18:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(20:0/18:0) + H2O")
    expect_equal(l$reaction_product, "PC(20:0/0:0) + FA(18:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 => CHEBI:58168 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:58168 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:15803"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(20:0/18:0) + H2O <= PC(20:0/0:0) + FA(18:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(20:0/18:0) + H2O")
    expect_equal(l$reaction_product, "PC(20:0/0:0) + FA(18:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <= CHEBI:58168 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:58168 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:15804"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(20:0/18:0) + H2O <=> PC(20:0/0:0) + FA(18:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(20:0/18:0) + H2O")
    expect_equal(l$reaction_product, "PC(20:0/0:0) + FA(18:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <=> CHEBI:58168 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:58168 + CHEBI:28868 + CHEBI:15378")
    
    ## pa_to_cdpdg
    pa <- "PA(14:0/16:0)"
    substrates <- list(PA = pa)
    
    reaction <- "RHEA:16229"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PA(14:0/16:0) + CTP + H+ = CDP-DG(14:0/16:0) + PPi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PA(14:0/16:0) + CTP + H+")
    expect_equal(l$reaction_product, "CDP-DG(14:0/16:0) + PPi")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58608 + CHEBI:37563 + CHEBI:15378 = CHEBI:58332 + CHEBI:33019")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58608 + CHEBI:37563 + CHEBI:15378")
    expect_equal(l$reaction_product_chebi, "CHEBI:58332 + CHEBI:33019")
    
    reaction <- "RHEA:16230"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PA(14:0/16:0) + CTP + H+ => CDP-DG(14:0/16:0) + PPi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PA(14:0/16:0) + CTP + H+")
    expect_equal(l$reaction_product, "CDP-DG(14:0/16:0) + PPi")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58608 + CHEBI:37563 + CHEBI:15378 => CHEBI:58332 + CHEBI:33019")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58608 + CHEBI:37563 + CHEBI:15378")
    expect_equal(l$reaction_product_chebi, "CHEBI:58332 + CHEBI:33019")
    
    reaction <- "RHEA:16231"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PA(14:0/16:0) + CTP + H+ <= CDP-DG(14:0/16:0) + PPi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PA(14:0/16:0) + CTP + H+")
    expect_equal(l$reaction_product, "CDP-DG(14:0/16:0) + PPi")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58608 + CHEBI:37563 + CHEBI:15378 <= CHEBI:58332 + CHEBI:33019")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58608 + CHEBI:37563 + CHEBI:15378")
    expect_equal(l$reaction_product_chebi, "CHEBI:58332 + CHEBI:33019")
    
    reaction <- "RHEA:16232"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PA(14:0/16:0) + CTP + H+ <=> CDP-DG(14:0/16:0) + PPi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PA(14:0/16:0) + CTP + H+")
    expect_equal(l$reaction_product, "CDP-DG(14:0/16:0) + PPi")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58608 + CHEBI:37563 + CHEBI:15378 <=> CHEBI:58332 + CHEBI:33019")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58608 + CHEBI:37563 + CHEBI:15378")
    expect_equal(l$reaction_product_chebi, "CHEBI:58332 + CHEBI:33019")
    
    ## lpep_to_pep
    lpep <- "PE(P-16:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPEP = lpep, AcylCoA = acylcoa)
    
    reaction <- "RHEA:16245"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/0:0) + CoA(18:0) = PE(P-16:0/18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/0:0) + CoA(18:0)")
    expect_equal(l$reaction_product, "PE(P-16:0/18:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77288 + CHEBI:58342 = CHEBI:77290 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77288 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:77290 + CHEBI:57287")
    
    reaction <- "RHEA:16246"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/0:0) + CoA(18:0) => PE(P-16:0/18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/0:0) + CoA(18:0)")
    expect_equal(l$reaction_product, "PE(P-16:0/18:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77288 + CHEBI:58342 => CHEBI:77290 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77288 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:77290 + CHEBI:57287")
    
    reaction <- "RHEA:16247"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/0:0) + CoA(18:0) <= PE(P-16:0/18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/0:0) + CoA(18:0)")
    expect_equal(l$reaction_product, "PE(P-16:0/18:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77288 + CHEBI:58342 <= CHEBI:77290 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77288 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:77290 + CHEBI:57287")
    
    reaction <- "RHEA:16248"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/0:0) + CoA(18:0) <=> PE(P-16:0/18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/0:0) + CoA(18:0)")
    expect_equal(l$reaction_product, "PE(P-16:0/18:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77288 + CHEBI:58342 <=> CHEBI:77290 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77288 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:77290 + CHEBI:57287")
    
    ## sn2mg_to_dg
    sn2mg <- "MG(0:0/14:0/0:0)"
    acylcoa <- "CoA(16:0)"
    substrates <- list(sn2MG = sn2mg, AcylCoA = acylcoa)
    
    reaction <- "RHEA:16741"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(0:0/14:0/0:0) + CoA(16:0) = DG(16:0/14:0/0:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(0:0/14:0/0:0) + CoA(16:0)")
    expect_equal(l$reaction_product, "DG(16:0/14:0/0:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 = CHEBI:49172 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17389 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:49172 + CHEBI:57287")
    
    reaction <- "RHEA:16742"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(0:0/14:0/0:0) + CoA(16:0) => DG(16:0/14:0/0:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(0:0/14:0/0:0) + CoA(16:0)")
    expect_equal(l$reaction_product, "DG(16:0/14:0/0:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 => CHEBI:49172 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17389 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:49172 + CHEBI:57287")
    
    reaction <- "RHEA:16743"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(0:0/14:0/0:0) + CoA(16:0) <= DG(16:0/14:0/0:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(0:0/14:0/0:0) + CoA(16:0)")
    expect_equal(l$reaction_product, "DG(16:0/14:0/0:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 <= CHEBI:49172 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17389 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:49172 + CHEBI:57287")
    
    reaction <- "RHEA:16744"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(0:0/14:0/0:0) + CoA(16:0) <=> DG(16:0/14:0/0:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(0:0/14:0/0:0) + CoA(16:0)")
    expect_equal(l$reaction_product, "DG(16:0/14:0/0:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 <=> CHEBI:49172 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17389 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:49172 + CHEBI:57287")
    
    reaction <- "RHEA:32947"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(0:0/14:0/0:0) + CoA(16:0) = DG(16:0/14:0/0:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(0:0/14:0/0:0) + CoA(16:0)")
    expect_equal(l$reaction_product, "DG(16:0/14:0/0:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 = CHEBI:17815 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17389 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:57287")
    
    reaction <- "RHEA:32948"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(0:0/14:0/0:0) + CoA(16:0) => DG(16:0/14:0/0:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(0:0/14:0/0:0) + CoA(16:0)")
    expect_equal(l$reaction_product, "DG(16:0/14:0/0:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 => CHEBI:17815 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17389 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:57287")
    
    reaction <- "RHEA:32949"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(0:0/14:0/0:0) + CoA(16:0) <= DG(16:0/14:0/0:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(0:0/14:0/0:0) + CoA(16:0)")
    expect_equal(l$reaction_product, "DG(16:0/14:0/0:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 <= CHEBI:17815 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17389 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:57287")
    
    reaction <- "RHEA:32950"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(0:0/14:0/0:0) + CoA(16:0) <=> DG(16:0/14:0/0:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(0:0/14:0/0:0) + CoA(16:0)")
    expect_equal(l$reaction_product, "DG(16:0/14:0/0:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 <=> CHEBI:17815 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17389 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:57287")
    
    ## lpep_to_fal
    sn1lpep <- "PE(P-16:0/0:0)"
    substrates = list(sn1LPEP = sn1lpep)
    
    reaction <- "RHEA:16905"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/0:0) + H2O = FAL(16:0) + Glycerophosphoethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "FAL(16:0) + Glycerophosphoethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 = CHEBI:73359 + CHEBI:143890")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77288 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:73359 + CHEBI:143890")
    
    reaction <- "RHEA:16906"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/0:0) + H2O => FAL(16:0) + Glycerophosphoethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "FAL(16:0) + Glycerophosphoethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 => CHEBI:73359 + CHEBI:143890")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77288 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:73359 + CHEBI:143890")
    
    reaction <- "RHEA:16907"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/0:0) + H2O <= FAL(16:0) + Glycerophosphoethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "FAL(16:0) + Glycerophosphoethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 <= CHEBI:73359 + CHEBI:143890")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77288 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:73359 + CHEBI:143890")
    
    reaction <- "RHEA:16908"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/0:0) + H2O <=> FAL(16:0) + Glycerophosphoethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "FAL(16:0) + Glycerophosphoethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 <=> CHEBI:73359 + CHEBI:143890")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77288 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:73359 + CHEBI:143890")
    
    ## nae_to_fa
    nae <- "NAE(18:0)"
    substrates <- list(NAE = nae)
    
    reaction <- "RHEA:17505"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + NAE(18:0) = FA(18:0) + Ethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + NAE(18:0)")
    expect_equal(l$reaction_product, "FA(18:0) + Ethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:15897 = CHEBI:57560 + CHEBI:57603")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:15897")
    expect_equal(l$reaction_product_chebi, "CHEBI:57560 + CHEBI:57603")
    
    reaction <- "RHEA:17506"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + NAE(18:0) => FA(18:0) + Ethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + NAE(18:0)")
    expect_equal(l$reaction_product, "FA(18:0) + Ethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:15897 => CHEBI:57560 + CHEBI:57603")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:15897")
    expect_equal(l$reaction_product_chebi, "CHEBI:57560 + CHEBI:57603")
    
    reaction <- "RHEA:17507"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + NAE(18:0) <= FA(18:0) + Ethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + NAE(18:0)")
    expect_equal(l$reaction_product, "FA(18:0) + Ethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:15897 <= CHEBI:57560 + CHEBI:57603")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:15897")
    expect_equal(l$reaction_product_chebi, "CHEBI:57560 + CHEBI:57603")
    
    reaction <- "RHEA:17508"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + NAE(18:0) <=> FA(18:0) + Ethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + NAE(18:0)")
    expect_equal(l$reaction_product, "FA(18:0) + Ethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:15897 <=> CHEBI:57560 + CHEBI:57603")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:15897")
    expect_equal(l$reaction_product_chebi, "CHEBI:57560 + CHEBI:57603")
    
    reaction <- "RHEA:39995"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "H2O + NAE(18:0) = FA(18:0) + Ethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + NAE(18:0)")
    expect_equal(l$reaction_product, "FA(18:0) + Ethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:52640 = CHEBI:28868 + CHEBI:57603")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:52640")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:57603")
    
    reaction <- "RHEA:39996"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + NAE(18:0) => FA(18:0) + Ethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + NAE(18:0)")
    expect_equal(l$reaction_product, "FA(18:0) + Ethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:52640 => CHEBI:28868 + CHEBI:57603")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:52640")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:57603")
    
    reaction <- "RHEA:39997"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + NAE(18:0) <= FA(18:0) + Ethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + NAE(18:0)")
    expect_equal(l$reaction_product, "FA(18:0) + Ethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:52640 <= CHEBI:28868 + CHEBI:57603")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:52640")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:57603")
    
    reaction <- "RHEA:39998"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + NAE(18:0) <=> FA(18:0) + Ethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + NAE(18:0)")
    expect_equal(l$reaction_product, "FA(18:0) + Ethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:52640 <=> CHEBI:28868 + CHEBI:57603")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:52640")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:57603")
    
    ## coa_to_acyldhap
    acylcoa <- "CoA(18:0)"
    substrates <- list(AcylCoA = acylcoa)
    
    reaction <- "RHEA:17657"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CoA(18:0) + Dihydroxyacetone-P = DHAP(18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CoA(18:0) + Dihydroxyacetone-P")
    expect_equal(l$reaction_product, "DHAP(18:0) + CoA") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57642 = CHEBI:57534 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58342 + CHEBI:57642")
    expect_equal(l$reaction_product_chebi, "CHEBI:57534 + CHEBI:57287")
    
    reaction <- "RHEA:17658"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CoA(18:0) + Dihydroxyacetone-P => DHAP(18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CoA(18:0) + Dihydroxyacetone-P")
    expect_equal(l$reaction_product, "DHAP(18:0) + CoA") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57642 => CHEBI:57534 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58342 + CHEBI:57642")
    expect_equal(l$reaction_product_chebi, "CHEBI:57534 + CHEBI:57287")
    
    reaction <- "RHEA:17659"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CoA(18:0) + Dihydroxyacetone-P <= DHAP(18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CoA(18:0) + Dihydroxyacetone-P")
    expect_equal(l$reaction_product, "DHAP(18:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57642 <= CHEBI:57534 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58342 + CHEBI:57642")
    expect_equal(l$reaction_product_chebi, "CHEBI:57534 + CHEBI:57287")
    
    reaction <- "RHEA:17660"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CoA(18:0) + Dihydroxyacetone-P <=> DHAP(18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CoA(18:0) + Dihydroxyacetone-P")
    expect_equal(l$reaction_product, "DHAP(18:0) + CoA") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57642 <=> CHEBI:57534 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58342 + CHEBI:57642")
    expect_equal(l$reaction_product_chebi, "CHEBI:57534 + CHEBI:57287")
    
    ## coa_to_ce
    acylcoa <- "CoA(18:0)"
    substrates <- list(AcylCoA = acylcoa)
    
    reaction <- "RHEA:17729"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "Cholesterol + CoA(18:0) = CE(18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "Cholesterol + CoA(18:0)")
    expect_equal(l$reaction_product, "CE(18:0) + CoA") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:16113 + CHEBI:58342 = CHEBI:17002 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:16113 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:17002 + CHEBI:57287")
    
    reaction <- "RHEA:17730"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "Cholesterol + CoA(18:0) => CE(18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "Cholesterol + CoA(18:0)")
    expect_equal(l$reaction_product, "CE(18:0) + CoA") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:16113 + CHEBI:58342 => CHEBI:17002 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:16113 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:17002 + CHEBI:57287")
    
    reaction <- "RHEA:17731"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "Cholesterol + CoA(18:0) <= CE(18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "Cholesterol + CoA(18:0)")
    expect_equal(l$reaction_product, "CE(18:0) + CoA") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:16113 + CHEBI:58342 <= CHEBI:17002 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:16113 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:17002 + CHEBI:57287")
    
    reaction <- "RHEA:17732"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "Cholesterol + CoA(18:0) <=> CE(18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "Cholesterol + CoA(18:0)")
    expect_equal(l$reaction_product, "CE(18:0) + CoA") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:16113 + CHEBI:58342 <=> CHEBI:17002 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:16113 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:17002 + CHEBI:57287")
    
    ## cer_to_cerp
    cer <- "Cer(d16:1(3OH,4OH)(15Me)/18:0)"
    substrates <- list(Cer = cer)
    
    reaction <- "RHEA:17929"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction, 
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "Cer(d16:1(3OH,4OH)(15Me)/18:0) + ATP = ADP + CerP(d16:1(3OH,4OH)(15Me)/18:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "Cer(d16:1(3OH,4OH)(15Me)/18:0) + ATP")
    expect_equal(l$reaction_product, "ADP + CerP(d16:1(3OH,4OH)(15Me)/18:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:52639 + CHEBI:30616 = CHEBI:456216 + CHEBI:57674 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:52639 + CHEBI:30616")
    expect_equal(l$reaction_product_chebi, "CHEBI:456216 + CHEBI:57674 + CHEBI:15378")
    
    reaction <- "RHEA:17930"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "Cer(d16:1(3OH,4OH)(15Me)/18:0) + ATP => ADP + CerP(d16:1(3OH,4OH)(15Me)/18:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "Cer(d16:1(3OH,4OH)(15Me)/18:0) + ATP")
    expect_equal(l$reaction_product, "ADP + CerP(d16:1(3OH,4OH)(15Me)/18:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:52639 + CHEBI:30616 => CHEBI:456216 + CHEBI:57674 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:52639 + CHEBI:30616")
    expect_equal(l$reaction_product_chebi, "CHEBI:456216 + CHEBI:57674 + CHEBI:15378")
    
    reaction <- "RHEA:17931"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "Cer(d16:1(3OH,4OH)(15Me)/18:0) + ATP <= ADP + CerP(d16:1(3OH,4OH)(15Me)/18:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "Cer(d16:1(3OH,4OH)(15Me)/18:0) + ATP")
    expect_equal(l$reaction_product, "ADP + CerP(d16:1(3OH,4OH)(15Me)/18:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:52639 + CHEBI:30616 <= CHEBI:456216 + CHEBI:57674 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:52639 + CHEBI:30616")
    expect_equal(l$reaction_product_chebi, "CHEBI:456216 + CHEBI:57674 + CHEBI:15378")
    
    reaction <- "RHEA:17932"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "Cer(d16:1(3OH,4OH)(15Me)/18:0) + ATP <=> ADP + CerP(d16:1(3OH,4OH)(15Me)/18:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "Cer(d16:1(3OH,4OH)(15Me)/18:0) + ATP")
    expect_equal(l$reaction_product, "ADP + CerP(d16:1(3OH,4OH)(15Me)/18:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:52639 + CHEBI:30616 <=> CHEBI:456216 + CHEBI:57674 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:52639 + CHEBI:30616")
    expect_equal(l$reaction_product_chebi, "CHEBI:456216 + CHEBI:57674 + CHEBI:15378")

    ## pi_to_sn1lpi
    pi <- "PI(16:0/18:1(9Z))"
    substrates <- list(PI = pi)
    
    reaction <- "RHEA:18001"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PI(16:0/18:1(9Z)) + H2O = PI(16:0/0:0) + FA(18:1(9Z)) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PI(16:0/18:1(9Z)) + H2O")
    expect_equal(l$reaction_product, "PI(16:0/0:0) + FA(18:1(9Z)) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 = CHEBI:64771 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57880 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64771 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:18002"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PI(16:0/18:1(9Z)) + H2O => PI(16:0/0:0) + FA(18:1(9Z)) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PI(16:0/18:1(9Z)) + H2O")
    expect_equal(l$reaction_product, "PI(16:0/0:0) + FA(18:1(9Z)) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 => CHEBI:64771 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57880 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64771 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:18003"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PI(16:0/18:1(9Z)) + H2O <= PI(16:0/0:0) + FA(18:1(9Z)) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PI(16:0/18:1(9Z)) + H2O")
    expect_equal(l$reaction_product, "PI(16:0/0:0) + FA(18:1(9Z)) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 <= CHEBI:64771 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57880 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64771 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:18004"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PI(16:0/18:1(9Z)) + H2O <=> PI(16:0/0:0) + FA(18:1(9Z)) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PI(16:0/18:1(9Z)) + H2O")
    expect_equal(l$reaction_product, "PI(16:0/0:0) + FA(18:1(9Z)) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 <=> CHEBI:64771 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57880 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64771 + CHEBI:28868 + CHEBI:15378")
    
    ## pc_to_sn2lpc
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    
    reaction <- "RHEA:18689"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(20:0/18:0) + H2O = PC(0:0/18:0) + FA(20:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(20:0/18:0) + H2O")
    expect_equal(l$reaction_product, "PC(0:0/18:0) + FA(20:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 = CHEBI:57875 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:57875 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:18690"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(20:0/18:0) + H2O => PC(0:0/18:0) + FA(20:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(20:0/18:0) + H2O")
    expect_equal(l$reaction_product, "PC(0:0/18:0) + FA(20:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 => CHEBI:57875 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:57875 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:18691"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(20:0/18:0) + H2O <= PC(0:0/18:0) + FA(20:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(20:0/18:0) + H2O")
    expect_equal(l$reaction_product, "PC(0:0/18:0) + FA(20:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <= CHEBI:57875 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:57875 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:18692"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(20:0/18:0) + H2O <=> PC(0:0/18:0) + FA(20:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(20:0/18:0) + H2O")
    expect_equal(l$reaction_product, "PC(0:0/18:0) + FA(20:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <=> CHEBI:57875 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:57875 + CHEBI:28868 + CHEBI:15378")
    
    ## cer_to_sm
    cer <- "Cer(d16:1(3OH,4OH)(15Me)/18:0)"
    substrates <- list(Cer = cer)
    
    reaction <-  "RHEA:18765"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC + Cer(d16:1(3OH,4OH)(15Me)/18:0) = DG + SM(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC + Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_product, "DG + SM(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:52639 = CHEBI:17815 + CHEBI:17636")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:52639")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:17636")
    
    reaction <-  "RHEA:18766"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC + Cer(d16:1(3OH,4OH)(15Me)/18:0) => DG + SM(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC + Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_product, "DG + SM(d16:1(3OH,4OH)(15Me)/18:0)") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:52639 => CHEBI:17815 + CHEBI:17636")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:52639")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:17636")
    
    reaction <-  "RHEA:18767"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC + Cer(d16:1(3OH,4OH)(15Me)/18:0) <= DG + SM(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC + Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_product, "DG + SM(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:52639 <= CHEBI:17815 + CHEBI:17636")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:52639")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:17636")
    
    reaction <-  "RHEA:18768"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC + Cer(d16:1(3OH,4OH)(15Me)/18:0) <=> DG + SM(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC + Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_product, "DG + SM(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:52639 <=> CHEBI:17815 + CHEBI:17636")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:52639")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:17636")
    
    ## lpa_to_pa
    lpa <- "PA(18:0/0:0)"
    acylcoa <- "CoA(14:0)"
    substrates <- list(sn1LPA = lpa, AcylCoA = acylcoa)
    
    reaction <- "RHEA:19709"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PA(18:0/0:0) + CoA(14:0) = PA(18:0/14:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PA(18:0/0:0) + CoA(14:0)")
    expect_equal(l$reaction_product, "PA(18:0/14:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57970 + CHEBI:58342 = CHEBI:58608 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57970 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:58608 + CHEBI:57287")
    
    reaction <- "RHEA:19710"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PA(18:0/0:0) + CoA(14:0) => PA(18:0/14:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PA(18:0/0:0) + CoA(14:0)")
    expect_equal(l$reaction_product, "PA(18:0/14:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57970 + CHEBI:58342 => CHEBI:58608 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57970 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:58608 + CHEBI:57287")
    
    reaction <- "RHEA:19711"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PA(18:0/0:0) + CoA(14:0) <= PA(18:0/14:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PA(18:0/0:0) + CoA(14:0)")
    expect_equal(l$reaction_product, "PA(18:0/14:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57970 + CHEBI:58342 <= CHEBI:58608 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57970 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:58608 + CHEBI:57287")
    
    reaction <- "RHEA:19712"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PA(18:0/0:0) + CoA(14:0) <=> PA(18:0/14:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PA(18:0/0:0) + CoA(14:0)")
    expect_equal(l$reaction_product, "PA(18:0/14:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57970 + CHEBI:58342 <=> CHEBI:58608 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57970 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:58608 + CHEBI:57287")
    
    ## ps_to_pe
    ps <- "PS(14:0/14:0)"
    substrates = list(PS = ps)
    
    reaction <- "RHEA:20828"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PS(14:0/14:0) + H+ = PE(14:0/14:0) + CO2")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PS(14:0/14:0) + H+")
    expect_equal(l$reaction_product, "PE(14:0/14:0) + CO2")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57262 + CHEBI:15378 = CHEBI:64612 + CHEBI:16526")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57262 + CHEBI:15378")
    expect_equal(l$reaction_product_chebi, "CHEBI:64612 + CHEBI:16526")
    
    reaction <- "RHEA:20829"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PS(14:0/14:0) + H+ => PE(14:0/14:0) + CO2")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PS(14:0/14:0) + H+")
    expect_equal(l$reaction_product, "PE(14:0/14:0) + CO2")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57262 + CHEBI:15378 => CHEBI:64612 + CHEBI:16526")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57262 + CHEBI:15378")
    expect_equal(l$reaction_product_chebi, "CHEBI:64612 + CHEBI:16526")
    
    reaction <- "RHEA:20830"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PS(14:0/14:0) + H+ <= PE(14:0/14:0) + CO2")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PS(14:0/14:0) + H+")
    expect_equal(l$reaction_product, "PE(14:0/14:0) + CO2")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57262 + CHEBI:15378 <= CHEBI:64612 + CHEBI:16526")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57262 + CHEBI:15378")
    expect_equal(l$reaction_product_chebi, "CHEBI:64612 + CHEBI:16526")
    
    reaction <- "RHEA:20831"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PS(14:0/14:0) + H+ <=> PE(14:0/14:0) + CO2")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PS(14:0/14:0) + H+")
    expect_equal(l$reaction_product, "PE(14:0/14:0) + CO2")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57262 + CHEBI:15378 <=> CHEBI:64612 + CHEBI:16526")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57262 + CHEBI:15378")
    expect_equal(l$reaction_product_chebi, "CHEBI:64612 + CHEBI:16526")
    
    ## peo_to_pep
    peo <- "PE(O-16:0/14:0)"
    substrates <- list(PEO = peo)
    
    reaction <- "RHEA:22956"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2 = PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2")
    expect_equal(l$reaction_product, "PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    expect_equal(l$reaction_formula_chebi, "CHEBI:75028 + CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 = CHEBI:77290 + CHEBI:29034 + 2 CHEBI:15377")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:75028 + CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379")
    expect_equal(l$reaction_product_chebi, "CHEBI:77290 + CHEBI:29034 + 2 CHEBI:15377")
    
    reaction <- "RHEA:22957"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2 => PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2")
    expect_equal(l$reaction_product, "PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    expect_equal(l$reaction_formula_chebi, "CHEBI:75028 + CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 => CHEBI:77290 + CHEBI:29034 + 2 CHEBI:15377")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:75028 + CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379")
    expect_equal(l$reaction_product_chebi, "CHEBI:77290 + CHEBI:29034 + 2 CHEBI:15377")
    
    reaction <- "RHEA:22958"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2 <= PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2")
    expect_equal(l$reaction_product, "PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    expect_equal(l$reaction_formula_chebi, "CHEBI:75028 + CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 <= CHEBI:77290 + CHEBI:29034 + 2 CHEBI:15377")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:75028 + CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379")
    expect_equal(l$reaction_product_chebi, "CHEBI:77290 + CHEBI:29034 + 2 CHEBI:15377")
    
    reaction <- "RHEA:22959"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2 <=> PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2")
    expect_equal(l$reaction_product, "PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    expect_equal(l$reaction_formula_chebi, "CHEBI:75028 + CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 <=> CHEBI:77290 + CHEBI:29034 + 2 CHEBI:15377")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:75028 + CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379")
    expect_equal(l$reaction_product_chebi, "CHEBI:77290 + CHEBI:29034 + 2 CHEBI:15377")
    
    ## lpco_to_pco
    sn1lpco <- "PC(O-16:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPCO = sn1lpco, AcylCoA = acylcoa)
    
    reaction <- "RHEA:23992"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(O-16:0/0:0) + CoA(18:0) = PC(O-16:0/18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(O-16:0/0:0) + CoA(18:0)")
    expect_equal(l$reaction_product, "PC(O-16:0/18:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:30909 + CHEBI:58342 = CHEBI:36702 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:30909 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:36702 + CHEBI:57287")
    
    reaction <- "RHEA:23993"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(O-16:0/0:0) + CoA(18:0) => PC(O-16:0/18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(O-16:0/0:0) + CoA(18:0)")
    expect_equal(l$reaction_product, "PC(O-16:0/18:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:30909 + CHEBI:58342 => CHEBI:36702 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:30909 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:36702 + CHEBI:57287")
    
    reaction <- "RHEA:23994"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(O-16:0/0:0) + CoA(18:0) <= PC(O-16:0/18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(O-16:0/0:0) + CoA(18:0)")
    expect_equal(l$reaction_product, "PC(O-16:0/18:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:30909 + CHEBI:58342 <= CHEBI:36702 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:30909 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:36702 + CHEBI:57287")
    
    reaction <- "RHEA:23995"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(O-16:0/0:0) + CoA(18:0) <=> PC(O-16:0/18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(O-16:0/0:0) + CoA(18:0)")
    expect_equal(l$reaction_product, "PC(O-16:0/18:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:30909 + CHEBI:58342 <=> CHEBI:36702 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:30909 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:36702 + CHEBI:57287")
    
    ## pa_to_dg
    pa <- "PA(14:0/16:0)"
    substrates = list(PA = pa)
    
    reaction <- "RHEA:27429"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PA(14:0/16:0) + H2O = DG(14:0/16:0/0:0) + Pi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PA(14:0/16:0) + H2O")
    expect_equal(l$reaction_product, "DG(14:0/16:0/0:0) + Pi")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58608 + CHEBI:15377 = CHEBI:17815 + CHEBI:43474")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58608 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:43474")
    
    reaction <- "RHEA:27430"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PA(14:0/16:0) + H2O => DG(14:0/16:0/0:0) + Pi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PA(14:0/16:0) + H2O")
    expect_equal(l$reaction_product, "DG(14:0/16:0/0:0) + Pi")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58608 + CHEBI:15377 => CHEBI:17815 + CHEBI:43474")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58608 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:43474")
    
    reaction <- "RHEA:27431"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PA(14:0/16:0) + H2O <= DG(14:0/16:0/0:0) + Pi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PA(14:0/16:0) + H2O")
    expect_equal(l$reaction_product, "DG(14:0/16:0/0:0) + Pi")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58608 + CHEBI:15377 <= CHEBI:17815 + CHEBI:43474")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58608 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:43474")
    
    reaction <- "RHEA:27432"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PA(14:0/16:0) + H2O <=> DG(14:0/16:0/0:0) + Pi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PA(14:0/16:0) + H2O")
    expect_equal(l$reaction_product, "DG(14:0/16:0/0:0) + Pi")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58608 + CHEBI:15377 <=> CHEBI:17815 + CHEBI:43474")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58608 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:43474")
    
    ## pe_to_ps
    pe <- "PE(14:0/16:0)"
    substrates <- list(PE = pe)
    
    reaction <- "RHEA:27606" 
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(14:0/16:0) + L-Serine = PS(14:0/16:0) + Ethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(14:0/16:0) + L-Serine")
    expect_equal(l$reaction_product, "PS(14:0/16:0) + Ethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64612 + CHEBI:33384 = CHEBI:57262 + CHEBI:57603")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64612 + CHEBI:33384")
    expect_equal(l$reaction_product_chebi, "CHEBI:57262 + CHEBI:57603")
    
    reaction <- "RHEA:27607"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(14:0/16:0) + L-Serine => PS(14:0/16:0) + Ethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(14:0/16:0) + L-Serine")
    expect_equal(l$reaction_product, "PS(14:0/16:0) + Ethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64612 + CHEBI:33384 => CHEBI:57262 + CHEBI:57603")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64612 + CHEBI:33384")
    expect_equal(l$reaction_product_chebi, "CHEBI:57262 + CHEBI:57603")
    
    reaction <- "RHEA:27608"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(14:0/16:0) + L-Serine <= PS(14:0/16:0) + Ethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(14:0/16:0) + L-Serine")
    expect_equal(l$reaction_product, "PS(14:0/16:0) + Ethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64612 + CHEBI:33384 <= CHEBI:57262 + CHEBI:57603")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64612 + CHEBI:33384")
    expect_equal(l$reaction_product_chebi, "CHEBI:57262 + CHEBI:57603")
    
    reaction <- "RHEA:27609"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(14:0/16:0) + L-Serine <=> PS(14:0/16:0) + Ethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(14:0/16:0) + L-Serine")
    expect_equal(l$reaction_product, "PS(14:0/16:0) + Ethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64612 + CHEBI:33384 <=> CHEBI:57262 + CHEBI:57603")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64612 + CHEBI:33384")
    expect_equal(l$reaction_product_chebi, "CHEBI:57262 + CHEBI:57603")
    
    ## sn2mg_to_fa
    sn2mg <- "MG(0:0/14:0/0:0)"
    substrates <- list(sn2MG = sn2mg)
    
    reaction <- "RHEA:32871"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + MG(0:0/14:0/0:0) = FA(14:0) + Glycerol + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + MG(0:0/14:0/0:0)")
    expect_equal(l$reaction_product, "FA(14:0) + Glycerol + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:17389 = CHEBI:29067 + CHEBI:17754 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:17389")
    expect_equal(l$reaction_product_chebi, "CHEBI:29067 + CHEBI:17754 + CHEBI:15378")
    
    reaction <- "RHEA:32872"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + MG(0:0/14:0/0:0) => FA(14:0) + Glycerol + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + MG(0:0/14:0/0:0)")
    expect_equal(l$reaction_product, "FA(14:0) + Glycerol + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:17389 => CHEBI:29067 + CHEBI:17754 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:17389")
    expect_equal(l$reaction_product_chebi, "CHEBI:29067 + CHEBI:17754 + CHEBI:15378")
    
    reaction <- "RHEA:32873"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + MG(0:0/14:0/0:0) <= FA(14:0) + Glycerol + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + MG(0:0/14:0/0:0)")
    expect_equal(l$reaction_product, "FA(14:0) + Glycerol + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:17389 <= CHEBI:29067 + CHEBI:17754 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:17389")
    expect_equal(l$reaction_product_chebi, "CHEBI:29067 + CHEBI:17754 + CHEBI:15378")
    
    reaction <- "RHEA:32874"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + MG(0:0/14:0/0:0) <=> FA(14:0) + Glycerol + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + MG(0:0/14:0/0:0)")
    expect_equal(l$reaction_product, "FA(14:0) + Glycerol + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:17389 <=> CHEBI:29067 + CHEBI:17754 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:17389")
    expect_equal(l$reaction_product_chebi, "CHEBI:29067 + CHEBI:17754 + CHEBI:15378")
    
    ## pg_to_cl
    pg <- "PG(18:4(6Z,9Z,12Z,15Z)/14:0)"
    cdpdg <- "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)"
    substrates = list(PG = pg, CDPDG = cdpdg)
    
    reaction <- "RHEA:32931"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PG(18:4(6Z,9Z,12Z,15Z)/14:0) + CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0) = CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0]) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PG(18:4(6Z,9Z,12Z,15Z)/14:0) + CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(l$reaction_product, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0]) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64716 + CHEBI:58332 = CHEBI:62237 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64716 + CHEBI:58332")
    expect_equal(l$reaction_product_chebi, "CHEBI:62237 + CHEBI:60377 + CHEBI:15378")
    
    reaction <- "RHEA:32932"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PG(18:4(6Z,9Z,12Z,15Z)/14:0) + CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0) => CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0]) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PG(18:4(6Z,9Z,12Z,15Z)/14:0) + CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(l$reaction_product, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0]) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64716 + CHEBI:58332 => CHEBI:62237 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64716 + CHEBI:58332")
    expect_equal(l$reaction_product_chebi, "CHEBI:62237 + CHEBI:60377 + CHEBI:15378")
    
    reaction <- "RHEA:32933"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PG(18:4(6Z,9Z,12Z,15Z)/14:0) + CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0) <= CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0]) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PG(18:4(6Z,9Z,12Z,15Z)/14:0) + CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(l$reaction_product, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0]) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64716 + CHEBI:58332 <= CHEBI:62237 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64716 + CHEBI:58332")
    expect_equal(l$reaction_product_chebi, "CHEBI:62237 + CHEBI:60377 + CHEBI:15378")
    
    reaction <- "RHEA:32934"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PG(18:4(6Z,9Z,12Z,15Z)/14:0) + CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0) <=> CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0]) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PG(18:4(6Z,9Z,12Z,15Z)/14:0) + CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(l$reaction_product, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0]) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64716 + CHEBI:58332 <=> CHEBI:62237 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64716 + CHEBI:58332")
    expect_equal(l$reaction_product_chebi, "CHEBI:62237 + CHEBI:60377 + CHEBI:15378")
    
    ## cl_to_lcl
    cl <- "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])"
    substrates <- list(CL = cl)
    
    reaction <- "RHEA:32935"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + H2O = CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + FA(18:4(6Z,9Z,12Z,15Z)) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + H2O")
    expect_equal(l$reaction_product, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + FA(18:4(6Z,9Z,12Z,15Z)) + H+") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:62237 + CHEBI:15377 = CHEBI:64743 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:62237 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64743 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:32936"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + H2O => CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + FA(18:4(6Z,9Z,12Z,15Z)) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + H2O")
    expect_equal(l$reaction_product, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + FA(18:4(6Z,9Z,12Z,15Z)) + H+") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:62237 + CHEBI:15377 => CHEBI:64743 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:62237 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64743 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:32937"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + H2O <= CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + FA(18:4(6Z,9Z,12Z,15Z)) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + H2O")
    expect_equal(l$reaction_product, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + FA(18:4(6Z,9Z,12Z,15Z)) + H+") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:62237 + CHEBI:15377 <= CHEBI:64743 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:62237 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64743 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:32938"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + H2O <=> CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + FA(18:4(6Z,9Z,12Z,15Z)) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + H2O")
    expect_equal(l$reaction_product, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + FA(18:4(6Z,9Z,12Z,15Z)) + H+") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:62237 + CHEBI:15377 <=> CHEBI:64743 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:62237 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64743 + CHEBI:28868 + CHEBI:15378")
    
    ## dg_to_pc
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    
    reaction <- "RHEA:32939"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + CDP-Choline = PC(18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + CDP-Choline")
    expect_equal(l$reaction_product, "PC(18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58779 = CHEBI:57643 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:58779")
    expect_equal(l$reaction_product_chebi, "CHEBI:57643 + CHEBI:60377 + CHEBI:15378")
    
    reaction <- "RHEA:32940"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + CDP-Choline => PC(18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + CDP-Choline")
    expect_equal(l$reaction_product, "PC(18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58779 => CHEBI:57643 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:58779")
    expect_equal(l$reaction_product_chebi, "CHEBI:57643 + CHEBI:60377 + CHEBI:15378")
    
    reaction <- "RHEA:32941"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + CDP-Choline <= PC(18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + CDP-Choline")
    expect_equal(l$reaction_product, "PC(18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58779 <= CHEBI:57643 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:58779")
    expect_equal(l$reaction_product_chebi, "CHEBI:57643 + CHEBI:60377 + CHEBI:15378")
    
    reaction <- "RHEA:32942"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + CDP-Choline <=> PC(18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + CDP-Choline")
    expect_equal(l$reaction_product, "PC(18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58779 <=> CHEBI:57643 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:58779")
    expect_equal(l$reaction_product_chebi, "CHEBI:57643 + CHEBI:60377 + CHEBI:15378")
    
    ## dg_to_pe
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    
    reaction <- "RHEA:32943"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + CDP-Ethanolamine = PE(18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + CDP-Ethanolamine")
    expect_equal(l$reaction_product, "PE(18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:57876 = CHEBI:64612 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:57876")
    expect_equal(l$reaction_product_chebi, "CHEBI:64612 + CHEBI:60377 + CHEBI:15378")
    
    reaction <- "RHEA:32944"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + CDP-Ethanolamine => PE(18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + CDP-Ethanolamine")
    expect_equal(l$reaction_product, "PE(18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:57876 => CHEBI:64612 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:57876")
    expect_equal(l$reaction_product_chebi, "CHEBI:64612 + CHEBI:60377 + CHEBI:15378")
    
    reaction <- "RHEA:32945"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + CDP-Ethanolamine <= PE(18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + CDP-Ethanolamine")
    expect_equal(l$reaction_product, "PE(18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:57876 <= CHEBI:64612 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:57876")
    expect_equal(l$reaction_product_chebi, "CHEBI:64612 + CHEBI:60377 + CHEBI:15378")
    
    reaction <- "RHEA:32946"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + CDP-Ethanolamine <=> PE(18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + CDP-Ethanolamine")
    expect_equal(l$reaction_product, "PE(18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:57876 <=> CHEBI:64612 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:57876")
    expect_equal(l$reaction_product_chebi, "CHEBI:64612 + CHEBI:60377 + CHEBI:15378")
    
    ## sn1lpe_to_fa
    pe <- "PE(14:0/0:0)"
    substrates <- list(sn1LPE = pe)
    
    reaction <- "RHEA:32967"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(14:0/0:0) + H2O = FA(14:0) + H+ + Glycerophosphoethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(14:0/0:0) + H2O")
    expect_equal(l$reaction_product, "FA(14:0) + H+ + Glycerophosphoethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64381 + CHEBI:15377 = CHEBI:28868 + CHEBI:15378 + CHEBI:143890")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64381 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:15378 + CHEBI:143890")
    
    reaction <- "RHEA:32968"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(14:0/0:0) + H2O => FA(14:0) + H+ + Glycerophosphoethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(14:0/0:0) + H2O")
    expect_equal(l$reaction_product, "FA(14:0) + H+ + Glycerophosphoethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64381 + CHEBI:15377 => CHEBI:28868 + CHEBI:15378 + CHEBI:143890")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64381 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:15378 + CHEBI:143890")
    
    reaction <- "RHEA:32969"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(14:0/0:0) + H2O <= FA(14:0) + H+ + Glycerophosphoethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(14:0/0:0) + H2O")
    expect_equal(l$reaction_product, "FA(14:0) + H+ + Glycerophosphoethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64381 + CHEBI:15377 <= CHEBI:28868 + CHEBI:15378 + CHEBI:143890")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64381 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:15378 + CHEBI:143890")
    
    reaction <- "RHEA:32970"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(14:0/0:0) + H2O <=> FA(14:0) + H+ + Glycerophosphoethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(14:0/0:0) + H2O")
    expect_equal(l$reaction_product, "FA(14:0) + H+ + Glycerophosphoethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64381 + CHEBI:15377 <=> CHEBI:28868 + CHEBI:15378 + CHEBI:143890")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64381 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:15378 + CHEBI:143890")
    
    ## sn1lpe_to_pe
    pe <- "PE(14:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPE = pe, AcylCoA = acylcoa)
    
    reaction <- "RHEA:32995"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(14:0/0:0) + CoA(18:0) = PE(14:0/18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(14:0/0:0) + CoA(18:0)")
    expect_equal(l$reaction_product, "PE(14:0/18:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64381 + CHEBI:58342 = CHEBI:64612 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64381 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:64612 + CHEBI:57287")
    
    reaction <- "RHEA:32996"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(14:0/0:0) + CoA(18:0) => PE(14:0/18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(14:0/0:0) + CoA(18:0)")
    expect_equal(l$reaction_product, "PE(14:0/18:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64381 + CHEBI:58342 => CHEBI:64612 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64381 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:64612 + CHEBI:57287")
    
    reaction <- "RHEA:32997"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(14:0/0:0) + CoA(18:0) <= PE(14:0/18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(14:0/0:0) + CoA(18:0)")
    expect_equal(l$reaction_product, "PE(14:0/18:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64381 + CHEBI:58342 <= CHEBI:64612 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64381 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:64612 + CHEBI:57287")
    
    reaction <- "RHEA:32998"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(14:0/0:0) + CoA(18:0) <=> PE(14:0/18:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(14:0/0:0) + CoA(18:0)")
    expect_equal(l$reaction_product, "PE(14:0/18:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64381 + CHEBI:58342 <=> CHEBI:64612 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64381 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:64612 + CHEBI:57287")
    
    ## nape_to_nae
    nape <- "NAPE(14:0/16:0/18:0)"
    substrates = list(NAPE = nape)
    
    reaction <- "RHEA:33159"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + NAPE(14:0/16:0/18:0) = PA(14:0/16:0) + NAE(18:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + NAPE(14:0/16:0/18:0)")
    expect_equal(l$reaction_product, "PA(14:0/16:0) + NAE(18:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 = CHEBI:58608 + CHEBI:52640 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:62537")
    expect_equal(l$reaction_product_chebi, "CHEBI:58608 + CHEBI:52640 + CHEBI:15378")
    
    reaction <- "RHEA:33160"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + NAPE(14:0/16:0/18:0) => PA(14:0/16:0) + NAE(18:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + NAPE(14:0/16:0/18:0)")
    expect_equal(l$reaction_product, "PA(14:0/16:0) + NAE(18:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 => CHEBI:58608 + CHEBI:52640 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:62537")
    expect_equal(l$reaction_product_chebi, "CHEBI:58608 + CHEBI:52640 + CHEBI:15378")
    
    reaction <- "RHEA:33161"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + NAPE(14:0/16:0/18:0) <= PA(14:0/16:0) + NAE(18:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + NAPE(14:0/16:0/18:0)")
    expect_equal(l$reaction_product, "PA(14:0/16:0) + NAE(18:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 <= CHEBI:58608 + CHEBI:52640 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:62537")
    expect_equal(l$reaction_product_chebi, "CHEBI:58608 + CHEBI:52640 + CHEBI:15378")
    
    reaction <- "RHEA:33162"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + NAPE(14:0/16:0/18:0) <=> PA(14:0/16:0) + NAE(18:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + NAPE(14:0/16:0/18:0)")
    expect_equal(l$reaction_product, "PA(14:0/16:0) + NAE(18:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 <=> CHEBI:58608 + CHEBI:52640 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:62537")
    expect_equal(l$reaction_product_chebi, "CHEBI:58608 + CHEBI:52640 + CHEBI:15378")
    
    ## sn1lpi_to_pi
    sn1lpi <- "PI(16:0/0:0)"
    acylcoa <- "CoA(18:1(9Z))"
    substrates = list(sn1LPI = sn1lpi, AcylCoA = acylcoa)
    
    reaction <- "RHEA:33195"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PI(16:0/0:0) + CoA(18:1(9Z)) = PI(16:0/18:1(9Z)) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PI(16:0/0:0) + CoA(18:1(9Z))")
    expect_equal(l$reaction_product, "PI(16:0/18:1(9Z)) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64771 + CHEBI:58342 = CHEBI:57880 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64771 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:57880 + CHEBI:57287")
    
    reaction <- "RHEA:33196"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PI(16:0/0:0) + CoA(18:1(9Z)) => PI(16:0/18:1(9Z)) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PI(16:0/0:0) + CoA(18:1(9Z))")
    expect_equal(l$reaction_product, "PI(16:0/18:1(9Z)) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64771 + CHEBI:58342 => CHEBI:57880 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64771 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:57880 + CHEBI:57287")
    
    reaction <- "RHEA:33197"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PI(16:0/0:0) + CoA(18:1(9Z)) <= PI(16:0/18:1(9Z)) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PI(16:0/0:0) + CoA(18:1(9Z))")
    expect_equal(l$reaction_product, "PI(16:0/18:1(9Z)) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64771 + CHEBI:58342 <= CHEBI:57880 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64771 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:57880 + CHEBI:57287")
    
    reaction <- "RHEA:33198"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PI(16:0/0:0) + CoA(18:1(9Z)) <=> PI(16:0/18:1(9Z)) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PI(16:0/0:0) + CoA(18:1(9Z))")
    expect_equal(l$reaction_product, "PI(16:0/18:1(9Z)) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64771 + CHEBI:58342 <=> CHEBI:57880 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64771 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:57880 + CHEBI:57287")
    
    ## sn1lpg_to_pg
    sn1lpg <- "PG(16:0/0:0)"
    acylcoa <- "CoA(14:0)"
    substrates = list(sn1LPG = sn1lpg, AcylCoA = acylcoa)
    
    reaction <- "RHEA:33203"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PG(16:0/0:0) + CoA(14:0) = PG(16:0/14:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PG(16:0/0:0) + CoA(14:0)")
    expect_equal(l$reaction_product, "PG(16:0/14:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64840 + CHEBI:58342 = CHEBI:64716 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64840 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:64716 + CHEBI:57287")
    
    reaction <- "RHEA:33204"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PG(16:0/0:0) + CoA(14:0) => PG(16:0/14:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PG(16:0/0:0) + CoA(14:0)")
    expect_equal(l$reaction_product, "PG(16:0/14:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64840 + CHEBI:58342 => CHEBI:64716 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64840 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:64716 + CHEBI:57287")
    
    reaction <- "RHEA:33205"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PG(16:0/0:0) + CoA(14:0) <= PG(16:0/14:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PG(16:0/0:0) + CoA(14:0)")
    expect_equal(l$reaction_product, "PG(16:0/14:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64840 + CHEBI:58342 <= CHEBI:64716 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64840 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:64716 + CHEBI:57287")
    
    reaction <- "RHEA:33206"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PG(16:0/0:0) + CoA(14:0) <=> PG(16:0/14:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PG(16:0/0:0) + CoA(14:0)")
    expect_equal(l$reaction_product, "PG(16:0/14:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64840 + CHEBI:58342 <=> CHEBI:64716 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64840 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:64716 + CHEBI:57287")
    
    
    ## tg_to_dg
    tg <- "TG(18:0/16:0/14:0)"
    substrates = list(TG = tg)
    
    reaction <- "RHEA:33271"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, c("", ""))
    expect_equal(l$reaction_formula, 
        c("TG(18:0/16:0/14:0) + H2O = DG(14:0/16:0/0:0) + FA(18:0) + H+", 
          "TG(18:0/16:0/14:0) + H2O = DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(l$reaction_isReversible, c("", ""))
    expect_equal(l$reaction_geneAssociation, c("", ""))
    expect_equal(l$reaction_pathway, c("", ""))
    expect_equal(l$reaction_substrate, 
        c("TG(18:0/16:0/14:0) + H2O", "TG(18:0/16:0/14:0) + H2O"))
    expect_equal(l$reaction_product, 
        c("DG(14:0/16:0/0:0) + FA(18:0) + H+", "DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(l$reaction_formula_chebi, 
        c("CHEBI:64615 + CHEBI:15377 = CHEBI:17815 + CHEBI:28868 + CHEBI:15378", 
          "CHEBI:64615 + CHEBI:15377 = CHEBI:17815 + CHEBI:28868 + CHEBI:15378"))
    expect_equal(l$reaction_substrate_chebi, 
        c("CHEBI:64615 + CHEBI:15377", "CHEBI:64615 + CHEBI:15377"))
    expect_equal(l$reaction_product_chebi, 
        c("CHEBI:17815 + CHEBI:28868 + CHEBI:15378", 
          "CHEBI:17815 + CHEBI:28868 + CHEBI:15378"))
    
    reaction <- "RHEA:33272"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, c("", ""))
    expect_equal(l$reaction_formula, 
        c("TG(18:0/16:0/14:0) + H2O => DG(14:0/16:0/0:0) + FA(18:0) + H+", 
          "TG(18:0/16:0/14:0) + H2O => DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(l$reaction_isReversible, c("", ""))
    expect_equal(l$reaction_geneAssociation, c("", ""))
    expect_equal(l$reaction_pathway, c("", ""))
    expect_equal(l$reaction_substrate, 
        c("TG(18:0/16:0/14:0) + H2O", "TG(18:0/16:0/14:0) + H2O"))
    expect_equal(l$reaction_product, 
        c("DG(14:0/16:0/0:0) + FA(18:0) + H+", "DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(l$reaction_formula_chebi, 
        c("CHEBI:64615 + CHEBI:15377 => CHEBI:17815 + CHEBI:28868 + CHEBI:15378", 
          "CHEBI:64615 + CHEBI:15377 => CHEBI:17815 + CHEBI:28868 + CHEBI:15378"))
    expect_equal(l$reaction_substrate_chebi, 
        c("CHEBI:64615 + CHEBI:15377", "CHEBI:64615 + CHEBI:15377"))
    expect_equal(l$reaction_product_chebi, 
        c("CHEBI:17815 + CHEBI:28868 + CHEBI:15378", 
          "CHEBI:17815 + CHEBI:28868 + CHEBI:15378"))
    
    reaction <- "RHEA:33273"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, c("", ""))
    expect_equal(l$reaction_formula, 
        c("TG(18:0/16:0/14:0) + H2O <= DG(14:0/16:0/0:0) + FA(18:0) + H+", 
          "TG(18:0/16:0/14:0) + H2O <= DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(l$reaction_isReversible, c("", ""))
    expect_equal(l$reaction_geneAssociation, c("", ""))
    expect_equal(l$reaction_pathway, c("", ""))
    expect_equal(l$reaction_substrate, 
        c("TG(18:0/16:0/14:0) + H2O", "TG(18:0/16:0/14:0) + H2O"))
    expect_equal(l$reaction_product, 
        c("DG(14:0/16:0/0:0) + FA(18:0) + H+", "DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(l$reaction_formula_chebi, 
        c("CHEBI:64615 + CHEBI:15377 <= CHEBI:17815 + CHEBI:28868 + CHEBI:15378", 
          "CHEBI:64615 + CHEBI:15377 <= CHEBI:17815 + CHEBI:28868 + CHEBI:15378"))
    expect_equal(l$reaction_substrate_chebi, 
        c("CHEBI:64615 + CHEBI:15377", "CHEBI:64615 + CHEBI:15377"))
    expect_equal(l$reaction_product_chebi, 
        c("CHEBI:17815 + CHEBI:28868 + CHEBI:15378", 
          "CHEBI:17815 + CHEBI:28868 + CHEBI:15378"))
    
    reaction <- "RHEA:33274"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, c("", ""))
    expect_equal(l$reaction_formula, 
        c("TG(18:0/16:0/14:0) + H2O <=> DG(14:0/16:0/0:0) + FA(18:0) + H+", 
          "TG(18:0/16:0/14:0) + H2O <=> DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(l$reaction_isReversible, c("", ""))
    expect_equal(l$reaction_geneAssociation, c("", ""))
    expect_equal(l$reaction_pathway, c("", ""))
    expect_equal(l$reaction_substrate, 
        c("TG(18:0/16:0/14:0) + H2O", "TG(18:0/16:0/14:0) + H2O"))
    expect_equal(l$reaction_product, 
        c("DG(14:0/16:0/0:0) + FA(18:0) + H+", "DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(l$reaction_formula_chebi, 
        c("CHEBI:64615 + CHEBI:15377 <=> CHEBI:17815 + CHEBI:28868 + CHEBI:15378", 
          "CHEBI:64615 + CHEBI:15377 <=> CHEBI:17815 + CHEBI:28868 + CHEBI:15378"))
    expect_equal(l$reaction_substrate_chebi, 
        c("CHEBI:64615 + CHEBI:15377", "CHEBI:64615 + CHEBI:15377"))
    expect_equal(l$reaction_product_chebi, 
        c("CHEBI:17815 + CHEBI:28868 + CHEBI:15378", 
          "CHEBI:17815 + CHEBI:28868 + CHEBI:15378"))
    
    reaction <- "RHEA:44864"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, c("", ""))
    expect_equal(l$reaction_formula, 
        c("TG(18:0/16:0/14:0) + H2O = DG(14:0/16:0/0:0) + FA(18:0) + H+", 
          "TG(18:0/16:0/14:0) + H2O = DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(l$reaction_isReversible, c("", ""))
    expect_equal(l$reaction_geneAssociation, c("", ""))
    expect_equal(l$reaction_pathway, c("", ""))
    expect_equal(l$reaction_substrate, 
        c("TG(18:0/16:0/14:0) + H2O", "TG(18:0/16:0/14:0) + H2O"))
    expect_equal(l$reaction_product, 
        c("DG(14:0/16:0/0:0) + FA(18:0) + H+", "DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(l$reaction_formula_chebi, 
        c("CHEBI:64615 + CHEBI:15377 = CHEBI:49172 + CHEBI:28868 + CHEBI:15378", 
          "CHEBI:64615 + CHEBI:15377 = CHEBI:49172 + CHEBI:28868 + CHEBI:15378"))
    expect_equal(l$reaction_substrate_chebi, 
        c("CHEBI:64615 + CHEBI:15377", "CHEBI:64615 + CHEBI:15377"))
    expect_equal(l$reaction_product_chebi, 
        c("CHEBI:49172 + CHEBI:28868 + CHEBI:15378", 
          "CHEBI:49172 + CHEBI:28868 + CHEBI:15378"))
    
    reaction <- "RHEA:44865"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, c("", ""))
    expect_equal(l$reaction_formula, 
        c("TG(18:0/16:0/14:0) + H2O => DG(14:0/16:0/0:0) + FA(18:0) + H+", 
          "TG(18:0/16:0/14:0) + H2O => DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(l$reaction_isReversible, c("", ""))
    expect_equal(l$reaction_geneAssociation, c("", ""))
    expect_equal(l$reaction_pathway, c("", ""))
    expect_equal(l$reaction_substrate, 
        c("TG(18:0/16:0/14:0) + H2O", "TG(18:0/16:0/14:0) + H2O"))
    expect_equal(l$reaction_product, 
        c("DG(14:0/16:0/0:0) + FA(18:0) + H+", "DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(l$reaction_formula_chebi, 
        c("CHEBI:64615 + CHEBI:15377 => CHEBI:49172 + CHEBI:28868 + CHEBI:15378", 
          "CHEBI:64615 + CHEBI:15377 => CHEBI:49172 + CHEBI:28868 + CHEBI:15378"))
    expect_equal(l$reaction_substrate_chebi, 
        c("CHEBI:64615 + CHEBI:15377", "CHEBI:64615 + CHEBI:15377"))
    expect_equal(l$reaction_product_chebi, 
        c("CHEBI:49172 + CHEBI:28868 + CHEBI:15378", 
          "CHEBI:49172 + CHEBI:28868 + CHEBI:15378"))
    
    reaction <- "RHEA:44866"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, c("", ""))
    expect_equal(l$reaction_formula, 
        c("TG(18:0/16:0/14:0) + H2O <= DG(14:0/16:0/0:0) + FA(18:0) + H+", 
          "TG(18:0/16:0/14:0) + H2O <= DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(l$reaction_isReversible, c("", ""))
    expect_equal(l$reaction_geneAssociation, c("", ""))
    expect_equal(l$reaction_pathway, c("", ""))
    expect_equal(l$reaction_substrate, 
        c("TG(18:0/16:0/14:0) + H2O", "TG(18:0/16:0/14:0) + H2O"))
    expect_equal(l$reaction_product, 
        c("DG(14:0/16:0/0:0) + FA(18:0) + H+", "DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(l$reaction_formula_chebi, 
        c("CHEBI:64615 + CHEBI:15377 <= CHEBI:49172 + CHEBI:28868 + CHEBI:15378", 
          "CHEBI:64615 + CHEBI:15377 <= CHEBI:49172 + CHEBI:28868 + CHEBI:15378"))
    expect_equal(l$reaction_substrate_chebi, 
        c("CHEBI:64615 + CHEBI:15377", "CHEBI:64615 + CHEBI:15377"))
    expect_equal(l$reaction_product_chebi, 
        c("CHEBI:49172 + CHEBI:28868 + CHEBI:15378", 
          "CHEBI:49172 + CHEBI:28868 + CHEBI:15378"))
    
    reaction <- "RHEA:44867"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, c("", ""))
    expect_equal(l$reaction_formula, 
        c("TG(18:0/16:0/14:0) + H2O <=> DG(14:0/16:0/0:0) + FA(18:0) + H+", 
          "TG(18:0/16:0/14:0) + H2O <=> DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(l$reaction_isReversible, c("", ""))
    expect_equal(l$reaction_geneAssociation, c("", ""))
    expect_equal(l$reaction_pathway, c("", ""))
    expect_equal(l$reaction_substrate, 
        c("TG(18:0/16:0/14:0) + H2O", "TG(18:0/16:0/14:0) + H2O"))
    expect_equal(l$reaction_product, 
        c("DG(14:0/16:0/0:0) + FA(18:0) + H+", "DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(l$reaction_formula_chebi, 
        c("CHEBI:64615 + CHEBI:15377 <=> CHEBI:49172 + CHEBI:28868 + CHEBI:15378", 
          "CHEBI:64615 + CHEBI:15377 <=> CHEBI:49172 + CHEBI:28868 + CHEBI:15378"))
    expect_equal(l$reaction_substrate_chebi, 
        c("CHEBI:64615 + CHEBI:15377", "CHEBI:64615 + CHEBI:15377"))
    expect_equal(l$reaction_product_chebi, 
        c("CHEBI:49172 + CHEBI:28868 + CHEBI:15378", 
          "CHEBI:49172 + CHEBI:28868 + CHEBI:15378"))
    
    ## dg_to_sn2mg
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    
    reaction <- "RHEA:33275"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O = MG(0:0/16:0/0:0) + FA(18:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "MG(0:0/16:0/0:0) + FA(18:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:15377 = CHEBI:17389 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:17389 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:33276"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O => MG(0:0/16:0/0:0) + FA(18:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "MG(0:0/16:0/0:0) + FA(18:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:15377 => CHEBI:17389 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:17389 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:33277"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O <= MG(0:0/16:0/0:0) + FA(18:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "MG(0:0/16:0/0:0) + FA(18:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:15377 <= CHEBI:17389 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:17389 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:33278"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O <=> MG(0:0/16:0/0:0) + FA(18:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "MG(0:0/16:0/0:0) + FA(18:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:15377 <=> CHEBI:17389 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:17389 + CHEBI:28868 + CHEBI:15378")
    
    ## cerp_to_cer
    cerp <- "CerP(d16:1(3OH,4OH)(15Me)/18:0)"
    substrates <- list(CerP = cerp)
    
    reaction <- "RHEA:33743"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "H2O + CerP(d16:1(3OH,4OH)(15Me)/18:0) = Pi + Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + CerP(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_product, "Pi + Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:57674 = CHEBI:43474 + CHEBI:52639")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:57674")
    expect_equal(l$reaction_product_chebi, "CHEBI:43474 + CHEBI:52639")
    
    reaction <- "RHEA:33744"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "H2O + CerP(d16:1(3OH,4OH)(15Me)/18:0) => Pi + Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + CerP(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_product, "Pi + Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:57674 => CHEBI:43474 + CHEBI:52639")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:57674")
    expect_equal(l$reaction_product_chebi, "CHEBI:43474 + CHEBI:52639")
    
    reaction <- "RHEA:33745"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "H2O + CerP(d16:1(3OH,4OH)(15Me)/18:0) <= Pi + Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + CerP(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_product, "Pi + Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:57674 <= CHEBI:43474 + CHEBI:52639")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:57674")
    expect_equal(l$reaction_product_chebi, "CHEBI:43474 + CHEBI:52639")
    
    reaction <- "RHEA:33746"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "H2O + CerP(d16:1(3OH,4OH)(15Me)/18:0) <=> Pi + Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + CerP(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_product, "Pi + Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:57674 <=> CHEBI:43474 + CHEBI:52639")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:57674")
    expect_equal(l$reaction_product_chebi, "CHEBI:43474 + CHEBI:52639")
    
    ## sn1mg_to_lpa
    sn1mg <- "MG(14:0/0:0/0:0)"
    substrates <- list(sn1MG = sn1mg)
    
    reaction <- "RHEA:33747"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(14:0/0:0/0:0) + ATP = PA(14:0/0:0) + ADP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(14:0/0:0/0:0) + ATP")
    expect_equal(l$reaction_product, "PA(14:0/0:0) + ADP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64683 + CHEBI:30616 = CHEBI:57970 + CHEBI:456216 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64683 + CHEBI:30616")
    expect_equal(l$reaction_product_chebi, "CHEBI:57970 + CHEBI:456216 + CHEBI:15378")
    
    reaction <- "RHEA:33748"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(14:0/0:0/0:0) + ATP => PA(14:0/0:0) + ADP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(14:0/0:0/0:0) + ATP")
    expect_equal(l$reaction_product, "PA(14:0/0:0) + ADP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64683 + CHEBI:30616 => CHEBI:57970 + CHEBI:456216 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64683 + CHEBI:30616")
    expect_equal(l$reaction_product_chebi, "CHEBI:57970 + CHEBI:456216 + CHEBI:15378")
    
    reaction <- "RHEA:33749"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(14:0/0:0/0:0) + ATP <= PA(14:0/0:0) + ADP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(14:0/0:0/0:0) + ATP")
    expect_equal(l$reaction_product, "PA(14:0/0:0) + ADP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64683 + CHEBI:30616 <= CHEBI:57970 + CHEBI:456216 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64683 + CHEBI:30616")
    expect_equal(l$reaction_product_chebi, "CHEBI:57970 + CHEBI:456216 + CHEBI:15378")
    
    reaction <- "RHEA:33750"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(14:0/0:0/0:0) + ATP <=> PA(14:0/0:0) + ADP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(14:0/0:0/0:0) + ATP")
    expect_equal(l$reaction_product, "PA(14:0/0:0) + ADP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64683 + CHEBI:30616 <=> CHEBI:57970 + CHEBI:456216 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64683 + CHEBI:30616")
    expect_equal(l$reaction_product_chebi, "CHEBI:57970 + CHEBI:456216 + CHEBI:15378")
    
    ## pgp_to_pg
    pgp <- "PGP(16:0/14:0)"
    substrates <- list(PGP = pgp)
    
    reaction <- "RHEA:33751"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PGP(16:0/14:0) + H2O = PG(16:0/14:0) + Pi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PGP(16:0/14:0) + H2O")
    expect_equal(l$reaction_product, "PG(16:0/14:0) + Pi")
    expect_equal(l$reaction_formula_chebi, "CHEBI:60110 + CHEBI:15377 = CHEBI:64716 + CHEBI:43474")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:60110 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64716 + CHEBI:43474")
    
    reaction <- "RHEA:33752"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PGP(16:0/14:0) + H2O => PG(16:0/14:0) + Pi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PGP(16:0/14:0) + H2O")
    expect_equal(l$reaction_product, "PG(16:0/14:0) + Pi")
    expect_equal(l$reaction_formula_chebi, "CHEBI:60110 + CHEBI:15377 => CHEBI:64716 + CHEBI:43474")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:60110 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64716 + CHEBI:43474")
    
    reaction <- "RHEA:33753"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PGP(16:0/14:0) + H2O <= PG(16:0/14:0) + Pi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PGP(16:0/14:0) + H2O")
    expect_equal(l$reaction_product, "PG(16:0/14:0) + Pi")
    expect_equal(l$reaction_formula_chebi, "CHEBI:60110 + CHEBI:15377 <= CHEBI:64716 + CHEBI:43474")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:60110 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64716 + CHEBI:43474")
    
    reaction <- "RHEA:33754"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PGP(16:0/14:0) + H2O <=> PG(16:0/14:0) + Pi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PGP(16:0/14:0) + H2O")
    expect_equal(l$reaction_product, "PG(16:0/14:0) + Pi")
    expect_equal(l$reaction_formula_chebi, "CHEBI:60110 + CHEBI:15377 <=> CHEBI:64716 + CHEBI:43474")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:60110 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64716 + CHEBI:43474")
    
    ## sn1mg_to_fa
    sn1mg <- "MG(14:0/0:0/0:0)"
    substrates <- list(sn1MG = sn1mg)
    
    reaction <- "RHEA:34019"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(14:0/0:0/0:0) + H2O = FA(14:0) + Glycerol + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(14:0/0:0/0:0) + H2O")
    expect_equal(l$reaction_product, "FA(14:0) + Glycerol + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:35759 + CHEBI:15377 = CHEBI:28868 + CHEBI:17754 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:35759 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:17754 + CHEBI:15378")
    
    reaction <- "RHEA:34020"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(14:0/0:0/0:0) + H2O => FA(14:0) + Glycerol + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(14:0/0:0/0:0) + H2O")
    expect_equal(l$reaction_product, "FA(14:0) + Glycerol + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:35759 + CHEBI:15377 => CHEBI:28868 + CHEBI:17754 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:35759 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:17754 + CHEBI:15378")
    
    reaction <- "RHEA:34021"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(14:0/0:0/0:0) + H2O <= FA(14:0) + Glycerol + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(14:0/0:0/0:0) + H2O")
    expect_equal(l$reaction_product, "FA(14:0) + Glycerol + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:35759 + CHEBI:15377 <= CHEBI:28868 + CHEBI:17754 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:35759 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:17754 + CHEBI:15378")
    
    reaction <- "RHEA:34022"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(14:0/0:0/0:0) + H2O <=> FA(14:0) + Glycerol + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(14:0/0:0/0:0) + H2O")
    expect_equal(l$reaction_product, "FA(14:0) + Glycerol + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:35759 + CHEBI:15377 <=> CHEBI:28868 + CHEBI:17754 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:35759 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:17754 + CHEBI:15378")
    
    ## dg_to_sn1mg
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    
    reaction <- "RHEA:35663"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O = MG(18:0/0:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "MG(18:0/0:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:15377 = CHEBI:64683 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64683 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:35664"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O => MG(18:0/0:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "MG(18:0/0:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:15377 => CHEBI:64683 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64683 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:35665"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O <= MG(18:0/0:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "MG(18:0/0:0/0:0) + FA(16:0) + H+") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:15377 <= CHEBI:64683 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64683 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:35666"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O <=> MG(18:0/0:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "MG(18:0/0:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:17815 + CHEBI:15377 <=> CHEBI:64683 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:17815 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64683 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:44712"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O = MG(18:0/0:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "MG(18:0/0:0/0:0) + FA(16:0) + H+") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:49172 + CHEBI:15377 = CHEBI:35759 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:49172 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:35759 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:44713"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O => MG(18:0/0:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "MG(18:0/0:0/0:0) + FA(16:0) + H+") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:49172 + CHEBI:15377 => CHEBI:35759 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:49172 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:35759 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:44714"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O <= MG(18:0/0:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "MG(18:0/0:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:49172 + CHEBI:15377 <= CHEBI:35759 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:49172 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:35759 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:44715"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O <=> MG(18:0/0:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(18:0/16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "MG(18:0/0:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:49172 + CHEBI:15377 <=> CHEBI:35759 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:49172 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:35759 + CHEBI:28868 + CHEBI:15378")
    
    ## lcl_to_cl
    lcl <- "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])" ###################
    acylcoa <- "CoA(18:4(6Z,9Z,12Z,15Z))"
    substrates <- list(LCL = lcl, AcylCoA = acylcoa)
    
    reaction <- "RHEA:35839"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z)) = CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(l$reaction_product, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64743 + CHEBI:58342 = CHEBI:62237 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64743 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:62237 + CHEBI:57287")
    
    reaction <- "RHEA:35840"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z)) => CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(l$reaction_product, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64743 + CHEBI:58342 => CHEBI:62237 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64743 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:62237 + CHEBI:57287")
    
    reaction <- "RHEA:35841"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z)) <= CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(l$reaction_product, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64743 + CHEBI:58342 <= CHEBI:62237 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64743 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:62237 + CHEBI:57287")
    
    reaction <- "RHEA:35842"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z)) <=> CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(l$reaction_product, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64743 + CHEBI:58342 <=> CHEBI:62237 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64743 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:62237 + CHEBI:57287")
    
    ## lpco_to_mgo
    sn1lpco <- "PC(O-16:0/0:0)"
    substrates <- list(sn1LPCO = sn1lpco)
    
    reaction <- "RHEA:36083"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(O-16:0/0:0) + H2O = MG(O-16:0/0:0/0:0) + H+ + Phosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(O-16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "MG(O-16:0/0:0/0:0) + H+ + Phosphocholine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 = CHEBI:15850 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:30909 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:15850 + CHEBI:15378 + CHEBI:295975")
    
    reaction <- "RHEA:36084"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(O-16:0/0:0) + H2O => MG(O-16:0/0:0/0:0) + H+ + Phosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(O-16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "MG(O-16:0/0:0/0:0) + H+ + Phosphocholine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 => CHEBI:15850 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:30909 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:15850 + CHEBI:15378 + CHEBI:295975")
    
    reaction <- "RHEA:36085"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(O-16:0/0:0) + H2O <= MG(O-16:0/0:0/0:0) + H+ + Phosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(O-16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "MG(O-16:0/0:0/0:0) + H+ + Phosphocholine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 <= CHEBI:15850 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:30909 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:15850 + CHEBI:15378 + CHEBI:295975")
    
    reaction <- "RHEA:36086"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(O-16:0/0:0) + H2O <=> MG(O-16:0/0:0/0:0) + H+ + Phosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(O-16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "MG(O-16:0/0:0/0:0) + H+ + Phosphocholine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 <=> CHEBI:15850 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:30909 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:15850 + CHEBI:15378 + CHEBI:295975")
    
    ## acyldhap_to_alkyldhap
    acyldhap <- "DHAP(18:0)"
    fao <- "FAO(16:0)"
    substrates <- list(AcylDHAP = acyldhap, FAO = fao)
    
    reaction <- "RHEA:36171"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "DHAP(18:0) + FAO(16:0) = DHAP(O-16:0) + FA(18:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DHAP(18:0) + FAO(16:0)")
    expect_equal(l$reaction_product, "DHAP(O-16:0) + FA(18:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57534 + CHEBI:17135 = CHEBI:73315 + CHEBI:57560 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57534 + CHEBI:17135")
    expect_equal(l$reaction_product_chebi, "CHEBI:73315 + CHEBI:57560 + CHEBI:15378")
    
    reaction <- "RHEA:36172"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "DHAP(18:0) + FAO(16:0) => DHAP(O-16:0) + FA(18:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DHAP(18:0) + FAO(16:0)")
    expect_equal(l$reaction_product, "DHAP(O-16:0) + FA(18:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57534 + CHEBI:17135 => CHEBI:73315 + CHEBI:57560 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57534 + CHEBI:17135")
    expect_equal(l$reaction_product_chebi, "CHEBI:73315 + CHEBI:57560 + CHEBI:15378")
    
    reaction <- "RHEA:36173"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "DHAP(18:0) + FAO(16:0) <= DHAP(O-16:0) + FA(18:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DHAP(18:0) + FAO(16:0)")
    expect_equal(l$reaction_product, "DHAP(O-16:0) + FA(18:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57534 + CHEBI:17135 <= CHEBI:73315 + CHEBI:57560 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57534 + CHEBI:17135")
    expect_equal(l$reaction_product_chebi, "CHEBI:73315 + CHEBI:57560 + CHEBI:15378")
    
    reaction <- "RHEA:36174"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "DHAP(18:0) + FAO(16:0) <=> DHAP(O-16:0) + FA(18:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DHAP(18:0) + FAO(16:0)")
    expect_equal(l$reaction_product, "DHAP(O-16:0) + FA(18:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57534 + CHEBI:17135 <=> CHEBI:73315 + CHEBI:57560 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57534 + CHEBI:17135")
    expect_equal(l$reaction_product_chebi, "CHEBI:73315 + CHEBI:57560 + CHEBI:15378")
    
    ## alkyldhap_to_lpao
    alkyldhap <- "DHAP(O-18:0)"
    substrates = list(AlkylDHAP = alkyldhap)
    
    reaction <- "RHEA:36175"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "DHAP(O-18:0) + H+ + NADPH = PA(O-18:0/0:0) + NADP+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DHAP(O-18:0) + H+ + NADPH")
    expect_equal(l$reaction_product, "PA(O-18:0/0:0) + NADP+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:73315 + CHEBI:15378 + CHEBI:57783 = CHEBI:58014 + CHEBI:58349")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:73315 + CHEBI:15378 + CHEBI:57783")
    expect_equal(l$reaction_product_chebi, "CHEBI:58014 + CHEBI:58349")
    
    reaction <- "RHEA:36176"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "DHAP(O-18:0) + H+ + NADPH => PA(O-18:0/0:0) + NADP+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DHAP(O-18:0) + H+ + NADPH")
    expect_equal(l$reaction_product, "PA(O-18:0/0:0) + NADP+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:73315 + CHEBI:15378 + CHEBI:57783 => CHEBI:58014 + CHEBI:58349")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:73315 + CHEBI:15378 + CHEBI:57783")
    expect_equal(l$reaction_product_chebi, "CHEBI:58014 + CHEBI:58349")
    
    reaction <- "RHEA:36177"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "DHAP(O-18:0) + H+ + NADPH <= PA(O-18:0/0:0) + NADP+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DHAP(O-18:0) + H+ + NADPH")
    expect_equal(l$reaction_product, "PA(O-18:0/0:0) + NADP+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:73315 + CHEBI:15378 + CHEBI:57783 <= CHEBI:58014 + CHEBI:58349")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:73315 + CHEBI:15378 + CHEBI:57783")
    expect_equal(l$reaction_product_chebi, "CHEBI:58014 + CHEBI:58349")
    
    reaction <- "RHEA:36178"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, "DHAP(O-18:0) + H+ + NADPH <=> PA(O-18:0/0:0) + NADP+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DHAP(O-18:0) + H+ + NADPH")
    expect_equal(l$reaction_product, "PA(O-18:0/0:0) + NADP+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:73315 + CHEBI:15378 + CHEBI:57783 <=> CHEBI:58014 + CHEBI:58349")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:73315 + CHEBI:15378 + CHEBI:57783")
    expect_equal(l$reaction_product_chebi, "CHEBI:58014 + CHEBI:58349")
    
    ## dgo_to_pco
    dgo <- "DG(O-18:0/16:0/0:0)"
    substrates <- list(DGO = dgo)
    
    reaction <- "RHEA:36179"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(O-18:0/16:0/0:0) + CDP-Choline = PC(O-18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(O-18:0/16:0/0:0) + CDP-Choline")
    expect_equal(l$reaction_product, "PC(O-18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:52595 + CHEBI:58779 = CHEBI:36702 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:52595 + CHEBI:58779")
    expect_equal(l$reaction_product_chebi, "CHEBI:36702 + CHEBI:60377 + CHEBI:15378")
    
    reaction <- "RHEA:36180"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(O-18:0/16:0/0:0) + CDP-Choline => PC(O-18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(O-18:0/16:0/0:0) + CDP-Choline")
    expect_equal(l$reaction_product, "PC(O-18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:52595 + CHEBI:58779 => CHEBI:36702 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:52595 + CHEBI:58779")
    expect_equal(l$reaction_product_chebi, "CHEBI:36702 + CHEBI:60377 + CHEBI:15378")
    
    reaction <- "RHEA:36181"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(O-18:0/16:0/0:0) + CDP-Choline <= PC(O-18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(O-18:0/16:0/0:0) + CDP-Choline")
    expect_equal(l$reaction_product, "PC(O-18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:52595 + CHEBI:58779 <= CHEBI:36702 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:52595 + CHEBI:58779")
    expect_equal(l$reaction_product_chebi, "CHEBI:36702 + CHEBI:60377 + CHEBI:15378")
    
    reaction <- "RHEA:36182"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(O-18:0/16:0/0:0) + CDP-Choline <=> PC(O-18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(O-18:0/16:0/0:0) + CDP-Choline")
    expect_equal(l$reaction_product, "PC(O-18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:52595 + CHEBI:58779 <=> CHEBI:36702 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:52595 + CHEBI:58779")
    expect_equal(l$reaction_product_chebi, "CHEBI:36702 + CHEBI:60377 + CHEBI:15378")
    
    ## dgo_to_peo
    dgo <- "DG(O-18:0/16:0/0:0)"
    substrates <- list(DGO = dgo)
    
    reaction <- "RHEA:36187"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(O-18:0/16:0/0:0) + CDP-Ethanolamine = PE(O-18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(O-18:0/16:0/0:0) + CDP-Ethanolamine")
    expect_equal(l$reaction_product, "PE(O-18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:52595 + CHEBI:57876 = CHEBI:60520 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:52595 + CHEBI:57876")
    expect_equal(l$reaction_product_chebi, "CHEBI:60520 + CHEBI:60377 + CHEBI:15378")
    
    reaction <- "RHEA:36188"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(O-18:0/16:0/0:0) + CDP-Ethanolamine => PE(O-18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(O-18:0/16:0/0:0) + CDP-Ethanolamine")
    expect_equal(l$reaction_product, "PE(O-18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:52595 + CHEBI:57876 => CHEBI:60520 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:52595 + CHEBI:57876")
    expect_equal(l$reaction_product_chebi, "CHEBI:60520 + CHEBI:60377 + CHEBI:15378")
    
    reaction <- "RHEA:36189"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(O-18:0/16:0/0:0) + CDP-Ethanolamine <= PE(O-18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(O-18:0/16:0/0:0) + CDP-Ethanolamine")
    expect_equal(l$reaction_product, "PE(O-18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:52595 + CHEBI:57876 <= CHEBI:60520 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:52595 + CHEBI:57876")
    expect_equal(l$reaction_product_chebi, "CHEBI:60520 + CHEBI:60377 + CHEBI:15378")
    
    reaction <- "RHEA:36190"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "DG(O-18:0/16:0/0:0) + CDP-Ethanolamine <=> PE(O-18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "DG(O-18:0/16:0/0:0) + CDP-Ethanolamine")
    expect_equal(l$reaction_product, "PE(O-18:0/16:0) + CMP + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:52595 + CHEBI:57876 <=> CHEBI:60520 + CHEBI:60377 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:52595 + CHEBI:57876")
    expect_equal(l$reaction_product_chebi, "CHEBI:60520 + CHEBI:60377 + CHEBI:15378")
    
    ## pep_to_lpep
    pep <- "PE(P-16:0/14:0)"
    substrates <- list(PEP = pep)
    
    reaction <- "RHEA:36195"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/14:0) + H2O = PE(P-16:0/0:0) + FA(14:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/14:0) + H2O")
    expect_equal(l$reaction_product, "PE(P-16:0/0:0) + FA(14:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77290 + CHEBI:15377 = CHEBI:77288 + CHEBI:29067 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77290 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:77288 + CHEBI:29067 + CHEBI:15378")
    
    reaction <- "RHEA:36196"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/14:0) + H2O => PE(P-16:0/0:0) + FA(14:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/14:0) + H2O")
    expect_equal(l$reaction_product, "PE(P-16:0/0:0) + FA(14:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77290 + CHEBI:15377 => CHEBI:77288 + CHEBI:29067 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77290 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:77288 + CHEBI:29067 + CHEBI:15378")
    
    reaction <- "RHEA:36197"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/14:0) + H2O <= PE(P-16:0/0:0) + FA(14:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/14:0) + H2O")
    expect_equal(l$reaction_product, "PE(P-16:0/0:0) + FA(14:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77290 + CHEBI:15377 <= CHEBI:77288 + CHEBI:29067 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77290 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:77288 + CHEBI:29067 + CHEBI:15378")
    
    reaction <- "RHEA:36198"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/14:0) + H2O <=> PE(P-16:0/0:0) + FA(14:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/14:0) + H2O")
    expect_equal(l$reaction_product, "PE(P-16:0/0:0) + FA(14:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77290 + CHEBI:15377 <=> CHEBI:77288 + CHEBI:29067 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77290 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:77288 + CHEBI:29067 + CHEBI:15378")
    
    ## lpep_to_mgp
    sn1lpep <- "PE(P-16:0/0:0)"
    substrates <- list(sn1LPEP = sn1lpep)
    
    reaction <- "RHEA:36199"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/0:0) + H2O = MG(P-16:0/0:0/0:0) + H+ + Phosphoethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "MG(P-16:0/0:0/0:0) + H+ + Phosphoethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 = CHEBI:77297 + CHEBI:15378 + CHEBI:58190")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77288 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:77297 + CHEBI:15378 + CHEBI:58190")
    
    reaction <- "RHEA:36200"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/0:0) + H2O => MG(P-16:0/0:0/0:0) + H+ + Phosphoethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "MG(P-16:0/0:0/0:0) + H+ + Phosphoethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 => CHEBI:77297 + CHEBI:15378 + CHEBI:58190")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77288 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:77297 + CHEBI:15378 + CHEBI:58190")
    
    reaction <- "RHEA:36201"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/0:0) + H2O <= MG(P-16:0/0:0/0:0) + H+ + Phosphoethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "MG(P-16:0/0:0/0:0) + H+ + Phosphoethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 <= CHEBI:77297 + CHEBI:15378 + CHEBI:58190")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77288 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:77297 + CHEBI:15378 + CHEBI:58190")
    
    reaction <- "RHEA:36202"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/0:0) + H2O <=> MG(P-16:0/0:0/0:0) + H+ + Phosphoethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "MG(P-16:0/0:0/0:0) + H+ + Phosphoethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 <=> CHEBI:77297 + CHEBI:15378 + CHEBI:58190")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77288 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:77297 + CHEBI:15378 + CHEBI:58190")
    
    ## lpep_to_lpap
    sn1lpep <- "PE(P-16:0/0:0)"
    substrates <- list(sn1LPEP = sn1lpep)
    
    reaction <- "RHEA:36203"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/0:0) + H2O = PA(P-16:0/0:0) + Ethanolamine + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "PA(P-16:0/0:0) + Ethanolamine + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 = CHEBI:77283 + CHEBI:57603 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77288 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:77283 + CHEBI:57603 + CHEBI:15378")
    
    reaction <- "RHEA:36204"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/0:0) + H2O => PA(P-16:0/0:0) + Ethanolamine + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "PA(P-16:0/0:0) + Ethanolamine + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 => CHEBI:77283 + CHEBI:57603 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77288 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:77283 + CHEBI:57603 + CHEBI:15378")
    
    reaction <- "RHEA:36205"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/0:0) + H2O <= PA(P-16:0/0:0) + Ethanolamine + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "PA(P-16:0/0:0) + Ethanolamine + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 <= CHEBI:77283 + CHEBI:57603 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77288 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:77283 + CHEBI:57603 + CHEBI:15378")
    
    reaction <- "RHEA:36206"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/0:0) + H2O <=> PA(P-16:0/0:0) + Ethanolamine + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "PA(P-16:0/0:0) + Ethanolamine + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 <=> CHEBI:77283 + CHEBI:57603 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77288 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:77283 + CHEBI:57603 + CHEBI:15378")
    
    ## pco_to_lpco
    pco <- "PC(O-16:0/14:0)"
    substrates <- list(PCO = pco)
    
    reaction <- "RHEA:36231"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(O-16:0/14:0) + H2O = PC(O-16:0/0:0) + FA(14:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(O-16:0/14:0) + H2O")
    expect_equal(l$reaction_product, "PC(O-16:0/0:0) + FA(14:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:36702 + CHEBI:15377 = CHEBI:30909 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:36702 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:30909 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:36232"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(O-16:0/14:0) + H2O => PC(O-16:0/0:0) + FA(14:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(O-16:0/14:0) + H2O")
    expect_equal(l$reaction_product, "PC(O-16:0/0:0) + FA(14:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:36702 + CHEBI:15377 => CHEBI:30909 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:36702 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:30909 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:36233"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(O-16:0/14:0) + H2O <= PC(O-16:0/0:0) + FA(14:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(O-16:0/14:0) + H2O")
    expect_equal(l$reaction_product, "PC(O-16:0/0:0) + FA(14:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:36702 + CHEBI:15377 <= CHEBI:30909 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:36702 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:30909 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:36234"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(O-16:0/14:0) + H2O <=> PC(O-16:0/0:0) + FA(14:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(O-16:0/14:0) + H2O")
    expect_equal(l$reaction_product, "PC(O-16:0/0:0) + FA(14:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:36702 + CHEBI:15377 <=> CHEBI:30909 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:36702 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:30909 + CHEBI:28868 + CHEBI:15378")
    
    ## lpao_to_pao
    lpao <- "PA(O-18:0/0:0)"
    acylcoa <- "CoA(14:0)"
    substrates <- list(sn1LPAO = lpao, AcylCoA = acylcoa)
    
    reaction <- "RHEA:36235"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PA(O-18:0/0:0) + CoA(14:0) = PA(O-18:0/14:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PA(O-18:0/0:0) + CoA(14:0)")
    expect_equal(l$reaction_product, "PA(O-18:0/14:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58014 + CHEBI:58342 = CHEBI:73332 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58014 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:73332 + CHEBI:57287")
    
    reaction <- "RHEA:36236"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PA(O-18:0/0:0) + CoA(14:0) => PA(O-18:0/14:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PA(O-18:0/0:0) + CoA(14:0)")
    expect_equal(l$reaction_product, "PA(O-18:0/14:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58014 + CHEBI:58342 => CHEBI:73332 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58014 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:73332 + CHEBI:57287")
    
    reaction <- "RHEA:36237"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PA(O-18:0/0:0) + CoA(14:0) <= PA(O-18:0/14:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PA(O-18:0/0:0) + CoA(14:0)")
    expect_equal(l$reaction_product, "PA(O-18:0/14:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58014 + CHEBI:58342 <= CHEBI:73332 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58014 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:73332 + CHEBI:57287")
    
    reaction <- "RHEA:36238"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PA(O-18:0/0:0) + CoA(14:0) <=> PA(O-18:0/14:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PA(O-18:0/0:0) + CoA(14:0)")
    expect_equal(l$reaction_product, "PA(O-18:0/14:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58014 + CHEBI:58342 <=> CHEBI:73332 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58014 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:73332 + CHEBI:57287")
    
    ## pao_to_dgo
    pao <- "PA(O-14:0/16:0)"
    substrates = list(PAO = pao)
    
    reaction <- "RHEA:36239"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PA(O-14:0/16:0) + H2O = DG(O-14:0/16:0/0:0) + Pi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PA(O-14:0/16:0) + H2O")
    expect_equal(l$reaction_product, "DG(O-14:0/16:0/0:0) + Pi")
    expect_equal(l$reaction_formula_chebi, "CHEBI:73332 + CHEBI:15377 = CHEBI:52595 + CHEBI:43474")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:73332 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:52595 + CHEBI:43474")
    
    reaction <- "RHEA:36240"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PA(O-14:0/16:0) + H2O => DG(O-14:0/16:0/0:0) + Pi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PA(O-14:0/16:0) + H2O")
    expect_equal(l$reaction_product, "DG(O-14:0/16:0/0:0) + Pi")
    expect_equal(l$reaction_formula_chebi, "CHEBI:73332 + CHEBI:15377 => CHEBI:52595 + CHEBI:43474")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:73332 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:52595 + CHEBI:43474")
    
    reaction <- "RHEA:36241"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PA(O-14:0/16:0) + H2O <= DG(O-14:0/16:0/0:0) + Pi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PA(O-14:0/16:0) + H2O")
    expect_equal(l$reaction_product, "DG(O-14:0/16:0/0:0) + Pi")
    expect_equal(l$reaction_formula_chebi, "CHEBI:73332 + CHEBI:15377 <= CHEBI:52595 + CHEBI:43474")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:73332 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:52595 + CHEBI:43474")
    
    reaction <- "RHEA:36242"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PA(O-14:0/16:0) + H2O <=> DG(O-14:0/16:0/0:0) + Pi")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PA(O-14:0/16:0) + H2O")
    expect_equal(l$reaction_product, "DG(O-14:0/16:0/0:0) + Pi")
    expect_equal(l$reaction_formula_chebi, "CHEBI:73332 + CHEBI:15377 <=> CHEBI:52595 + CHEBI:43474")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:73332 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:52595 + CHEBI:43474")
    
    ## sn1mg_to_dg
    sn1mg <- "MG(14:0/0:0/0:0)"
    acylcoa <- "CoA(16:0)"
    substrates <- list(sn1MG = sn1mg, AcylCoA = acylcoa)
    
    reaction <- "RHEA:38463"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) = DG(14:0/16:0/0:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(14:0/0:0/0:0) + CoA(16:0)")
    expect_equal(l$reaction_product, "DG(14:0/16:0/0:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64683 + CHEBI:58342 = CHEBI:17815 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64683 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:57287")
    
    reaction <- "RHEA:38464"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) => DG(14:0/16:0/0:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(14:0/0:0/0:0) + CoA(16:0)")
    expect_equal(l$reaction_product, "DG(14:0/16:0/0:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64683 + CHEBI:58342 => CHEBI:17815 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64683 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:57287")
    
    reaction <- "RHEA:38465"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) <= DG(14:0/16:0/0:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(14:0/0:0/0:0) + CoA(16:0)")
    expect_equal(l$reaction_product, "DG(14:0/16:0/0:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64683 + CHEBI:58342 <= CHEBI:17815 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64683 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:57287")
    
    reaction <- "RHEA:38466"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) <=> DG(14:0/16:0/0:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(14:0/0:0/0:0) + CoA(16:0)")
    expect_equal(l$reaction_product, "DG(14:0/16:0/0:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64683 + CHEBI:58342 <=> CHEBI:17815 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64683 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:57287")
    
    reaction <- "RHEA:39943"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) = DG(14:0/16:0/0:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(14:0/0:0/0:0) + CoA(16:0)")
    expect_equal(l$reaction_product, "DG(14:0/16:0/0:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:35759 + CHEBI:58342 = CHEBI:49172 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:35759 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:49172 + CHEBI:57287")
    
    reaction <- "RHEA:39944"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) => DG(14:0/16:0/0:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(14:0/0:0/0:0) + CoA(16:0)")
    expect_equal(l$reaction_product, "DG(14:0/16:0/0:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:35759 + CHEBI:58342 => CHEBI:49172 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:35759 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:49172 + CHEBI:57287")
    
    reaction <- "RHEA:39945"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula,
        "MG(14:0/0:0/0:0) + CoA(16:0) <= DG(14:0/16:0/0:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(14:0/0:0/0:0) + CoA(16:0)")
    expect_equal(l$reaction_product, "DG(14:0/16:0/0:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:35759 + CHEBI:58342 <= CHEBI:49172 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:35759 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:49172 + CHEBI:57287")
    
    reaction <- "RHEA:39946"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) <=> DG(14:0/16:0/0:0) + CoA")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(14:0/0:0/0:0) + CoA(16:0)")
    expect_equal(l$reaction_product, "DG(14:0/16:0/0:0) + CoA")
    expect_equal(l$reaction_formula_chebi, "CHEBI:35759 + CHEBI:58342 <=> CHEBI:49172 + CHEBI:57287")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:35759 + CHEBI:58342")
    expect_equal(l$reaction_product_chebi, "CHEBI:49172 + CHEBI:57287")
    
    ## lpco_to_lpao
    sn1lpco <- "PC(O-16:0/0:0)"
    substrates <- list(sn1LPCO = sn1lpco)
    
    reaction <- "RHEA:39927"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(O-16:0/0:0) + H2O = PA(O-16:0/0:0) + Choline + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(O-16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "PA(O-16:0/0:0) + Choline + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 = CHEBI:58014 + CHEBI:15354 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:30909 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:58014 + CHEBI:15354 + CHEBI:15378")
    
    reaction <- "RHEA:39928"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(O-16:0/0:0) + H2O => PA(O-16:0/0:0) + Choline + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(O-16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "PA(O-16:0/0:0) + Choline + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 => CHEBI:58014 + CHEBI:15354 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:30909 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:58014 + CHEBI:15354 + CHEBI:15378")
    
    reaction <- "RHEA:39929"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(O-16:0/0:0) + H2O <= PA(O-16:0/0:0) + Choline + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(O-16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "PA(O-16:0/0:0) + Choline + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 <= CHEBI:58014 + CHEBI:15354 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:30909 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:58014 + CHEBI:15354 + CHEBI:15378")
    
    reaction <- "RHEA:39930"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(O-16:0/0:0) + H2O <=> PA(O-16:0/0:0) + Choline + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(O-16:0/0:0) + H2O")
    expect_equal(l$reaction_product, "PA(O-16:0/0:0) + Choline + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 <=> CHEBI:58014 + CHEBI:15354 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:30909 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:58014 + CHEBI:15354 + CHEBI:15378")
    
    ## pi_to_dg
    pi <- "PI(16:0/18:1(9Z))"
    substrates <- list(PI = pi)
    
    reaction <- "RHEA:43484"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PI(16:0/18:1(9Z)) + H2O = myo-Inositol-1-P + DG(16:0/18:1(9Z)) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PI(16:0/18:1(9Z)) + H2O")
    expect_equal(l$reaction_product, "myo-Inositol-1-P + DG(16:0/18:1(9Z)) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 = CHEBI:58433 + CHEBI:17815 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57880 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:58433 + CHEBI:17815 + CHEBI:15378")
    
    reaction <- "RHEA:43485"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PI(16:0/18:1(9Z)) + H2O => myo-Inositol-1-P + DG(16:0/18:1(9Z)) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PI(16:0/18:1(9Z)) + H2O")
    expect_equal(l$reaction_product, "myo-Inositol-1-P + DG(16:0/18:1(9Z)) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 => CHEBI:58433 + CHEBI:17815 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57880 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:58433 + CHEBI:17815 + CHEBI:15378")
    
    reaction <- "RHEA:43486"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PI(16:0/18:1(9Z)) + H2O <= myo-Inositol-1-P + DG(16:0/18:1(9Z)) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PI(16:0/18:1(9Z)) + H2O")
    expect_equal(l$reaction_product, "myo-Inositol-1-P + DG(16:0/18:1(9Z)) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 <= CHEBI:58433 + CHEBI:17815 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57880 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:58433 + CHEBI:17815 + CHEBI:15378")
    
    reaction <- "RHEA:43487"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PI(16:0/18:1(9Z)) + H2O <=> myo-Inositol-1-P + DG(16:0/18:1(9Z)) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PI(16:0/18:1(9Z)) + H2O")
    expect_equal(l$reaction_product, "myo-Inositol-1-P + DG(16:0/18:1(9Z)) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 <=> CHEBI:58433 + CHEBI:17815 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57880 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:58433 + CHEBI:17815 + CHEBI:15378")
    
    ## pe_to_sn2lpe
    pe <- "PE(14:0/16:0)"
    substrates = list(PE = pe)
    
    reaction <- "RHEA:44408"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(14:0/16:0) + H2O = PE(0:0/14:0) + FA(16:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(14:0/16:0) + H2O")
    expect_equal(l$reaction_product, "PE(0:0/14:0) + FA(16:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 = CHEBI:65213 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64612 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:65213 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:44409"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(14:0/16:0) + H2O => PE(0:0/14:0) + FA(16:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(14:0/16:0) + H2O")
    expect_equal(l$reaction_product, "PE(0:0/14:0) + FA(16:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 => CHEBI:65213 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64612 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:65213 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:44410"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(14:0/16:0) + H2O <= PE(0:0/14:0) + FA(16:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(14:0/16:0) + H2O")
    expect_equal(l$reaction_product, "PE(0:0/14:0) + FA(16:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 <= CHEBI:65213 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64612 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:65213 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:44411"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(14:0/16:0) + H2O <=> PE(0:0/14:0) + FA(16:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(14:0/16:0) + H2O")
    expect_equal(l$reaction_product, "PE(0:0/14:0) + FA(16:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 <=> CHEBI:65213 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64612 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:65213 + CHEBI:28868 + CHEBI:15378")
    
    ## pg_to_sn1lpg
    pg <- "PG(14:0/16:0)"
    substrates = list(PG = pg)
    
    reaction <- "RHEA:44416"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PG(14:0/16:0) + H2O = PG(14:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PG(14:0/16:0) + H2O")
    expect_equal(l$reaction_product, "PG(14:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64716 + CHEBI:15377 = CHEBI:64840 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64716 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64840 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:44417"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PG(14:0/16:0) + H2O => PG(14:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PG(14:0/16:0) + H2O")
    expect_equal(l$reaction_product, "PG(14:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64716 + CHEBI:15377 => CHEBI:64840 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64716 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64840 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:44418"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PG(14:0/16:0) + H2O <= PG(14:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PG(14:0/16:0) + H2O")
    expect_equal(l$reaction_product, "PG(14:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64716 + CHEBI:15377 <= CHEBI:64840 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64716 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64840 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:44419"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PG(14:0/16:0) + H2O <=> PG(14:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PG(14:0/16:0) + H2O")
    expect_equal(l$reaction_product, "PG(14:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64716 + CHEBI:15377 <=> CHEBI:64840 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64716 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64840 + CHEBI:28868 + CHEBI:15378")
    
    ## pe_to_sn1lpe
    pe <- "PE(14:0/16:0)"
    substrates = list(PE = pe)
    
    reaction <- "RHEA:44604"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(14:0/16:0) + H2O = PE(14:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(14:0/16:0) + H2O")
    expect_equal(l$reaction_product, "PE(14:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 = CHEBI:64381 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64612 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64381 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:44605"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(14:0/16:0) + H2O => PE(14:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(14:0/16:0) + H2O")
    expect_equal(l$reaction_product, "PE(14:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 => CHEBI:64381 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64612 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64381 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:44606"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(14:0/16:0) + H2O <= PE(14:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(14:0/16:0) + H2O")
    expect_equal(l$reaction_product, "PE(14:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 <= CHEBI:64381 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64612 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64381 + CHEBI:28868 + CHEBI:15378")
    
    reaction <- "RHEA:44607"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(14:0/16:0) + H2O <=> PE(14:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(14:0/16:0) + H2O")
    expect_equal(l$reaction_product, "PE(14:0/0:0) + FA(16:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 <=> CHEBI:64381 + CHEBI:28868 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64612 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:64381 + CHEBI:28868 + CHEBI:15378")
    
    ## dhcer_to_dhsm
    dhcer <- "Cer(d16:0(3OH,4OH)(15Me)/12:0)"
    substrates <- list(DhCer = dhcer)
    
    reaction <- "RHEA:44620"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC + Cer(d16:0(3OH,4OH)(15Me)/12:0) = DG + SM(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC + Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(l$reaction_product, "DG + SM(d16:0(3OH,4OH)(15Me)/12:0)")  
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:31488 = CHEBI:17815 + CHEBI:67090")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:31488")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:67090")
    
    reaction <- "RHEA:44621"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC + Cer(d16:0(3OH,4OH)(15Me)/12:0) => DG + SM(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC + Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(l$reaction_product, "DG + SM(d16:0(3OH,4OH)(15Me)/12:0)")  
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:31488 => CHEBI:17815 + CHEBI:67090")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:31488")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:67090")
    
    reaction <- "RHEA:44622"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC + Cer(d16:0(3OH,4OH)(15Me)/12:0) <= DG + SM(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC + Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(l$reaction_product, "DG + SM(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:31488 <= CHEBI:17815 + CHEBI:67090")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:31488")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:67090")
    
    reaction <- "RHEA:44623"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC + Cer(d16:0(3OH,4OH)(15Me)/12:0) <=> DG + SM(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC + Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(l$reaction_product, "DG + SM(d16:0(3OH,4OH)(15Me)/12:0)") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:31488 <=> CHEBI:17815 + CHEBI:67090")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:31488")
    expect_equal(l$reaction_product_chebi, "CHEBI:17815 + CHEBI:67090")
    
    ## sn2lpc_to_fa
    sn2lpc <- "PC(0:0/14:0)"
    substrates <- list(sn2LPC = sn2lpc)
    
    reaction <- "RHEA:44696"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(0:0/14:0) + H2O = FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(0:0/14:0) + H2O")
    expect_equal(l$reaction_product, "FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57875 + CHEBI:15377 = CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57875 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    
    reaction <- "RHEA:44697"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(0:0/14:0) + H2O => FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(0:0/14:0) + H2O")
    expect_equal(l$reaction_product, "FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57875 + CHEBI:15377 => CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57875 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    
    reaction <- "RHEA:44698"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(0:0/14:0) + H2O <= FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(0:0/14:0) + H2O")
    expect_equal(l$reaction_product, "FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57875 + CHEBI:15377 <= CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57875 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    
    reaction <- "RHEA:44699"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(0:0/14:0) + H2O <=> FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(0:0/14:0) + H2O")
    expect_equal(l$reaction_product, "FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57875 + CHEBI:15377 <=> CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57875 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    
    ## pc_to_ps
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    
    reaction <- "RHEA:45088"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(20:0/18:0) + L-Serine = PS(20:0/18:0) + Choline")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(20:0/18:0) + L-Serine")
    expect_equal(l$reaction_product, "PS(20:0/18:0) + Choline")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:33384 = CHEBI:57262 + CHEBI:15354")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:33384")
    expect_equal(l$reaction_product_chebi, "CHEBI:57262 + CHEBI:15354")
    
    reaction <- "RHEA:45089"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(20:0/18:0) + L-Serine => PS(20:0/18:0) + Choline")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(20:0/18:0) + L-Serine")
    expect_equal(l$reaction_product, "PS(20:0/18:0) + Choline")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:33384 => CHEBI:57262 + CHEBI:15354")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:33384")
    expect_equal(l$reaction_product_chebi, "CHEBI:57262 + CHEBI:15354")
    
    reaction <- "RHEA:45090"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(20:0/18:0) + L-Serine <= PS(20:0/18:0) + Choline")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(20:0/18:0) + L-Serine")
    expect_equal(l$reaction_product, "PS(20:0/18:0) + Choline")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:33384 <= CHEBI:57262 + CHEBI:15354")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:33384")
    expect_equal(l$reaction_product_chebi, "CHEBI:57262 + CHEBI:15354")
    
    reaction <- "RHEA:45091"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(20:0/18:0) + L-Serine <=> PS(20:0/18:0) + Choline")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(20:0/18:0) + L-Serine")
    expect_equal(l$reaction_product, "PS(20:0/18:0) + Choline")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:33384 <=> CHEBI:57262 + CHEBI:15354")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:33384")
    expect_equal(l$reaction_product_chebi, "CHEBI:57262 + CHEBI:15354")
    
    ## pe_to_nape_sn1
    pe <- "PE(14:0/16:0)"
    pc <- "PC(18:0/20:0)"
    substrates <- list(PE = pe, PC = pc)
    
    reaction <- "RHEA:45188"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(18:0/20:0) + PE(14:0/16:0) = PC(0:0/20:0) + H+ + NAPE(14:0/16:0/18:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(18:0/20:0) + PE(14:0/16:0)")
    expect_equal(l$reaction_product, "PC(0:0/20:0) + H+ + NAPE(14:0/16:0/18:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 = CHEBI:57875 + CHEBI:15378 + CHEBI:62537")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:64612")
    expect_equal(l$reaction_product_chebi, "CHEBI:57875 + CHEBI:15378 + CHEBI:62537")
    
    reaction <- "RHEA:45189"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(18:0/20:0) + PE(14:0/16:0) => PC(0:0/20:0) + H+ + NAPE(14:0/16:0/18:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(18:0/20:0) + PE(14:0/16:0)")
    expect_equal(l$reaction_product, "PC(0:0/20:0) + H+ + NAPE(14:0/16:0/18:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 => CHEBI:57875 + CHEBI:15378 + CHEBI:62537")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:64612")
    expect_equal(l$reaction_product_chebi, "CHEBI:57875 + CHEBI:15378 + CHEBI:62537")
    
    reaction <- "RHEA:45190"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(18:0/20:0) + PE(14:0/16:0) <= PC(0:0/20:0) + H+ + NAPE(14:0/16:0/18:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(18:0/20:0) + PE(14:0/16:0)")
    expect_equal(l$reaction_product, "PC(0:0/20:0) + H+ + NAPE(14:0/16:0/18:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 <= CHEBI:57875 + CHEBI:15378 + CHEBI:62537")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:64612")
    expect_equal(l$reaction_product_chebi, "CHEBI:57875 + CHEBI:15378 + CHEBI:62537")
    
    reaction <- "RHEA:45191"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(18:0/20:0) + PE(14:0/16:0) <=> PC(0:0/20:0) + H+ + NAPE(14:0/16:0/18:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(18:0/20:0) + PE(14:0/16:0)")
    expect_equal(l$reaction_product, "PC(0:0/20:0) + H+ + NAPE(14:0/16:0/18:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 <=> CHEBI:57875 + CHEBI:15378 + CHEBI:62537")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:64612")
    expect_equal(l$reaction_product_chebi, "CHEBI:57875 + CHEBI:15378 + CHEBI:62537")
    
    ## pe_to_nape_sn2
    pe <- "PE(14:0/16:0)"
    pc <- "PC(18:0/20:0)"
    substrates <- list(PE = pe, PC = pc)
    
    reaction <- "RHEA:45192"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(18:0/20:0) + PE(14:0/16:0) = PC(18:0/0:0) + H+ + NAPE(14:0/16:0/20:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(18:0/20:0) + PE(14:0/16:0)")
    expect_equal(l$reaction_product, "PC(18:0/0:0) + H+ + NAPE(14:0/16:0/20:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 = CHEBI:58168 + CHEBI:15378 + CHEBI:62537")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:64612")
    expect_equal(l$reaction_product_chebi, "CHEBI:58168 + CHEBI:15378 + CHEBI:62537")
    
    reaction <- "RHEA:45193"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(18:0/20:0) + PE(14:0/16:0) => PC(18:0/0:0) + H+ + NAPE(14:0/16:0/20:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(18:0/20:0) + PE(14:0/16:0)")
    expect_equal(l$reaction_product, "PC(18:0/0:0) + H+ + NAPE(14:0/16:0/20:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 => CHEBI:58168 + CHEBI:15378 + CHEBI:62537")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:64612")
    expect_equal(l$reaction_product_chebi, "CHEBI:58168 + CHEBI:15378 + CHEBI:62537")
    
    reaction <- "RHEA:45194"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(18:0/20:0) + PE(14:0/16:0) <= PC(18:0/0:0) + H+ + NAPE(14:0/16:0/20:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(18:0/20:0) + PE(14:0/16:0)")
    expect_equal(l$reaction_product, "PC(18:0/0:0) + H+ + NAPE(14:0/16:0/20:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 <= CHEBI:58168 + CHEBI:15378 + CHEBI:62537")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:64612")
    expect_equal(l$reaction_product_chebi, "CHEBI:58168 + CHEBI:15378 + CHEBI:62537")
    
    reaction <- "RHEA:45195"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PC(18:0/20:0) + PE(14:0/16:0) <=> PC(18:0/0:0) + H+ + NAPE(14:0/16:0/20:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PC(18:0/20:0) + PE(14:0/16:0)")
    expect_equal(l$reaction_product, "PC(18:0/0:0) + H+ + NAPE(14:0/16:0/20:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 <=> CHEBI:58168 + CHEBI:15378 + CHEBI:62537")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:57643 + CHEBI:64612")
    expect_equal(l$reaction_product_chebi, "CHEBI:58168 + CHEBI:15378 + CHEBI:62537")
    
    ## dhsm_to_dhcer
    dhsm <- "SM(d16:0(3OH,4OH)(15Me)/12:0)"
    substrates <- list(DhSM = dhsm)
    
    reaction <- "RHEA:45300"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "SM(d16:0(3OH,4OH)(15Me)/12:0) + H2O = Cer(d16:0(3OH,4OH)(15Me)/12:0) + H+ + Phosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "SM(d16:0(3OH,4OH)(15Me)/12:0) + H2O")
    expect_equal(l$reaction_product, "Cer(d16:0(3OH,4OH)(15Me)/12:0) + H+ + Phosphocholine")        
    expect_equal(l$reaction_formula_chebi, "CHEBI:64583 + CHEBI:15377 = CHEBI:83273 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64583 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:83273 + CHEBI:15378 + CHEBI:295975")
    
    reaction <- "RHEA:45301"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "SM(d16:0(3OH,4OH)(15Me)/12:0) + H2O => Cer(d16:0(3OH,4OH)(15Me)/12:0) + H+ + Phosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "SM(d16:0(3OH,4OH)(15Me)/12:0) + H2O")
    expect_equal(l$reaction_product, "Cer(d16:0(3OH,4OH)(15Me)/12:0) + H+ + Phosphocholine")    
    expect_equal(l$reaction_formula_chebi, "CHEBI:64583 + CHEBI:15377 => CHEBI:83273 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64583 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:83273 + CHEBI:15378 + CHEBI:295975")
    
    reaction <- "RHEA:45302"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "SM(d16:0(3OH,4OH)(15Me)/12:0) + H2O <= Cer(d16:0(3OH,4OH)(15Me)/12:0) + H+ + Phosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "SM(d16:0(3OH,4OH)(15Me)/12:0) + H2O")
    expect_equal(l$reaction_product, "Cer(d16:0(3OH,4OH)(15Me)/12:0) + H+ + Phosphocholine")      
    expect_equal(l$reaction_formula_chebi, "CHEBI:64583 + CHEBI:15377 <= CHEBI:83273 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64583 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:83273 + CHEBI:15378 + CHEBI:295975")
    
    reaction <- "RHEA:45303"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "SM(d16:0(3OH,4OH)(15Me)/12:0) + H2O <=> Cer(d16:0(3OH,4OH)(15Me)/12:0) + H+ + Phosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "SM(d16:0(3OH,4OH)(15Me)/12:0) + H2O")
    expect_equal(l$reaction_product, "Cer(d16:0(3OH,4OH)(15Me)/12:0) + H+ + Phosphocholine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64583 + CHEBI:15377 <=> CHEBI:83273 + CHEBI:15378 + CHEBI:295975")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64583 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:83273 + CHEBI:15378 + CHEBI:295975")
    
    ## lnape_to_gpnae
    lnape <- "NAPE(14:0/0:0/0:0)"
    substrates <- list(LNAPE = lnape)
    
    reaction <- "RHEA:45420"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + NAPE(14:0/0:0/0:0) = FA(14:0) + H+ + GPNAE(0:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + NAPE(14:0/0:0/0:0)")
    expect_equal(l$reaction_product, "FA(14:0) + H+ + GPNAE(0:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:85216 = CHEBI:28868 + CHEBI:15378 + CHEBI:85225")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:85216")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:15378 + CHEBI:85225")
    
    reaction <- "RHEA:45421"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + NAPE(14:0/0:0/0:0) => FA(14:0) + H+ + GPNAE(0:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + NAPE(14:0/0:0/0:0)")
    expect_equal(l$reaction_product, "FA(14:0) + H+ + GPNAE(0:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:85216 => CHEBI:28868 + CHEBI:15378 + CHEBI:85225")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:85216")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:15378 + CHEBI:85225")
    
    reaction <- "RHEA:45422"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + NAPE(14:0/0:0/0:0) <= FA(14:0) + H+ + GPNAE(0:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + NAPE(14:0/0:0/0:0)")
    expect_equal(l$reaction_product, "FA(14:0) + H+ + GPNAE(0:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:85216 <= CHEBI:28868 + CHEBI:15378 + CHEBI:85225")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:85216")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:15378 + CHEBI:85225")
    
    reaction <- "RHEA:45423"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + NAPE(14:0/0:0/0:0) <=> FA(14:0) + H+ + GPNAE(0:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + NAPE(14:0/0:0/0:0)")
    expect_equal(l$reaction_product, "FA(14:0) + H+ + GPNAE(0:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:85216 <=> CHEBI:28868 + CHEBI:15378 + CHEBI:85225")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:85216")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:15378 + CHEBI:85225")
    
    ## nape_to_lnape
    nape <- "NAPE(14:0/16:0/18:0)"
    substrates = list(NAPE = nape)
    
    reaction <- "RHEA:45460"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + NAPE(14:0/16:0/18:0) = FA(16:0) + H+ + NAPE(14:0/0:0/18:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + NAPE(14:0/16:0/18:0)")
    expect_equal(l$reaction_product, "FA(16:0) + H+ + NAPE(14:0/0:0/18:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 = CHEBI:28868 + CHEBI:15378 + CHEBI:85216")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:62537")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:15378 + CHEBI:85216")
    
    reaction <- "RHEA:45461"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + NAPE(14:0/16:0/18:0) => FA(16:0) + H+ + NAPE(14:0/0:0/18:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + NAPE(14:0/16:0/18:0)")
    expect_equal(l$reaction_product, "FA(16:0) + H+ + NAPE(14:0/0:0/18:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 => CHEBI:28868 + CHEBI:15378 + CHEBI:85216")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:62537")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:15378 + CHEBI:85216")
    
    reaction <- "RHEA:45462"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + NAPE(14:0/16:0/18:0) <= FA(16:0) + H+ + NAPE(14:0/0:0/18:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + NAPE(14:0/16:0/18:0)")
    expect_equal(l$reaction_product, "FA(16:0) + H+ + NAPE(14:0/0:0/18:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 <= CHEBI:28868 + CHEBI:15378 + CHEBI:85216")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:62537")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:15378 + CHEBI:85216")
    
    reaction <- "RHEA:45463"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + NAPE(14:0/16:0/18:0) <=> FA(16:0) + H+ + NAPE(14:0/0:0/18:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + NAPE(14:0/16:0/18:0)")
    expect_equal(l$reaction_product, "FA(16:0) + H+ + NAPE(14:0/0:0/18:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 <=> CHEBI:28868 + CHEBI:15378 + CHEBI:85216")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:62537")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:15378 + CHEBI:85216")
    
    ## sm_to_cer
    sm <- "SM(16:0(3OH,4OH,15Me)/12:0)"
    substrates <- list(SM = sm)
    
    reaction <- "RHEA:45644"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + SM(16:0(3OH,4OH,15Me)/12:0) = H+ + Cer(16:0(3OH,4OH,15Me)/12:0) + Phosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + SM(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l$reaction_product, "H+ + Cer(16:0(3OH,4OH,15Me)/12:0) + Phosphocholine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:78646 = CHEBI:15378 + CHEBI:72959 + CHEBI:295975")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:78646")
    expect_equal(l$reaction_product_chebi, "CHEBI:15378 + CHEBI:72959 + CHEBI:295975")
    
    reaction <- "RHEA:45645"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + SM(16:0(3OH,4OH,15Me)/12:0) => H+ + Cer(16:0(3OH,4OH,15Me)/12:0) + Phosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + SM(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l$reaction_product, "H+ + Cer(16:0(3OH,4OH,15Me)/12:0) + Phosphocholine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:78646 => CHEBI:15378 + CHEBI:72959 + CHEBI:295975")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:78646")
    expect_equal(l$reaction_product_chebi, "CHEBI:15378 + CHEBI:72959 + CHEBI:295975")
    
    reaction <- "RHEA:45646"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + SM(16:0(3OH,4OH,15Me)/12:0) <= H+ + Cer(16:0(3OH,4OH,15Me)/12:0) + Phosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + SM(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l$reaction_product, "H+ + Cer(16:0(3OH,4OH,15Me)/12:0) + Phosphocholine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:78646 <= CHEBI:15378 + CHEBI:72959 + CHEBI:295975")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:78646")
    expect_equal(l$reaction_product_chebi, "CHEBI:15378 + CHEBI:72959 + CHEBI:295975")
    
    reaction <- "RHEA:45647"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + SM(16:0(3OH,4OH,15Me)/12:0) <=> H+ + Cer(16:0(3OH,4OH,15Me)/12:0) + Phosphocholine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + SM(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(l$reaction_product, "H+ + Cer(16:0(3OH,4OH,15Me)/12:0) + Phosphocholine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:78646 <=> CHEBI:15378 + CHEBI:72959 + CHEBI:295975")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:78646")
    expect_equal(l$reaction_product_chebi, "CHEBI:15378 + CHEBI:72959 + CHEBI:295975")
    
    ## dhcer_to_cer
    dhcer <- "Cer(d16:0(3OH,4OH)(15Me)/12:0)"
    substrates <- list(DhCer = dhcer)
    
    reaction <- "RHEA:46544"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "Cer(d16:0(3OH,4OH)(15Me)/12:0) + 2 Fe2+-cytochrome b5 + 2 H+ + O2 = 2 Fe3+-cytochrome b5 + Cer(d16:1(3OH,4OH)(15Me)/12:0) + 2 H2O") ## was: Cer(d16:1(4E)(3OH,4OH)(15Me)/12:0)
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "Cer(d16:0(3OH,4OH)(15Me)/12:0) + 2 Fe2+-cytochrome b5 + 2 H+ + O2")
    expect_equal(l$reaction_product, "2 Fe3+-cytochrome b5 + Cer(d16:1(3OH,4OH)(15Me)/12:0) + 2 H2O") ## was: Cer(d16:1(4E)(3OH,4OH)(15Me)/12:0)
    expect_equal(l$reaction_formula_chebi, "CHEBI:31488 + 2 CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 = 2 CHEBI:57540 + CHEBI:52639 + 2 CHEBI:15377")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:31488 + 2 CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379")
    expect_equal(l$reaction_product_chebi, "2 CHEBI:57540 + CHEBI:52639 + 2 CHEBI:15377")
    
    reaction <- "RHEA:46545"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "Cer(d16:0(3OH,4OH)(15Me)/12:0) + 2 Fe2+-cytochrome b5 + 2 H+ + O2 => 2 Fe3+-cytochrome b5 + Cer(d16:1(3OH,4OH)(15Me)/12:0) + 2 H2O") ## was: Cer(d16:1(4E)(3OH,4OH)(15Me)/12:0)
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "Cer(d16:0(3OH,4OH)(15Me)/12:0) + 2 Fe2+-cytochrome b5 + 2 H+ + O2")
    expect_equal(l$reaction_product, "2 Fe3+-cytochrome b5 + Cer(d16:1(3OH,4OH)(15Me)/12:0) + 2 H2O") ## was: Cer(d16:1(4E)(3OH,4OH)(15Me)/12:0)
    expect_equal(l$reaction_formula_chebi, "CHEBI:31488 + 2 CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 => 2 CHEBI:57540 + CHEBI:52639 + 2 CHEBI:15377")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:31488 + 2 CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379")
    expect_equal(l$reaction_product_chebi, "2 CHEBI:57540 + CHEBI:52639 + 2 CHEBI:15377")
    
    reaction <- "RHEA:46546"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "Cer(d16:0(3OH,4OH)(15Me)/12:0) + 2 Fe2+-cytochrome b5 + 2 H+ + O2 <= 2 Fe3+-cytochrome b5 + Cer(d16:1(3OH,4OH)(15Me)/12:0) + 2 H2O") ## was: Cer(d16:1(4E)(3OH,4OH)(15Me)/12:0)
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "Cer(d16:0(3OH,4OH)(15Me)/12:0) + 2 Fe2+-cytochrome b5 + 2 H+ + O2")
    expect_equal(l$reaction_product, "2 Fe3+-cytochrome b5 + Cer(d16:1(3OH,4OH)(15Me)/12:0) + 2 H2O") ## was: Cer(d16:1(4E)(3OH,4OH)(15Me)/12:0)
    expect_equal(l$reaction_formula_chebi, "CHEBI:31488 + 2 CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 <= 2 CHEBI:57540 + CHEBI:52639 + 2 CHEBI:15377")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:31488 + 2 CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379")
    expect_equal(l$reaction_product_chebi, "2 CHEBI:57540 + CHEBI:52639 + 2 CHEBI:15377")
    
    reaction <- "RHEA:46547"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "Cer(d16:0(3OH,4OH)(15Me)/12:0) + 2 Fe2+-cytochrome b5 + 2 H+ + O2 <=> 2 Fe3+-cytochrome b5 + Cer(d16:1(3OH,4OH)(15Me)/12:0) + 2 H2O") ## was: Cer(d16:1(4E)(3OH,4OH)(15Me)/12:0)
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "Cer(d16:0(3OH,4OH)(15Me)/12:0) + 2 Fe2+-cytochrome b5 + 2 H+ + O2")
    expect_equal(l$reaction_product, "2 Fe3+-cytochrome b5 + Cer(d16:1(3OH,4OH)(15Me)/12:0) + 2 H2O") ## was: Cer(d16:1(4E)(3OH,4OH)(15Me)/12:0)
    expect_equal(l$reaction_formula_chebi, "CHEBI:31488 + 2 CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 <=> 2 CHEBI:57540 + CHEBI:52639 + 2 CHEBI:15377")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:31488 + 2 CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379")
    expect_equal(l$reaction_product_chebi, "2 CHEBI:57540 + CHEBI:52639 + 2 CHEBI:15377")
    
    ## coa_to_FAO
    acylcoa <- "CoA(18:0)"
    substrates <- list(AcylCoA = acylcoa)
    
    reaction <- "RHEA:52716"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CoA(18:0) + 2 H+ + 2 NADPH = FAO(18:0) + CoA + 2 NADP+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CoA(18:0) + 2 H+ + 2 NADPH")
    expect_equal(l$reaction_product, "FAO(18:0) + CoA + 2 NADP+") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:83139 + 2 CHEBI:15378 + 2 CHEBI:57783 = CHEBI:77396 + CHEBI:57287 + 2 CHEBI:58349")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:83139 + 2 CHEBI:15378 + 2 CHEBI:57783")
    expect_equal(l$reaction_product_chebi, "CHEBI:77396 + CHEBI:57287 + 2 CHEBI:58349")
    
    reaction <- "RHEA:52717"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CoA(18:0) + 2 H+ + 2 NADPH => FAO(18:0) + CoA + 2 NADP+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CoA(18:0) + 2 H+ + 2 NADPH")
    expect_equal(l$reaction_product, "FAO(18:0) + CoA + 2 NADP+") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:83139 + 2 CHEBI:15378 + 2 CHEBI:57783 => CHEBI:77396 + CHEBI:57287 + 2 CHEBI:58349")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:83139 + 2 CHEBI:15378 + 2 CHEBI:57783")
    expect_equal(l$reaction_product_chebi, "CHEBI:77396 + CHEBI:57287 + 2 CHEBI:58349")
    
    reaction <- "RHEA:52718"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CoA(18:0) + 2 H+ + 2 NADPH <= FAO(18:0) + CoA + 2 NADP+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CoA(18:0) + 2 H+ + 2 NADPH")
    expect_equal(l$reaction_product, "FAO(18:0) + CoA + 2 NADP+") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:83139 + 2 CHEBI:15378 + 2 CHEBI:57783 <= CHEBI:77396 + CHEBI:57287 + 2 CHEBI:58349")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:83139 + 2 CHEBI:15378 + 2 CHEBI:57783")
    expect_equal(l$reaction_product_chebi, "CHEBI:77396 + CHEBI:57287 + 2 CHEBI:58349")
    
    reaction <- "RHEA:52719"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CoA(18:0) + 2 H+ + 2 NADPH <=> FAO(18:0) + CoA + 2 NADP+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CoA(18:0) + 2 H+ + 2 NADPH")
    expect_equal(l$reaction_product, "FAO(18:0) + CoA + 2 NADP+") 
    expect_equal(l$reaction_formula_chebi, "CHEBI:83139 + 2 CHEBI:15378 + 2 CHEBI:57783 <=> CHEBI:77396 + CHEBI:57287 + 2 CHEBI:58349")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:83139 + 2 CHEBI:15378 + 2 CHEBI:57783")
    expect_equal(l$reaction_product_chebi, "CHEBI:77396 + CHEBI:57287 + 2 CHEBI:58349")
    
    ## sphinga_to_dhcer
    acylcoa <- "CoA(12:0)"
    sph <- "SPH(d16:0(1OH,3OH)(15Me))"
    substrates <- list(AcylCoA = acylcoa, SPH = sph)
    
    reaction <- "RHEA:53424"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me)) = Cer(d16:0(1OH,3OH)(15Me)/12:0) + CoA + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me))")
    expect_equal(l$reaction_product, "Cer(d16:0(1OH,3OH)(15Me)/12:0) + CoA + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77636 + CHEBI:84410 = CHEBI:83273 + CHEBI:57287 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77636 + CHEBI:84410")
    expect_equal(l$reaction_product_chebi, "CHEBI:83273 + CHEBI:57287 + CHEBI:15378")
    
    reaction <- "RHEA:53425"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me)) => Cer(d16:0(1OH,3OH)(15Me)/12:0) + CoA + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me))")
    expect_equal(l$reaction_product, "Cer(d16:0(1OH,3OH)(15Me)/12:0) + CoA + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77636 + CHEBI:84410 => CHEBI:83273 + CHEBI:57287 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77636 + CHEBI:84410")
    expect_equal(l$reaction_product_chebi, "CHEBI:83273 + CHEBI:57287 + CHEBI:15378")
    
    reaction <- "RHEA:53426"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me)) <= Cer(d16:0(1OH,3OH)(15Me)/12:0) + CoA + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me))")
    expect_equal(l$reaction_product, "Cer(d16:0(1OH,3OH)(15Me)/12:0) + CoA + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77636 + CHEBI:84410 <= CHEBI:83273 + CHEBI:57287 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77636 + CHEBI:84410")
    expect_equal(l$reaction_product_chebi, "CHEBI:83273 + CHEBI:57287 + CHEBI:15378")
    
    reaction <- "RHEA:53427"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me)) <=> Cer(d16:0(1OH,3OH)(15Me)/12:0) + CoA + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me))")
    expect_equal(l$reaction_product, "Cer(d16:0(1OH,3OH)(15Me)/12:0) + CoA + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77636 + CHEBI:84410 <=> CHEBI:83273 + CHEBI:57287 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77636 + CHEBI:84410")
    expect_equal(l$reaction_product_chebi, "CHEBI:83273 + CHEBI:57287 + CHEBI:15378")
    
    ## pep_to_napep_sn1
    pep <- "PE(P-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    substrates <- list(PEP = pep, PC = pc)
    
    reaction <- "RHEA:63596"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/14:0) + PC(20:0/18:0) = PC(0:0/18:0) + H+ + NAPE(P-16:0/14:0/20:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/14:0) + PC(20:0/18:0)")
    expect_equal(l$reaction_product, "PC(0:0/18:0) + H+ + NAPE(P-16:0/14:0/20:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77290 + CHEBI:57643 = CHEBI:57875 + CHEBI:15378 + CHEBI:140451")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77290 + CHEBI:57643")
    expect_equal(l$reaction_product_chebi, "CHEBI:57875 + CHEBI:15378 + CHEBI:140451")
    
    reaction <- "RHEA:63597"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/14:0) + PC(20:0/18:0) => PC(0:0/18:0) + H+ + NAPE(P-16:0/14:0/20:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/14:0) + PC(20:0/18:0)")
    expect_equal(l$reaction_product, "PC(0:0/18:0) + H+ + NAPE(P-16:0/14:0/20:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77290 + CHEBI:57643 => CHEBI:57875 + CHEBI:15378 + CHEBI:140451")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77290 + CHEBI:57643")
    expect_equal(l$reaction_product_chebi, "CHEBI:57875 + CHEBI:15378 + CHEBI:140451")
    
    reaction <- "RHEA:63598"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/14:0) + PC(20:0/18:0) <= PC(0:0/18:0) + H+ + NAPE(P-16:0/14:0/20:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/14:0) + PC(20:0/18:0)")
    expect_equal(l$reaction_product, "PC(0:0/18:0) + H+ + NAPE(P-16:0/14:0/20:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77290 + CHEBI:57643 <= CHEBI:57875 + CHEBI:15378 + CHEBI:140451")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77290 + CHEBI:57643")
    expect_equal(l$reaction_product_chebi, "CHEBI:57875 + CHEBI:15378 + CHEBI:140451")
    
    reaction <- "RHEA:63599"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/14:0) + PC(20:0/18:0) <=> PC(0:0/18:0) + H+ + NAPE(P-16:0/14:0/20:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/14:0) + PC(20:0/18:0)")
    expect_equal(l$reaction_product, "PC(0:0/18:0) + H+ + NAPE(P-16:0/14:0/20:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77290 + CHEBI:57643 <=> CHEBI:57875 + CHEBI:15378 + CHEBI:140451")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77290 + CHEBI:57643")
    expect_equal(l$reaction_product_chebi, "CHEBI:57875 + CHEBI:15378 + CHEBI:140451")
    
    ## pe_to_dg
    pe <- "PE(14:0/16:0)"
    substrates <- list(PE = pe)
    
    reaction <- "RHEA:78951"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + PE(14:0/16:0) = P-Ethanolamine + DG(14:0/16:0/0:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + PE(14:0/16:0)")
    expect_equal(l$reaction_product, "P-Ethanolamine + DG(14:0/16:0/0:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:64612 = CHEBI:58190 + CHEBI:17815 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:64612")
    expect_equal(l$reaction_product_chebi, "CHEBI:58190 + CHEBI:17815 + CHEBI:15378")
    
    reaction <- "RHEA:78952"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + PE(14:0/16:0) => P-Ethanolamine + DG(14:0/16:0/0:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + PE(14:0/16:0)")
    expect_equal(l$reaction_product, "P-Ethanolamine + DG(14:0/16:0/0:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:64612 => CHEBI:58190 + CHEBI:17815 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:64612")
    expect_equal(l$reaction_product_chebi, "CHEBI:58190 + CHEBI:17815 + CHEBI:15378")
    
    reaction <- "RHEA:78953"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + PE(14:0/16:0) <= P-Ethanolamine + DG(14:0/16:0/0:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + PE(14:0/16:0)")
    expect_equal(l$reaction_product, "P-Ethanolamine + DG(14:0/16:0/0:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:64612 <= CHEBI:58190 + CHEBI:17815 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:64612")
    expect_equal(l$reaction_product_chebi, "CHEBI:58190 + CHEBI:17815 + CHEBI:15378")
    
    reaction <- "RHEA:78954"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "H2O + PE(14:0/16:0) <=> P-Ethanolamine + DG(14:0/16:0/0:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "H2O + PE(14:0/16:0)")
    expect_equal(l$reaction_product, "P-Ethanolamine + DG(14:0/16:0/0:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:15377 + CHEBI:64612 <=> CHEBI:58190 + CHEBI:17815 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:15377 + CHEBI:64612")
    expect_equal(l$reaction_product_chebi, "CHEBI:58190 + CHEBI:17815 + CHEBI:15378")
    
    ## sn2lpe_to_fa
    pe <- "PE(0:0/14:0)"
    substrates <- list(sn2LPE = pe)
    
    reaction <- "sn2lpe_to_fa"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(0:0/14:0) + H2O <=> FA(14:0) + H+ + Glycerophosphoethanolamine")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(0:0/14:0) + H2O")
    expect_equal(l$reaction_product, "FA(14:0) + H+ + Glycerophosphoethanolamine")
    expect_equal(l$reaction_formula_chebi, "CHEBI:65213 + CHEBI:15377 <=> CHEBI:28868 + CHEBI:15378 + CHEBI:143890")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:65213 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:28868 + CHEBI:15378 + CHEBI:143890")
    
    ## lpeo_to_peo
    lpeo <- "PE(O-16:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPEO = lpeo, AcylCoA = acylcoa)
    
    reaction <- "lpeo_to_peo"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "CoA(18:0) + PE(O-16:0/0:0) <=> CoA + PE(O-16:0/18:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "CoA(18:0) + PE(O-16:0/0:0)")
    expect_equal(l$reaction_product, "CoA + PE(O-16:0/18:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:58342 +  <=> CHEBI:57287 + CHEBI:75028") ####################
    expect_equal(l$reaction_substrate_chebi, "CHEBI:58342 + ")
    expect_equal(l$reaction_product_chebi, "CHEBI:57287 + CHEBI:75028")
    
    ## sn2mg_to_sn1mg
    sn2mg <- "MG(0:0/14:0/0:0)"
    substrates = list(sn2MG = sn2mg)
    
    reaction <- "sn2mg_to_sn1mg"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "MG(0:0/14:0/0:0) <=> MG(14:0/0:0/0:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "MG(0:0/14:0/0:0)")
    expect_equal(l$reaction_product, "MG(14:0/0:0/0:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:35759 <=> CHEBI:17389")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:35759")
    expect_equal(l$reaction_product_chebi, "CHEBI:17389")

    ## nape_to_pnae
    nape <- "NAPE(14:0/16:0/18:0)"
    substrates <- list(NAPE = nape)
    
    reaction <- "nape_to_pnae"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "NAPE(14:0/16:0/18:0) + H2O <=> PNAE(18:0) + DG(14:0/16:0/0:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "NAPE(14:0/16:0/18:0) + H2O")
    expect_equal(l$reaction_product, "PNAE(18:0) + DG(14:0/16:0/0:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:62537 + CHEBI:15377 <=> CHEBI:145538 + CHEBI:17815")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:62537 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:145538 + CHEBI:17815")
    
    ## napeo_to_nae
    napeo <- "NAPE(O-18:0/16:0/14:0)"
    substrates = list(NAPEO = napeo)
    
    reaction <- "napeo_to_nae"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "NAPE(O-18:0/16:0/14:0) + H2O <=> NAE(14:0) + PA(O-18:0/16:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "NAPE(O-18:0/16:0/14:0) + H2O")
    expect_equal(l$reaction_product, "NAE(14:0) + PA(O-18:0/16:0)")
    expect_equal(l$reaction_formula_chebi, " + CHEBI:15377 <=> CHEBI:52640 + CHEBI:73332") ##########################
    expect_equal(l$reaction_substrate_chebi, " + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:52640 + CHEBI:73332")
    
    ## pe_to_pa
    pe <- "PE(14:0/16:0)"
    substrates <- list(PE = pe)
    
    reaction <- "pe_to_pa"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(14:0/16:0) + H2O <=> PA(14:0/16:0) + Ethanolamine + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(14:0/16:0) + H2O")
    expect_equal(l$reaction_product, "PA(14:0/16:0) + Ethanolamine + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 <=> CHEBI:58608 + CHEBI:57603 + CHEBI:15378")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:64612 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, "CHEBI:58608 + CHEBI:57603 + CHEBI:15378")
    
    ## peo_to_lpeo
    peo <- "PE(O-16:0/14:0)"
    substrates <- list(PEO = peo)
    
    reaction <- "peo_to_lpeo"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(O-16:0/14:0) + H2O <=> PE(O-16:0/0:0) + FA(14:0) + H+")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(O-16:0/14:0) + H2O")
    expect_equal(l$reaction_product, "PE(O-16:0/0:0) + FA(14:0) + H+")
    expect_equal(l$reaction_formula_chebi, "CHEBI:75028 + CHEBI:15377 <=>  + CHEBI:28868 + CHEBI:15378") ########################
    expect_equal(l$reaction_substrate_chebi, "CHEBI:75028 + CHEBI:15377")
    expect_equal(l$reaction_product_chebi, " + CHEBI:28868 + CHEBI:15378")
    
    ## peo_to_napeo_sn1
    peo <- "PE(O-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    substrates <- list(PEO = peo, PC = pc) ################
    
    reaction <- "peo_to_napeo_sn1"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(O-16:0/14:0) + PC(20:0/18:0) <=> NAPE(O-16:0/14:0/20:0) + PC(0:0/18:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(O-16:0/14:0) + PC(20:0/18:0)")
    expect_equal(l$reaction_product, "NAPE(O-16:0/14:0/20:0) + PC(0:0/18:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:75028 + CHEBI:57643 <=>  + CHEBI:58168") ########################
    expect_equal(l$reaction_substrate_chebi, "CHEBI:75028 + CHEBI:57643")
    expect_equal(l$reaction_product_chebi, " + CHEBI:58168")
    
    ## peo_to_napeo_sn2
    peo <- "PE(O-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    substrates <- list(PEO = peo, PC = pc) ################
    
    reaction <- "peo_to_napeo_sn2"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(O-16:0/14:0) + PC(20:0/18:0) <=> NAPE(O-16:0/14:0/18:0) + PC(20:0/0:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(O-16:0/14:0) + PC(20:0/18:0)")
    expect_equal(l$reaction_product, "NAPE(O-16:0/14:0/18:0) + PC(20:0/0:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:75028 + CHEBI:57643 <=>  + CHEBI:58168") ########################
    expect_equal(l$reaction_substrate_chebi, "CHEBI:75028 + CHEBI:57643")
    expect_equal(l$reaction_product_chebi, " + CHEBI:58168")
    
    ## pep_to_napep_sn2
    pep <- "PE(P-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    substrates <- list(PEP = pep, PC = pc)
    
    reaction <- "pep_to_napep_sn2"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    l <- .create_list_with_template(df_reaction = df_reaction,
        template = template, reaction = reaction)
    expect_equal(l$reaction_name, "")
    expect_equal(l$reaction_formula, 
        "PE(P-16:0/14:0) + PC(20:0/18:0) <=> NAPE(P-16:0/14:0/18:0) + PC(20:0/0:0)")
    expect_equal(l$reaction_isReversible, "")
    expect_equal(l$reaction_geneAssociation, "")
    expect_equal(l$reaction_pathway, "")
    expect_equal(l$reaction_substrate, "PE(P-16:0/14:0) + PC(20:0/18:0)")
    expect_equal(l$reaction_product, "NAPE(P-16:0/14:0/18:0) + PC(20:0/0:0)")
    expect_equal(l$reaction_formula_chebi, "CHEBI:77290 + CHEBI:57643 <=> CHEBI:140451 + CHEBI:58168")
    expect_equal(l$reaction_substrate_chebi, "CHEBI:77290 + CHEBI:57643")
    expect_equal(l$reaction_product_chebi, "CHEBI:140451 + CHEBI:58168")
    
})

