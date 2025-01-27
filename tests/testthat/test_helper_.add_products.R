## function .create_reaction
test_that(".add_products works", {
    
    ## acyldhap_to_alkyldhap
    acyldhap <- "DHAP(18:0)"
    fao <- "FAO(16:0)"
    substrates <- list(AcylDHAP = acyldhap, FAO = fao)
    .names <- c("AcylDHAP", "FAO", "FA", "AlkylDHAP")
    .values <- c("DHAP(18:0)", "FAO(16:0)", "FA(18:0)", "DHAP(O-16:0)")
    
    reaction <- "RHEA:36171"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36172"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36173"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36174"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, 
        reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
                
    ## alkyldhap_to_lpao
    alkyldhap <- "DHAP(O-18:0)"
    substrates = list(AlkylDHAP = alkyldhap)
    .names <- c("AlkylDHAP", "sn1LPAO")
    .values <- c("DHAP(O-18:0)", "PA(O-18:0/0:0)")
    
    reaction <- "RHEA:36175"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, 
                        reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36176"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, 
                        reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36177"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, 
                        reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36178"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, 
                        reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## cerp_to_cer
    cerp <- "CerP(d16:1(3OH,4OH)(15Me)/18:0)"
    substrates <- list(CerP = cerp)
    .names <- c("CerP", "Cer")
    .values <- c("CerP(d16:1(3OH,4OH)(15Me)/18:0)", "Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    
    reaction <- "RHEA:33743"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33744"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33745"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33746"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## cdpdg_to_pgp
    cdpdg <- "CDP-DG(12:0(11Me)/14:0)"
    substrates <- list(CDPDG = cdpdg)
    .names <- c("CDPDG", "PGP")
    .values <- c("CDP-DG(12:0(11Me)/14:0)", "PGP(12:0(11Me)/14:0)")
    
    reaction <- "RHEA:12593"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:12594"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:12595"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:12596"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## cdpdg_to_pi
    cdgdg <- "CDP-DG(12:0(11Me)/14:0)"
    substrates <- list(CDPDG = cdpdg)
    .names <- c("CDPDG", "PI")
    .values <- c("CDP-DG(12:0(11Me)/14:0)", "PI(12:0(11Me)/14:0)")
    
    reaction <- "RHEA:11580"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:11581"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:11582"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:11583"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## cer_to_cerp
    cer <- "Cer(d16:1(3OH,4OH)(15Me)/18:0)"
    .names <- c("Cer", "CerP")
    .values <- c("Cer(d16:1(3OH,4OH)(15Me)/18:0)", "CerP(d16:1(3OH,4OH)(15Me)/18:0)")
    substrates <- list(Cer = cer)
    
    reaction <- "RHEA:17929"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:17930"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:17931"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:17932"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## cer_to_glccer
    cer <- "Cer(d16:1(3OH,4OH)(15Me)/18:0)"
    substrates <- list(Cer = cer)
    .names <- c("Cer", "GlcCer")
    .values <- c("Cer(d16:1(3OH,4OH)(15Me)/18:0)", "GlcCer(d16:1(3OH,4OH)(15Me)/18:0)")
    
    reaction <- "RHEA:12088"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:12089"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:12090"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:12091"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## cer_to_sm
    cer <- "Cer(d16:1(3OH,4OH)(15Me)/18:0)"
    substrates <- list(Cer = cer)
    .names <- c("Cer", "SM")
    .values <- c("Cer(d16:1(3OH,4OH)(15Me)/18:0)", "SM(d16:1(3OH,4OH)(15Me)/18:0)")
    
    reaction <- "RHEA:18765"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:18766"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:18767"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:18768"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## cl_to_lcl
    cl <- "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])"
    substrates <- list(CL = cl)
    .names <- c("CL", "LCL", "FA")
    .values <- c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])", 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])",
        "FA(18:4(6Z,9Z,12Z,15Z))")
    
    reaction <- "RHEA:32935"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32936"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32937"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32938"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## coa_to_acyldhap
    acylcoa <- "CoA(18:0)"
    substrates <- list(AcylCoA = acylcoa)
    .names <- c("AcylCoA", "AcylDHAP")
    .values <- c("CoA(18:0)", "DHAP(18:0)")
    
    reaction <- "RHEA:17657"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:17658"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:17659"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:17660"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## coa_to_ce
    acylcoa <- "CoA(18:0)"
    substrates <- list(AcylCoA = acylcoa)
    .names <- c("AcylCoA", "CE")
    .values <- c("CoA(18:0)", "CE(18:0)")
    
    reaction <- "RHEA:17729"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:17730"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:17731"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:17732"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## coa_to_FAO
    acylcoa <- "CoA(18:0)"
    substrates <- list(AcylCoA = acylcoa)
    .names <- c("AcylCoA", "FAO")
    .values <- c("CoA(18:0)", "FAO(18:0)")
    
    reaction <- "RHEA:52716"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:52717"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:52718"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:52719"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## coa_to_lpa
    acylcoa <- "CoA(18:0)"
    substrates <- list(AcylCoA = acylcoa)
    .names <- c("AcylCoA", "sn1LPA")
    .values <- c("CoA(18:0)", "PA(18:0/0:0)")
    
    reaction <- "RHEA:15325"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:15326"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:15327"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:15328"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## dg_to_sn1mg
    dg <- "DG(18:0/16:0/0:0)" 
    substrates <- list(DG = dg)
    .names <- c("DG", "sn1MG", "FA")
    .values <- c("DG(18:0/16:0/0:0)", "MG(18:0/0:0/0:0)", "FA(16:0)")
    
    reaction <- "RHEA:44712"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:44713"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:44714"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:44715"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:35663"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:35664"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:35665"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:35666"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## dg_to_sn2mg
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    .names <- c("DG", "sn2MG", "FA")
    .values <- c("DG(18:0/16:0/0:0)", "MG(0:0/16:0/0:0)", "FA(18:0)")
    
    reaction <- "RHEA:33275"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33276"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33277"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33278"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## dg_to_pa
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    .names <- c("DG", "PA")
    .values <- c("DG(18:0/16:0/0:0)", "PA(18:0/16:0)")
    
    reaction <- "RHEA:10272"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:10273"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:10274"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:10275"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## dg_to_pc
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    .names <- c("DG", "PC")
    .values <- c("DG(18:0/16:0/0:0)", "PC(18:0/16:0)")
    
    reaction <- "RHEA:32939"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32940"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32941"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32942"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## dg_to_pe
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    .names <- c("DG", "PE")
    .values <- c("DG(18:0/16:0/0:0)", "PE(18:0/16:0)")
    
    reaction <- "RHEA:32943"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32944"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32945"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32946"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## dg_to_tg
    dg <- "DG(18:0/16:0/0:0)"
    acylcoa <- "CoA(14:0)"
    substrates <- list(DG = dg, AcylCoA = acylcoa)
    .names <- c("DG", "AcylCoA", "TG")
    .values <- c("DG(18:0/16:0/0:0)", "CoA(14:0)", "TG(18:0/16:0/14:0)")
    
    reaction <- "RHEA:10868"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:10869"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:10870"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:10871"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## dgo_to_pco
    dgo <- "DG(O-18:0/16:0/0:0)"
    substrates <- list(DGO = dgo)
    .names <- c("DGO", "PCO")
    .values <- c("DG(O-18:0/16:0/0:0)", "PC(O-18:0/16:0)")
    
    reaction <- "RHEA:36179"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36180"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36181"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36182"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## dgo_to_peo
    dgo <- "DG(O-18:0/16:0/0:0)"
    substrates <- list(DGO = dgo)
    .names <- c("DGO", "PEO")
    .values <- c("DG(O-18:0/16:0/0:0)", "PE(O-18:0/16:0)")
    
    reaction <- "RHEA:36187"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36188"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36189"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36190"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## dhcer_to_cer
    dhcer <- "Cer(d16:0(3OH,4OH)(15Me)/12:0)"
    substrates <- list(DhCer = dhcer)
    .names <- c("DhCer", "Cer")
    .values <- c("Cer(d16:0(3OH,4OH)(15Me)/12:0)", "Cer(d16:1(3OH,4OH)(15Me)/12:0)") ## was: "Cer(d16:1(4E)(3OH,4OH)(15Me)/12:0)"
    
    reaction <- "RHEA:46544"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:46545"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:46546"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:46547"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## dhcer_to_dhsm
    dhcer <- "Cer(d16:0(3OH,4OH)(15Me)/12:0)"
    substrates <- list(DhCer = dhcer)
    .names <- c("DhCer", "DhSM")
    .values <- c("Cer(d16:0(3OH,4OH)(15Me)/12:0)", "SM(d16:0(3OH,4OH)(15Me)/12:0)")
    
    reaction <- "RHEA:44620"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:44621"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:44622"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:44623"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## dhsm_to_dhcer
    dhsm <- "SM(d16:0(3OH,4OH)(15Me)/12:0)" ###################################################
    substrates <- list(DhSM = dhsm)
    .names <- c("DhSM", "DhCer")
    .values <- c("SM(d16:0(3OH,4OH)(15Me)/12:0)", "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    
    reaction <- "RHEA:45300"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45301"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45302"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45303"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## fa_to_coa
    fa <- "FA(18:0)"
    substrates = list(FA = fa)
    .names <- c("FA", "AcylCoA")
    .values <- c("FA(18:0)", "CoA(18:0)")
    
    reaction <- "RHEA:15421"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:15422"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:15423"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:15424"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:38883"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:38884"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:38885"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:38886"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lcl_to_cl
    lcl <- "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])" ###################
    acylcoa <- "CoA(18:4(6Z,9Z,12Z,15Z))"
    substrates <- list(LCL = lcl, AcylCoA = acylcoa)
    .names <- c("LCL", "AcylCoA", "CL")
    .values <- c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])", 
        "CoA(18:4(6Z,9Z,12Z,15Z))",
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    
    reaction <- "RHEA:35839"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:35840"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:35841"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:35842"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lnape_to_gpnae
    lnape <- "NAPE(14:0/0:0/0:0)"
    substrates <- list(LNAPE = lnape)
    .names <- c("LNAPE", "GPNAE", "FA")
    .values <- c("NAPE(14:0/0:0/0:0)", "GPNAE(0:0)", "FA(14:0)")
    
    reaction <- "RHEA:45420"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45421"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45422"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45423"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lpa_to_pa
    lpa <- "PA(18:0/0:0)"
    acylcoa <- "CoA(14:0)"
    substrates <- list(sn1LPA = lpa, AcylCoA = acylcoa)
    .names <- c("sn1LPA", "AcylCoA", "PA")
    .values <- c("PA(18:0/0:0)", "CoA(14:0)", "PA(18:0/14:0)")
    
    reaction <- "RHEA:19709"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:19710"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:19711"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:19712"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lpao_to_pao
    lpao <- "PA(O-18:0/0:0)"
    acylcoa <- "CoA(14:0)"
    substrates <- list(sn1LPAO = lpao, AcylCoA = acylcoa)
    .names <- c("sn1LPAO", "AcylCoA", "PAO")
    .values <- c("PA(O-18:0/0:0)", "CoA(14:0)", "PA(O-18:0/14:0)")
    
    reaction <- "RHEA:36235"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36236"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36237"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36238"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn1lpc_to_fa
    sn1lpc <- "PC(14:0/0:0)"
    substrates <- list(sn1LPC = sn1lpc)
    .names <- c("sn1LPC", "FA")
    .values <- c("PC(14:0/0:0)", "FA(14:0)")
    
    reaction <- "RHEA:15177"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:15178"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:15179"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:15180"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)

    ## sn2lpc_to_fa
    sn2lpc <- "PC(0:0/14:0)"
    substrates <- list(sn2LPC = sn2lpc)
    .names <- c("sn2LPC", "FA")
    .values <- c("PC(0:0/14:0)", "FA(14:0)")
    
    reaction <- "RHEA:44696"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:44697"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:44698"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:44699"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn1lpc_to_pc
    sn1lpc <- "PC(14:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPC = sn1lpc, AcylCoA = acylcoa)
    .names <- c("sn1LPC", "AcylCoA", "PC")
    .values <- c("PC(14:0/0:0)", "CoA(18:0)", "PC(14:0/18:0)")
    
    reaction <- "RHEA:12937"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:12938"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:12939"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:12940"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn1lpe_to_fa
    pe <- "PE(14:0/0:0)"
    substrates <- list(sn1LPE = pe)
    .names <- c("sn1LPE", "FA")
    .values <- c("PE(14:0/0:0)", "FA(14:0)")
    
    reaction <- "RHEA:32967"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32968"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32969"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32970"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn2lpe_to_fa
    pe <- "PE(0:0/14:0)"
    substrates <- list(sn2LPE = pe)
    .names <- c("sn2LPE", "FA")
    .values <- c("PE(0:0/14:0)", "FA(14:0)")
    
    reaction <- "sn2lpe_to_fa"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn1lpe_to_pe
    pe <- "PE(14:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPE = pe, AcylCoA = acylcoa)
    .names <- c("sn1LPE", "AcylCoA", "PE")
    .values <- c("PE(14:0/0:0)", "CoA(18:0)", "PE(14:0/18:0)")
    
    reaction <- "RHEA:32995"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32996"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32997"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32998"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn1lpg_to_pg
    sn1lpg <- "PG(16:0/0:0)"
    acylcoa <- "CoA(14:0)"
    substrates <- list(sn1LPG = sn1lpg, AcylCoA = acylcoa)
    .names <- c("sn1LPG", "AcylCoA", "PG")
    .values <- c("PG(16:0/0:0)", "CoA(14:0)", "PG(16:0/14:0)")
    
    reaction <- "RHEA:33203"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33204"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33205"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33206"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn1lpi_to_pi
    sn1lpi <- "PI(16:0/0:0)"
    acylcoa <- "CoA(18:1(9Z))"
    substrates = list(sn1LPI = sn1lpi, AcylCoA = acylcoa)
    .names <- c("sn1LPI", "AcylCoA", "PI")
    .values <- c("PI(16:0/0:0)", "CoA(18:1(9Z))", "PI(16:0/18:1(9Z))")
    
    reaction <- "RHEA:33195"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33196"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33197"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33198"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lpeo_to_peo
    lpeo <- "PE(O-16:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPEO = lpeo, AcylCoA = acylcoa)
    .names <- c("sn1LPEO", "AcylCoA", "PEO")
    .values <- c("PE(O-16:0/0:0)", "CoA(18:0)", "PE(O-16:0/18:0)")
    
    reaction <- "lpeo_to_peo"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lpep_to_pep
    lpep <- "PE(P-16:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPEP = lpep, AcylCoA = acylcoa)
    .names <- c("sn1LPEP", "AcylCoA", "PEP")
    .values <- c("PE(P-16:0/0:0)", "CoA(18:0)", "PE(P-16:0/18:0)")
    
    reaction <- "RHEA:16245"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:16246"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:16247"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:16248"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn1mg_to_dg
    sn1mg <- "MG(14:0/0:0/0:0)"
    acylcoa <- "CoA(16:0)"
    substrates <- list(sn1MG = sn1mg, AcylCoA = acylcoa)
    .names <- c("sn1MG", "AcylCoA", "DG")
    .values <- c("MG(14:0/0:0/0:0)", "CoA(16:0)", "DG(14:0/16:0/0:0)")
    
    reaction <- "RHEA:38463"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:38464"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:38465"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:38466"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:39943"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:39944"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:39945"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:39946"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn2mg_to_dg
    sn2mg <- "MG(0:0/14:0/0:0)"
    acylcoa <- "CoA(16:0)"
    substrates <- list(sn2MG = sn2mg, AcylCoA = acylcoa)
    .names <- c("sn2MG", "AcylCoA", "DG")
    .values <- c("MG(0:0/14:0/0:0)", "CoA(16:0)", "DG(16:0/14:0/0:0)")
    
    reaction <- "RHEA:32947"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32948"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32949"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32950"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:16741"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:16742"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:16743"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:16744"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn1mg_to_fa
    sn1mg <- "MG(14:0/0:0/0:0)"
    substrates <- list(sn1MG = sn1mg)
    .names <- c("sn1MG", "FA")
    .values <- c("MG(14:0/0:0/0:0)", "FA(14:0)")
    
    reaction <- "RHEA:34019"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:34020"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:34021"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:34022"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn2mg_to_fa
    sn2mg <- "MG(0:0/14:0/0:0)"
    substrates <- list(sn2MG = sn2mg)
    .names <- c("sn2MG", "FA")
    .values <- c("MG(0:0/14:0/0:0)", "FA(14:0)")
    
    reaction <- "RHEA:32871"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32872"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32873"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32874"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)

    ## sn1mg_to_lpa
    sn1mg <- "MG(14:0/0:0/0:0)"
    substrates <- list(sn1MG = sn1mg)
    .names <- c("sn1MG", "sn1LPA")
    .values <- c("MG(14:0/0:0/0:0)", "PA(14:0/0:0)")
    
    reaction <- "RHEA:33747"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33748"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33749"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33750"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn2mg_to_sn1mg
    sn2mg <- "MG(0:0/14:0/0:0)"
    substrates = list(sn2MG = sn2mg)
    .names <- c("sn2MG", "sn1MG")
    .values <- c("MG(0:0/14:0/0:0)", "MG(14:0/0:0/0:0)")
    
    reaction <- "sn2mg_to_sn1mg"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## nae_to_fa
    nae <- "NAE(18:0)"
    substrates <- list(NAE = nae)
    .names <- c("NAE", "FA")
    .values <- c("NAE(18:0)", "FA(18:0)")
    
    reaction <- "RHEA:17505"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:17506"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:17507"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:17508"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:39995"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:39996"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:39997"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:39998"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## nape_to_lnape
    nape <- "NAPE(14:0/16:0/18:0)"
    substrates = list(NAPE = nape)
    .names <- c("NAPE", "LNAPE", "FA")
    .values <- c("NAPE(14:0/16:0/18:0)", "NAPE(14:0/0:0/18:0)", "FA(16:0)")
    
    reaction <- "RHEA:45460"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45461"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45462"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45463"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## nape_to_nae
    nape <- "NAPE(14:0/16:0/18:0)"
    substrates = list(NAPE = nape)
    .names <- c("NAPE", "NAE", "PA")
    .values <- c("NAPE(14:0/16:0/18:0)", "NAE(18:0)", "PA(14:0/16:0)")
    
    reaction <- "RHEA:33159"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33160"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33161"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33162"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## nape_to_pnae
    nape <- "NAPE(14:0/16:0/18:0)"
    substrates <- list(NAPE = nape)
    .names <- c("NAPE", "PNAE", "DG")
    .values <- c("NAPE(14:0/16:0/18:0)", "PNAE(18:0)", "DG(14:0/16:0/0:0)")
    
    reaction <- "nape_to_pnae"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## napeo_to_nae
    napeo <- "NAPE(O-18:0/16:0/14:0)"
    substrates = list(NAPEO = napeo)
    .names <- c("NAPEO", "NAE", "PAO")
    .values <- c("NAPE(O-18:0/16:0/14:0)", "NAE(14:0)", "PA(O-18:0/16:0)")
    
    reaction <- "napeo_to_nae"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pa_to_cdpdg
    pa <- "PA(14:0/16:0)"
    substrates <- list(PA = pa)
    .names <- c("PA", "CDPDG")
    .values <- c("PA(14:0/16:0)", "CDP-DG(14:0/16:0)")
    
    reaction <- "RHEA:16229"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:16230"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:16231"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:16232"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pa_to_dg
    pa <- "PA(14:0/16:0)"
    substrates = list(PA = pa)
    .names <- c("PA", "DG")
    .values <- c("PA(14:0/16:0)", "DG(14:0/16:0/0:0)")
    
    reaction <- "RHEA:27429"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:27430"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:27431"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:27432"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    
    ## pao_to_dgo
    pao <- "PA(O-14:0/16:0)"
    substrates = list(PAO = pao)
    .names <- c("PAO", "DGO")
    .values <- c("PA(O-14:0/16:0)", "DG(O-14:0/16:0/0:0)")
    
    reaction <- "RHEA:36239"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36240"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36241"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36242"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pc_to_dg
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    .names <- c("PC", "DG")
    .values <- c("PC(20:0/18:0)", "DG(20:0/18:0/0:0)")
    
    reaction <- "RHEA:10604"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:10605"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:10606"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:10607"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pc_to_sn1lpc
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    .names <- c("PC", "sn1LPC", "FA")
    .values <- c("PC(20:0/18:0)", "PC(20:0/0:0)", "FA(18:0)")
    
    reaction <- "RHEA:15801"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:15802"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:15803"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:15804"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pc_to_sn2lpc
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    .names <- c("PC", "sn2LPC", "FA")
    .values <- c("PC(20:0/18:0)", "PC(0:0/18:0)", "FA(20:0)")
    
    reaction <- "RHEA:18689"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:18690"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:18691"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:18692"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pc_to_pa
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    .names <- c("PC", "PA")
    .values <- c("PC(20:0/18:0)", "PA(20:0/18:0)")
    
    reaction <- "RHEA:14445"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:14446"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:14447"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:14448"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pc_to_ps
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    .names <- c("PC", "PS")
    .values <- c("PC(20:0/18:0)", "PS(20:0/18:0)")
    
    reaction <- "RHEA:45088"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45089"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45090"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45091"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pco_to_lpco
    pco <- "PC(O-16:0/14:0)"
    substrates <- list(PCO = pco)
    .names <- c("PCO", "sn1LPCO", "FA")
    .values <- c("PC(O-16:0/14:0)", "PC(O-16:0/0:0)", "FA(14:0)")
    
    reaction <- "RHEA:36231"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36232"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36233"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36234"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lpco_to_lpao
    sn1lpco <- "PC(O-16:0/0:0)"
    substrates <- list(sn1LPCO = sn1lpco)
    .names <- c("sn1LPCO", "sn1LPAO")
    .values <- c("PC(O-16:0/0:0)", "PA(O-16:0/0:0)")
    
    reaction <- "RHEA:39927"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:39928"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:39929"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:39930"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lpco_to_mgo
    sn1lpco <- "PC(O-16:0/0:0)"
    substrates <- list(sn1LPCO = sn1lpco)
    .names <- c("sn1LPCO", "sn1MGO")
    .values <- c("PC(O-16:0/0:0)", "MG(O-16:0/0:0/0:0)")
    
    reaction <- "RHEA:36083"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36084"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36085"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36086"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lpco_to_pco
    sn1lpco <- "PC(O-16:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPCO = sn1lpco, AcylCoA = acylcoa)
    .names <- c("sn1LPCO", "AcylCoA", "PCO")
    .values <- c("PC(O-16:0/0:0)", "CoA(18:0)", "PC(O-16:0/18:0)")
    
    reaction <- "RHEA:23992"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:23993"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:23994"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:23995"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pe_to_dg
    pe <- "PE(14:0/16:0)"
    substrates <- list(PE = pe)
    .names <- c("PE", "DG")
    .values <- c("PE(14:0/16:0)", "DG(14:0/16:0/0:0)")
    
    reaction <- "RHEA:78951"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:78952"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:78953"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:78954"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pe_to_sn1lpe
    pe <- "PE(14:0/16:0)"
    substrates = list(PE = pe)
    .names <- c("PE", "sn1LPE", "FA")
    .values <- c("PE(14:0/16:0)", "PE(14:0/0:0)", "FA(16:0)")
    
    reaction <- "RHEA:44604"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:44605"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:44606"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:44607"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pe_to_sn2lpe
    pe <- "PE(14:0/16:0)"
    substrates = list(PE = pe)
    .names <- c("PE", "sn2LPE", "FA")
    .values <- c("PE(14:0/16:0)", "PE(0:0/14:0)", "FA(16:0)")
    
    reaction <- "RHEA:44408"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:44409"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:44410"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:44411"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pe_to_nape_sn1
    pe <- "PE(14:0/16:0)"
    pc <- "PC(18:0/20:0)"
    substrates <- list(PE = pe, PC = pc)
    .names <- c("PE", "PC", "sn2LPC", "NAPE")
    .values <- c("PE(14:0/16:0)", "PC(18:0/20:0)", "PC(0:0/20:0)", "NAPE(14:0/16:0/18:0)")
    
    reaction <- "RHEA:45188"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45189"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45190"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45191"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45191"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pe_to_nape_sn2
    pe <- "PE(14:0/16:0)"
    pc <- "PC(18:0/20:0)"
    substrates <- list(PE = pe, PC = pc)
    .names <- c("PE", "PC", "sn1LPC", "NAPE")
    .values <- c("PE(14:0/16:0)", "PC(18:0/20:0)", "PC(18:0/0:0)", "NAPE(14:0/16:0/20:0)")
    
    reaction <- "RHEA:45192"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45193"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45194"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45195"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pe_to_pa
    pe <- "PE(14:0/16:0)"
    substrates <- list(PE = pe)
    .names <- c("PE", "PA")
    .values <- c("PE(14:0/16:0)", "PA(14:0/16:0)")
    
    reaction <- "pe_to_pa"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pe_to_ps
    pe <- "PE(14:0/16:0)"
    substrates <- list(PE = pe)
    .names <- c("PE", "PS")
    .values <- c("PE(14:0/16:0)", "PS(14:0/16:0)")
    
    reaction <- "RHEA:27606"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:27607"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:27608"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:27609"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## peo_to_lpeo
    peo <- "PE(O-16:0/14:0)"
    substrates <- list(PEO = peo)
    .names <- c("PEO", "sn1LPEO", "FA")
    .values <- c("PE(O-16:0/14:0)", "PE(O-16:0/0:0)", "FA(14:0)")
    
    reaction <- "peo_to_lpeo"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## peo_to_napeo_sn1
    peo <- "PE(O-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    substrates <- list(PEO = peo, PC = pc)
    .names <- c("PEO", "PC", "sn2LPC", "NAPEO") 
    .values <- c("PE(O-16:0/14:0)", "PC(20:0/18:0)", "PC(0:0/18:0)", "NAPE(O-16:0/14:0/20:0)")
    
    reaction <- "peo_to_napeo_sn1"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## peo_to_napeo_sn2
    peo <- "PE(O-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    substrates <- list(PEO = peo, PC = pc)
    .names <- c("PEO", "PC", "sn1LPC", "NAPEO") ### sn2LPC??
    .values <- c("PE(O-16:0/14:0)", "PC(20:0/18:0)", "PC(20:0/0:0)", "NAPE(O-16:0/14:0/18:0)")
    
    reaction <- "peo_to_napeo_sn2"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## peo_to_pep
    peo <- "PE(O-16:0/14:0)"
    substrates <- list(PEO = peo)
    .names <- c("PEO", "PEP")
    .values <- c("PE(O-16:0/14:0)", "PE(P-16:0/14:0)")
    
    reaction <- "RHEA:22956"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:22957"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:22958"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:22959"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pep_to_lpep
    pep <- "PE(P-16:0/14:0)"
    substrates <- list(PEP = pep)
    .names <- c("PEP", "sn1LPEP", "FA")
    .values <- c("PE(P-16:0/14:0)", "PE(P-16:0/0:0)", "FA(14:0)")
    
    reaction <- "RHEA:36195"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36196"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36197"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36198"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)

    ## lpep_to_fal
    sn1lpep <- "PE(P-16:0/0:0)"
    substrates = list(sn1LPEP = sn1lpep)
    .names <- c("sn1LPEP", "FAL")
    .values <- c("PE(P-16:0/0:0)", "FAL(16:0)")
    
    reaction <- "RHEA:16905"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:16906"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:16907"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:16908"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lpep_to_lpap
    sn1lpep <- "PE(P-16:0/0:0)"
    substrates <- list(sn1LPEP = sn1lpep)
    .names <- c("sn1LPEP", "sn1LPAP")
    .values <- c("PE(P-16:0/0:0)", "PA(P-16:0/0:0)")
    
    reaction <- "RHEA:36203"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36204"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36205"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36206"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lpep_to_mgp
    sn1lpep <- "PE(P-16:0/0:0)"
    substrates <- list(sn1LPEP = sn1lpep)
    .names <- c("sn1LPEP", "sn1MGP")
    .values <- c("PE(P-16:0/0:0)", "MG(P-16:0/0:0/0:0)")
    
    reaction <- "RHEA:36199"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36200"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36201"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:36202"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pep_to_napep_sn1
    pep <- "PE(P-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    substrates <- list(PEP = pep, PC = pc)
    .names <- c("PEP", "PC", "sn2LPC", "NAPEP")
    .values <- c("PE(P-16:0/14:0)", "PC(20:0/18:0)", "PC(0:0/18:0)", "NAPE(P-16:0/14:0/20:0)")
    
    reaction <- "RHEA:63596"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:63597"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:63598"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:63599"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pep_to_napep_sn2
    pep <- "PE(P-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    substrates <- list(PEP = pep, PC = pc)
    .names <- c("PEP", "PC", "sn1LPC", "NAPEP")
    .values <- c("PE(P-16:0/14:0)", "PC(20:0/18:0)", "PC(20:0/0:0)", "NAPE(P-16:0/14:0/18:0)")
    
    reaction <- "pep_to_napep_sn2"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pg_to_cl
    pg <- "PG(18:4(6Z,9Z,12Z,15Z)/14:0)"
    cdpdg <- "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)"
    substrates = list(PG = pg, CDPDG = cdpdg)
    .names <- c("PG", "CDPDG", "PGs2", "CDPDGs2", "CL")
    .values <- c("PG(18:4(6Z,9Z,12Z,15Z)/14:0)", 
        "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)",
        "18:4(6Z,9Z,12Z,15Z)/14:0",
        "18:4(6Z,9Z,12Z,15Z)/14:0",
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])") #########################
    
    reaction <- "RHEA:32931"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:32932"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pgp_to_pg
    pgp <- "PGP(16:0/14:0)"
    substrates <- list(PGP = pgp)
    .names <- c("PGP", "PG")
    .values <- c("PGP(16:0/14:0)", "PG(16:0/14:0)")
    
    reaction <- "RHEA:33751"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33752"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33753"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33754"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pi_to_dg
    pi <- "PI(16:0/18:1(9Z))"
    substrates <- list(PI = pi)
    .names <- c("PI", "DG")
    .values <- c("PI(16:0/18:1(9Z))", "DG(16:0/18:1(9Z))")
    
    reaction <- "RHEA:43484"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:43485"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:43486"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:43487"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pi_to_sn1lpi
    pi <- "PI(16:0/18:1(9Z))"
    substrates <- list(PI = pi)
    .names <- c("PI", "sn1LPI", "FA")
    .values <- c("PI(16:0/18:1(9Z))", "PI(16:0/0:0)", "FA(18:1(9Z))")
    
    reaction <- "RHEA:18001"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:18002"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:18003"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:18004"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## ps_to_pe
    ps <- "PS(14:0/14:0)"
    substrates = list(PS = ps)
    .names <- c("PS", "PE")
    .values <- c("PS(14:0/14:0)", "PE(14:0/14:0)")
    
    reaction <- "RHEA:20828"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:20829"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:20830"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:20831"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sm_to_cer
    sm <- "SM(d16:1(3OH,4OH)(15Me)/12:0)"
    substrates <- list(SM = sm)
    .names <- c("SM", "Cer")
    .values <- c("SM(d16:1(3OH,4OH)(15Me)/12:0)", "Cer(d16:1(3OH,4OH)(15Me)/12:0)")
    
    reaction <- "RHEA:45644"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45645"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45646"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:45647"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sphinga_to_dhcer
    acylcoa <- "CoA(12:0)"
    sph <- "SPH(d16:0(3OH,4OH)(15Me))"
    substrates <- list(AcylCoA = acylcoa, SPH = sph)
    .names <- c("AcylCoA", "SPH", "DhCer")
    .values <- c("CoA(12:0)", "SPH(d16:0(3OH,4OH)(15Me))", "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    
    reaction <- "RHEA:53424"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:53425"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:53426"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:53427"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## tg_to_dg
    tg <- "TG(18:0/16:0/14:0)"
    substrates = list(TG = tg)
    .names <- c("TG", "sn1Loss_DG", "sn1Loss_FA", "sn3Loss_DG", "sn3Loss_FA")
    .values <- c("TG(18:0/16:0/14:0)", "DG(14:0/16:0/0:0)", "FA(18:0)", "DG(18:0/16:0/0:0)", "FA(14:0)")
    
    reaction <- "RHEA:33271"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33272"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33273"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:33274"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:44864"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:44865"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:44866"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    reaction <- "RHEA:44867"
    template <- .create_template(reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    df <- .add_products(substrates = df_substrates, reaction = reaction)
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)

})

