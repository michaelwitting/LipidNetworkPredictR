## function .create_reaction
test_that(".add_reaction works", {
    
    ## acyldhap_to_alkyldhap
    acyldhap <- "DHAP(18:0)"
    fao <- "FAO(16:0)"
    substrates <- list(AcylDHAP = acyldhap, FAO = fao)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("AcylDHAP", "FAO", "FA", "AlkylDHAP")
    .values <- c("DHAP(18:0)", "FAO(16:0)", "FA(18:0)", "DHAP(O-16:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36171")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36172")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36173")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36174")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
                
    ## alkyldhap_to_lpao
    alkyldhap <- "DHAP(O-18:0)"
    substrates = list(AlkylDHAP = alkyldhap)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("AlkylDHAP", "sn1LPAO")
    .values <- c("DHAP(O-18:0)", "PA(O-18:0/0:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36175")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36176")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36177")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36178")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## cerp_to_cer
    cerp <- "CerP(16:0(3OH,4OH,15Me)/18:0)"
    substrates <- list(CerP = cerp)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("CerP", "Cer")
    .values <- c("CerP(16:0(3OH,4OH,15Me)/18:0)", "Cer(16:0(3OH,4OH,15Me)/18:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "cerp_to_cer")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## cdpdg_to_pgp
    cdpdg <- "CDP-DG(12:0(11Me)/14:0)"
    substrates <- list(CDPDG = cdpdg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("CDPDG", "PGP")
    .values <- c("CDP-DG(12:0(11Me)/14:0)", "PGP(12:0(11Me)/14:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:12593")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:12594")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:12595")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:12596")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## cdpdg_to_pi
    cdgdg <- "CDP-DG(12:0(11Me)/14:0)"
    substrates <- list(CDPDG = cdpdg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("CDPDG", "PI")
    .values <- c("CDP-DG(12:0(11Me)/14:0)", "PI(12:0(11Me)/14:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:11580")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:11581")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:11582")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:11583")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## cer_to_cerp
    cer <- "Cer(16:0(3OH,4OH,15Me)/18:0)"
    substrates <- list(Cer = cer)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("Cer", "CerP")
    .values <- c("Cer(16:0(3OH,4OH,15Me)/18:0)", "CerP(16:0(3OH,4OH,15Me)/18:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:17929")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:17930")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)   
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:17931")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)  
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:17932")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)   
    
    ## cer_to_glccer
    cer <- "Cer(16:0(3OH,4OH,15Me)/18:0)"
    substrates <- list(Cer = cer)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("Cer", "GlcCer")
    .values <- c("Cer(16:0(3OH,4OH,15Me)/18:0)", "GlcCer(16:0(3OH,4OH,15Me)/18:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:12088")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:12089")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:12090")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:12091")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## cer_to_sm
    cer <- "Cer(16:0(3OH,4OH,15Me)/18:0)"
    substrates <- list(Cer = cer)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("Cer", "SM")
    .values <- c("Cer(16:0(3OH,4OH,15Me)/18:0)", "SM(16:0(3OH,4OH,15Me)/18:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:18765")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:18766")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:18767")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:18768")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## cl_to_lcl
    cl <- "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])"
    substrates <- list(CL = cl)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("CL", "LCL", "FA")
    .values <- c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])", 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])",
        "FA(18:4(6Z,9Z,12Z,15Z))")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32935")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32936")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32937")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32938")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## coa_to_acyldhap
    acylcoa <- "CoA(18:0)"
    substrates <- list(AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("AcylCoA", "AcylDHAP")
    .values <- c("CoA(18:0)", "DHAP(18:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:17657")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:17658")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:17659")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:17660")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## coa_to_FAO
    acylcoa <- "CoA(18:0)"
    substrates <- list(AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("AcylCoA", "FAO")
    .values <- c("CoA(18:0)", "FAO(18:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:52716")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:52717")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:52718")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:52719")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## coa_to_lpa
    acylcoa <- "CoA(18:0)"
    substrates <- list(AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("AcylCoA", "sn1LPA")
    .values <- c("CoA(18:0)", "PA(18:0/0:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:15325")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:15326")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:15327")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:15328")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## dg_to_sn1mg
    dg <- "DG(18:0/16:0/0:0)" 
    substrates <- list(DG = dg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("DG", "sn1MG", "FA")
    .values <- c("DG(18:0/16:0/0:0)", "MG(18:0/0:0/0:0)", "FA(16:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44712")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44713")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44714")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44715")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:35663")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:35664")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:35665")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:35666")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## dg_to_sn2mg
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("DG", "sn2MG", "FA")
    .values <- c("DG(18:0/16:0/0:0)", "MG(0:0/16:0/0:0)", "FA(18:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:33275")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:33276")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:33277")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:33278")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## dg_to_pa
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("DG", "PA")
    .values <- c("DG(18:0/16:0/0:0)", "PA(18:0/16:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:10272")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:10273")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:10274")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:10275")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## dg_to_pc
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("DG", "PC")
    .values <- c("DG(18:0/16:0/0:0)", "PC(18:0/16:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32939")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32940")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32941")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32942")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## dg_to_pe
    dg <- "DG(18:0/16:0/0:0)"
    substrates <- list(DG = dg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("DG", "PE")
    .values <- c("DG(18:0/16:0/0:0)", "PE(18:0/16:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32943")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32944")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32945")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32946")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## dg_to_tg
    dg <- "DG(18:0/16:0/0:0)"
    acylcoa <- "CoA(14:0)"
    substrates <- list(DG = dg, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("DG", "AcylCoA", "TG")
    .values <- c("DG(18:0/16:0/0:0)", "CoA(14:0)", "TG(18:0/16:0/14:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:10868")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:10869")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:10870")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:10871")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## dgo_to_pco
    dgo <- "DG(O-18:0/16:0/0:0)"
    substrates <- list(DGO = dgo)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("DGO", "PCO")
    .values <- c("DG(O-18:0/16:0/0:0)", "PC(O-18:0/16:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36179")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36180")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36181")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36182")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## dgo_to_peo
    dgo <- "DG(O-18:0/16:0/0:0)"
    substrates <- list(DGO = dgo)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("DGO", "PEO")
    .values <- c("DG(O-18:0/16:0/0:0)", "PE(O-18:0/16:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36187")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36188")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36189")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36190")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## dhcer_to_cer
    dhcer <- "Cer(d16:0(3OH,4OH)(15Me)/12:0)"
    substrates <- list(DhCer = dhcer)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("DhCer", "Cer")
    .values <- c("Cer(d16:0(3OH,4OH)(15Me)/12:0)", "Cer(d16:1(4E)(3OH,4OH)(15Me)/12:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "dhcer_to_cer")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## dhcer_to_dhsm
    dhcer <- "Cer(16:1(3OH,4OH,15Me)/12:0)"
    substrates <- list(DhCer = dhcer)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("DhCer", "DhSM")
    .values <- c("Cer(16:1(3OH,4OH,15Me)/12:0)", "SM(16:1(3OH,4OH,15Me)/12:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44620")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)  
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44621")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)  
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44622")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)  
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44623")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)  
    
    ## dhsm_to_dhcer
    dhsm <- "SM(16:1(3OH,4OH,15Me)/12:0)" ###################################################
    substrates <- list(DhSM = dhsm)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("DhSM", "DhCer")
    .values <- c("SM(16:1(3OH,4OH,15Me)/12:0)", "Cer(16:1(3OH,4OH,15Me)/12:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45300")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45301")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45302")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45303")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## fa_to_coa
    fa <- "FA(18:0)"
    substrates = list(FA = fa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("FA", "AcylCoA")
    .values <- c("FA(18:0)", "CoA(18:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:15421")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:15422")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:15423")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:15424")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:38883")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:38884")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:38885")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:38886")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lcl_to_cl
    lcl <- "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])" ###################
    acylcoa <- "CoA(18:4(6Z,9Z,12Z,15Z))"
    substrates <- list(LCL = lcl, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("LCL", "AcylCoA", "CL")
    .values <- c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])", 
        "CoA(18:4(6Z,9Z,12Z,15Z))",
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:35839")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:35840")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:35841")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:35842")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lnape_to_gpnae
    lnape <- "NAPE(14:0/0:0/0:0)"
    substrates <- list(LNAPE = lnape)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("LNAPE", "GPNAE", "FA")
    .values <- c("NAPE(14:0/0:0/0:0)", "GPNAE(0:0)", "FA(14:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45420")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45421")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45422")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45423")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lpa_to_pa
    lpa <- "PA(18:0/0:0)"
    acylcoa <- "CoA(14:0)"
    substrates <- list(sn1LPA = lpa, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn1LPA", "AcylCoA", "PA")
    .values <- c("PA(18:0/0:0)", "CoA(14:0)", "PA(18:0/14:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:19709")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:19710")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:19711")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:19712")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lpao_to_pao
    lpao <- "PA(O-18:0/0:0)"
    acylcoa <- "CoA(14:0)"
    substrates <- list(sn1LPAO = lpao, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn1LPAO", "AcylCoA", "PAO")
    .values <- c("PA(O-18:0/0:0)", "CoA(14:0)", "PA(O-18:0/14:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36235")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36236")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36237")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36238")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn1lpc_to_fa
    sn1lpc <- "PC(14:0/0:0)"
    substrates <- list(sn1LPC = sn1lpc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn1LPC", "FA")
    .values <- c("PC(14:0/0:0)", "FA(14:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:15177")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:15178")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:15179")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:15180")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)

    ## sn2lpc_to_fa
    sn2lpc <- "PC(0:0/14:0)"
    substrates <- list(sn2LPC = sn2lpc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn2LPC", "FA")
    .values <- c("PC(0:0/14:0)", "FA(14:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44696")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44697")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44698")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44699")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn1lpc_to_pc
    sn1lpc <- "PC(14:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPC = sn1lpc, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn1LPC", "AcylCoA", "PC")
    .values <- c("PC(14:0/0:0)", "CoA(18:0)", "PC(14:0/18:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:12937")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:12938")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:12939")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:12940")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn1lpe_to_fa
    pe <- "PE(14:0/0:0)"
    substrates <- list(sn1LPE = pe)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn1LPE", "FA")
    .values <- c("PE(14:0/0:0)", "FA(14:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32967")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32968")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32969")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32970")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn2lpe_to_fa
    pe <- "PE(0:0/14:0)"
    substrates <- list(sn2LPE = pe)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn2LPE", "FA")
    .values <- c("PE(0:0/14:0)", "FA(14:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "sn2lpe_to_fa")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn1lpe_to_pe
    pe <- "PE(14:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPE = pe, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn1LPE", "AcylCoA", "PE")
    .values <- c("PE(14:0/0:0)", "CoA(18:0)", "PE(14:0/18:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32995")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32996")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32997")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32998")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn1lpi_to_pi
    sn1lpi <- "PI(16:0/0:0)"
    acylcoa <- "CoA(18:1(9Z))"
    substrates = list(sn1LPI = sn1lpi, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn1LPI", "AcylCoA", "PI")
    .values <- c("PI(16:0/0:0)", "CoA(18:1(9Z))", "PI(16:0/18:1(9Z))")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:33195")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:33196")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:33197")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:33198")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lpeo_to_peo
    lpeo <- "PE(O-16:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPEO = lpeo, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn1LPEO", "AcylCoA", "PEO")
    .values <- c("PE(O-16:0/0:0)", "CoA(18:0)", "PE(O-16:0/18:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "lpeo_to_peo")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lpep_to_pep
    lpep <- "PE(P-16:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPEP = lpep, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn1LPEP", "AcylCoA", "PEP")
    .values <- c("PE(P-16:0/0:0)", "CoA(18:0)", "PE(P-16:0/18:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:16245")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:16246")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:16247")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:16248")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn1mg_to_dg
    sn1mg <- "MG(14:0/0:0/0:0)"
    acylcoa <- "CoA(16:0)"
    substrates <- list(sn1MG = sn1mg, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn1MG", "AcylCoA", "DG")
    .values <- c("MG(14:0/0:0/0:0)", "CoA(16:0)", "DG(14:0/16:0/0:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:38463")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:38464")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:38465")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:38466")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:39943")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:39944")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:39945")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:39946")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn2mg_to_dg
    sn2mg <- "MG(0:0/14:0/0:0)"
    acylcoa <- "CoA(16:0)"
    substrates <- list(sn2MG = sn2mg, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn2MG", "AcylCoA", "DG")
    .values <- c("MG(0:0/14:0/0:0)", "CoA(16:0)", "DG(16:0/14:0/0:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32947")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32948")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32949")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32950")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:16741")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:16742")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:16743")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:16744")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)

    ## sn1mg_to_fa
    sn1mg <- "MG(14:0/0:0/0:0)"
    substrates <- list(sn1MG = sn1mg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn1MG", "FA")
    .values <- c("MG(14:0/0:0/0:0)", "FA(14:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:34019")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:34020")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:34021")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:34022")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn2mg_to_fa
    sn2mg <- "MG(0:0/14:0/0:0)"
    substrates <- list(sn2MG = sn2mg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn2MG", "FA")
    .values <- c("MG(0:0/14:0/0:0)", "FA(14:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32871")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32872")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32873")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32874")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn1mg_to_lpa
    sn1mg <- "MG(14:0/0:0/0:0)"
    substrates <- list(sn1MG = sn1mg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn1MG", "sn1LPA")
    .values <- c("MG(14:0/0:0/0:0)", "PA(14:0/0:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:33747")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:33748")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:33749")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:33750")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sn2mg_to_sn1mg
    sn2mg <- "MG(0:0/14:0/0:0)"
    substrates = list(sn2MG = sn2mg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn2MG", "sn1MG")
    .values <- c("MG(0:0/14:0/0:0)", "MG(14:0/0:0/0:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "sn2mg_to_sn1mg")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## nae_to_fa
    nae <- "NAE(18:0)"
    substrates <- list(NAE = nae)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("NAE", "FA")
    .values <- c("NAE(18:0)", "FA(18:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:17505")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:17506")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:17507")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:17508")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:39995")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:39996")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:39997")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:39998")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## nape_to_lnape
    nape <- "NAPE(14:0/16:0/18:0)"
    substrates = list(NAPE = nape)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("NAPE", "LNAPE", "FA")
    .values <- c("NAPE(14:0/16:0/18:0)", "NAPE(14:0/0:0/18:0)", "FA(16:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "nape_to_lnape")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## nape_to_nae
    nape <- "NAPE(14:0/16:0/18:0)"
    substrates = list(NAPE = nape)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("NAPE", "NAE", "PA")
    .values <- c("NAPE(14:0/16:0/18:0)", "NAE(18:0)", "PA(14:0/16:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "nape_to_nae")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## nape_to_pnae
    nape <- "NAPE(14:0/16:0/18:0)"
    substrates <- list(NAPE = nape)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("NAPE", "PNAE", "DG")
    .values <- c("NAPE(14:0/16:0/18:0)", "PNAE(18:0)", "DG(14:0/16:0/0:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "nape_to_pnae")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## napeo_to_nae
    napeo <- "NAPE(O-18:0/16:0/14:0)"
    substrates = list(NAPEO = napeo)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("NAPEO", "NAE", "PAO")
    .values <- c("NAPE(O-18:0/16:0/14:0)", "NAE(14:0)", "PA(O-18:0/16:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "napeo_to_nae")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pa_to_cdpdg
    pa <- "PA(14:0/16:0)"
    substrates <- list(PA = pa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PA", "CDPDG")
    .values <- c("PA(14:0/16:0)", "CDP-DG(14:0/16:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:16229")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:16230")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:16231")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:16232")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pa_to_dg
    pa <- "PA(14:0/16:0)"
    substrates = list(PA = pa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PA", "DG")
    .values <- c("PA(14:0/16:0)", "DG(14:0/16:0/0:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:27429")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:27430")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:27431")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:27432")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pao_to_dgo
    pa0 <- "PA(O-14:0/16:0)"
    substrates = list(PAO = pa0)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PAO", "DGO")
    .values <- c("PA(O-14:0/16:0)", "DG(O-14:0/16:0/0:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36239")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36240")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36241")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36242")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pc_to_dg
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PC", "DG")
    .values <- c("PC(20:0/18:0)", "DG(20:0/18:0/0:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:10604")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:10605")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:10606")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:10607")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)

    ## pc_to_sn1lpc
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PC", "sn1LPC", "FA")
    .values <- c("PC(20:0/18:0)", "PC(20:0/0:0)", "FA(18:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:15801")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:15802")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:15803")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:15804")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pc_to_sn2lpc
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PC", "sn2LPC", "FA")
    .values <- c("PC(20:0/18:0)", "PC(0:0/18:0)", "FA(20:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:18689")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:18690")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:18691")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:18692")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pc_to_pa
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PC", "PA")
    .values <- c("PC(20:0/18:0)", "PA(20:0/18:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:14445")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:14446")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:14447")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:14448")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pc_to_ps
    pc <- "PC(20:0/18:0)"
    substrates <- list(PC = pc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PC", "PS")
    .values <- c("PC(20:0/18:0)", "PS(20:0/18:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45088")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45089")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45090")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45091")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)

    ## pco_to_lpco
    pco <- "PC(O-16:0/14:0)"
    substrates <- list(PCO = pco)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PCO", "sn1LPCO", "FA")
    .values <- c("PC(O-16:0/14:0)", "PC(O-16:0/0:0)", "FA(14:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36231")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36232")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36233")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36234")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lpco_to_lpao
    sn1lpco <- "PC(O-16:0/0:0)"
    substrates <- list(sn1LPCO = sn1lpco)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn1LPCO", "sn1LPAO")
    .values <- c("PC(O-16:0/0:0)", "PA(O-16:0/0:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:39927")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:39928")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:39929")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:39930")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lpco_to_mgo
    sn1lpco <- "PC(O-16:0/0:0)"
    substrates <- list(sn1LPCO = sn1lpco)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn1LPCO", "sn1MGO")
    .values <- c("PC(O-16:0/0:0)", "MG(O-16:0/0:0/0:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36083")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36084")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36085")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36086")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lpco_to_pco
    sn1lpco <- "PC(O-16:0/0:0)"
    acylcoa <- "CoA(18:0)"
    substrates <- list(sn1LPCO = sn1lpco, AcylCoA = acylcoa)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn1LPCO", "AcylCoA", "PCO")
    .values <- c("PC(O-16:0/0:0)", "CoA(18:0)", "PC(O-16:0/18:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:23992")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:23993")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:23994")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:23995")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pe_to_dg
    pe <- "PE(14:0/16:0)"
    substrates <- list(PE = pe)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PE", "DG")
    .values <- c("PE(14:0/16:0)", "DG(14:0/16:0/0:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:78951")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:78952")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:78953")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:78954")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pe_to_sn1lpe
    pe <- "PE(14:0/16:0)"
    substrates = list(PE = pe)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PE", "sn1LPE", "FA")
    .values <- c("PE(14:0/16:0)", "PE(14:0/0:0)", "FA(16:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44604")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44605")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44606")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44607")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pe_to_sn2lpe
    pe <- "PE(14:0/16:0)"
    substrates = list(PE = pe)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PE", "sn2LPE", "FA")
    .values <- c("PE(14:0/16:0)", "PE(0:0/14:0)", "FA(16:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44408")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44409")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44410")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44411")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pe_to_nape_sn1
    pe <- "PE(14:0/16:0)"
    pc <- "PC(18:0/20:0)"
    substrates <- list(PE = pe, PC = pc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PE", "PC", "sn2LPC", "NAPE")
    .values <- c("PE(14:0/16:0)", "PC(18:0/20:0)", "PC(0:0/20:0)", "NAPE(14:0/16:0/18:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45188")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45189")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45190")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45191")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pe_to_nape_sn2
    pe <- "PE(14:0/16:0)"
    pc <- "PC(18:0/20:0)"
    substrates <- list(PE = pe, PC = pc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PE", "PC", "sn1LPC", "NAPE")
    .values <- c("PE(14:0/16:0)", "PC(18:0/20:0)", "PC(18:0/0:0)", "NAPE(14:0/16:0/20:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45192")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45193")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45194")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45195")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pe_to_pa
    pe <- "PE(14:0/16:0)"
    substrates <- list(PE = pe)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PE", "PA")
    .values <- c("PE(14:0/16:0)", "PA(14:0/16:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "pe_to_pa")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pe_to_ps
    pe <- "PE(14:0/16:0)"
    substrates <- list(PE = pe)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PE", "PS")
    .values <- c("PE(14:0/16:0)", "PS(14:0/16:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:27606")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:27607")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:27608")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:27609")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## peo_to_lpeo
    peo <- "PE(O-16:0/14:0)"
    substrates <- list(PEO = peo)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PEO", "sn1LPEO", "FA")
    .values <- c("PE(O-16:0/14:0)", "PE(O-16:0/0:0)", "FA(14:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "peo_to_lpeo")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## peo_to_napeo_sn1
    peo <- "PE(O-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    substrates <- list(PEO = peo, PC = pc) ################
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PEO", "PC", "sn2LPC", "NAPEO") 
    .values <- c("PE(O-16:0/14:0)", "PC(20:0/18:0)", "PC(0:0/18:0)", "NAPE(O-16:0/14:0/20:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "peo_to_napeo_sn1")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## peo_to_napeo_sn2
    peo <- "PE(O-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    substrates <- list(PEO = peo, PC = pc) ################
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PEO", "PC", "sn1LPC", "NAPEO") ### sn2LPC??
    .values <- c("PE(O-16:0/14:0)", "PC(20:0/18:0)", "PC(20:0/0:0)", "NAPE(O-16:0/14:0/18:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "peo_to_napeo_sn2")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## peo_to_pep
    peo <- "PE(O-16:0/14:0)"
    substrates <- list(PEO = peo)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PEO", "PEP")
    .values <- c("PE(O-16:0/14:0)", "PE(P-16:0/14:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:22956")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:22957")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:22958")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:22959")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pep_to_lpep
    pep <- "PE(P-16:0/14:0)"
    substrates <- list(PEP = pep)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PEP", "sn1LPEP", "FA")
    .values <- c("PE(P-16:0/14:0)", "PE(P-16:0/0:0)", "FA(14:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36195")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36196")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36197")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36198")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)                     

    ## lpep_to_fal
    sn1lpep <- "PE(P-16:0/0:0)"
    substrates = list(sn1LPEP = sn1lpep)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn1LPEP", "FAL")
    .values <- c("PE(P-16:0/0:0)", "FAL(16:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:16905")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:16906")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:16907")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:16908")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lpep_to_lpap
    sn1lpep <- "PE(P-16:0/0:0)"
    substrates <- list(sn1LPEP = sn1lpep)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn1LPEP", "sn1LPAP")
    .values <- c("PE(P-16:0/0:0)", "PA(P-16:0/0:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36203")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36204")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36205")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36206")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## lpep_to_mgp
    sn1lpep <- "PE(P-16:0/0:0)"
    substrates <- list(sn1LPEP = sn1lpep)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("sn1LPEP", "sn1MGP")
    .values <- c("PE(P-16:0/0:0)", "MG(P-16:0/0:0/0:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36199")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36200")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36201")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:36202")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pep_to_napep_sn1
    pep <- "PE(P-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    substrates <- list(PEP = pep, PC = pc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PEP", "PC", "sn2LPC", "NAPEP")
    .values <- c("PE(P-16:0/14:0)", "PC(20:0/18:0)", "PC(0:0/18:0)", "NAPE(P-16:0/14:0/20:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "pep_to_napep_sn1")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pep_to_napep_sn2
    pep <- "PE(P-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    substrates <- list(PEP = pep, PC = pc)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PEP", "PC", "sn1LPC", "NAPEP")
    .values <- c("PE(P-16:0/14:0)", "PC(20:0/18:0)", "PC(20:0/0:0)", "NAPE(P-16:0/14:0/18:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "pep_to_napep_sn2")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pg_to_cl
    pg <- "PG(18:4(6Z,9Z,12Z,15Z)/14:0)"
    cdpdg <- "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)"
    substrates = list(PG = pg, CDPDG = cdpdg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PG", "CDPDG", "PGs2", "CDPDGs2", "CL")
    .values <- c("PG(18:4(6Z,9Z,12Z,15Z)/14:0)", 
        "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)",
        "18:4(6Z,9Z,12Z,15Z)/14:0",
        "18:4(6Z,9Z,12Z,15Z)/14:0",
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])") #########################
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32931")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32932")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32933")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:32934")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pgp_to_pg
    pgp <- "PGP(16:0/14:0)"
    substrates <- list(PGP = pgp)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PGP", "PG")
    .values <- c("PGP(16:0/14:0)", "PG(16:0/14:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:33751")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:33752")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:33753")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:33754")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pi_to_dg
    pi <- "PI(16:0/18:1(9Z))"
    substrates <- list(PI = pi)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PI", "DG")
    .values <- c("PI(16:0/18:1(9Z))", "DG(16:0/18:1(9Z))")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:43484")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:43485")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:43486")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:43487")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## pi_to_sn1lpi
    pi <- "PI(16:0/18:1(9Z))"
    substrates <- list(PI = pi)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PI", "sn1LPI", "FA")
    .values <- c("PI(16:0/18:1(9Z))", "PI(16:0/0:0)", "FA(18:1(9Z))")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:18001")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:18002")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:18003")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:18004")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## ps_to_pe
    ps <- "PS(14:0/14:0)"
    substrates = list(PS = ps)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("PS", "PE")
    .values <- c("PS(14:0/14:0)", "PE(14:0/14:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:20828")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:20829")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:20830")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:20831")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)                               
    
    ## sm_to_cer
    sm <- "SM(16:0(3OH,4OH,15Me)/12:0)"
    substrates <- list(SM = sm)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("SM", "Cer")
    .values <- c("SM(16:0(3OH,4OH,15Me)/12:0)", "Cer(16:0(3OH,4OH,15Me)/12:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45644")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45645")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45646")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:45647")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    
    ## sphinga_to_dhcer
    acylcoa <- "CoA(12:0)" 
    substrates <- list(AcylCoA = acylcoa) #######################################
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("AcylCoA", "DhCer")
    .values <- c("CoA(12:0)", "Cer(16:0(3OH,4OH,15Me)/12:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:53424")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:53425")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:53426")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:53427")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values) 
    
    ## tg_to_dg
    tg <- "TG(18:0/16:0/14:0)"
    substrates = list(TG = tg)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates)
    .names <- c("TG", "sn1Loss_DG", "sn1Loss_FA", "sn3Loss_DG", "sn3Loss_FA")
    .values <- c("TG(18:0/16:0/14:0)", "DG(14:0/16:0/0:0)", "FA(18:0)", "DG(18:0/16:0/0:0)", "FA(14:0)")
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:33271")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:33272") 
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:33273")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:33274")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values) 
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44864")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44865")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44866")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values)
    df <- .add_products(substrates = df_substrates, 
        reaction = "RHEA:44867")
    expect_equal(names(df), .names)
    expect_equal(as.character(unlist(df)), .values) 

})
