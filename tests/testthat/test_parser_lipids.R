## function isolate_radyls
test_that("isolate_radyls works.", {
    lipids <- c("PC(18:0/20:4(7E,9E,11Z,14Z))", "PC(16:0/18:1(9Z))")
    lipids_isolated <- isolate_radyls(lipids)
    expect_equal(lipids_isolated[[1]][1], "18:0")
    expect_equal(lipids_isolated[[1]][2], "20:4(7E,9E,11Z,14Z)")
    expect_equal(lipids_isolated[[2]][1], "16:0")
    expect_equal(lipids_isolated[[2]][2], "18:1(9Z)")
    
    ## glycerolipids (single lipid)
    expect_equal(isolate_radyls("MG(16:0/0:0/0:0)"), 
        list(c("16:0", "0:0", "0:0")))
    expect_equal(isolate_radyls("DG(16:0/18:1(9Z)/0:0)"), 
        list(c("16:0", "18:1(9Z)", "0:0")))
    expect_equal(isolate_radyls("TG(16:0/18:1(9Z)/18:0)"), 
        list(c("16:0", "18:1(9Z)", "18:0")))
    
    ## glycerolipids (multiple lipids)
    expect_equal(isolate_radyls(
        c("MG(16:0/0:0/0:0)", "DG(16:0/18:1(9Z)/0:0)", 
            "TG(16:0/18:1(9Z)/18:0)")), 
        list(c("16:0", "0:0", "0:0"), c("16:0", "18:1(9Z)", "0:0"),
            c("16:0", "18:1(9Z)", "18:0")))
    
    ## glycerphospholipids (single lipid)
    expect_equal(isolate_radyls("PC(16:0/18:1(9Z))"), 
        list(c("16:0", "18:1(9Z)")))
    expect_equal(isolate_radyls("PC(16:0/18:3(6Z,9Z,12Z))"), 
        list(c("16:0", "18:3(6Z,9Z,12Z)")))
    expect_equal(isolate_radyls("PC(O-16:0/18:1(9Z))"), 
        list(c("O-16:0", "18:1(9Z)")))
    expect_equal(isolate_radyls("PC(P-16:0/18:1(9Z))"), 
        list(c("P-16:0", "18:1(9Z)")))
    expect_equal(isolate_radyls("PC(P-16:0/16:0(15Me))"), 
        list(c("P-16:0", "16:0(15Me)")))
    
    ## glycerphospholipids (multiple lipids)
    expect_equal(isolate_radyls(
        c("PC(16:0/18:1(9Z))", "PC(16:0/18:3(6Z,9Z,12Z))", 
            "PC(O-16:0/18:1(9Z))", "PC(P-16:0/18:1(9Z))",
            "PC(P-16:0/16:0(15Me))")), 
        list(c("16:0", "18:1(9Z)"), c("16:0", "18:3(6Z,9Z,12Z)"),
            c("O-16:0", "18:1(9Z)"), c("P-16:0", "18:1(9Z)"),
            c("P-16:0", "16:0(15Me)")))
    
    ## oxidized lipid(single lipid)
    expect_equal(isolate_radyls("PC(18:0/20:4(7E,9E,11Z,14Z)(5OH[S],6OH[R])"), 
        list(c("18:0", "20:4(7E,9E,11Z,14Z)(5OH[S],6OH[R])")))
    
    ## sphingolipids (single lipid)
    expect_equal(isolate_radyls("Cer(d18:1/20:0)"), list(c("20:0")))
    expect_equal(isolate_radyls("Cer(d18:1/20:0(2OH))"), list(c("20:0(2OH)")))
    
    ## coenzyme A (single lipid)
    expect_equal(isolate_radyls("CoA(16:1(2E))"), list(c("16:1(2E)")))
})

## function .isolate_radyls_helper
test_that(".isolate_radyls_helper works.", {
    lipid_1 <- "PC(18:0/20:4(7E,9E,11Z,14Z))"
    lipid_1_isolated <- .isolate_radyls_helper(lipid_1)
    lipid_2 <- "PC(16:0/18:1(9Z))"
    lipid_2_isolated <- .isolate_radyls_helper(lipid_2)
    expect_equal(lipid_1_isolated[1], "18:0")
    expect_equal(lipid_1_isolated[2], "20:4(7E,9E,11Z,14Z)")
    expect_equal(lipid_2_isolated[1], "16:0")
    expect_equal(lipid_2_isolated[2], "18:1(9Z)")
})