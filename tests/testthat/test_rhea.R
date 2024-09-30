## function getRheaIDsFromProteinID

test_that("getRheaIDsFromProteinID works.", {
    
    expect_equal(getRheaIDsFromProteinID("Q920L6"), 
        c("RHEA:32727", "RHEA:35315", "RHEA:39675", "RHEA:60140", "RHEA:39167",
            "RHEA:36511", "RHEA:36503", "RHEA:36523"))
    expect_equal(getRheaIDsFromProteinID("Q9H5J4"),
        c("RHEA:32727", "RHEA:35315", "RHEA:39675", "RHEA:60140", "RHEA:39167",
            "RHEA:36511", "RHEA:36503", "RHEA:36523"))
    expect_error(getRheaIDsFromProteinID("Q92"), 
         "Failed to retrieve data from UniProt.")
})