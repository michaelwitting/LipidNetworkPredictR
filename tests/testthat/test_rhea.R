## function getRheaIDsFromProteinID

test_that("getRheaIDsFromProteinID works.", {
    
    # expect_equal(getRheaIDsFromProteinID("Q920L6"), 
    #     c("RHEA:32727", "RHEA:32728", "RHEA:35315", "RHEA:35316", "RHEA:39675",
    #         "RHEA:39676", "RHEA:60140", "RHEA:60141", "RHEA:39167", 
    #         "RHEA:39168", "RHEA:36511", "RHEA:36512", "RHEA:36503",
    #         "RHEA:36504", "RHEA:36523", "RHEA:36524"))
    # expect_equal(getRheaIDsFromProteinID("Q9H5J4"), 
    #     c("RHEA:32727", "RHEA:35315", "RHEA:35316", "RHEA:39675", "RHEA:39676",
    #         "RHEA:60140", "RHEA:60141", "RHEA:39167", "RHEA:39168",
    #         "RHEA:36511", "RHEA:36512", "RHEA:36503", "RHEA:36504",
    #         "RHEA:36523", "RHEA:36524"))
    # expect_error(getRheaIDsFromProteinID("Q92"), 
    #     "Protein ID not found.")
})