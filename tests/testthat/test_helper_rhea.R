## function .getRheaIDsFromProteinID
test_that(".getRheaIDsFromProteinID works", {
    if (Sys.info()["sysname"] != "Darwin") {
        expect_equal(getRheaIDsFromProteinID("Q920L6"),
            c("RHEA:32727", "RHEA:35315", "RHEA:39675", "RHEA:60140", "RHEA:39167",
                "RHEA:36511", "RHEA:36503", "RHEA:36523"))
    } 
    expect_error(getRheaIDsFromProteinID(NULL),
        "Failed to retrieve data from UniProt. Check the UniProt ID and try again.")
    expect_error(getRheaIDsFromProteinID("foo"),
        "Failed to retrieve data from UniProt. Check the UniProt ID and try again.")
})
