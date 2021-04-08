r_GPStrGen <- function(x, file_name) {
  
  # get path of perl script
  tools_path <-  system.file("extdata/lipidmapstools/bin", file = "GPStrGen.pl", package = "wormLipidPredictR")
  
  # write clipboard
  write.table(x, "clipboard.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

  cmd <- paste0("perl ",
                tools_path,
                " -ChainAbbrevMode Arbitrary ",
                "-mode AbbrevFileName ",
                "-o clipboard.txt")

  # create command
  shell(cmd)

}
