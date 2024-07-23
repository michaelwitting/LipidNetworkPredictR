#' @name r_GPStrGen
#' 
#' @title r_GPStrGen
#' 
#' @param x \code{matrix} or \code{data.frame}
#' @param file_name \code{character(1)}
#'
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#' 
#' @importFrom utils write.table
r_GPStrGen <- function(x, file_name) {
  
  # get path of perl script
  tools_path <-  system.file("extdata/lipidmapstools/bin", file = "GPStrGen.pl", package = "LipidNetworkPredictR")
  
  # write clipboard
  utils::write.table(x, "clipboard.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

  cmd <- paste0("perl ",
                tools_path,
                " -ChainAbbrevMode Arbitrary ",
                "-mode AbbrevFileName ",
                "-o clipboard.txt")

  # create command
  shell(cmd)
}
