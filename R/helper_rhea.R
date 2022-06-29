#' @name getRheaIDsFromProteinID
#' 
#' @title Get RHEA IDs from protein ID
#' 
#' @description 
#' The function \code{getRheaIDsFromProteinID} returns associated 
#' RHEA IDs for a given UNIPROT ID. The function will accept IDs in the 
#' form of UNIPROT IDs corresponding to the 'Entry' column in Uniprot, 
#' e.g. "Q920L6" for fatty acid synthase from \emph{Rattus norvegicus} or
#' "Q9H5J4" for fatty acid synthase from \emph{Homo sapiens}.
#' 
#' @details
#' The function will return an error if the protein ID is not found in 
#' UNIPROT. The error might be linked to internet connection problems or 
#' because the protein ID is incorrect.
#' 
#' @param protein_id \code{character(1)}, UNIPROT ID, e.g. \code{"Q9H5J4"}
#' 
#' @return vector
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @importFrom methods is
#' @importFrom utils URLencode read.csv
#' @importFrom httr GET
#' 
#' @examples
#' protein_id <- "Q920L6"
#' \dontrun{LipidNetworkPredictR:::getRheaIDsFromProteinID(protein_id)}
getRheaIDsFromProteinID <- function(protein_id) {
    
    ## check the protein_id argument for validity
    if (!(methods::is(protein_id, "character") & length(protein_id) == 1))
        stop("'protein_id' has to be of type 'character' and length 1")

    ## define the information to retrieve from UNIPROT
    columns <- c("id,rhea-id")
    uniprot_url <- "http://www.uniprot.org/uniprot/"

    ## check if internet connection exists 
    final_url <- paste0(uniprot_url, protein_id,".xml")
    request <- tryCatch(httr::GET(final_url), error = function(e) NULL)
    if (length(request) == 0)
        stop("Internet connection problem occured.")

    ## retrieve the information
    protein_url <- paste0("?query=accession:", protein_id, 
                                            "&format=tab&columns=", columns)
    request_url <- paste0(uniprot_url, protein_url)
    request_url <- utils::URLencode(request_url)

    ## prepare the information to be returned
    protein_data <- tryCatch(
        utils::read.csv(request_url, header = TRUE, sep = "\t"), 
        error = function(e) NULL)
    if (is.null(protein_data))
        stop("Protein ID not found.")
    rhea_id <- protein_data[, "Rhea.ID"]
    rhea_id <- strsplit(rhea_id, split = ";| ")[[1]]
    rhea_id <- rhea_id[grep(x = rhea_id, pattern = "RHEA:")]

    ## return the RHEA ids
    rhea_id
}
