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
#' @importFrom httr GET status_code content
#' @importFrom jsonlite fromJSON
#' 
#' @examples
#' protein_id <- "Q920L6"
#' \dontrun{LipidNetworkPredictR:::getRheaIDsFromProteinID(protein_id)}
getRheaIDsFromProteinID <- function(uniprot_id) {
    
    ## construct the UniProt URL for the API request
    url <- paste0("https://rest.uniprot.org/uniprotkb/", uniprot_id, ".json")
    
    ## make the HTTP GET request
    response <- httr::GET(url)
    
    ## check if the request was successful and stop if it was not succesful
    if (httr::status_code(response) != 200) 
        stop("Failed to retrieve data from UniProt. Check the UniProt ID and try again.")
        
    ## Parse the JSON content
    data <- httr::content(response, "text", encoding = "UTF-8")
    parsed_data <- jsonlite::fromJSON(data, flatten = TRUE)
        
    ## obtain RHEA ids
    rhea_ids <- lapply(parsed_data$comments$reaction.reactionCrossReferences,
        function(parsed_data_i) {
            grep(x = parsed_data_i[, "id"], pattern = "RHEA", value = TRUE)
        }) |>
        unlist()

    ## return the RHEA ids
        rhea_ids
}