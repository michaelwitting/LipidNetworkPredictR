#' @name .create_reaction
#' 
#' @title Create reaction for a given reaction of lipid metabolism
#' 
#' @description 
#' The function \code{.create_reaction} creates from a given template reaction
#' the specific reaction given a set of substrates.
#' 
#' The function returns a list of the products of the reaction (in the first
#' entry of the returned list), together with the \code{template} with the 
#' updated reactions (in the second entry of the returned list). 
#' 
#' @details 
#' If \code{template == NULL}, the function will load a template specific
#' to the \code{reaction} with mininum information on the reaction.
#' 
#' @param substrates \code{list}, with entry names equal to substrates
#' @param template \code{data.frame}
#' @param reaction \code{character(1)}
#' @param constraints \code{character}, with length \code{length(substrates)}
#' @param negate \code{logical}, with length \code{length(substrates)}
#' 
#' @export
#' @return list of length 2
#' 
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#'    and Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @examples 
#' FA <- c("FA(14:0(12Me))", "FA(16:0(14Me))", "FA(15:1(9Z)(14Me))",        
#'     "FA(17:0(16Me))", "FA(12:0(11Me))", "FA(13:0(12Me))", "FA(14:0(13Me))",
#'     "FA(15:0(14Me))", "FA(16:0(15Me))", "FA(12:0)", "FA(14:0)")
#'     
#' template <- list(reaction_name = "CoA biosynthesis", 
#'     reaction_formula = "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA",
#'     reaction_RHEA = "RHEA:15421",
#'     reaction_isReversible = "",
#'     reaction_geneAssociation = "",
#'     reaction_pathway = "Acyl-CoA synthetase")
#'  
#' ## run create_reaction (template = NULL)
#' .create_reaction(substrates = list(FA = FA), template = NULL,
#'     reaction = "RHEA:15421")
#'     
#' ## run create_reaction (template = template)
#' .create_reaction(substrates = list(FA = FA), template = template,
#'     reaction = "RHEA:15421")
.create_reaction <- function(substrates, template = NULL, 
    reaction = "RHEA:15421", 
    constraints = character(length(substrates)), 
    negate = logical(length(substrates))) {
    
    ## create template
    template <- .create_template(template = template, reaction = reaction)
    
    ## create data.frame of substrates
    ## check the validity of the colnames, only continue with the valid
    ## colnames
    cols <- .check_colnames_substrates_combinations(substrates = substrates, 
        reaction = reaction)
    substrates <- substrates[cols]
    df_substrates <- .create_substrates_combinations(substrates = substrates,
        constraints = constraints, negate = negate)
    .check_colnames_substrates_combinations(substrates = df_substrates,
        reaction = reaction)
    
    ## add products to data.frame
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    
    ## make new data.frame with reaction template
    df <- .create_df_with_template(df_reaction = df_reaction, 
        template = template, reaction = reaction)
    
    ## return results
    .create_list_reactants_with_template(df_reaction = df_reaction, 
        template = df)
}


#' @name create_reactions
#' 
#' @title Create reaction network for given reactions of lipid metabolism
#' 
#' @description 
#' The function \code{create_reactions} will create reactions following a 
#' given reaction order.
#' 
#' \code{create_reactions} will return a list containing one entry for 
#' each reaction type (RHEA id). Each reaction entry is again a list:
#' the first entry will contain a list of the products of the reaction,
#' the second entry will contain the filled \code{template} with the 
#' updated reactions.
#' 
#' @details 
#' The \code{data.frame} \code{reactions} has to contain the columns 
#' \code{order} and \code{RHEA}. It may in addition contain additional
#' columns such as \code{reaction} or \code{isReversible}. These columns will
#' be ignored by \code{create_reactions} but might be of use by the user.
#' 
#' The column \code{RHEA} will contain the ids that are used for 
#' matching the reaction type.
#' 
#' @param substrates list containing the substrates for the reaction, 
#' e.g. the initial substrate for the first reaction
#' @param reactions data.frame containing the reaction order (in the column
#' \code{"order"}) and the RHEA id (in the column \code{"RHEA"}).
#' 
#' @export
#' 
#' @return list containing the reactions
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @importFrom utils stack
#' 
#' @examples 
#' FA <- c("FA(14:0(12Me))", "FA(16:0(14Me))", "FA(15:1(9Z)(14Me))",        
#'     "FA(17:0(16Me))", "FA(12:0(11Me))", "FA(13:0(12Me))", "FA(14:0(13Me))",
#'     "FA(15:0(14Me))", "FA(16:0(15Me))", "FA(12:0)", "FA(14:0)")
#' 
#' ## create data.frame with reactions and reaction order
#' reactions <- rbind(
#'     c(1, "RHEA:15421", "M_ATP + M_CoA + M_FA <=> M_PPi + M_AMP + M_AcylCoA", FALSE),
#'     c(2, "RHEA:15325", "M_Glycerol-3-P + M_AcylCoA <=> M_CoA + M_LPA", FALSE),
#'     c(3, "RHEA:19709", "M_AcylCoA + M_LPA <=> M_CoA + M_PA", FALSE),
#'     c(4, "RHEA:27429", "M_H2O + M_PA <=> M_Pi + M_1,2-DG", FALSE)
#' )
#' reactions <- data.frame(order = reactions[, 1], RHEA = reactions[, 2],
#'     reactions = reactions[, 3], directed = reactions[, 4])
#' reactions$order <- as.numeric(reactions$order)
#' reactions$directed <- as.logical(reactions$directed)
#' 
#' ## run the function
#' create_reactions(substrates = list(FA = FA), reactions = reactions)
create_reactions <- function(substrates, reactions) {
    
    ## perform some checks for the objects
    if (!is.list(substrates)) stop("'substrates' is not a list")
    if (!is.data.frame(reactions))
        stop("'reactions' has to be a data.frame")
    if (!all(c("order", "RHEA", "reactions") %in% colnames(reactions)))
        stop("'reactions' has to contain the columns 'order', 'RHEA', and 'reactions'")
    
    ## order the reactions according to the reaction order and check if the order 
    ## does not contain any gaps
    reactions <- reactions[order(reactions$order), ]
    if (!all(reactions$order == seq_len(nrow(reactions))))
        stop("The column 'order' in 'reactions' should not contain any gaps.")
    
    ## create a list that will store the reaction results
    reaction_l <- list()
    
    for (reaction_i in reactions$order) {
        
        ## obtain the RHEA id (key) for the entry
        rhea <- reactions$RHEA[reaction_i]
        
        ## create the reaction using the given substrates from the pool of
        ## lipids (substrates is a list and will comprise the different 
        ## lipid species)
        reaction_l[[reaction_i]] <- .create_reaction(substrates = substrates, 
            reaction = rhea)
        
        ## append the latest products that might be the substrates and append
        ## to the list
        substrates <- append(substrates, reaction_l[[reaction_i]][[1]])
        
        ## combine/join the entries with identical entry names and make the
        ## entries unique (substrates is a list and represents the pool of 
        ## lipids)
        substrates <- with(utils::stack(substrates), split(values, ind)) |>
            lapply(FUN = unique)
    }
    reaction_l
}


#' @name create_reaction_adjacency_matrix
#' 
#' @title Create adjacency matrix from reactions
#' 
#' @description 
#' The function \code{create_reaction_adjacency_matrix} creates an adjacency
#' matrix connecting substrates and products that are linked by reactions.
#' 
#' @details
#' The function \code{create_reaction_adjacency_matrix} accepts the output of 
#' the \code{create_reactions} function.
#' 
#' The adjacency matrix can be used in subsequent analysis for network analysis,
#' e.g. by converting the adjacency matrix to a graph via
#' \code{igraph::graph_from_adjacency_matrix}.
#' 
#' @param reaction_l list as obtained from \code{create_reactions}
#' 
#' @export
#' 
#' @return matrix
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @examples 
#' FA <- c("FA(12:0)", "FA(14:0)", "FA(16:0)")
#' 
#' ## create data.frame with reactions and reaction order
#' reactions <- rbind(
#'     c(1, "RHEA:15421", "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA", FALSE),
#'     c(2, "RHEA:15325", "M_Glycerol-3-P + M_AcylCoA = M_CoA + M_LPA", FALSE),
#'     c(3, "RHEA:19709", "M_AcylCoA + M_LPA = M_CoA + M_PA", FALSE),
#'     c(4, "RHEA:27429", "M_H2O + M_PA = M_Pi + M_1,2-DG", FALSE)
#' )
#' reactions <- data.frame(order = reactions[, 1], RHEA = reactions[, 2],
#'     reactions = reactions[, 3], directed = reactions[, 4])
#' reactions$order <- as.numeric(reactions$order)
#' reactions$directed <- as.logical(reactions$directed)
#' 
#' ## create the list with reactions
#' reaction_l <- create_reactions(substrates = list(FA = FA), reactions = reactions)
#' 
#' ## create the adjacency matrix
#' create_reaction_adjacency_matrix(reaction_l)
create_reaction_adjacency_matrix <- function(reaction_l) {
    
    if (!is.list(reaction_l)) 
        stop("'reaction_l' has to be a list")
    
    ## obtain reaction_substrate: this information is stored in the second
    ## entry of the list in the column "reaction_substrate"
    substrates <- lapply(reaction_l, function(reaction_l_i)
        reaction_l_i[[2]][["reaction_substrate"]]) |>
        unlist()
    
    ## obtain reaction_product: this information is stored in the second
    ## entry of the list in the column "reaction_product"
    products <- lapply(reaction_l, function(reaction_l_i)
        reaction_l_i[[2]][["reaction_product"]]) |>
        unlist()
    
    ## create a matrix with rownames/colnames for all the substrates and 
    ## products
    substrates_s <- strsplit(substrates, split = " [+] ")
    products_s <- strsplit(products, split = " [+] ")
    substrates_products_names <- c(unlist(substrates_s), unlist(products_s)) |>
        unique() |> 
        sort()
    adj <- matrix(0, nrow = length(substrates_products_names),
           ncol = length(substrates_products_names),
           dimnames = list(substrates_products_names, substrates_products_names))
    
    ## fill the adjacency matrix with the information from the reactions
    for (i in seq_along(substrates_s))
        adj[substrates_s[[i]], products_s[[i]]] <- 1
    
    ## return the adjacency matrix
    adj
}
