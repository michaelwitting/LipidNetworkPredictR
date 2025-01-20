#' @name obtain_ChEBI_of_feature
#' 
#' @title Obtain ChEBI ids of features
#' 
#' @description
#' Lipid species are attributed to classes. These classes can be identified via 
#' the ChEBI ids in \code{reactions_l}. The function  obtains the corresponding
#' ChEBI ids for the features and returns a \code{data.frame} that contains 
#' the mappings between lipid species and lipid classes.
#' 
#' Mappings between the metabolites and ChEBI ids will be done via the data 
#' structure of \code{reaction_l}, as returned by \code{create_reactions}. 
#' 
#' @details
#' The function \code{obtain_ChEBI_of_feature} accepts the output of 
#' the \code{create_reactions} function. From \code{create_reactions}, the
#' ChEBI ids will be inferred.
#' 
#' @param reaction_l list as obtained from \code{create_reactions}
#' 
#' @export
#' 
#' @return data.frame
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @importFrom stringi stri_trim
#'
#' @examples 
#' FA <- c("FA(12:0)", "FA(14:0)", "FA(16:0)")
#' 
#' ## create data.frame with reactions and reaction order
#' reactions <- rbind(
#'     c(1, "RHEA:15421", "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA", FALSE),
#'     c(2, "RHEA:15325", "M_Glycerol-3-P + M_AcylCoA = M_CoA + M_LPA", FALSE),
#'     c(3, "RHEA:19709", "M_AcylCoA + M_LPA = M_CoA + M_PA", FALSE)
#' )
#' reactions <- data.frame(order = reactions[, 1], reaction_RHEA = reactions[, 2],
#'     reaction_formula = reactions[, 3], directed = reactions[, 4])
#' reactions$order <- as.numeric(reactions$order)
#' reactions$directed <- as.logical(reactions$directed)
#' 
#' ## create the list with reactions
#' reaction_l <- create_reactions(substrates = list(FA = FA), reactions = reactions)
#' 
#' ## run the function
#' obtain_ChEBI_of_feature(reaction_l = reaction_l)
obtain_ChEBI_of_feature <- function(reaction_l) {
    network <- lapply(seq_len(length(reaction_l)), function(i) {
        feature <- reaction_l[[i]][[2]][["reaction_formula"]] |>
            ## split only + when it is not precided by H, Fe2 or Fe3
            strsplit(split = "(?<!H)(?<!Fe2)(?<!Fe3)\\+|[=]|[=>]|[<=]|[<=>]", 
                perl = TRUE) |>
            unlist() |>
            strsplit(split = "^2 ") |>
            unlist() |>
            stringi::stri_trim()
        chebi <- reaction_l[[i]][[2]][["reaction_formula_chebi"]] |>
            ## split only + when it is not precided by H, Fe2 or Fe3
            strsplit(split = "(?<!H)(?<!Fe2)(?<!Fe3)\\+|[=]|[=>]|[<=]|[<=>]", 
                perl = TRUE) |>
            unlist() |>
            strsplit(split = "^2 ") |>
            unlist() |>
            stringi::stri_trim()
        data.frame(target = feature, source = chebi)
    }) |>
        do.call(what = "rbind")
    
    ## remove empty rows
    network <- network[
        !(network[["target"]] == "" & network[["source"]] == ""), ]
    
    ## order according (first priority) to column target and (with second
    ## priority) to column source, then remove the duplicate values
    network <- network[order(network[["target"]], network[["source"]]), ] 
    
    ## remove duplicated entries
    target_source <- paste(network[["target"]], network[["source"]])
    network <- network[!duplicated(target_source), ]
    
    ## return the object
    network
}


#' @name run_decouple
#' 
#' @title Run decouple
#' 
#' @description 
#' The function \code{run_decouple} runs \code{decouple} on \code{scores}
#' of lipid features. The scores can be e.g. logFC- or t-values. Mappings 
#' between the metabolites and ChEBI ids will be done via the data structure
#' of \code{reaction_l}, as returned by \code{create_reactions}. The function
#' \code{run_decouple} will run the following models from \code{decoupleR}:
#' \code{ulm}, \code{mlm}, \code{wsum}, and \code{wmean}.  
#' 
#' @details
#' The function \code{run_decouple} accepts the output of 
#' the \code{create_reactions} function. From \code{create_reactions}, the
#' ChEBI ids will be inferred.
#'
#' The function returns different type of values as results depending on the
#' model: 
#' \itemize{
#'     \item{ulm: t-values of the fitted linear model}
#'     \item{mlm: t-values of the fitted linear model}
#'     \item{wsum: scores of each target are multiplied
#'         by its associated weight and then summed}
#'     \item{norm_wsum: z-scores from a null distribution obtained by permutation 
#'         of random target features and calculation of weighted summed scores}
#'     \item{corr_wsum: corrected estimates of wsum scores obtained by 
#'         multiplying wsum scores by -log10(obtained empirical p-values)}
#'     \item{wmean: scores of each target are multiplied by its associated 
#'         weight and then averaged}
#'     \item{norm_wmean: z-scores from a null distribution obtained by 
#'         permutation of random target features and calculation of weighted
#'         averaged scores}
#'     \item{corr_wmean: corrected estimates of wmean scores obtained by 
#'         multiplying wmean scores by the -log10(obtained empirical p-values)}
#' }
#' 
#' The values correspond to the (de)regulation values for the lipid classes
#' represented as ChEBI ids.
#' 
#' @param reaction_l list as obtained from \code{create_reactions}
#' @param scores data.frame containing statistics for each feature
#' @param col_feature character of length 1, specifying the column in 
#' \code{scores} that contains the feature names
#' @param col_score character of length 1, specifying the column in 
#' \code{scores} that contains the scores
#' @param ... further arguments passed to \code{args} in \code{decouple} for 
#' the models \code{mlm}, \code{ulm}, \code{wsum}, \code{wmean}
#' 
#' @export
#' 
#' @return tibble
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @importFrom stringi stri_trim
#' @importFrom dplyr filter
#' @importFrom rlang sym
#' @importFrom tibble add_column
#' @importFrom decoupleR decouple run_mlm run_ulm run_wsum run_wmean
#'
#' @examples 
#' FA <- c("FA(12:0)", "FA(14:0)", "FA(16:0)")
#' 
#' ## create data.frame with reactions and reaction order
#' reactions <- rbind(
#'     c(1, "RHEA:15421", "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA", FALSE),
#'     c(2, "RHEA:15325", "M_Glycerol-3-P + M_AcylCoA = M_CoA + M_LPA", FALSE),
#'     c(3, "RHEA:19709", "M_AcylCoA + M_LPA = M_CoA + M_PA", FALSE)
#' )
#' reactions <- data.frame(order = reactions[, 1], reaction_RHEA = reactions[, 2],
#'     reaction_formula = reactions[, 3], directed = reactions[, 4])
#' reactions$order <- as.numeric(reactions$order)
#' reactions$directed <- as.logical(reactions$directed)
#' 
#' ## create the list with reactions
#' reaction_l <- create_reactions(substrates = list(FA = FA), reactions = reactions)
#'
#' ## create scores (simulate scores via rnorm with sd = 2)
#' scores <- rbind(
#'     c("AMP", 0.94446911),
#'     c("ATP", 0.78505230),
#'     c("CoA", 2.90333236),
#'     c("CoA(12:0)", 0.13577818),
#'     c("CoA(14:0)", 0.64465176),
#'     c("CoA(16:0)", 0.26128554),
#'     c("FA(12:0)", -1.98870069),
#'     c("FA(14:0)", -2.17736671),
#'     c("FA(16:0)", -2.41113650),
#'     c("Glycerol-3-P", -0.84185713 ),
#'     c("PA(12:0/0:0)", -2.12108629),
#'     c("PA(12:0/12:0)", 1.60395612),
#'     c("PA(12:0/14:0)", 0.07693279),
#'     c("PA(12:0/16:0)", 2.46639841),
#'     c("PA(14:0/0:0)", -1.73630919),
#'     c("PA(14:0/12:0)", 1.32713790),
#'     c("PA(14:0/14:0)", 3.78234917),
#'     c("PA(14:0/16:0)", 1.88168721),
#'     c("PA(16:0/0:0)", -0.84712629),
#'     c("PA(16:0/12:0)", 4.16900272),
#'     c("PA(16:0/14:0)", 2.07352961),
#'     c("PA(16:0/16:0)", 4.48556057), 
#'     c("PPi", -0.28139118)
#' )
#' scores <- data.frame(feature = scores[, 1], score = scores[, 2])
#' scores$feature <- as.character(scores$feature)
#' scores$score <- as.numeric(scores$score)
#' 
#' ## run the function
#' run_decouple(reaction_l = reaction_l, scores = scores)
run_decouple <- function(reaction_l, scores, col_feature = "feature", 
    col_score = "score", ...) {
    
    if (!is.list(reaction_l)) 
        stop("'reaction_l' has to be a list")
    
    if (!is.data.frame(scores))
        stop("'scores' has to be a data.frame")
    
    ## obtain the markers and prepare the object for decoupleR
    args_fct <- list(...)
    
    ## assign rownames and only keep the column "score"
    rownames(scores) <- scores[[col_feature]]
    scores <- scores[, col_score, drop = FALSE]
    
    ## lipid species are attributed to classes, these can be identified via 
    ## the ChEBI ids in reactions_l, obtain the corresponding CheBI ids for the 
    ## features, create a data.frame that contains the mappings between lipid
    ## species and lipid classes;
    ## create mapping between feature and ChEBI ids and add columns likelihood 
    ## and mor
    network <- obtain_ChEBI_of_feature(reaction_l = reaction_l) |>
        dplyr::filter(!duplicated(!!rlang::sym("target"))) |>
        tibble::add_column(likelihood = 1, mor = 1)
    
    ## obtain the argument names of the specific run_* functions and match
    ## with the ... arguments, these lists will be passed to args in decouple
    mlm_list <- args_fct[names(args_fct) %in%  names(formals("run_mlm"))]
    ulm_list <- args_fct[names(args_fct) %in%  names(formals("run_ulm"))]
    wsum_list <- args_fct[names(args_fct) %in%  names(formals("run_wsum"))]
    wmean_list <- args_fct[names(args_fct) %in%  names(formals("run_wmean"))]

    ## run the decoupleR functions and return
    decoupleR::decouple(mat = scores,
       network = network, .source = "source", .target = "target",
       minsize = 2, consensus_score = FALSE,
       statistics = c("mlm", "ulm", "wsum", "wmean"),
       args = list(
           mlm = mlm_list,
           ulm = ulm_list,
           wsum = wsum_list,
           wmean = wmean_list)
       )
}

#' @name plot_scores
#' 
#' @title Plot the deregulation scores of ChEBI identifiers
#' 
#' @description
#' The function \code{plot_scores} creates a barplot  that 
#' visualizes the score values returned by 
#' \code{run_decouple}. The function will return a \code{gtable} object that 
#' shows a barplot for each model returned by \code{run_decouple}. 
#' 
#' @details
#' The function requires the argument \code{scores} that can be created via 
#' \code{run_decouple}.
#' 
#' @param scores tibble, output from \code{run_decouple}
#' 
#' @return gtable
#' 
#' @export
#' 
#' @importFrom dplyr filter
#' @importFrom rlang sym
#' @importFrom ggplot2 ggplot aes geom_bar theme_minimal theme 
#' @importFrom ggplot2 element_text ylab
#' @importFrom gridExtra grid.arrange
#' 
#' @examples
#' FA <- c("FA(12:0)", "FA(14:0)", "FA(16:0)")
#' 
#' ## create data.frame with reactions and reaction order
#' reactions <- rbind(
#'     c(1, "RHEA:15421", "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA", FALSE),
#'     c(2, "RHEA:15325", "M_Glycerol-3-P + M_AcylCoA = M_CoA + M_LPA", FALSE),
#'     c(3, "RHEA:19709", "M_AcylCoA + M_LPA = M_CoA + M_PA", FALSE)
#' )
#' reactions <- data.frame(order = reactions[, 1], reaction_RHEA = reactions[, 2],
#'     reaction_formula = reactions[, 3], directed = reactions[, 4])
#' reactions$order <- as.numeric(reactions$order)
#' reactions$directed <- as.logical(reactions$directed)
#' 
#' ## create the list with reactions
#' reaction_l <- create_reactions(substrates = list(FA = FA), reactions = reactions)
#'
#' ## create scores (simulate scores via rnorm with sd = 2)
#' scores <- rbind(
#'     c("AMP", 0.94446911),
#'     c("ATP", 0.78505230),
#'     c("CoA", 2.90333236),
#'     c("CoA(12:0)", 0.13577818),
#'     c("CoA(14:0)", 0.64465176),
#'     c("CoA(16:0)", 0.26128554),
#'     c("FA(12:0)", -1.98870069),
#'     c("FA(14:0)", -2.17736671),
#'     c("FA(16:0)", -2.41113650),
#'     c("Glycerol-3-P", -0.84185713 ),
#'     c("PA(12:0/0:0)", -2.12108629),
#'     c("PA(12:0/12:0)", 1.60395612),
#'     c("PA(12:0/14:0)", 0.07693279),
#'     c("PA(12:0/16:0)", 2.46639841),
#'     c("PA(14:0/0:0)", -1.73630919),
#'     c("PA(14:0/12:0)", 1.32713790),
#'     c("PA(14:0/14:0)", 3.78234917),
#'     c("PA(14:0/16:0)", 1.88168721),
#'     c("PA(16:0/0:0)", -0.84712629),
#'     c("PA(16:0/12:0)", 4.16900272),
#'     c("PA(16:0/14:0)", 2.07352961),
#'     c("PA(16:0/16:0)", 4.48556057), 
#'     c("PPi", -0.28139118)
#' )
#' scores <- data.frame(feature = scores[, 1], score = scores[, 2])
#' scores$feature <- as.character(scores$feature)
#' scores$score <- as.numeric(scores$score)
#' 
#' ## run run_decouple
#' scores <- run_decouple(reaction_l = reaction_l, scores = scores)
#' 
#' ## run the function
#' plot_scores(scores)
plot_scores <- function(scores) {
    
    ## check scores object
    if (!"run_id" %in% colnames(scores)) 
        stop("column 'run_id' not in 'scores'")
    if (!"statistic" %in% colnames(scores)) 
        stop("column 'statistic' not in 'scores'")
    if (!"source" %in% colnames(scores)) 
        stop("column 'source' not in 'scores'")
    if (!"condition" %in% colnames(scores)) 
        stop("column 'condition' not in 'scores'")
    if (!"score" %in% colnames(scores)) 
        stop("column 'object' not in 'scores'")
    if (!"p_value" %in% colnames(scores)) 
        stop("column 'p_value' not in 'scores'")
    
    ## factorize the source column (used for x-axis)
    scores[["source"]] <- factor(scores[["source"]], 
        levels = unique(scores[["source"]]))
    
    ## iterate through unique values of scores[["statistic"]] and store the 
    ## ggplot object in the list g_l
    .statistic <- scores[["statistic"]] |>
        unique()
    
    ## create an empty list to store results in
    g_l <- list()
    
    for (.statistic_i in .statistic) {
        
        ## do the actual plotting
        
        g <- scores |>
            dplyr::filter(get("statistic") == .statistic_i) |> 
            ggplot2::ggplot(mapping = ggplot2::aes(y = !!rlang::sym("score"), 
                x = !!rlang::sym("source"))) +
            ggplot2::geom_bar(
                ggplot2::aes(fill = source, col = source), 
                position = "dodge", stat = "identity", color = "black")
        
        ## adjust the y-label depending on .statistic_i
        g <- g + 
            ggplot2::ylab(paste(.statistic_i, "score"))
        
        ## adjust the theme
        g <- g + ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(size = 5, 
                angle = 90, hjust = 1)) +
            ggplot2::theme(legend.position = "none")
        
        ## assign to g to an entry in g_l
        g_l[[.statistic_i]] <- g
    }
    
    ## assemble the final plot
    do.call("grid.arrange", args = list(grobs = g_l))
}

