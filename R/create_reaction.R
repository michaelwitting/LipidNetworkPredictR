#'
#' @param substrates list
#' @param template 
#' @param constraints
create_reaction <- function(substrates, template = NULL, reaction = "fa_to_coa", 
    constraints = character(length(substrates)), 
    negate = logical(length(substrates))) {
    
    
    ## create template
    template <- .create_template(template = template, reaction = reaction)
    
    ## create data.frame of substrates
    df_substrates <- .create_substrates_combinations(substrates = substrates, 
        constraints = constraints, negate = negate)
    .check_colnames_substrates_combinations(df_substrates, reaction = reaction)
    
    ## add products to data.frame
    df_reactions <- .add_products(substrates = df_substrates, reaction = reaction)
    .check_colnames_substrates_products_df(df_reactions, reaction = reaction)
    
    
    ## make new data.frame with reaction template
    df <- .create_df_with_template(df_reactions = df_reactions, 
        template = template)
    
    ## return results
    
    list(unique(df_reactions$PAs), df)
    
}