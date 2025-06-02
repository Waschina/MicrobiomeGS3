#' Get table of metabolites
#'
#' @description Produces a table with the model's metabolites (ID, name,
#' chemical formula, charge).
#'
#' @param mod Model of class `ModelOrg`
#'
#' @return data.table with the columns: (`met.id`) Metabolit IDs, (`met.name`)
#' Metabolite names, (`formula`) Chemical formulae of metabolites, (`charge`)
#' Charge of metabolites.
#'
#' @export
get_metabolite_table <- function(mod) {

  dt <- data.table(met.id   = mod@met_id,
                   met.name = mod@met_name,
                   formula  = mod@met_attr$chemicalFormula,
                   charge   = mod@met_attr$charge)

  return(dt)
}


#' Get table of reactions
#'
#' @description Produces a table with the model's reactions (ID, name, equation).
#'
#' @param mod Model of class `ModelOrg`
#'
#' @return data.table with the columns: (`react.id`) Reactions IDs, (`react.name`)
#' Reaction names, (`equation`) the reactions' equations. If models were built
#' with gapseq the additional columns `EC`, `gs.origin`, `keggID`, and `metacycID`
#' are also part of the output table.
#'
#' @export
get_reaction_table <- function(mod) {

  dt <- data.table(react.id   = mod@react_id,
                   react.name = mod@react_name,
                   equation   = printReaction(mod, mod@react_id))

  # in case of gapseq models
  if("gs.origin" %in% colnames(mod@react_attr)) {
    dt$EC        <- mod@react_attr$ec
    dt$gs.origin <- mod@react_attr$gs.origin
    dt$keggID    <- mod@react_attr$keggID
    dt$metacycID <- mod@react_attr$biocycID
  }

  return(dt)

}
