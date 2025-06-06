% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_analysis.R
\name{get_metabolite_table}
\alias{get_metabolite_table}
\title{Get table of metabolites}
\usage{
get_metabolite_table(mod)
}
\arguments{
\item{mod}{Model of class `ModelOrg`}
}
\value{
data.table with the columns: (`met.id`) Metabolit IDs, (`met.name`)
Metabolite names, (`formula`) Chemical formulae of metabolites, (`charge`)
Charge of metabolites.
}
\description{
Produces a table with the model's metabolites (ID, name,
chemical formula, charge).
}
