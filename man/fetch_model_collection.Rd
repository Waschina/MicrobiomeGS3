% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch_model_collection.R
\name{fetch_model_collection}
\alias{fetch_model_collection}
\title{Fetch a collection of gapseq models}
\usage{
fetch_model_collection(
  model.dir,
  IDs = NULL,
  file.type = "model",
  subset = NULL,
  entries = NULL,
  multi.thread = TRUE,
  ncores = NULL
)
}
\arguments{
\item{model.dir}{character. Path to the directory that contains the models'
RDS files.}

\item{IDs}{character vector that specifies which models should be retrieved
from the collection. If `NULL` all models from the collection will be
retrieved.}

\item{file.type}{character. Select which gapseq output file should be read.
One of 'model' (gap-filled model), 'draft' (draft network), 'reactions',
'pathways', 'transporters', 'medium' (predicted growth medium).}

\item{subset}{integer. For testing purposes you can choose the maximum number
of model files to be read.}

\item{entries}{In case 'file.type' is "pathways" or "reactions", the argument
can be used to limit the output to specific pathways or reactions. If `NULL`,
all entries are returned.}

\item{multi.thread}{logical. Indicating if parallel processing of models is
used.}

\item{ncores}{integer. Number of CPUs that are used in case of parallel
processing. If NULL, the number of available CPUs is detected.}
}
\value{
A named (IDs) list of with elements of class `ModelOrg` (in case of
file.type is 'model' or 'draft') or elements of class `data.table` otherwise.
}
\description{
Fetches a collection of genome-scale metabolic models
reconstructed by gapseq.
}
