#' @title Predicting auxotrophies
#'
#' @description Performs a FBA on the input model and on the model, where the
#' lower bound (i.e. influx) of the compounds is set each (individually) to 0,
#' and growth rates are compared. Works currently only for gapseq models.
#'
#' @param mod Model of class `ModelOrg` or a named list of these objects.
#' @param compounds character. IDs of the model compound whose lower bound
#' should be set to 0 in order to test its essentiality. No suffix for the
#' external compartment is required.
#' @param min.growth double. Minimum required growth of the original model.
#' @param min.growth.fraction double. Minimum growth fraction of the model
#' without the compound relative to the original model that is required to call
#' the model prototrophic for the compound.
#' @param open.bounds logigal. If TRUE, the lower bounds of all Exchange
#' reactions are set to be open (i.e., `-COBRAR_SETTINGS("MAXIMUM")`), except
#' for the exchange reaction of the metabolite of interest.
#' @param multi.thread logical. Indicating if parallel processing of models is
#' used.
#' @param ncores integer. Number of CPUs that are used in case of parallel
#' processing. If NULL, the number of available CPUs is detected.
#'
#' @return Named numeric vector with the same length and order as `compounds`.
#' Entry of 1 indicates prototrophy and 0 auxotrophy.
#'
#' @import cobrar
#' @import parallel
#'
#' @export
predict_auxotrophies <- function(mod, compounds = NULL, min.growth = 1e-6,
                                 min.growth.fraction = 1e-4,
                                 open.bounds = TRUE,
                                 multi.thread = TRUE,
                                 ncores = NULL) {

  single_mode <- FALSE

  if(is.list(mod)) {
    if(!all(unlist(lapply(mod, class)) == "ModelOrg"))
      stop("Not all models in list are of type 'ModelOrg'")
    if(is.null(names(mod)))
      stop("Model list is not a named list.")
    if(any(is.na(names(mod))))
      stop("Not all items in the list of models are named.")
  } else {
    if(class(mod) != "ModelOrg")
      stop("mod is not of class 'ModelOrg'")
    mod <- list(mod)
    names(mod)[1] <- mod[[1]]@mod_id
    single_mode <- TRUE
  }

  # parallel processing?
  n.cores <- ifelse(multi.thread, detectCores()-1, 1)
  if(!is.null(ncores))
    n.cores <- ncores
  n.cores <- min(c(n.cores, length(mod)))
  cl <- makeCluster(max(c(1,n.cores)))

  if(is.null(compounds) || compounds[1] == "amino acids") {
    compounds <- c(Ala = "cpd00035",
                   Val = "cpd00156",
                   Met = "cpd00060",
                   Leu = "cpd00107",
                   Ile = "cpd00322",
                   Pro = "cpd00129",
                   Trp = "cpd00065",
                   Phe = "cpd00066",
                   Lys = "cpd00039",
                   Arg = "cpd00051",
                   His = "cpd00119",
                   Tyr = "cpd00069",
                   Thr = "cpd00161",
                   Glu = "cpd00023",
                   Gln = "cpd00053",
                   Gly = "cpd00033",
                   Ser = "cpd00054",
                   Cys = "cpd00084",
                   Asp = "cpd00041",
                   Asn = "cpd00132"
    )
  }

  compounds <- gsub("^EX_|_.0$", "", compounds)
  if(is.null(names(compounds)))
    names(compounds) <- compounds

  clusterExport(cl, c("compounds","min.growth","min.growth.fraction"),
                envir = environment())

  if(open.bounds) {
    mod <- lapply(mod, function(modi) {
      ind_ex <- react_pos(modi, findExchReact(modi)$react_id)
      modi@lowbnd[ind_ex] <- -COBRAR_SETTINGS("MAXIMUM")
      return(modi)
    })
  }

  auxores <- parLapply(cl, mod, fun = worker_auxo_pred)
  stopCluster(cl)

  if(single_mode == FALSE)
    return(auxores)

  return(auxores[[1]])
}

#' Get growth rate / Value of objective function
#'
#' @description Uses a simple FBA to predict the models optimal value for the
#' objective functions (most commonly the growth rate).
#'
#' @param mod Model of type `ModelOrg`
#' @param multi.thread logical. Indicating if parallel processing of models is
#' used.
#' @param ncores integer. Number of CPUs that are used in case of parallel
#' processing. If NULL, the number of available CPUs is detected.
#'
#' @return Numeric. Value of optimal value of objective function. (e.g. growth
#' rate)
#'
#' @export
get_growth <- function(mod, multi.thread = TRUE, ncores = NULL) {

  single_mode <- FALSE

  if(is.list(mod)) {
    if(!all(unlist(lapply(mod, class)) == "ModelOrg"))
      stop("Not all models in list are of type 'ModelOrg'")
    if(is.null(names(mod)))
      stop("Model list is not a named list.")
    if(any(is.na(names(mod))))
      stop("Not all items in the list of models are named.")
  } else {
    if(class(mod) != "ModelOrg")
      stop("mod is not of class 'ModelOrg'")
    mod <- list(mod)
    names(mod)[1] <- mod[[1]]@mod_id
    single_mode <- TRUE
  }

  # parallel processing?
  n.cores <- ifelse(multi.thread, detectCores()-1, 1)
  if(!is.null(ncores))
    n.cores <- ncores
  n.cores <- min(c(n.cores, length(mod)))
  cl <- makeCluster(max(c(1,n.cores)))

  sol <- parLapply(cl, mod, fun = worker_growth_pred)
  stopCluster(cl)

  if(single_mode == FALSE)
    return(sol)

  return(sol[[1]])
}

#' Get flux distribution
#'
#' @description Predict all fluxes of the constrained model.
#'
#' @param mod Model of class `ModelOrg` or a named list of these objects.
#' @param exclude.unused Exclude reactions with zero fluxes from output table.
#' @param algorithm character. Algorithm to use to calculate flux distribution.
#' One of: "fba", "pfba", "pfbaHeuristic".
#' @param multi.thread logical. Indicating if parallel processing of models is
#' used.
#' @param ncores integer. Number of CPUs that are used in case of parallel
#' processing. If NULL, the number of available CPUs is detected.
#'
#' @return A data.table with columns: `rxn` - reaction ID, `flux` - predicted flux, and additional reaction attributes.
#'
#' @import parallel
#'
#' @export
get_flux_distribution <- function(mod, exclude.unused = F, algorithm = "pfbaHeuristic",
                                  multi.thread = TRUE, ncores = NULL) {
  single_mode <- FALSE

  if(is.list(mod)) {
    if(!all(unlist(lapply(mod, class)) == "ModelOrg"))
      stop("Not all models in list are of type 'ModelOrg'")
    if(is.null(names(mod)))
      stop("Model list is not a named list.")
    if(any(is.na(names(mod))))
      stop("Not all items in the list of models are named.")
  } else {
    if(class(mod) != "ModelOrg")
      stop("mod is not of class 'ModelOrg'")
    mod <- list(mod)
    names(mod)[1] <- mod[[1]]@mod_id
    single_mode <- TRUE
  }

  # parallel processing?
  n.cores <- ifelse(multi.thread, detectCores()-1, 1)
  if(!is.null(ncores))
    n.cores <- ncores
  n.cores <- min(c(n.cores, length(mod)))
  cl <- makeCluster(max(c(1,n.cores)))
  fbaalg <- algorithm
  clusterExport(cl, c("fbaalg","exclude.unused"), envir=environment())

  flxlst <- parLapply(cl, mod, fun = worker_fba)
  stopCluster(cl)

  if(single_mode == FALSE)
    return(flxlst)

  return(flxlst[[1]])
}


#' Get produced metabolites
#'
#' @description Predicts metabolite production.
#'
#' @param mod Model object of class `ModelOrg`
#' @param algorithm character. Algorithm to use to calculate flux distribution.
#' One of: "fba", "pfba", "pfbaHeuristic".
#' @param multi.thread logical. Indicating if parallel processing of models is
#' used.
#' @param ncores integer. Number of CPUs that are used in case of parallel
#' processing. If NULL, the number of available CPUs is detected.
#'
#' @return data.table with predicted production rate and flux ranges. Columns: `ex` - exchange reaction ID,
#' `name` - exchange reaction name, `flux` - predicted production flux, `fva.min` - minimum flux according to FVA, `fva.max` - maximum flux according to FVA.
#'
#' @details Ranges (`fva.min` & `fva.max`) are based on a FVA-derived method. This method attempts to prevent cases where a
#' nutrient is taken up from the environment, transformed to another compound that is produced, but without any
#' contribution to the organism's growth rate. When providing additional
#' constraints, please note that these constraints are not (yet) applied for the
#' flux variability analysis.
#'
#' @import parallel
#'
#' @export
get_produced_metabolites <- function(mod,  algorithm = "pfbaHeuristic",
                                     multi.thread = TRUE, ncores = NULL) {
  single_mode <- FALSE

  if(is.list(mod)) {
    if(!all(unlist(lapply(mod, class)) == "ModelOrg"))
      stop("Not all models in list are of type 'ModelOrg'")
    if(is.null(names(mod)))
      stop("Model list is not a named list.")
    if(any(is.na(names(mod))))
      stop("Not all items in the list of models are named.")
  } else {
    if(class(mod) != "ModelOrg")
      stop("mod is not of class 'ModelOrg'")
    mod <- list(mod)
    names(mod)[1] <- mod[[1]]@mod_id
    single_mode <- TRUE
  }

  # parallel processing?
  n.cores <- ifelse(multi.thread, detectCores()-1, 1)
  if(!is.null(ncores))
    n.cores <- ncores
  n.cores <- min(c(n.cores, length(mod)))
  cl <- makeCluster(max(c(1,n.cores)))
  fbaalg <- algorithm
  clusterExport(cl, c("fbaalg"), envir=environment())

  gpmets <- parLapply(cl, mod, fun = worker_get_prod_mets)
  stopCluster(cl)

  if(single_mode == FALSE)
    return(gpmets)

  return(gpmets[[1]])
}


#' Get utilized metabolites
#'
#' @description Predicts metabolite consumption.
#'
#' @param mod Model object of class `ModelOrg`
#' @param algorithm character. Algorithm to use to calculate flux distribution.
#' One of: "fba", "pfba", "pfbaHeuristic".
#' @param multi.thread logical. Indicating if parallel processing of models is
#' used.
#' @param ncores integer. Number of CPUs that are used in case of parallel
#' processing. If NULL, the number of available CPUs is detected.
#'
#' @return data.table with predicted consumption rate and flux ranges. Columns: `ex` - exchange reaction ID,
#' `rxn.name` - exchange reaction name, `l` - lower bound of reaction flux, `u` - upper bound, `mtf.flux` - predicted utilization flux,
#' `flux.at.limit` - character specifying if the uptake is at it's constraint limit. "*" if at limit and "" if not at limit.
#' When providing additional constraints, please note that these constraints are
#' not (yet) applied for the flux variability analysis.
#'
#' @import parallel
#'
#' @export
get_utilized_metabolites <- function(mod, algorithm = "pfbaHeuristic",
                                     multi.thread = TRUE, ncores = NULL) {
  single_mode <- FALSE

  if(is.list(mod)) {
    if(!all(unlist(lapply(mod, class)) == "ModelOrg"))
      stop("Not all models in list are of type 'ModelOrg'")
    if(is.null(names(mod)))
      stop("Model list is not a named list.")
    if(any(is.na(names(mod))))
      stop("Not all items in the list of models are named.")
  } else {
    if(class(mod) != "ModelOrg")
      stop("mod is not of class 'ModelOrg'")
    mod <- list(mod)
    names(mod)[1] <- mod[[1]]@mod_id
    single_mode <- TRUE
  }

  # parallel processing?
  n.cores <- ifelse(multi.thread, detectCores()-1, 1)
  if(!is.null(ncores))
    n.cores <- ncores
  n.cores <- min(c(n.cores, length(mod)))
  cl <- makeCluster(max(c(1,n.cores)))
  fbaalg <- algorithm
  clusterExport(cl, c("fbaalg"), envir=environment())

  gpmets <- parLapply(cl, mod, fun = worker_get_util_mets)
  stopCluster(cl)

  if(single_mode == FALSE)
    return(gpmets)

  return(gpmets[[1]])

}

#' Get table of exchange reactions
#'
#' @description Get a table with all exchange reactions, their predicted fluxes,
#' and lower+upper bounds.
#'
#' @param mod Model object of class `ModelOrg`
#' @param algorithm character. Algorithm to use to calculate flux distribution.
#'
#' @return data.table with the following columns: `ID` - identifiers of the
#' exchange reactions, `name` - name of the exchange reaction/compound, `flux` -
#' predicted flux (negative for uptake, positive for production), `lb` - lower
#' bound for flux, `ub` - upper bound for flux.
#'
#' @export
get_exchange_reactions <- function(mod, algorithm = "pfbaHeuristic") {

  if(algorithm == "fba") {
    sol <- fba(mod)
  } else if(algorithm == "pfbaHeuristic") {
    sol <- pfbaHeuristic(mod)
  } else if(algorithm == "pfba") {
    sol <- pfba(mod)
  }

  exr <- getExchanges(mod, sol)

  exr$lb <- mod@lowbnd[react_pos(mod, exr$ID)]
  exr$ub <- mod@uppbnd[react_pos(mod, exr$ID)]

  return(as.data.table(exr))
}

#' Prediction of maximum production rates
#'
#' @description Gets the maximum production capacity of metabolite \code{met}.
#' Implementation: A new reaction with the outflow of the metabolite of interest
#' is introduced and it's flux maximized
#' using FBA.
#'
#' @param mod Model object of class `ModelOrg`
#' @param met Character vector of metabolite IDs, whose production capacity is
#' predicted.
#'
#' @return A named numeric vector with individual maximum production rates.
#' Names correspond to metabolite IDs.
#'
#' @export
get_metabolite_production_capacity <- function(mod, met) {
  met_tmp <- met[met %in% mod@met_id]

  if(length(met_tmp) == 0) {
    stop("None of the metabolite IDs is part of the model. Returning empty vector.")
    return(numeric(0L))
  }

  if(length(met_tmp) != length(met)) {
    lost_mets <- met[!(met %in% met_tmp)]
    warning(paste("Follwing metabolite IDs are not part of the model:\n",paste(lost_mets, collapse = ", ")))
  }

  met <- met_tmp

  res <- c()
  for(moi in met) {
    mod.tmp <- addReact(mod,
                        id = "OUTFLOW_TMP",
                        met = moi,
                        Scoef = -1)

    mod.tmp <- changeObjFunc(mod.tmp, react = "OUTFLOW_TMP")

    sol <- fba(mod.tmp)
    res <- c(res, sol@obj)
  }

  names(res) <- met

  return(res)
}

#------------------------------------------------------------------------------#
# Parallel worker functions
#------------------------------------------------------------------------------#

#' @import cobrar
worker_auxo_pred <- function(x) {
  if("cobrarCPLEX" %in% rownames(utils::installed.packages())) {
    require(cobrarCPLEX)
    COBRAR_SETTINGS("SOLVER","cplex")
    okcode   <- c(1,2)
  } else {
    COBRAR_SETTINGS("SOLVER","glpk")
    okcode   <- c(2,5)
  }
  # COBRAR_SETTINGS("SOLVER","glpk")
  # okcode   <- c(2,5)
  # init output
  auxo_out <- rep(NA_real_, length(compounds))
  names(auxo_out) <- names(compounds)

  # get orig growth rate
  m0_growth <- fba(x)@obj
  if(m0_growth < min.growth) {
    warning(paste0("Model ('",x@mod_id,"') has a too low or zero growth rate."))
    return(auxo_out)
  }

  # checking auxotrophies
  for(i in 1:length(compounds)) {
    imet <- compounds[i]
    ex_id <- paste0("EX_",imet,"_e0")
    if(ex_id %in% x@react_id) {
      mod_tmp <- changeBounds(x, ex_id, lb = 0)
      sol_tmp <- fba(mod_tmp)
      m1_growth <- sol_tmp@obj
      auxo_out[i] <- m1_growth / m0_growth
    } else {
      auxo_out[i] <- 1
    }
  }

  rm(mod_tmp, sol_tmp, imet, m1_growth)

  auxo_out <- ifelse(auxo_out >= min.growth.fraction, 1, 0)

  return(auxo_out)
}

#' @import cobrar
worker_growth_pred <- function(x) {
  if("cobrarCPLEX" %in% rownames(utils::installed.packages())) {
    require(cobrarCPLEX)
    COBRAR_SETTINGS("SOLVER","cplex")
    okcode   <- c(1,2)
  } else {
    COBRAR_SETTINGS("SOLVER","glpk")
    okcode   <- c(2,5)
  }

  sol <- fba(x)

  if(sol@stat %notin% okcode)
    warning("LP solution infeasible.")

  out <- sol@obj
  rm(sol)

  return(out)
}

#' @import cobrar
#' @import data.table
worker_fba <- function(x) {
  if("cobrarCPLEX" %in% rownames(utils::installed.packages())) {
    require(cobrarCPLEX)
    COBRAR_SETTINGS("SOLVER","cplex")
    okcode   <- c(1,2)
  } else {
    COBRAR_SETTINGS("SOLVER","glpk")
    okcode   <- c(2,5)
  }

  if(fbaalg == "fba") {
    sol <- fba(x)
  } else if(fbaalg == "pfbaHeuristic") {
    sol <- pfbaHeuristic(x)
  } else if(fbaalg == "pfba") {
    sol <- pfba(x)
  }

  dt <- data.table(rxn = x@react_id, flux = sol@fluxes)

  dt <- cbind(dt, data.table(x@react_attr))

  colnames(dt)[duplicated(colnames(dt))] <- paste0(colnames(dt)[duplicated(colnames(dt))],"_2")

  if(exclude.unused)
    dt <- dt[flux != 0]

  dt$equation <- printReaction(x, react = dt$rxn)
  dt[, annotation := NULL]
  dt[, CVTerms := NULL]

  return(dt)
}

#' @import cobrar
#' @import data.table
worker_get_prod_mets <- function(x) {
  if("cobrarCPLEX" %in% rownames(utils::installed.packages())) {
    require(cobrarCPLEX)
    COBRAR_SETTINGS("SOLVER","cplex")
    okcode   <- c(1,2)
  } else {
    COBRAR_SETTINGS("SOLVER","glpk")
    okcode   <- c(2,5)
  }

  # get pfba solution
  if(fbaalg == "fba") {
    sol <- fba(x)
  } else if(fbaalg == "pfbaHeuristic") {
    sol <- pfbaHeuristic(x)
  } else if(fbaalg == "pfba") {
    sol <- pfba(x)
  }

  dt  <- data.table(ex = x@react_id,
                    flux = sol@fluxes)
  dt.tmp <- copy(dt[grepl("^EX_", ex)])

  # this following two lines are there to prevent the case that a nutrient (e.g. L-Lactate)
  # from the environment is taken up, and thus enables the production of D-Lactate.
  dt.tmp[flux > 0, flux := 0]
  model.tmp <- changeBounds(x, react = dt.tmp$ex, lb = dt.tmp$flux)


  # get FVA solution
  sol.fv <- fva(model.tmp, react = x@react_id[grep("^EX_", x@react_id)])

  dt.out <- data.table(ex = sol.fv$react,
                       name = x@react_name[grep("^EX_", x@react_id)],
                       fva.min = sol.fv$min.flux,
                       fva.max = sol.fv$max.flux)
  dt.out <- dt.out[(fva.max>1e-6 & fva.min >= 0) | (fva.max > 1)]
  dt.out <- merge(dt, dt.out, by = "ex")
  dt.out <- dt.out[,.(ex, name, flux, fva.min, fva.max)]

  return(dt.out[order(-fva.max)])
}

#' @import cobrar
#' @import data.table
worker_get_util_mets <- function(x) {
  if("cobrarCPLEX" %in% rownames(utils::installed.packages())) {
    require(cobrarCPLEX)
    COBRAR_SETTINGS("SOLVER","cplex")
    okcode   <- c(1,2)
  } else {
    COBRAR_SETTINGS("SOLVER","glpk")
    okcode   <- c(2,5)
  }

  # get pfba solution
  if(fbaalg == "fba") {
    sol <- fba(x)
  } else if(fbaalg == "pfbaHeuristic") {
    sol <- pfbaHeuristic(x)
  } else if(fbaalg == "pfba") {
    sol <- pfba(x)
  }

  dt  <- data.table(ex = x@react_id,
                    flux = sol@fluxes,
                    lb = x@lowbnd)
  dt.tmp <- copy(dt[grepl("^EX_", ex)])
  dt.tmp[flux < 0, flux := 0]
  model.tmp <- changeBounds(x, react = dt.tmp$ex, ub = dt.tmp$flux)
  #model.tmp <- x

  # get FV solution
  sol.fv <- fva(model.tmp, react = x@react_id[grep("^EX_", x@react_id)])

  dt.out <- data.table(ex = x@react_id[grep("^EX_", x@react_id)],
                       name = x@react_name[grep("^EX_", x@react_id)],
                       fva.min = sol.fv$min.flux,
                       fva.max = sol.fv$max.flux)
  dt.out <- dt.out[(fva.min < -1e-4 & fva.max <= 0) | (fva.min < -1)]


  dt.out <- merge(dt.out, dt, by = "ex")
  dt.out[, at.limit := ifelse(flux <= lb*0.999, "*","")]

  dt.out <- dt.out[,.(ex, name, lb, flux, fva.min, fva.max, at.limit)]

  return(dt.out[order(fva.min)])
}
