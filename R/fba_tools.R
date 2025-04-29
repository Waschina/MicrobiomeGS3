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

#------------------------------------------------------------------------------#
# Parallel worker function
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
