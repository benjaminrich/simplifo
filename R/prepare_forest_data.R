    md <- lo <- hi <- covar <- covval <- ymin <- ymax <- NULL
#' Prepare data for forest plotting
#'
#' @param model A function that expects two arguments:
#'  (1) a vector of model parameters;
#'  (2) a vector of covariates.
#' The function should return a `list` of parameters that will be shown in
#' the forest plot. See example for details.
#' @param params A `data.frame`. The first row consists of the final
#' parameter estimates from the model, while the remaining rows represent a
#' sample from the uncertainty distribution which can be obtained either by
#' bootstrap resampling or by sampling from a multivariate normal distribution
#' using the estimated variance-covariance matrix of the parameter vector.
#' @param covariate.values A named `list`. The names correspond to the
#' names of the covariate that are expected by the `model` function, and
#' each element is a vector of covariate values to include in the foreset plot.
#' @param reference.values A named `list`. The names correspond to the
#' names of the covariate that are expected by the `model` function, and
#' each element is a single covariate values which represents the reference
#' value for that covariate.
#' @return A `data.table` containing the data needed for the forest plot
#' in a convenient format.
#' @seealso [forest_plot()].
#' @examples
#' covariate.values <- list(
#'     sex = c("Male", "Female"),
#'     wt = c(40, 70, 120),
#'     crcl = c(50, 100, 150))
#' 
#' reference.values <- list(
#'     sex = "Male",
#'     wt = 70,      # kg
#'     crcl = 100)   # mL/min
#' 
#' model <- function(param, covar) {
#'     with(c(as.list(param), as.list(covar)), {
#'         cli <- cl *
#'             exp(cl_sex * (sex == "Female")) *
#'             ((wt/70)^cl_wt) *
#'             ((crcl/90)^cl_crcl)
#' 
#'         vci <- v * ((wt/70)^v_wt)
#' 
#'         amt <- 10 * wt  # 10 units/kg
#' 
#'         tmax <- linpk:::Tmax.oral1cpt(cli, vci, ka)
#'         t.obs <- sort(unique(c(tmax, seq(0, 8, 0.1))))  # 0-8h
#' 
#'         y <- linpk::pkprofile(t.obs, cl=cli, vc=vci, ka=ka, dose=list(amt=amt))
#'         s <- linpk::secondary(y)
#' 
#'         list(cl=cli, vc=vci, cmax=s$Cmax, auc=s$AUC)
#'     })
#' }
#' 
#' theta <- data.table::fread(text="
#'        cl         v        ka     cl_wt    cl_sex   cl_crcl      v_wt 
#' 10.303500 12.201500  4.243550  0.583708 -0.255737  0.891369  0.857272 
#' ")
#' 
#' 
#' # Sample from the variance-covariance matrix of theta
#' covmat <- read.table(header=TRUE, text="
#'             theta1       theta2       theta3       theta4       theta5       theta6       theta7
#' theta1  0.09795280  0.010994600  1.33304e-03 -0.015920800 -1.07415e-02  0.008765400  0.001717700
#' theta2  0.01099460  0.089091200 -1.93659e-04  0.002583560 -8.04885e-04  0.006237480  0.007393160
#' theta3  0.00133304 -0.000193659  5.01075e-03 -0.001949220 -4.14722e-05  0.001420430 -0.001423240
#' theta4 -0.01592080  0.002583560 -1.94922e-03  0.011589700  2.87094e-03 -0.000247892  0.000493268
#' theta5 -0.01074150 -0.000804885 -4.14722e-05  0.002870940  2.45663e-03 -0.000131245  0.000453010
#' theta6  0.00876540  0.006237480  1.42043e-03 -0.000247892 -1.31245e-04  0.014281200  0.000416080
#' theta7  0.00171770  0.007393160 -1.42324e-03  0.000493268  4.53010e-04  0.000416080  0.015352000
#' ")
#' 
#' params <- data.table::as.data.table(mvtnorm::rmvnorm(1000, as.numeric(theta), as.matrix(covmat)))
#' names(params) <- names(theta)
#' params <- rbind(theta, params) # First row is the final point estimates
#' 
#' fodat <- prepare_forest_data(model, params, covariate.values, reference.values)
#' 
#' @export
prepare_forest_data <- function(model, params, covariate.values, reference.values) {

    covar <- covval <- param <- NULL

    params <- data.table::setDT(params)

    expand.modelframe <- function(..., rv, covcol="covar") {
        args <- list(...)
        df <- lapply(args, function(x) x[[1]])
        df[names(rv)] <- rv
        res <- lapply(seq_along(rv), function(i) {
            df[[covcol]] <- names(rv)[i]
            df[[names(rv)[i]]] <- args[[names(rv)[i]]]
            data.table::as.data.table(df)
        })
        rbindlist(res)
    }

    mf <- do.call(expand.modelframe, c(covariate.values, list(rv=reference.values)))

    yref <- data.table::as.data.table(reference.values)
    yref <- yref[, model(params[1], .BY), by=yref]
    yref <- yref[, names(reference.values) := NULL]
    yref <- data.table::melt(yref, id.vars=NULL, measure.vars=names(yref), variable.name="param", value.name="refval")

    f <- function(covar) {
        y <- params[, model(.SD, covar), by=1:nrow(params)][, nrow := NULL]
        qy <- y[-1, lapply(.SD, stats::quantile, probs=c(0.5, 0.025, 0.975))]
        dat <- data.table::transpose(rbind(y[1], qy))
        dat <- stats::setNames(dat, c("est", "md", "lo", "hi"))
        dat[, param := names(y)]
    }
    temp <- mf[, f(.BY), by=mf]
    fodat <- merge(temp, yref)

    data.table::setkey(fodat, param)

    g <- function(x, y) {
        as.character(do.call(switch, c(list(as.character(x)), y)))
    }
    fodat[, covval := g(covar, .SD), by=1:nrow(fodat)]

    v1 <- unique(fodat$covval)
    v2 <- suppressWarnings(as.numeric(as.character(v1)))
    v3 <- c(as.character(sort(v2[!is.na(v2)])), as.character(v1[is.na(v2)]))
    fodat[, covval := factor(covval, levels=unique(unlist(covariate.values)))]

    fodat
}


