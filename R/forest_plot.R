#' @import data.table
#' @import ggplot2
NULL

#' Forest plot
#'
#' Create a forest plot from prepared data.
#'
#' @details `forest_plot_main()` returns just the main plot, and
#' `forest_plot_table()` returns just the table portion of the plot. Sometimes
#' it can be useful to have these separately.
#'
#' @param fodat A `data.frame` in the format produced by [prepare_forest_data()].
#' @param param The parameter to plot.
#' @param covars The covariates to plot.
#' @param lab Label for the x-axis.
#' @param lim Limits for the x-axis.
#' @param scale Can be either "relative" or "absolute". Relative scale means
#' that values will be normalized by their reference value, while absolute
#' scale means values remain unchanged.
#' @param refline Value of reference line.
#' @param refrange Limits for the reference range.
#' @param logscale Use log-scale for the x-axis?
#' @param widths A numeric vector of length 2 specifying the relative widths of
#' the main plot and the table containing the numeric values. By default, the
#' split is 70:30 for the main plot and the table.
#' @return An object of class `forest_plot`. This is a list of 2 components,
#' `main_plot` and `table_plot`, which are both `ggplot` objects. It can be
#' plotted and printed (with the same result).
#' @seealso [prepare_forest_data()].
#' @examples
#' fodat <- sample_forest_data
#' 
#' # Use nice labels
#' covar.labels <- mappings::mapping(c(
#'         sex  = "Sex",
#'         wt   = "Body Weight (kg)",
#'         crcl = "Creatinine Clearance (mL/min)"))
#' 
#' param.labels = mappings::mapping(c(
#'         cmax = "C[max]",
#'         auc  = "AUC['0-8h']",
#'         cl   = "CL",
#'         vc   = "Vc"))
#' 
#' my_forest_plot <- function(.param, .covar=c("sex", "wt", "crcl")) {
#'   fodatplot <- data.table::setDT(fodat)[(param == .param) & (covar %in% .covar)]
#'   fodatplot[, covar := covar.labels(covar)]
#'   lab <- param.labels(unique(fodatplot$param))
#'   lab <- parse(text=paste("`Fold Change in`", lab, "`Relative to Reference`", sep="~"))
#'   forest_plot(fodatplot, scale="relative", lab=lab, logscale=TRUE)
#' }
#' 
#' my_forest_plot("cl")
#' 
#' my_forest_plot("vc", "wt")
#' 
#' my_forest_plot("cmax")
#' 
#' my_forest_plot("auc")
#' 
#' @export
forest_plot <- function(
    fodat,
    param    = unique(fodat$param)[1],
    covars   = NULL,
    lab      = NULL,
    lim      = NULL,
    scale    = c("relative", "absolute"),
    refline  = if (scale=="relative") 1.0  else NULL,
    refrange = if (scale=="relative") c(0.8, 1.25) else NULL,
    logscale = FALSE,
    widths   = c(0.7, 0.3))
{
    scale <- match.arg(scale)

    fodat1 <- data.table::as.data.table(fodat)

    # Keep only 1 parameter
    .param <- param
    fodat1 <- fodat1[param == .param]

    # Keep only specific covariates
    if (is.null(covars)) {
        # Only covariates that have more than 1 value by default
        covars <- table(fodat1$covar)
        covars <- covars[covars > 1]
        covars <- names(covars)
    }
    fodat1 <- fodat1[covar %in% covars]

    fodat1 <- droplevels(fodat1)

    if (is.null(lab)) {
        lab <- unique(fodat1$param)
        if (scale == "relative") {
            lab <- parse(text=paste("`Fold Change in`", lab, "`Relative to Reference`", sep="~"))
        }
    }

    main_plot  <- forest_plot_main(fodat=fodat1, lab=lab, lim=lim,
        refline=refline, refrange=refrange, logscale=logscale)

    table_plot <- forest_plot_table(fodat=fodat1)

    structure(
        list(
            main_plot  = main_plot,
            table_plot = table_plot
        ),
        class  = "forest_plot",
        widths = widths
    )
}

#' @export
plot.forest_plot <- function(x, ...) {
    egg::ggarrange(x$main_plot, x$table_plot, nrow=1, widths=attr(x, "widths"))
}

#' @export
print.forest_plot <- function(x, ...) {
    plot.forest_plot(x, ...)
}

#' @rdname forest_plot
#' @export
forest_plot_main <- function(fodat, lab=NULL, lim=lim, scale=c("relative", "absolute"), refline=NULL, refrange=NULL, logscale=FALSE) {

    scale <- match.arg(scale)

    est <- md <- lo <- hi <- covar <- covval <- ymin <- ymax <- x <- NULL


    yax <- unique(fodat[, .(breaks=mappings::cf(covar, covval), labels=covval)])


    main_plot <-
        ggplot2::ggplot(fodat, ggplot2::aes(y=mappings::cf(covar, covval))) +
        ggplot2::labs(x=lab, y=NULL) +
        ggplot2::facet_grid(covar ~ ., scales="free", space="free", switch="y") +
        ggplot2::scale_y_discrete(expand=ggplot2::expansion(add=1), breaks=yax$breaks, labels=yax$labels, limits=rev)

    if (!is.null(refrange)) {
        main_plot <- main_plot +
            ggplot2::geom_ribbon(data=data.frame(x=refrange, ymin= -Inf, ymax=Inf),
                ggplot2::aes(x=x, y=NA, ymin=ymin, ymax=ymax), color=NA, fill="gray90")
            #geom_rect(xmin=min(refrange), xmax=max(refrange), ymin= -Inf, ymax=Inf, color=NA, fill="gray90")
    }

    if (!is.null(refline)) {
        main_plot <- main_plot +
            ggplot2::geom_vline(xintercept=refline, color="gray40", linetype="dotted")
    }


    if (logscale) {
        main_plot <- main_plot +
            ggplot2::scale_x_log10(breaks=c(0.5, 0.66, 0.8, 1.0, 1.25, 1.5, 2.0))
    } else if (!is.null(lim)) {
        main_plot <- main_plot +
            ggplot2::coord_cartesian(xlim=lim) + 
            ggplot2::scale_x_continuous(breaks=seq(0, max(lim), 0.2))
    }

    if (scale == "relative") {
        main_plot <- main_plot +
            ggplot2::geom_pointrange(ggplot2::aes(x=est/refval, xmin=lo/refval, xmax=hi/refval), color="#3A70B6") +
            ggplot2::geom_point(ggplot2::aes(x=est/refval), color="#3A70B6", size=3)
    } else {
        main_plot <- main_plot +
            ggplot2::geom_pointrange(ggplot2::aes(x=est, xmin=lo, xmax=hi), color="#3A70B6") +
            ggplot2::geom_point(ggplot2::aes(x=est), color="#3A70B6", size=3)
    }

    main_plot <- main_plot +
        theme_minimal() +
        #ggplot2::theme_bw() +
        ggplot2::theme(
            panel.spacing=ggplot2::unit(0, "cm"),
            panel.grid=ggplot2::element_blank(),
            panel.grid.minor=ggplot2::element_blank(),
            legend.position="bottom",
            strip.placement.y="outside",
            strip.background.y=ggplot2::element_blank(),
            strip.text.y.left=ggplot2::element_text(
                angle=0,
                hjust=1,
                vjust=1,
                face="bold",
                size=ggplot2::rel(1.1)),
            axis.line.x=ggplot2::element_line(),
            axis.text.x=ggplot2::element_text(margin=ggplot2::margin(t=4, b=4)),
            axis.title.x=ggplot2::element_text(),
            axis.ticks.x=ggplot2::element_line(size=0.3),
            axis.ticks.length.x=ggplot2::unit(.1, "cm"),
            axis.ticks.length.x.top=ggplot2::unit(.1, "cm"),
            axis.ticks.length.x.bottom=ggplot2::unit(.1, "cm"))


    main_plot
}

#' @rdname forest_plot
#' @export
forest_plot_table <- function(fodat, scale=c("relative", "absolute")) {

    scale <- match.arg(scale)

    est <- md <- lo <- hi <- covar <- covval <- NULL

    format_number <- function(x, digits=2) {
        is.inf <- function(x) { x == Inf }
        x2 <- table1::round_pad(x, digits=digits)
        x2[is.na(x)] <- "-"
        x2[is.inf(x)] <- "\U{221E}"
        x2
    }
    format_est_ci <- function(x, xmin, xmax) {
        sprintf("%s [%s, %s]",
            format_number(x), format_number(xmin), format_number(xmax))
    }

    yax <- unique(fodat[, .(breaks=mappings::cf(covar, covval), labels=covval)])

    table_plot <-
        ggplot2::ggplot(fodat, ggplot2::aes(x=1, y=mappings::cf(covar, covval))) +
        ggplot2::labs(x=NULL, y=NULL, title="Estimate [95% CI]") +
        ggplot2::facet_grid(covar ~ ., scales="free", space="free", switch="y") +
        ggplot2::scale_y_discrete(expand=ggplot2::expansion(add=1), breaks=yax$breaks, labels=yax$labels, limits=rev)

    if (scale == "relative") {
        table_plot <- table_plot +
            ggplot2::geom_text(ggplot2::aes(label=format_est_ci(est/refval, lo/refval, hi/refval)))
    } else {
        table_plot <- table_plot +
            ggplot2::geom_text(ggplot2::aes(label=format_est_ci(est, lo, hi)))
    }

    table_plot <- table_plot +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.title=ggplot2::element_text(size=ggplot2::rel(1), hjust=0.5),
            strip.text=ggplot2::element_blank(),
            panel.spacing=ggplot2::unit(0, "cm"))
}

#--------------------------------------------------------------------------------

# (c) 2022 Benjamin Rich

