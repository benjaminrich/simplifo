#' Styling for forest plots
#'
#' This function provides the default styling for [forest_plot()], which can be
#' overridden.
#'
#' @param point.color Color for points
#' @param pointrange.color Color for point ranges
#' @param refline.color Color for reference line
#' @param refline.linetype Line type for reference line
#' @param refrange.color Color for border of reference range
#' @param refrange.fill Fill color for reference range
#' @return An object of class `forest_plot_styling`.
#' @seealso [forest_plot()].
#' @examples
#' forest_plot_styling()
#' @export
forest_plot_styling <- function(
    point.color      = "#3333cc", #"#3A70B6",
    pointrange.color = "#3333cc", #"#3A70B6"
    refline.color    = "gray40",
    refline.linetype = "dotted",
    refrange.color   = NA,
    refrange.fill    = "gray90"
) {

    l <- mget(names(formals()), sys.frame(sys.nframe()))

    structure(l, class="forest_plot_styling")
}

#' @export
print.forest_plot_styling <- function(x, ...) {
    cat0 <- function(...) cat(..., sep="")

    cat0("Forest Plot Styling:\n",
         "────────────────────\n")


    for (i in 1:length(x)) {
        cat0(names(x)[i], ": ", x[[i]], "\n")
    }
}

#forest_plot_styling()
