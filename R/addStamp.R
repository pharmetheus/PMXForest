#' Add a source and time stamp to a plot
#'
#' @description Adds a caption recording where and when a figure was produced,
#'   so a plot pasted into a report can be traced back to the code that made it.
#'
#' @details
#'   Inside a `knitr` / R Markdown chunk the caption reads
#'   `<working directory name>/<input file>/<chunk label>` followed by the
#'   creation time. Outside one, the file and chunk are unknown and simply
#'   omitted, so the caption is the working directory name and the time.
#'
#'   Only the name of the working directory is used, not its full path -
#'   `PMXForest-private`, not `/home/you/github/PMXForest-private`. That
#'   matches `PhRame::add_stamp()` and keeps local paths out of a figure that
#'   may be pasted into a report, at the cost of being ambiguous between two
#'   checkouts whose leaf directory has the same name. Pass `source` to say
#'   something more specific.
#'
#'   This is a native reimplementation of the stamp Pharmetheus' internal
#'   `PhRame::add_stamp()` applies, kept dependency-free (`ggplot2` only) so it
#'   can be used from the public packages. Only the "return the annotated plot"
#'   behaviour is reproduced - saving and printing are left to the caller, which
#'   is all any of the plotting functions here need.
#'
#' @param p A `ggplot` object.
#' @param size Text size of the caption. Defaults to `6`, matching
#'   `PhRame::add_stamp()`.
#' @param source Optional character string used instead of the derived
#'   directory / file / chunk path, for a caller that knows better.
#'
#' @return The `ggplot` object with the caption added.
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point()
#' addStamp(p)
#'
#' @seealso
#'   `vignette("Part2-walkthrough", package = "PMXForest")` for how this fits the whole workflow.
#'
#' @export
addStamp <- function(p, size = 6, source = NULL) {
  if (!inherits(p, "ggplot")) {
    stop("`p` must be a ggplot object.", call. = FALSE)
  }
  p +
    ggplot2::labs(caption = stampText(source)) +
    ggplot2::theme(plot.caption = ggplot2::element_text(size = size))
}


#' The stamp caption text
#'
#' Split out so the plotting functions can be tested against the text without
#' building a plot, and so the knitr lookup stays in one place.
#'
#' @param source Optional character string replacing the derived path.
#' @return A single string: a path-like source, a newline, and the time.
#' @keywords internal
#' @noRd
stampText <- function(source = NULL) {
  if (is.null(source)) {
    ## knitr is a Suggests: outside a knit session, or without it installed,
    ## the input file and chunk label are simply unknown. Drop the empty
    ## segments rather than emitting "dir//" as PhRame::add_stamp() does.
    inputFile <- NULL
    chunk <- NULL
    if (requireNamespace("knitr", quietly = TRUE)) {
      inputFile <- tryCatch(knitr::current_input(), error = function(e) NULL)
      chunk <- tryCatch(knitr::opts_current$get("label"), error = function(e) NULL)
    }
    parts <- c(basename(getwd()), inputFile, chunk)
    source <- paste(parts[nzchar(parts) & !is.na(parts)], collapse = "/")
  }
  paste0(source, "\nCreated:", Sys.time())
}
