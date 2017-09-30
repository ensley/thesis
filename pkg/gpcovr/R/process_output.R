#' Plot the trace of spline coefficient MCMC
#'
#' @param betas A matrix where each row represents a set of coefficients
#' @param n The number of trace plots to draw. beta_1 through beta_n will be
#'   displayed.
#'
#' @return Draws a plot. Returns \code{betas} invisibly.
#' @export
traceplot <- function(betas, n = 4)
{
  if (n > ncol(beta)) stop('n cannot be greater than the number of coefficients')
  bmcmc <- coda::mcmc(betas[ ,seq_len(n)])
  graphics::plot(bmcmc)
  invisible(betas)
}


#' Plot the log-log spectral density after estimation
#'
#' @param beta A vector of spline coefficients (typically the final result from
#'   \code{\link{estbeta}})
#' @param gpModel A \code{GPmodel} object
#' @param knots The knot locations
#' @param show_knots To hide the vertical dashed lines indicating where the
#'   knots are located, set this to FALSE. Defaults to TRUE
#' @param beta_true The "true" spline coefficients
#' @param beta_init The initial spline coefficient values
#' @param err_bounds A 2-column matrix containing lower and upper error bounds
#'
#' @details All of \code{beta_true}, \code{beta_init}, and \code{err_bounds} are \code{NULL} by default. If any are not known or not desired on the final plot, leave them as \code{NULL}. Otherwise, explicitly provide them.
#'
#' @return The \code{GPmodel} object, invisibly. Creates the plot as a side effect.
#' @export
loglogdensity_plot <- function(beta, gpModel, knots, show_knots = TRUE, beta_true = NULL, beta_init = NULL, err_bounds = NULL)
{
  seq_length <- if (!is.null(err_bounds)) nrow(err_bounds) else 500
  if (gpModel$family == 'matern') {
    x <- seq(-5, 5, length = seq_length)
  } else if (gpModel$family == 'dampedcos') {
    x <- seq(-5, 3, length = seq_length)
  }
  y <- gpModel$specdens(x)
  logx <- log(x)
  logy <- logx + log(y)
  w <- round(seq(1, length(x), length = 50))
  xx <- logx[w]
  yy <- logy[w]

  spl <- nsbasis(logx, knots)

  pts <- yy - log(normalize_dlogspline(spl, beta_true, knots))
  l1 <- dloglogspline(logx, beta, knots)
  l2 <- if (!is.null(beta_true)) dloglogspline(logx, beta_true, knots) else NULL
  l3 <- if (!is.null(beta_init)) dloglogspline(logx, beta_init, knots) else NULL

  if (!is.null(err_bounds)) {
    min_lower_bd <- min(err_bounds[ ,1])
    min_upper_bd <- min(err_bounds[ ,2])
    max_lower_bd <- max(err_bounds[ ,1])
    max_upper_bd <- max(err_bounds[ ,2])
  } else {
    min_lower_bd <- min_upper_bd <- max_lower_bd <- max_upper_bd <- NULL
  }

  plot_ymin <- min(min(pts), min_lower_bd, min_upper_bd, min(l1), min(l2), min(l3))
  plot_ymax <- max(max(pts), max_lower_bd, max_upper_bd, max(l1), max(l2), max(l3))

  graphics::plot(range(logx), c(plot_ymin, plot_ymax), type = 'n',
       xlab = expression(paste(log, ' ', omega)),
       ylab = expression(paste(log, ' f(', log, ' ', omega, ')')),
       main = 'log-log density')
  if (!is.null(err_bounds)) {
    graphics::polygon(c(logx, rev(logx)),
            c(err_bounds[ ,1], rev(err_bounds[ ,2])),
            col = 'grey90',
            border = NA)
  }
  graphics::points(xx, pts, pch = 20, cex = 1.5, col = 'grey50')
  graphics::lines(logx, l1, lwd = 2)
  if (!is.null(beta_true)) graphics::lines(logx, l2, col = 'blue', lwd = 2)
  if (!is.null(beta_init)) graphics::lines(logx, l3, lty = 2, lwd = 2, col = 'orange')
  if (show_knots) graphics::abline(v = knots, col = 'grey', lty = 3)
  graphics::legend('topleft',
         c('true points', 'closest fit', 'estimated', 'initial'),
         col = c('grey50', 'black', 'blue', 'orange'),
         lty = c(NA, 1, 1, 2),
         lwd = c(NA, 2, 2, 2),
         pch = c(16, NA, NA, NA))

  invisible(gpModel)
}
