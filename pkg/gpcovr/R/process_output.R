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


#' Plot the spectral density after estimation (log or log-log)
#'
#' @param beta A vector of spline coefficients (typically the final result from
#'   \code{\link{estbeta}})
#' @param knots The knot locations
#' @param show_knots To hide the vertical dashed lines indicating where the
#'   knots are located, set this to FALSE. Defaults to TRUE
#' @param show_points To show points along the true log-log spectral density,
#'   set this to TRUE. Defaults to FALSE.
#' @param gpModel A \code{GPmodel} object
#' @param beta_true The "true" spline coefficients
#' @param beta_init The initial spline coefficient values
#' @param err_bounds A 2-column matrix containing lower and upper error bounds
#'
#' @details All of \code{gpModel}, \code{beta_true}, \code{beta_init}, and \code{err_bounds} are
#'   \code{NULL} by default. If any are not known or not desired on the final
#'   plot, leave them as \code{NULL}. Otherwise, explicitly provide them.
#'
#'   If \code{show_points} is set to TRUE but \code{gpModel} is NULL,
#'   \code{show_points} will get forced to FALSE.
#'
#' @return The \code{GPmodel} object, invisibly. Creates the plot as a side
#'   effect.
#' @export
loglogdensity_plot <- function(beta, knots, show_knots = TRUE, show_points = FALSE, gpModel = NULL, beta_true = NULL, beta_init = NULL, err_bounds = NULL)
{
  seq_length <- if (!is.null(err_bounds)) nrow(err_bounds) else 500
  if (!is.null(gpModel)) {
    x <- exp(seq(gpModel$logx[1], gpModel$logx[2], length = seq_length))
    y <- gpModel$specdens(x)
    logx <- log(x)
    logy <- logx + log(y)
    w <- round(seq(1, length(x), length = 50))
    xx <- logx[w]
    yy <- logy[w]

    spl <- nsbasis(logx, knots)
    pts <- yy - log(normalize_dlogspline(spl, beta_true, knots))
  } else {
    show_points <- FALSE
    logx <- seq(-5, 5, length = seq_length)
    pts <- NULL
  }

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

  plot_ymin <- min(suppressWarnings(min(pts)), min_lower_bd, min_upper_bd, min(l1), suppressWarnings(min(l2)), suppressWarnings(min(l3)))
  plot_ymax <- max(suppressWarnings(max(pts)), max_lower_bd, max_upper_bd, max(l1), suppressWarnings(max(l2)), suppressWarnings(max(l3)))

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
  if (show_points) graphics::points(xx, pts, pch = 20, cex = 1.5, col = 'grey50')
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



#' @rdname loglogdensity_plot
#' @export
logdensity_plot <- function(beta, knots, show_knots = TRUE, show_points = FALSE, gpModel = NULL, beta_true = NULL, beta_init = NULL, err_bounds = NULL)
{
  seq_length <- if (!is.null(err_bounds)) nrow(err_bounds) else 500
  if (!is.null(gpModel)) {
    x <- exp(seq(gpModel$logx[1], gpModel$logx[2], length = seq_length))
    y <- gpModel$specdens(x)
    logx <- log(x)
    logy <- logx + log(y)
    w <- round(seq(1, length(x), length = 50))
    xx <- logx[w]
    yy <- logy[w]

    spl <- nsbasis(logx, knots)
    pts <- exp(yy) / normalize_dlogspline(spl, beta_true, knots)
  } else {
    show_points <- FALSE
    logx <- seq(-5, 5, length = seq_length)
    pts <- NULL
  }

  l1 <- dlogspline(logx, beta, knots)
  l2 <- if (!is.null(beta_true)) dlogspline(logx, beta_true, knots) else NULL
  l3 <- if (!is.null(beta_init)) dlogspline(logx, beta_init, knots) else NULL

  if (!is.null(err_bounds)) {
    min_lower_bd <- min(exp(err_bounds[ ,1]))
    min_upper_bd <- min(exp(err_bounds[ ,2]))
    max_lower_bd <- max(exp(err_bounds[ ,1]))
    max_upper_bd <- max(exp(err_bounds[ ,2]))
  } else {
    min_lower_bd <- min_upper_bd <- max_lower_bd <- max_upper_bd <- NULL
  }

  plot_ymin <- min(suppressWarnings(min(pts)), min_lower_bd, min_upper_bd, min(l1), suppressWarnings(min(l2)), suppressWarnings(min(l3)))
  plot_ymax <- max(suppressWarnings(max(pts)), max_lower_bd, max_upper_bd, max(l1), suppressWarnings(max(l2)), suppressWarnings(max(l3)))

  graphics::plot(range(logx), c(plot_ymin, plot_ymax), type = 'n',
                 xlab = expression(paste(log, ' ', omega)),
                 ylab = expression(paste('f(', log, ' ', omega, ')')),
                 main = 'log density')
  if (!is.null(err_bounds)) {
    graphics::polygon(c(logx, rev(logx)),
                      exp(c(err_bounds[ ,1], rev(err_bounds[ ,2]))),
                      col = 'grey90',
                      border = NA)
  }
  if (show_points) graphics::points(xx, pts, pch = 20, cex = 1.5, col = 'grey50')
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



#' Plot the covariance function after estimation
#'
#' @param h The lag values to plot the covariance function over
#' @param beta A vector of spline coefficients (typically the final result from
#'   \code{\link{estbeta}})
#' @param knots The knot locations
#' @param gpModel A \code{GPmodel} object
#' @param beta_true The "true" spline coefficients
#' @param beta_init The initial spline coefficient values
#' @param err_bounds A 2-column matrix containing lower and upper error bounds
#'
#' @details All of \code{gpModel}, \code{beta_true}, \code{beta_init}, and
#'   \code{err_bounds} are \code{NULL} by default. If any are not known or not
#'   desired on the final plot, leave them as \code{NULL}. Otherwise, explicitly
#'   provide them.
#'
#' @return The \code{GPmodel} object, invisibly. Creates the plot as a side
#'   effect.
#' @export
covar_plot <- function(h, beta, knots, gpModel = NULL, beta_true = NULL, beta_init = NULL, err_bounds = NULL)
{
  l1 <- gpumc::mc(h, exp(rlogspline(50000, beta, knots)))
  l2 <- if (!is.null(beta_true)) gpumc::mc(h, exp(rlogspline(50000, beta_true, knots))) else NULL
  l3 <- if (!is.null(gpModel)) RandomFields::RFcov(gpModel$model, h)

  if (!is.null(err_bounds)) {
    min_lower_bd <- min(exp(err_bounds[ ,1]))
    min_upper_bd <- min(exp(err_bounds[ ,2]))
    max_lower_bd <- max(exp(err_bounds[ ,1]))
    max_upper_bd <- max(exp(err_bounds[ ,2]))
  } else {
    min_lower_bd <- min_upper_bd <- max_lower_bd <- max_upper_bd <- NULL
  }

  plot_ymin <- min(min_lower_bd, min_upper_bd, min(l1), suppressWarnings(min(l2)), suppressWarnings(min(l3)))
  plot_ymax <- max(max_lower_bd, max_upper_bd, max(l1), suppressWarnings(max(l2)), suppressWarnings(max(l3)))

  graphics::plot(range(h), c(plot_ymin, plot_ymax), type = 'n',
       xlab = 'h',
       ylab = 'C(h)',
       main = 'Covariance function')

  if (!is.null(err_bounds)) graphics::polygon(c(h, rev(h)), c(err_bounds[ ,1], rev(err_bounds[ ,2])), col = 'grey90', border = NA)
  graphics::lines(h, l1, lwd = 2)
  if (!is.null(beta_true)) graphics::lines(h, l2, col = 'blue', lwd = 2)
  if (!is.null(gpModel)) graphics::lines(h, l3, lwd = 2)
  graphics::abline(h = 0, lty = 3)
  graphics::legend('topright', c('true', 'estimated'), col = c('black', 'blue'), lwd = 2)

  invisible(gpModel)
}
