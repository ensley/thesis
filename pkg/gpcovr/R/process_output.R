#' Plot the trace of spline coefficient MCMC
#'
#' @param betas A matrix where each row represents a set of coefficients
#' @param n The number of trace plots to draw. beta_1 through beta_n will be
#'   displayed.
#'
#' @return Draws a plot. Returns \code{betas} invisibly.
#' @export
trace_plot <- function(betas, n = 4)
{
  if (n > ncol(betas)) stop('n cannot be greater than the number of coefficients')
  bmcmc <- coda::mcmc(betas[ ,seq_len(n)])
  graphics::plot(bmcmc)
  invisible(betas)
}


#' Plot the spectral density after estimation (log or log-log)
#'
#' @param beta A vector of spline coefficients (typically the final result from
#'   \code{\link{estbeta_aao}} or \code{\link{estbeta_ooat}})
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
    pts <- if (!is.null(beta_true)) yy - log(normalize_dlogspline(spl, beta_true, knots)) else NULL
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
  graphics::lines(logx, l1, lwd = 2, col = 'blue')
  if (!is.null(beta_true)) graphics::lines(logx, l2, lwd = 2, lty = 2)
  if (!is.null(beta_init)) graphics::lines(logx, l3, lty = 2, lwd = 2, col = 'orange')
  if (show_knots) graphics::abline(v = knots, col = 'grey', lty = 3)
  graphics::legend('topleft',
         c('true points', 'closest fit', 'estimated', 'initial'),
         col = c('grey50', 'black', 'blue', 'orange'),
         lty = c(NA, 2, 1, 2),
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
  graphics::lines(logx, l1, lwd = 2, col = 'blue')
  if (!is.null(beta_true)) graphics::lines(logx, l2, lty = 2, lwd = 2)
  if (!is.null(beta_init)) graphics::lines(logx, l3, lty = 2, lwd = 2, col = 'orange')
  if (show_knots) graphics::abline(v = knots, col = 'grey', lty = 3)
  graphics::legend('topleft',
                   c('true points', 'closest fit', 'estimated', 'initial'),
                   col = c('grey50', 'black', 'blue', 'orange'),
                   lty = c(NA, 2, 1, 2),
                   lwd = c(NA, 2, 2, 2),
                   pch = c(16, NA, NA, NA))

  invisible(gpModel)
}



#' Plot the covariance function after estimation
#'
#' @param h The lag values to plot the covariance function over
#' @param beta A vector of spline coefficients (typically the final result from
#'   \code{\link{estbeta_aao}} or \code{\link{estbeta_ooat}})
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
    min_lower_bd <- min(err_bounds[ ,1])
    min_upper_bd <- min(err_bounds[ ,2])
    max_lower_bd <- max(err_bounds[ ,1])
    max_upper_bd <- max(err_bounds[ ,2])
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
  graphics::lines(h, l1, lwd = 2, col = 'blue')
  if (!is.null(beta_true)) graphics::lines(h, l2, lwd = 2, lty = 2)
  if (!is.null(gpModel)) graphics::lines(h, l3, col = 'orange', lwd = 2, lty = 2)
  graphics::abline(h = 0, lty = 3)
  graphics::legend('topright', c('estimated', 'optimal', 'true'), col = c('blue', 'black', 'orange'), lwd = 2, lty = c(1, 2, 2))

  invisible(gpModel)
}


#' Predict points from a GPsimulated object
#'
#' @param object A \code{GPsimulated} object
#' @param ... Not used
#' @param beta The spline coefficients
#' @param knots The knot locations
#'
#' @return A vector of predictions
#' @export
predict.GPsimulated <- function(object, ..., beta, knots) {
  grid <- object$locs[ ,1:2]
  is_pred <- which(object$locs$type == 'pred')
  pts_obs <- object$Y[-is_pred]
  dist_pred2obs <- as.matrix(stats::dist(grid))[is_pred,-is_pred]

  draws <- exp(rlogspline(50000, beta, knots))
  gamma1 <- matrix(gpumc::mc(dist_pred2obs, draws),
                   nrow = length(is_pred))
  G <- matrix(gpumc::mc(object$dist_obs, draws),
              nrow = length(pts_obs))

  while(methods::is(try(solve(G), silent = TRUE), 'try-error')) {
    cat('matrix inversion failed. drawing different spectral density samples.\n')
    draws <- exp(rlogspline(50000, beta, knots))
    gamma1 <- matrix(gpumc::mc(dist_pred2obs, draws),
                     nrow = length(is_pred))
    G <- matrix(gpumc::mc(object$dist_obs, draws),
                nrow = length(pts_obs))
  }

  drop(gamma1 %*% solve(G) %*% pts_obs)
}



#' Optimal prediction from GPsimulated object using known covariance function
#'
#' Returns predictions assuming the covariance model used to simulate the data
#' was known
#'
#' @param object A \code{GPsimulated} object
#'
#' @return A vector of predictions
#' @export
optimal_predict <- function(object) {
  grid <- object$locs[ ,1:2]
  is_pred <- which(object$locs$type == 'pred')
  pts_obs <- object$Y[-is_pred]
  dist_pred2obs <- as.matrix(stats::dist(grid))[is_pred,-is_pred]
  gamma1 <- matrix(RandomFields::RFcov(object$m$model, as.numeric(dist_pred2obs)),
                   nrow = length(is_pred))
  G <- matrix(RandomFields::RFcov(object$m$model, as.numeric(object$dist_obs)),
              nrow = length(pts_obs))

  drop(gamma1 %*% solve(G) %*% pts_obs)
}



#' Optimally predict points using an arbitrary covariance model
#'
#' To predict points from another model, such as the best-fitting Matern model,
#' create it using \code{RandomFields} and specify it here.
#'
#' @param object A \code{GPsimulated} object
#' @param model A covariance model made using the \code{RandomFields} package
#'
#' @return A vector of predictions
#' @export
othermodel_predict <- function(object, model) {
  object$m$model <- model
  optimal_predict(object)
}



#' Mean squared error
#'
#' @param pred Predicted values
#' @param actual Actual values
#'
#' @return The mean squared error: \code{mean((pred - actual)^2)}
#' @export
mse <- function(pred, actual) {
  mean((pred - actual)^2)
}



#' Predict points from a GPsimulated object and find the mean squared error
#'
#' If the predictions have already been found via
#' \code{\link{predict.GPsimulated}}, use \code{\link{mse}} instead.
#'
#' @param gpEst A \code{GPsimulated} object
#' @param beta The spline coefficients
#' @param knots The knot locations
#'
#' @return The mean squared error
#' @export
gp_mse <- function(gpEst, beta, knots) {
  is_pred <- which(gpEst$locs$type == 'pred')
  pts_pred <- gpEst$Y[is_pred]
  preds <- predict.GPsimulated(gpEst, beta = beta, knots = knots)
  mean((preds - pts_pred)^2)
}



#' Optimally predict points from a GPsimulated object and find the mean squared
#' error
#'
#' If the predictions have already been found via \code{\link{optimal_predict}},
#' use \code{\link{mse}} instead.
#'
#' @param gpEst A \code{GPsimulated} object
#'
#' @return The mean squared error
#' @export
gp_mse_optimal <- function(gpEst) {
  is_pred <- which(gpEst$locs$type == 'pred')
  pts_pred <- gpEst$Y[is_pred]
  preds <- optimal_predict(gpEst)
  mean((preds - pts_pred)^2)
}


#' Optimally predict points using an arbitrary covariance model and find the
#' mean squared error
#'
#' To predict points from another model, such as the best-fitting Matern model,
#' create it using \code{RandomFields} and specify it here.
#'
#' If the predictions have already been found via
#' \code{\link{predict.GPsimulated}}, use \code{\link{mse}} instead.
#'
#' @param gpEst A \code{GPsimulated} object
#' @param model A covariance model made using the \code{RandomFields} package
#'
#' @return The mean squared error
#' @export
gp_mse_othermodel <- function(gpEst, model) {
  gpEst$m$model <- model
  gp_mse_optimal(gpEst)
}



#' Fit a Matern model to a GPsimulated object
#'
#' Uses \code{\link[spBayes]{spLM}} to fit a Matern model to the data.
#'
#' @param gpEst A \code{GPsimulated} object
#' @param n_samps Number of MCMC samples to take
#' @param starting Starting parameter values
#' @param tuning Tuning parameter values
#' @param priors Priors and hyperparameters
#'
#' @details See \code{\link[spBayes]{spLM}} for details on the arguments. Some sensible defaults are provided but may need to be tweaked.
#'
#' @return A list containing the best-fitting Matern model, the result of \code{\link[spBayes]{spLM}}, and the parameters (transformed back into our parameterization)
#' @export
fit_matern <- function(gpEst, n_samps, starting = list('sigma.sq' = 1, 'tau.sq' = 0.1, 'phi' = 0.5, 'nu' = 1), tuning = list('sigma.sq' = 0.01, 'tau.sq' = 0.01, 'phi' = 0.005, 'nu' = 0.005), priors = list('sigma.sq.ig' = c(2, 1), 'tau.sq.ig' = c(2, 1), 'phi.unif' = c(1e-5, 20), 'nu.unif' = c(1e-5, 5))) {
  is_obs <- which(gpEst$locs$type == 'obs')
  coords <- as.matrix(gpEst$locs[is_obs,1:2])
  y <- gpEst$Y[is_obs]
  burn.in <- floor(0.75*n_samps)

  m_i <- spBayes::spLM(y ~ 1, coords = coords, starting = starting, tuning = tuning, priors = priors, cov.model = 'matern', n.samples = n_samps, n.report = 200)

  pars <- apply(m_i$p.theta.samples[burn.in:nrow(m_i$p.theta.samples), ], 2, stats::median)
  model <- RandomFields::RMwhittle(nu = pars['nu'],
                                   notinvnu = TRUE,
                                   var = pars['sigma.sq'],
                                   scale = 1/pars['phi']) +
    RandomFields::RMnugget(var = pars['tau.sq'])

  params <- c(pars['nu'],
              (2 * sqrt(pars['nu'])) / pars['phi'],
              pars['sigma.sq'],
              pars['tau.sq'])
  names(params) <- c('nu', 'rho', 'sigma', 'nugget')

  list(model = model, fitobj = m_i, params = params)
}
