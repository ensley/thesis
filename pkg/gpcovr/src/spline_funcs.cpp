// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <algorithm>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_spline.h>
using namespace Rcpp;

// --------- these are all helper functions and structs relating to the
// --------- integration done in normalize_dlogspline().


struct integrand_params { RcppGSL::Vector knots; RcppGSL::Vector coefs; };


double d_one(int k, double x, RcppGSL::Vector knots) {
  int K = knots.size();
  double knot_k = knots[k];
  double knot_K = knots[K-1];
  double xk = x - knot_k;
  double xK = x - knot_K;
  return (std::max(xk*xk*xk, 0.0) - std::max(xK*xK*xK, 0.0)) / (knot_K - knot_k);
}


double H_one(int i, double x, RcppGSL::Vector knots) {
  int K = knots.size();
  if(i == 0) {
    return 1.0;
  } else if(i == 1) {
    return x;
  } else {
    return d_one(i-2, x, knots) - d_one(K-2, x, knots);
  }
}


RcppGSL::Vector nsbasis_one(double x, RcppGSL::Vector knots) {
  int nk = knots.size();
  RcppGSL::Vector out(nk);
  for(int i = 0; i < nk; ++i) {
    out[i] = H_one(i, x, knots);
  }
  return out;
}


double integrand(double x, void *p) {
  struct integrand_params *params = (struct integrand_params *) p;
  RcppGSL::Vector knots = (params->knots);
  RcppGSL::Vector coefs = (params->coefs);
  double y;
  RcppGSL::Vector basis = nsbasis_one(x, knots);
  int n = basis.size();
  gsl_blas_ddot(basis, coefs, &y);
  return exp(y);
}


// ----------------------------------------------------------------------------


// ----------- MAIN FUNCTIONS



RcppGSL::Vector d(int k, RcppGSL::Vector X, RcppGSL::Vector knots) {
  int K = knots.size(), xn = X.size();
  RcppGSL::Vector out(xn); // result vector
  double knot_k = knots[k];
  double knot_K = knots[K-1];

  for(int i = 0; i < xn; ++i) {
    double xi = X[i];
    double xk = xi - knot_k;
    double xK = xi - knot_K;
    out[i] = (std::max(xk*xk*xk, 0.0) - std::max(xK*xK*xK, 0.0)) / (knot_K - knot_k);
  }
  return out;
}



RcppGSL::Vector H(int i, RcppGSL::Vector X, RcppGSL::Vector knots) {
  int K = knots.size(), n = X.size();

  if(i == 0) {
    RcppGSL::Vector y(n);
    gsl_vector_set_all(y, 1.0);
    return y;
  } else if(i == 1) {
    return X;
  } else {
    RcppGSL::Vector y = d(i-2, X, knots);
    RcppGSL::Vector d2 = d(K-2, X, knots);
    gsl_vector_sub(y, d2);
    return y;
  }
}


//' Create a natural cubic spline basis
//'
//' @param X the locations at which the spline functions will be evaluated. The
//' matrix will have \code{length(X)} rows.
//' @param knots the knot locations. The first and last knots are the boundary
//' knots, and the rest are internal knots. The matrix will have
//' \code{length(knots)} knots.
//' @return the basis matrix
//' @details The basis used by this function is described in \emph{The Elements
//' of Statistical Learning} by Hastie and Tibshirani, in the appendix of
//' Chapter 5 (2nd edition).
//' @export
// [[Rcpp::export]]
RcppGSL::Matrix nsbasis(RcppGSL::Vector X, RcppGSL::Vector knots) {
  int nx = X.size();
  int nk = knots.size();
  RcppGSL::Matrix out(nx, nk);
  for(int i = 0; i < nk; ++i) {
    RcppGSL::Vector v = H(i, X, knots);
    for(int j = 0; j < nx; ++j) {
      out(j,i) = v[j];
    }
  }
  return out;
}


//' Predict values from spline basis
//'
//' Given a natural spline basis and a set of coefficients, returns the
//' values along the fitted curve.
//'
//' @param basis An \code{n} by \code{k} spline basis matrix, where \code{n}
//' is the number of points the curve will be evaluated at and \code{k} is the
//' number of knots. Should be the output of \code{\link{nsbasis}}.
//' @param coefs A vector of length \code{k} representing the basis coefficients.
//' @return A vector of length \code{n} containing the fitted values.
//' @seealso \code{\link{dloglogspline}} for the same thing but shifted so that
//' it integrates to 1 when exponentiated.
//' @export
// [[Rcpp::export]]
RcppGSL::Vector predict_natspl(const RcppGSL::Matrix &basis, const RcppGSL::Vector &coefs) {
  int n = basis.nrow(), k = basis.ncol();
  if(coefs.size() != k) Rcpp::stop("dimensions do not match");
  RcppGSL::Vector y(n);
  gsl_blas_dgemv(CblasNoTrans, 1.0, basis, coefs, 1.0, y);
  return y;
}


//' Evaluate the un-normalized log spectral density
//'
//' @param basis An \code{n} by \code{k} spline basis matrix, where \code{n}
//' is the number of points the curve will be evaluated at and \code{k} is the
//' number of knots. Should be the output of \code{\link{nsbasis}}.
//' @param coefs A vector of length \code{k} representing the basis coefficients.
//' @param knots A vector of length \code{k} representing the knot locations.
//' @return A vector of length \code{n} containing the log density.
//' @seealso \code{\link{predict_natspl}} for the log of this function.
//' @export
// [[Rcpp::export]]
RcppGSL::Vector dlogspline_unnorm(RcppGSL::Matrix basis, RcppGSL::Vector coefs, RcppGSL::Vector knots) {
  int n = basis.nrow(), k = basis.ncol();
  if(coefs.size() != k) Rcpp::stop("dimensions do not match");
  RcppGSL::Vector y(n);
  gsl_blas_dgemv(CblasNoTrans, 1.0, basis, coefs, 1.0, y);
  for(int i=0; i < n; ++i) {
    y[i] = exp(y[i]);
  }
  return y;
}


//' Calculate the normalizing constant for the log spectral density
//'
//' The log density given by \code{\link{dlogspline_unnorm}} is not a valid
//' density because it does not generally integrate to 1. This function finds
//' the constant C such that 1/C times \code{dlogspine_unnorm} will integrate
//' to 1.
//'
//' @param basis An \code{n} by \code{k} spline basis matrix, where \code{n}
//' is the number of points the curve will be evaluated at and \code{k} is the
//' number of knots. Should be the output of \code{\link{nsbasis}}.
//' @param coefs A vector of length \code{k} representing the basis coefficients.
//' @param knots A vector of length \code{k} representing the knot locations.
//' @return the normalizing constant C
//' @seealso \code{\link{dlogspline_unnorm}}, \code{\link{dlogspline}}
//' @export
// [[Rcpp::export]]
double normalize_dlogspline(RcppGSL::Matrix basis, RcppGSL::Vector coefs, RcppGSL::Vector knots) {
  gsl_set_error_handler_off();
  gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);
  gsl_function F;
  struct integrand_params params = { knots, coefs };
  F.function = &integrand;
  F.params = &params;
  double result, error;

  int intresult = gsl_integration_qagi(&F, 0, 1e-7, 1000, workspace, &result, &error);
  if(ISNAN(intresult)) {
    Rcpp::stop("Integration error");
  }

  return result;
}


//' Evaluate the log spectral density
//'
//' @param x The values at which to evaluate the log density.
//' @param coefs A vector of length \code{k} representing the basis coefficients.
//' @param knots A vector of length \code{k} representing the knot locations.
//' @return A vector of log density values.
//' @export
// [[Rcpp::export]]
RcppGSL::Vector dlogspline(RcppGSL::Vector x, RcppGSL::Vector coefs, RcppGSL::Vector knots) {
  RcppGSL::Matrix basis = nsbasis(x, knots);
  RcppGSL::Vector y = dlogspline_unnorm(basis, coefs, knots);
  double scale = 1.0 / normalize_dlogspline(basis, coefs, knots);
  for(int i=0; i < y.size(); ++i) {
    y[i] = scale * y[i];
  }
  return y;
}


//' Evaluate the log of the log spectral density
//'
//' @param x The values at which to evaluate the log of the log density.
//' @param coefs A vector of length \code{k} representing the basis coefficients.
//' @param knots A vector of length \code{k} representing the knot locations.
//' @return A vector of log log density values.
//' @export
// [[Rcpp::export]]
RcppGSL::Vector dloglogspline(RcppGSL::Vector x, RcppGSL::Vector coefs, RcppGSL::Vector knots) {
  RcppGSL::Vector y = dlogspline(x, coefs, knots);
  int n = y.size();
  for(int i=0; i < n; ++i) {
    y[i] = log(y[i]);
  }
  return y;
}


RcppGSL::Vector get_slopes(RcppGSL::Vector coefs, RcppGSL::Vector knots) {
  double kmin, kmax;
  gsl_vector_minmax(knots, &kmin, &kmax);

  RcppGSL::Vector x(4);
  x[0] = kmin - 1;
  x[1] = kmin;
  x[2] = kmax;
  x[3] = kmax + 1;

  RcppGSL::Matrix basis = nsbasis(x, knots);
  RcppGSL::Vector y = predict_natspl(basis, coefs);
  RcppGSL::Vector out(2);
  out[0] = y[1] - y[0];
  out[1] = y[3] - y[2];
  return out;
}

//' Evaluate points beyond the knot boundaries
//'
//' Since we are using a natural spline basis, the fitted curve beyond the
//' boundary knots is guaranteed to be linear. To quickly find the slopes
//' of these linear tails, this function evaluates the fitted curve at four
//' points: the boundary knots, as well as at a point slightly
//' beyond each boundary knot.
//'
//' @param coefs A vector of length \code{k} representing the basis coefficients.
//' @param knots A vector of length \code{k} representing the knot locations.
//' @return A vector of length 4 containing the following values:
//' \enumerate{
//'   \item \eqn{f(kmin - 1)}
//'   \item \eqn{f(kmin)}
//'   \item \eqn{f(kmax)}
//'   \item \eqn{f(kmax + 1)}
//' }
//' where \code{kmin} and \code{kmax} are the minimum and maximum knots and
//' \code{f} is the log log density. Then the left tail slope is
//' \eqn{f(kmin) - f(kmin - 1)} and the right tail slope is
//' \eqn{f(kmax+ 1) - f(kmax)}.
//' @export
// [[Rcpp::export]]
RcppGSL::Vector get_slope_pts(RcppGSL::Vector coefs, RcppGSL::Vector knots) {
  double kmin, kmax;
  gsl_vector_minmax(knots, &kmin, &kmax);

  RcppGSL::Vector x(4);
  x[0] = kmin - 1;
  x[1] = kmin;
  x[2] = kmax;
  x[3] = kmax + 1;

  RcppGSL::Matrix basis = nsbasis(x, knots);
  RcppGSL::Vector y = predict_natspl(basis, coefs);

  double c = log(normalize_dlogspline(basis, coefs, knots));
  for(int i = 0; i < 4; ++i) {
    y[i] = y[i] - c;
  }
  return y;
}


//' Evaluate the CDF of the log spectral density
//'
//' This function should probably not be used directly. It relies on there
//' being a (relatively) large number of \code{x} values supplied between the
//' boundary knots. The values of \code{x} are used to create a grid over
//' which the trapezoid rule is applied to determine the CDF values. If very
//' few \code{x} values are given, this function's behavior will be unreliable
//' and probably inaccurate.
//'
//' @param x The values at which to evaluate the CDF.
//' @param coefs A vector of length \code{k} representing the basis coefficients.
//' @param knots A vector of length \code{k} representing the knot locations.
//' @return A vector of CDF values. These are always returned in ascending order.
//' This is another reason you probably don't want to use this function directly.
//' @export
// [[Rcpp::export]]
RcppGSL::Vector plogspline(RcppGSL::Vector x, RcppGSL::Vector coefs, RcppGSL::Vector knots) {
  // make sure the input values are sorted first
  gsl_sort_vector(x);
  int n = x.size();                         // length of x
  double kmin, kmax;                        // store the minimum and maximum
  gsl_vector_minmax(knots, &kmin, &kmax);   //     knot values
  RcppGSL::Vector y(n);                     // storage for output values
  // xmid will contain the values inside the knot boundaries
  std::vector<double> xmid;
  xmid.push_back(kmin);                     // set kmin as the first value

  // calculate the slopes
  RcppGSL::Vector points = get_slope_pts(coefs, knots);
  double a = points[1] - points[0];         // left tail slope (positive)
  double b = points[3] - points[2];         // right tail slope (negative)
  int first_md = 0;                         // index of first value greater than kmin

  // loop through x
  for(int i = 0; i < n; ++i) {
    if(x[i] < kmin) {
      // if x[i] is in the left tail, we can integrate exactly
      y[i] = 1/a * exp(points[1] - a * (kmin - x[i]));
      first_md++;                           // haven't gotten to xmin yet
    } else if(x[i] > kmax) {
      // if x[i] is in the right tail, we can integrate exactly
      y[i] = 1 + 1/b * exp(points[2] - b * (kmax - x[i]));
    } else {
      // otherwise, add it to xmid
      xmid.push_back(x[i]);
    }
  }

  int m = xmid.size();
  if(m == 1) {
    // if there were no points between kmin and kmax, y should be all filled up
    return y;
  }

  // otherwise, perform numerical integration
  RcppGSL::Vector xm(m);
  // not sure how else to convert from std::vector to RcppGSL::Vector ...
  for(int i = 0; i < m; ++i) {
    xm[i] = xmid[i];
  }

  RcppGSL::Vector pdf_vals = dlogspline(xm, coefs, knots);
  // trapezoid rule
  y[first_md] = 1/a * exp(points[1]) + (pdf_vals[1] + pdf_vals[0]) / 2 * (xm[1] - xm[0]);
  for(int i = 1; i < m-1; ++i) {
    y[first_md+i] = y[first_md+i-1] + (pdf_vals[i+1] + pdf_vals[i]) / 2 * (xm[i+1] - xm[i]);
  }

  return y;
}


//' Evaluate the inverse CDF of the log spectral density
//'
//' AKA the quantile function.
//'
//' @param p The quantiles at which to evaluate the quantile function
//' @param coefs A vector of length \code{k} representing the basis coefficients
//' @param knots A vector of length \code{k} representing the knot locations
//' @param grid_size The size of the grid that the integral will be evaluated on
//' @return A vector of inverse CDF values. These will always be returned in
//' increasing order, regardless of the original order of the quantiles.
//' @export
// [[Rcpp::export]]
RcppGSL::Vector qlogspline(RcppGSL::Vector p, RcppGSL::Vector coefs, RcppGSL::Vector knots, double grid_size = 1e-2) {
  int n = p.size();                       // number of quantiles
  double kmin, kmax;                      // store min and max knot locations
  RcppGSL::Vector y(n);                   // store the output
  // pmid will contain the quantiles that correspond to CDF values within
  // the knot boundaries
  std::vector<double> pmid;
  int first_md = 0;                       // index of the first quantile in pmid

  // make sure p is sorted first
  gsl_sort_vector(p);

  gsl_vector_minmax(knots, &kmin, &kmax); // kmin and kmax now have correct values

  int len = floor((kmax - kmin) / grid_size);

  RcppGSL::Vector points = get_slope_pts(coefs, knots);
  double a = points[1] - points[0];
  double b = points[3] - points[2];

  double sm_cutoff = 1/a * exp(points[1]);
  double lg_cutoff = 1 + 1/b * exp(points[2]);

  // Rprintf("points: (%g,%g) (%g,%g) (%g,%g) (%g,%g)\n", kmin-1, gsl_vector_get(points,0), kmin, gsl_vector_get(points,1), kmax, gsl_vector_get(points,2), kmax+1, gsl_vector_get(points,3));
  // Rprintf("slopes: %g, %g\n", a, b);
  // Rprintf("cutoffs: %g, %g\n", sm_cutoff, lg_cutoff);

  for(int i = 0; i < n; ++i) {
    if(p[i] < sm_cutoff) {
      y[i] = 1/a * log(a*p[i]) + kmin - points[1]/a;
      first_md++;
    } else if(p[i] > lg_cutoff) {
      y[i] = 1/b * log(b*(p[i]-1)) + kmax - points[2]/b;
    } else {
      pmid.push_back(p[i]);
    }
  }

  int m = pmid.size();
  RcppGSL::Vector pm(m);
  for(int i = 0; i < m; ++i) {
    pm[i] = pmid[i];
  }


  RcppGSL::Vector x_vals(len+1);
  for(int i = 0; i < len; ++i) {
    x_vals[i] = grid_size * i + kmin;
  }
  x_vals[len] = kmax;

  RcppGSL::Vector cdf_vals = plogspline(x_vals, coefs, knots);
  double dx_val[len+1], dcdf_val[len+1];
  for(int j = 0; j < len+1; ++j) {
    dx_val[j] = x_vals[j];
    dcdf_val[j] = cdf_vals[j];
  }

  // return cdf_vals;

  gsl_set_error_handler_off();
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, len+1);
  int err = gsl_spline_init(spline, dcdf_val, dx_val, len+1);
  if(ISNAN(err)) {
    Rprintf("error = %s\n", gsl_strerror(err));
    Rcpp::stop("error in interpolation");
  }

  for(int i = 0; i < m; ++i) {
    // y[first_md+i] = gsl_spline_eval(spline, pm[i], acc);
    double yi;
    int err = gsl_spline_eval_e(spline, (double)pm[i], acc, &yi);
    if(err != 0) {
      Rprintf("error = %s\n", gsl_strerror(err));
      Rprintf("\tinput = %g\n", (double)pm[i]);
      y[first_md+i] = 0.0;
    } else {
      y[first_md+i] = yi;
    }
  }

  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  return y;

}


//' Generate random draws according to the log spectral density
//'
//' @param n The number of draws
//' @param coefs A vector of length \code{k} representing the basis coefficients
//' @param knots A vector of length \code{k} representing the knot locations
//' @param grid_size The size of the grid that the integral will be evaluated on
//' @return A vector of length \code{n} of random draws from the log spectral density
//' defined by \code{coefs} and \code{knots}
//' @export
// [[Rcpp::export]]
RcppGSL::Vector rlogspline(int n, RcppGSL::Vector coefs, RcppGSL::Vector knots, double grid_size = 1e-2) {
  RcppGSL::Vector r(n);
  NumericVector rand = runif(n, 0, 1);
  for(int i = 0; i < n; ++i) {
    r[i] = rand[i];
  }
  return qlogspline(r, coefs, knots, grid_size);
}

