// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;


double d_prior_beta(arma::vec beta, double tau, bool l) {
  if (tau <= 0) {
    if (l) {
      return -1*datum::inf;
    } else {
      return 0;
    }
  }

  double res = R::dnorm(beta(0), 0., sqrt(tau), true) + R::dnorm(beta(1), beta(0), sqrt(tau), true);
  int n = beta.n_elem;

  for(int i = 2; i < n; ++i) {
    res += R::dnorm(beta(i), 2*beta(i-1) - beta(i-2), sqrt(tau), true);
  }

  if(l) {
    return res;
  } else {
    return exp(res);
  }
}


double d_prior_sigma(double sigma, double a, double b, bool l) {
  if (sigma <= 0) {
    if (l) {
      return -1*datum::inf;
    } else {
      return 0;
    }
  }

  double x = a * log(b) - lgamma(a) - (a+1) * log(sigma) - b / sigma;
  if (l) {
    return x;
  } else {
    return exp(x);
  }
}


double d_prior_tau(double tau, double a, double b, bool l) {
  if (tau <= 0) {
    if (l) {
      return -1*datum::inf;
    } else {
      return 0;
    }
  }

  double x = a * log(b) - lgamma(a) - (a+1) * log(tau) - b / tau;
  if (l) {
    return x;
  } else {
    return exp(x);
  }
}


double d_prior(arma::vec beta,
               double sigma,
               double tau,
               double s_a,
               double s_b,
               double t_a,
               double t_b,
               bool l) {
  double x = d_prior_beta(beta, tau, true) + d_prior_sigma(sigma, s_a, s_b, true) + d_prior_tau(tau, t_a, t_b, true);
  if (l) {
    return x;
  } else {
    return exp(x);
  }
}


double lik(arma::vec y, arma::vec beta, double sigma, arma::mat spline, bool l) {
  if (sigma <= 0) {
    if (l) {
      return -1*datum::inf;
    } else {
      return 0;
    }
  }

  arma::vec avg = spline * beta;
  int n = avg.n_elem;
  double p = 0;

  for(int i = 0; i < n; ++i) {
    p += R::dnorm(y(i), avg(i), sqrt(sigma), l);
  }
  return p;
}

double r_q_beta(double beta, double v) {
  return R::rnorm(beta, sqrt(v));
}

double r_q_sigma(double sigma, double v) {
  return R::rnorm(sigma, sqrt(v));
}

double r_q_tau(double tau, double v) {
  return R::rnorm(tau, sqrt(v));
}

double d_q_beta(double beta_star, double beta, double v, bool l) {
  return R::dnorm(beta_star, beta, sqrt(v), l);
}

double d_q_sigma(double sigma_star, double sigma, double v, bool l) {
  return R::dnorm(sigma_star, sigma, sqrt(v), l);
}

double d_q_tau(double tau_star, double tau, double v, bool l) {
  return R::dnorm(tau_star, tau, sqrt(v), l);
}

//' Bayesian spline fitting, in C++
//'
//' Don't use this directly. Use the wrapper function \code{\link{fitspline}}.
//' @param B B
//' @param y y
//' @param spl spl
//' @param burnin burning
//' @param initBeta initBeta
//' @param initSigma initSigma
//' @param initTau initTau
//' @param initBetaTune initBetaTune
//' @param initSigmaTune initSigmaTune
//' @param initTauTune initTauTune
//' @param s_a s_a
//' @param s_b s_b
//' @param t_a t_a
//' @param t_b t_b
//' @param c0 c0
//' @param c1 c1
//' @param k k
//' @param r r
//' @param oar oar
//' @param progress progress
//' @return A list of length 3. The first element is a matrix where each row is a
//'   random \code{beta} draw. The second and third are vectors where each
//'   element is a random draw of \code{sigma} and \code{tau} respectively.
//' @export
// [[Rcpp::export]]
List fitsplinecpp(int B,
                arma::vec y,
                arma::mat spl,
                int burnin,
                arma::vec initBeta,
                double initSigma,
                double initTau,
                arma::vec initBetaTune,
                double initSigmaTune,
                double initTauTune,
                double s_a,
                double s_b,
                double t_a,
                double t_b,
                double c0,
                double c1,
                double k,
                double r,
                double oar,
                bool progress) {
  int n = spl.n_cols;

  // storage
  arma::mat beta_samp(B, n);
  arma::vec sigma_samp(B);
  arma::vec tau_samp(B);

  arma::mat beta_accepts(B, n, arma::fill::zeros);
  arma::vec sigma_accepts(B, arma::fill::zeros);
  arma::vec tau_accepts(B, arma::fill::zeros);

  arma::mat v_beta(B, n);
  arma::vec v_sigma(B);
  arma::vec v_tau(B);

  // starting values
  beta_accepts.row(0).fill(1);
  sigma_accepts(0) = 1;
  tau_accepts(0) = 1;

  v_beta.row(0) = initBetaTune.t();
  v_sigma(0) = initSigmaTune;
  v_tau(0) = initTauTune;

  arma::vec beta = initBeta;
  double sigma = initSigma;
  double tau = initTau;

  // loop
  for(int i = 1; i < B; ++i) {
    // adaptive tuning
    double gamma1 = c0 / std::pow(i + k, c1);
    arma::rowvec lfr_beta(n);
    double lfr_sigma;
    double lfr_tau;



    if(i == 1) {
      lfr_beta = beta_accepts.row(i-1);
      lfr_sigma = sigma_accepts(i-1);
      lfr_tau = tau_accepts(i-1);
    } else if(i < r+1) {
      lfr_beta = arma::mean(beta_accepts.rows(0, i-1), 0);
      lfr_sigma = arma::mean(sigma_accepts.subvec(0, i-1));
      lfr_tau = arma::mean(tau_accepts.subvec(0, i-1));
    } else {
      lfr_beta = arma::mean(beta_accepts.rows(i-r, i-1), 0);
      lfr_sigma = arma::mean(sigma_accepts.subvec(i-r, i-1));
      lfr_tau = arma::mean(tau_accepts.subvec(i-r, i-1));
    }

    v_beta.row(i) = exp(log(v_beta.row(i-1)) + gamma1 * (lfr_beta - oar));
    v_sigma(i) = exp(log(v_sigma(i-1)) + gamma1 * (lfr_sigma - oar));
    v_tau(i) = exp(log(v_tau(i-1)) + gamma1 * (lfr_tau - oar));

    // update each of the betas one after another
    for(int j = 0; j < n; ++j) {
      arma::vec beta_star = beta;
      beta_star(j) = r_q_beta(beta(j), v_beta(i,j));

      // numerator
      double num = lik(y, beta_star, sigma, spl, true) +
        d_prior(beta_star, sigma, tau, s_a, s_b, t_a, t_b, true) +
        d_q_beta(beta(j), beta_star(j), v_beta(i,j), true);

      // denominator
      double denom = lik(y, beta, sigma, spl, true) +
        d_prior(beta, sigma, tau, s_a, s_b, t_a, t_b, true) +
        d_q_beta(beta_star(j), beta(j), v_beta(i,j), true);

      double logr = num - denom;

      if(log(R::runif(0., 1.)) < logr) {
        // accept the proposal
        beta(j) = beta_star(j);
        beta_accepts(i,j) = 1;
      }

    }


    // update sigma
    double sigma_star = r_q_sigma(sigma, v_sigma(i));

    // numerator
    double num = lik(y, beta, sigma_star, spl, true) +
      d_prior(beta, sigma_star, tau, s_a, s_b, t_a, t_b, true) +
      d_q_sigma(sigma, sigma_star, v_sigma(i), true);

    // denominator
    double denom = lik(y, beta, sigma, spl, true) +
      d_prior(beta, sigma, tau, s_a, s_b, t_a, t_b, true) +
      d_q_sigma(sigma_star, sigma, v_sigma(i), true);

    double logr = num - denom;

    if(log(R::runif(0., 1.)) < logr) {
      // accept the proposal
      sigma = sigma_star;
      sigma_accepts(i) = 1;
    }


    // update tau
    double tau_star = r_q_tau(tau, v_tau(i));

    // numerator
    num = d_prior(beta, sigma, tau_star, s_a, s_b, t_a, t_b, true) +
      d_q_tau(tau, tau_star, v_tau(i), true);

    // denominator
    denom = d_prior(beta, sigma, tau, s_a, s_b, t_a, t_b, true) +
      d_q_tau(tau_star, tau, v_tau(i), true);

    logr = num - denom;

    if(log(R::runif(0., 1.)) < logr) {
      // accept the proposal
      tau = tau_star;
      tau_accepts(i) = 1;
    }

    // save the current values
    beta_samp.row(i) = beta.t();
    sigma_samp(i) = sigma;
    tau_samp(i) = tau;

  }

  arma::mat beta_final = beta_samp.tail_rows(B - burnin);
  arma::vec sigma_final = sigma_samp.tail(B - burnin);
  arma::vec tau_final = tau_samp.tail(B - burnin);

  // return
  List ret;
  ret["b"] = beta_final;
  ret["s"] = sigma_final;
  ret["t"] = tau_final;

  return ret;

}
