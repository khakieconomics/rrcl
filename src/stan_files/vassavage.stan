functions {
  // A convenience function to avoid having to do the same matrix multiplication inside both functions below
  matrix tastes(real alpha, vector beta, vector tau, matrix Omega_L, matrix Z, int alpha_negative) {
    matrix[rows(Z), rows(beta) + 1] pars = rep_matrix(append_row(alpha, beta)', rows(Z)) + Z * diag_matrix(tau) * Omega_L;
    if(alpha_negative == 1) {
      pars[:,1] = -exp(pars[:,1]);
    }
    return(pars);
  }
  matrix raw_utilities(vector price, matrix X, vector xi, matrix tastes) {
    matrix[rows(tastes), rows(price) + 1] raw_utils = append_col(rep_vector(0.0, rows(tastes)),
                                                                tastes *
                                                                append_col(price, X)' +
                                                                rep_matrix(xi', rows(tastes)));
    return(raw_utils);
  }

  matrix probs(matrix raw_utils) {
    matrix[rows(raw_utils), cols(raw_utils)] pp;
    for(i in 1:rows(raw_utils)) {
      pp[i] = softmax(raw_utils[i]')';
    }
    return(pp);
  }
  // Performs numerical intergration across the random tastes constructed with tastes()
  vector shares(matrix probs) {
    vector[cols(probs)] share;
    for(j in 1:cols(probs)) {
      share[j] = mean(probs[:,j]);
    }
    return(share);
  }

  // Takes the gradient of market shares with respect to price
  vector dsdp(matrix probs, matrix tastes) {
    matrix[rows(probs), cols(probs)-1] derivatives;
    vector[cols(probs)-1] out;
    derivatives = rep_matrix(tastes[:,1], cols(probs)-1) .* probs[:,2:] .* (1 - probs[:,2:]);
    for(j in 1:cols(derivatives)) {
      out[j] = mean(derivatives[:,j]);
    }
    return(out);
  }

  // The best response price function given marginal prices and competitors' prices,
  // asusming each producer owns only a single product.
  vector price_fn(vector marginal_cost, vector price,
                  matrix X, vector xi, matrix tastes) {
    matrix[rows(tastes), rows(X) + 1] raw_utils = raw_utilities(price, X, xi, tastes);
    matrix[rows(tastes), rows(X) + 1] pp = probs(raw_utils);
    vector[rows(X)+1] shares_pred = shares(pp);
    vector[rows(X)] gradient = dsdp(pp, tastes);
    vector[rows(X)] out = marginal_cost - shares_pred[2:] ./ gradient;
    return(out);

  }

  // Simulates equilibrium in the market, ie. finds market clearing prices
  vector find_equilibrium_price(vector marginal_cost,
                                vector initial_price,
                                matrix X, vector xi, real alpha,
                                vector beta, vector tau, matrix Omega_L,
                                matrix Z, int alpha_negative) {
    real diff = 1;
    vector[rows(X)] old_price = initial_price;
    int count = 0;
    matrix[rows(Z), cols(X)+1] part_worths = tastes(alpha, beta, tau, Omega_L, Z, alpha_negative);
    while(diff > 0.01) {
      vector[rows(X)] new_price = price_fn(marginal_cost, old_price, X, xi, part_worths);
      count += 1;
      if(max(new_price) > 1e5 || min(new_price) < 0) {
        reject("Iteration diverging -- check that marginal costs are positive and alphas are negative");
      }
      if(count > 300) {
        reject("did not converge. Min new price = ", min(new_price), " max new price = ", max(new_price));
      }
      diff = sum(fabs(new_price - old_price));
      old_price = new_price;
    }
    return(old_price);
  }

}
data {
  int N; // rows of sales
  int sales[N]; // sales count for each product
  vector[N] price; // the price
  int T; // number of markets
  int market[T]; // index of market
  int outside_sales[T]; // outside sales in each market
  int start_index[T]; // which observation is the first in market?
  int end_index[T]; // which observation is the last?
  int P; // number of product attributes
  int P2; // number of instruments
  matrix[N, P] X; // product attributes
  matrix[N, P2] W; // instrument matrix
  int<lower = 0, upper = 1> alpha_negative;
  int NS; // number of sims
  matrix[NS, P+1] Z;
}
parameters {
  vector[N] xi; // structural shocks
  vector[P] beta; // means of preferences on product attributes
  real alpha; // mean on preference on price
  cholesky_factor_corr[P+1] L_Omega; //  Cholesky factor of correlation matrix of preferences across people
  vector<lower = 0>[P+1] tau;
  vector[P + P2] gamma;
  real<lower = 0> sigma_mc; // noise in implicit marginal cost
}
model {
  // priors
  xi~ student_t(4, 0, 2);
  beta[2:] ~ student_t(4, 0, .5);
  beta[1] ~ student_t(4, 0, 2);

  alpha ~ normal(0, 2);
  L_Omega ~ lkj_corr_cholesky(10);
  tau ~ normal(0, .5);
  gamma[1] ~ student_t(4, 0, 3);
  gamma[2:] ~ student_t(6, 0, 1);
  sigma_mc ~ student_t(4, 0, 1);
  // model for price and



  //

  {
    matrix[NS, P+1] ranef = tastes(alpha, beta, tau, L_Omega, Z, alpha_negative);
    for(t in 1:T) {
      matrix[NS, end_index[t] - start_index[t] + 2] raw_utils = raw_utilities(price[start_index[t]:end_index[t]],
                                                                              X[start_index[t]:end_index[t],:],
                                                                              xi[start_index[t]:end_index[t]],
                                                                              ranef);
      matrix[NS, end_index[t] - start_index[t] + 2] pp = probs(raw_utils);
      vector[end_index[t] - start_index[t] + 2] market_shares = shares(pp);
      vector[end_index[t] - start_index[t] + 1] derivs = dsdp(pp, ranef);
      int v[end_index[t] - start_index[t] + 2];
      v[1] =  outside_sales[t];
      v[2:] = sales[start_index[t]:end_index[t]];

      // sales model
      target += multinomial_lpmf(v | market_shares);

      // price model
      {
        vector[end_index[t] - start_index[t] + 1] mc_location = append_col(X[start_index[t]:end_index[t],:], W[start_index[t]:end_index[t],:]) * gamma;

        target += normal_lpdf(price[start_index[t]:end_index[t]] | mc_location - market_shares[2:] ./ derivs , sigma_mc);
      }


    }
  }
}
