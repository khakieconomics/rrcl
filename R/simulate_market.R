#' sim_markets
#'
#' Simulates market equilibrium prices for a number of markets given key structural parameters
#' under random coefficients logit and the assumption that each firm produces one product and
#' chooses only price to profit maximize.
#'
#' @param P number of choice attributes (these will be binary)
#' @param P2 number of cost-shifting excluded instruments
#' @param J number of products to be simulated-- each market will have different numbers of product options
#' @param T number of markets
#' @param NS the number of fake decision-makers to simulate. 500 by default
#' @param alpha alpha_i ~ normal(alpha, ...) for each decision-maker, with alpha_i either being the
#' coefficient on price for a given decisionmaker (when alpha_negative = 0) or -exp(alpha) being their
#' coefficient on price. 0.5 by default
#' @param alpha_negative 1 for alpha to be transformed to -exp(alpha), the default, 0 otherwise.
#' @param beta_mean the beta--"part worth" coefficients are drawn from beta ~ normal(beta_mean, sigma_beta)
#' @param sigma_beta see beta_mean
#' @param sigma_tau the marginal scale of the covariance matrix. tau ~ normal+(0, sigma_tau)
#' @param gamma_intercept the intercept in the marginal cost model. Increase this if you're getting negative
#' marginal costs. 4 by default.
#' @param sigma_gamma scale of the coefficients on the cost-shifters
#' @param market_size the number of participants in each market. 50k by default
#' @param xi_halfrange in the generative model, we assume xi is uniformly
#' distributed with zero mean. xi_halfrange is half the range (that is,
#' xi is uniform between -xi_halfrange, xi_halfrange). 3 by default.
#' @export
#' @import dplyr
#' @return a data frame containing shares

sim_markets <- function(P, P2, J, T, NS = 500, alpha = 0.5, alpha_negative = 1, beta_mean = 1,
                        sigma_beta = 1, sigma_tau = 0.5, gamma_intercept = 4, sigma_gamma = 0.3,
                        market_size = 50000, xi_halfrange = 3) {
  X <- matrix(abs(rbinom((P)*J, 1, .2)), J, P)
  beta <- rnorm(P, beta_mean, sigma_beta)
  gamma <- rnorm(P + P2, 0, sigma_gamma)
  gamma[1] <- gamma_intercept
  gamma <- gamma
  tau <- abs(rnorm(P+1, 0, sigma_tau))
  Omega <- cor(matrix(rnorm((P+2)*(P+1)), P+2, (P+1)))
  L <- t(chol(Omega))

  cost_shifting_instruments <- cbind(1, matrix(rbinom(J*(P2-1), 1, .2), J, P2-1)) %>% as_tibble(.name_repair = ~ c("intercept_col", paste0("demand_instruments", 0:(P2-2)))) %>%
    mutate(product_ids = 1:J)
  Z <- matrix(rnorm(NS*(P+1)), NS, P+1)
  preferences <- tastes(alpha, beta, tau, L, Z, alpha_negative = 1)

  tibble(market_ids = 1:T) %>%
    group_by(market_ids) %>%
    mutate(product_ids = list(sample(1:J, sample(3:20, 1)))) %>%
    unnest() %>%
    arrange(market_ids, product_ids) %>%
    left_join(bind_cols(product_ids = 1:J, as.data.frame(X)), by = "product_ids") %>%
    mutate(xi = runif(n(), -xi_halfrange, xi_halfrange)) %>%
    left_join(cost_shifting_instruments, by = "product_ids") %>%
    group_by(market_ids) %>%
    do({
      df <- as.data.frame(.)
      mc <- (df %>% select("intercept_col", dplyr::contains("demand_instruments"), dplyr::contains("V")) %>% as.matrix) %*% gamma
      mc <- mc + rnorm(nrow(df), 0, .3)
      if(min(mc) < 0 | max(mc) > 300) stop("Marginal costs must be positive and less than $300")
      initial_prices <- rep(0.1, nrow(df))
      X <- df %>% select(dplyr::contains("V")) %>% as.matrix
      prices <- find_equilibrium_price(marginal_cost = mc, initial_price = initial_prices, X = X, xi = df$xi, alpha = alpha, beta = beta, tau = tau, Omega_L = L,Z =  Z, alpha_negative = 1)
      utilities <- raw_utilities(price = prices, X = X, xi = df$xi, tastes = preferences)
      market_share <- shares(probs(utilities))
      sales <- rmultinom(1, size = market_size, prob = market_share)
      recorded_shares <- sales/sum(sales)
      bind_cols(df, tibble(prices = prices, costs = mc,shares = pmax(0.000001, recorded_shares[-1]), sales = sales[-1], outside_sales = market_size - sum(sales), outside_share = recorded_shares[1]))
    }) %>% ungroup -> dat_out

  list(market_data = dat_out, alpha = alpha, beta = beta, L = L, tau = tau, gamma = gamma)
}

#' gp
#'
#' get parameter point estimates from a Stan fit (optimizing, vb or sampling)
#'
#' @param f the fit
#' @param p the parameter name
#' @importFrom rstan get_posterior_mean
#' @export

gp <- function(f, p) {
  if(class(f) == "list") {
    f$par[grepl(p, names(f$par))]
  } else {
    yy <- get_posterior_mean(f, pars = p)
    yy[,ncol(yy)]
  }
}
