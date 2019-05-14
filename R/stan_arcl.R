#' stan_arcl
#'
#' Fits an aggregate random coefficients logit model in Stan
#' @param data your dataset
#' @param NS the number of individuals to be simulated in each market
#' @param sales string. The column name in data where sales data belong. \code{sales} by default
#' @param market_ids string. Column name for the identifier of the markets.  \code{market} by default
#' @param prices string. Column name for the price.  \code{price} by default
#' @param X_names string. Column names for the product attributes
#' @param W_names string. Column names for excluded instruments that affect marginal costs
#' @param alpha_negative should we constrain alpha to be negative? 1 for negatgive, 0 otherwise.
#' @param init passed to optimizing
#' @export
#' @importFrom dplyr "%>%" ungroup mutate group_by summarise first last select_
#' @importFrom rstan optimizing

stan_arcl <- function(data, NS = 100, sales = "sales", market_ids = "market_ids", prices = "prices", X_names, W_names, alpha_negative = 1, init = 0) {
  data %>%
    ungroup %>%
    mutate(row = 1:n()) %>%
    group_by(market_ids) %>%
    summarise(start_index = first(row),
              end_index = last(row),
              outside_sales = first(outside_sales)) %>%
    as.data.frame()-> Market_data


  data_list <- list(N = nrow(data),
                    NS = NS,
                    P = length(X_names),
                    P2 = length(W_names),
                    sales = as.numeric(data[[sales]]),
                    market = as.numeric(Market_data[,market_ids]),
                    outside_sales = as.numeric(Market_data[,"outside_sales"]),
                    start_index = as.numeric(Market_data[,"start_index"]),
                    end_index = as.numeric(Market_data[,"end_index"]),
                    T = max(Market_data[,market_ids]),
                    price = as.numeric(data[[prices]]),
                    X = data %>% select(X_names),
                    W = data %>% select(W_names),
                    Z = matrix(rnorm(NS * (length(X_names)+1)), NS, length(X_names)+1),
                    alpha_negative = alpha_negative)

  out <- rstan::optimizing(stanmodels$vassavage, data = data_list, init = init)
  return(out)
}
