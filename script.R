library(ggplot2)
library(dplyr)
library(gendist)
library(cowplot)


# 0. functions part
# Harell-Davis quantile estimator
hd <- function(data, q, clean_na=TRUE) {
  if (clean_na) {
    data <- data[!is.na(data)]
  }
  
  n <- length(data)
  a <- q * (n + 1)
  b <- (1 - q) * (n + 1)
  vec <- seq(along = data)
  weights <- pbeta(vec / n, a, b) - pbeta((vec - 1) / n, a, b)
  
  sorted_data <- sort(data)
  hd <- sum(weights * sorted_data)
  hd
}

# Shift function
shift_func <- function(x, y) {
  m <- matrix(0,9,4)
  
  for (d in 1:9) {
    q <- d/10
    m[d,1] <- q
    m[d,2] <- hd(x, q)
    m[d,3] <- hd(y, q)
    m[d,4] <- m[d,2] - m[d,3]
  }
  
  res <- data.frame(m)
  names(res) <- c('quantile',
                'group1','group2',
                'emp_diff')
  res
}



# global vars
N <- 1e4
set.seed(17)

# actual pipeline
# 1. generate distribution

gen_dists <- function() {
  # (1) U[11, 17] vs U[0, 900]
  unif1_min <- 11
  unif1_max <- 17
  
  unif2_min <- 0
  unif2_max <- 900
  unif1 <- runif(N, unif1_min, unif1_max)
  unif2 <- runif(N, unif2_min, unif2_max)

  res <- data.frame(unif1, unif2)
  
  # (2) Cauchy(0, 0.05) vs Beta(3, 3)
  cauchy_loc <- 0
  cauchy_scale <- 0.05
  
  beta_a <- 3
  beta_b <- 3
  
  cauchy_d <- rcauchy(N, cauchy_loc, cauchy_scale)
  beta_d <- rbeta(N, beta_a, beta_b)
  
  res <- cbind(res, cauchy_d, beta_d)
  
  # (3) Bimodal: 0.5 N(4, 2) + 0.5 N(20, 2) vs Bimodal: 0.5 N(7, 2) + 0.5 N(15, 2)
  mixt_d1 <- rmixt(N, 0.5, 
                   spec1="norm", arg1 = list(mean=4, sd=2),
                   spec2="norm", arg2 = list(mean=20, sd=2),
                   interval=c(-1e3,1e3))
  mixt_d2 <- rmixt(N, 0.5, 
                   spec1="norm", arg1 = list(mean=7, sd=2),
                   spec2="norm", arg2 = list(mean=15, sd=2),
                   interval=c(-1e3,1e3))
  
  
  
  res <- cbind(res, mixt_d1, mixt_d2)
  
  
  # (4) Unimodal right-skewed: Weibull(1, 1.5) vs Exp(2)
  weibull_d <- rweibull(N, 1, 1.5)
  exp_d <- rexp(N, 2)
  
  res <- cbind(res, weibull_d, exp_d)
  
  res
}

gen_qs <- function() {
  res <- data.frame(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
  
  
  # (1) U[11, 17] vs U[0, 900]
  unif1_min <- 11
  unif1_max <- 17
  
  unif2_min <- 0
  unif2_max <- 900
  qs <- c()
  for (d in 1:9) {
    q <- d / 10
    diff <- qunif(q, unif1_min, unif1_max) - qunif(q, unif2_min, unif2_max)
    qs <- append(qs, diff)
    
  }
  res <- cbind(res, qs)
  
  # (2) Cauchy(0, 0.5) vs Beta(3, 3)
  cauchy_loc <- 0
  cauchy_scale <- 0.05
  
  beta_a <- 3
  beta_b <- 3
  qs <- c()
  for (d in 1:9) {
    q <- d / 10
    diff <- qcauchy(q, cauchy_loc, cauchy_scale) - qbeta(q, beta_a, beta_b)
    qs <- append(qs, diff)
    
  }
  
  res <- cbind(res, qs)
  
  # (3) Bimodal: 0.5 N(4, 2) + 0.5 N(20, 2) vs Bimodal: 0.5 N(7, 2) + 0.5 N(15, 2)
  qs <- c()
  for (d in 1:9) {
    q <- d / 10
    q1 <- qmixt(q, 0.5, 
                  spec1="norm", arg1 = list(mean=4, sd=2),
                  spec2="norm", arg2 = list(mean=20, sd=2),
                  interval=c(-1e3,1e3))
    q2 <- qmixt(q, 0.5, 
               spec1="norm", arg1 = list(mean=7, sd=2),
               spec2="norm", arg2 = list(mean=15, sd=2),
               interval=c(-1e3,1e3))
    
    diff <- q1 - q2
    qs <- append(qs, diff)
    
  }
  res <- cbind(res, qs)
  
  # (4) Unimodal right-skewed: Weibull(1, 1.5) vs Exp(2)
  qs <- c()
  for (d in 1:9) {
    q <- d / 10
    diff <- qweibull(q, 1, 1.5) - qexp(q, 2)
    qs <- append(qs, diff)
    
  }
  res <- cbind(res, qs)
  
  names(res) <- 
    c(
    'quantile',
    'uni_uni',
    'cauchy_beta',
    'mixt_norm',
    'weibull_exp'
    )
  res
}

distributions <- gen_dists()
theo_quant_diffs <- gen_qs()


# 2. construct final output
# (1) U[11, 17] vs U[0, 900]
add_theo_quants <- function(df, theo_df, case) {
  #df <- cbind(df, theo_df[case])
  z <- theo_df[case]
  df$theo_diff <- unlist(unname(z))
  df$diff_in_diff <- abs(df$theo_diff - df$emp_diff)
  
  df
}

# (1)
shift_res <- shift_func(distributions$unif1, distributions$unif2)
uni_uni_out <- add_theo_quants(shift_res, theo_quant_diffs, "uni_uni")
plot_diff(uni_uni_out)

# (2)
shift_res <- shift_func(distributions$cauchy_d, distributions$beta_d)
cauchy_beta_out <- add_theo_quants(shift_res, theo_quant_diffs, "cauchy_beta")
plot_diff(cauchy_beta_out)

# (3)
shift_res <- shift_func(distributions$mixt_d1, distributions$mixt_d2)
mixt_norm_out <- add_theo_quants(shift_res, theo_quant_diffs, "mixt_norm")
plot_diff(mixt_norm_out)

# (4)
shift_res <- shift_func(distributions$weibull_d, distributions$exp_d)
weibull_exp_out <- add_theo_quants(shift_res, theo_quant_diffs, "weibull_exp")
plot_diff(weibull_exp_out)



# 3. plot all results

clrs <- c('Theoretical' = "#55c40b",
          'Empirical' = "#e27830",
          'Absolute error' = "#16537e"
)

plot_diff <- function(data) {
  p <- ggplot(data, aes_string(x='quantile')) +
    geom_line(aes(y=theo_diff, color='Theoretical'), linetype="solid", size=0.3) +
    geom_line(aes(y=emp_diff, color='Empirical'),  linetype="solid", size=0.3) +
    geom_line(aes(y=diff_in_diff, color='Absolute error'),  linetype="solid", size=1.5) +
    ylab("difference") +
    labs(color="Legend") +
    scale_color_manual(values=clrs)
  p
}

final_plot <- function(d1, d2, out) {
  df <- data.frame(d1, d2)
  dp1 <- ggplot(df, aes(x=d1)) + 
    geom_density(color="#000000", fill="#999999", alpha=0.55)
  
  dp2 <- ggplot(df, aes(x=d2)) + 
    geom_density(color="#000000", fill="#999999", alpha=0.55)
  
  diff_plot <- plot_diff(out)
  plot_grid(diff_plot, 
            dp1, 
            NULL, dp2, 
            ncol=2,
            rel_widths = c(2, 1))
}


final_plot(distributions$unif1, distributions$unif2, uni_uni_out)

final_plot(distributions$cauchy_d, distributions$beta_d, cauchy_beta_out)

final_plot(distributions$mixt_d1, distributions$mixt_d2, mixt_norm_out)

final_plot(distributions$weibull_d, distributions$exp_d, weibull_exp_out)

