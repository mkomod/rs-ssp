library("rstan")

n <- 100
omega <- 1
censoring_lvl <- 0.4

set.seed(1)
b <- c(1, 1, rep(0, 8))
p <- length(b)
X <- matrix(rnorm(n * p), nrow=n)
y <- runif(nrow(X))
Y <- log(1 - y) / - (exp(X %*% b) * omega)

delta  <- runif(n) > censoring_lvl   # 0: censored, 1: uncensored
Y[!delta] <- Y[!delta] * runif(sum(!delta))
ind_sorted_T <- order(Y)
ind_failure_T <- which(as.logical(delta[order(Y)]))
N_delta <- sum(delta)

# laplace prior
m.1 <- stan(file="l1.stan", 
	    data=list(N = n, N_delta = N_delta, P=p, ind_sorted_T=ind_sorted_T,
		      ind_failure_T = ind_failure_T, ind = N_delta:1, X=X, lambda=0.1), 
	    cores=3, iter=1e4, warmup=1e3)

traceplot(m.1)
print(m.1)
summary(m.1)
stan_ac(m.1)
stan_hist(m.1)
stan_plot(m.1)
stan_ess(m.1)
stan_rhat(m.1)
stan_diag(m.1)
stan_dens(m.1)
get_posterior_mean(m.1)
get_elapsed_time(m.1)

# horseshoe prior
m.2 <- stan(file="horseshoe.stan", control = list(adapt_delta = 0.99, max_treedepth=12),
	    data=list(N = n, N_delta = N_delta, P=p, ind_sorted_T=ind_sorted_T,
		      ind_failure_T = ind_failure_T, ind = N_delta:1, X=X, lambda=0.1), 
	    cores=3, iter=2e4, warmup=5e3, chains=2)

traceplot(m.2)
summary(m.2)


