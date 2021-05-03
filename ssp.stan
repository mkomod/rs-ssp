functions {
  real pl_lpdf(data matrix X, vector b, data int[] ind_sorted_T, 
        data int[] ind_failure_T, data int[] ind, data int N) 
  {
    vector[N] xb = X * b;
    real a = max(xb);

    real risk = 0.0;
    real tot = 0.0;
    int ind_last_failure = N + 1;
    int ind_curr_failure = 0;

    for (i in ind) {
      ind_curr_failure = ind_failure_T[i];
      risk = risk + sum(exp(xb[ind_sorted_T[ind_curr_failure:(ind_last_failure - 1)]] - a));
      tot += xb[ind_sorted_T[ind_curr_failure]] - (a + log(risk));
      ind_last_failure = ind_curr_failure;
    }
    return tot;
  }
}
data {
  int<lower=0> N;
  int<lower=0> N_delta;
  int<lower=0> P;
  int ind_sorted_T[N];
  int ind_failure_T[N_delta];
  int ind[N_delta];
  matrix[N, P] X;
  real<lower=0> lambda;
}
parameters {
  vector[P] b;
}
model {
  b ~ double_exponential(0, lambda);
  target += pl_lpdf(X | b, ind_sorted_T, ind_failure_T, ind, N);
}
