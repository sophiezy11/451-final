data {
  int<lower=0> N;
  array[N] real motif_score;
  array[N] real chip_signal;
  array[N] int<lower=0, upper=1> y;
}
parameters {
  real alpha;
  real beta_motif;
  real mu_0;
  real<lower=0> sigma_0;
  real mu_1;
  real<lower=0> sigma_1;
}

model {
    // ========== PRIORS ========== //
    // Logistic regression parameters
    alpha ~ normal(0, 2);
    beta_motif ~ normal(0, 2);
    
    // ChIP-seq signal parameters
    mu_0 ~ normal(0, 2);
    mu_1 ~ normal(0, 2);
    sigma_0 ~ normal(0, 1);
    sigma_1 ~ normal(0, 1);
    
    // ========== LIKELIHOOD ========== //
    // Binding labels (Bernoulli-logit)
    for (n in 1:N)
        y ~ bernoulli_logit(alpha + beta_motif * motif_score[n]);
    
    // ChIP-seq signals (single likelihood block)
    for (n in 1:N) {
        if (y[n] == 1) {
            chip_signal[n] ~ normal(mu_1, sigma_1);
        } else {
            chip_signal[n] ~ normal(mu_0, sigma_0);
        }
    }
}
