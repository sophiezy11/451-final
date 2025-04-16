data {
  int<lower=0> N;                // Number of data points
  array[N] real motif_score;     // Motif scores
  array[N] real chip_signal;     // ChIP-seq signal
  array[N] int<lower=0, upper=1> y;   // Binding labels (0 or 1)
}

parameters {
    real alpha;                     // intercept for the logistic regression
    real beta_motif;                // coefficient for motif score
    real mu_0;                       // mean of chip signal when binding = 0
    real<lower=0> sigma_0;          // standard deviation of chip signal when binding = 0
    real mu_1;                       // mean of chip signal when binding = 1
    real<lower=0> sigma_1;          // standard deviation of chip signal when binding = 1
}

model {
    // Logistic regression prior for binding probability
    for (n in 1:N) {
        y[n] ~ bernoulli_logit(alpha + beta_motif * motif_score[n]);
    }

    // Conditional Gaussian likelihood for chip-seq signal
    for (n in 1:N) {
        if (y[n] == 1) {
            chip_signal[n] ~ normal(mu_1, sigma_1);  // Binding = 1
        } else {
            chip_signal[n] ~ normal(mu_0, sigma_0);  // Binding = 0
        }
    }
}
