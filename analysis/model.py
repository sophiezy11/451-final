import cmdstanpy
import pandas as pd
import numpy as np

# Load your training dataset
df = pd.read_csv("../features/training_data.csv")

# Prepare the data dictionary for Stan
stan_data = {
    'N': len(df),  # Number of data points
    'motif_score': df["motif_score"].values,
    'chip_signal': df["chip_signal"].values,
    'y': df["label"].values
}

# Path to the Stan model
stan_model_path = 'model.stan'

# Compile the Stan model
stan_model = cmdstanpy.CmdStanModel(stan_file=stan_model_path)

# Fit the model using the data
fit = stan_model.sample(data=stan_data)

# Extract the posterior samples for model parameters
posterior_alpha = fit.stan_variable("alpha")
posterior_beta_motif = fit.stan_variable("beta_motif")
posterior_mu_0 = fit.stan_variable("mu_0")
posterior_sigma_0 = fit.stan_variable("sigma_0")
posterior_mu_1 = fit.stan_variable("mu_1")
posterior_sigma_1 = fit.stan_variable("sigma_1")

# Check the shape of the posterior samples (for debugging purposes)
print(f"Shape of posterior_alpha: {posterior_alpha.shape}")
print(f"Shape of posterior_beta_motif: {posterior_beta_motif.shape}")
print(f"Shape of posterior_mu_0: {posterior_mu_0.shape}")
print(f"Shape of posterior_sigma_0: {posterior_sigma_0.shape}")

# For each data point, calculate the posterior prediction for binding probability
# Using the posterior samples of alpha, beta_motif, mu_0, mu_1, sigma_0, sigma_1

# Get the mean values of the posterior for parameters (or use samples for more detailed analysis)
alpha_mean = np.mean(posterior_alpha)
beta_motif_mean = np.mean(posterior_beta_motif)

mu_0_mean = np.mean(posterior_mu_0)
sigma_0_mean = np.mean(posterior_sigma_0)
mu_1_mean = np.mean(posterior_mu_1)
sigma_1_mean = np.mean(posterior_sigma_1)

# Calculate the posterior probabilities (probability of binding given motif and chip signal)
# Logistic regression for the probability of binding
logit_values = alpha_mean + beta_motif_mean * df["motif_score"].values
binding_probs = 1 / (1 + np.exp(-logit_values))  # Sigmoid function to get probabilities

# Now calculate the likelihood for each data point (using conditional Gaussian)
chip_signal_log = np.log1p(df["chip_signal"].values)
likelihood_0 = (1 / (sigma_0_mean * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((chip_signal_log - mu_0_mean) / sigma_0_mean) ** 2)
likelihood_1 = (1 / (sigma_1_mean * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((chip_signal_log - mu_1_mean) / sigma_1_mean) ** 2)

# Combine prior and likelihood to calculate posterior (using Bayes' rule)
posterior_binding = (binding_probs * likelihood_1) / (binding_probs * likelihood_1 + (1 - binding_probs) * likelihood_0)

# Convert posterior predictions to a DataFrame
posterior_df = pd.DataFrame(posterior_binding, columns=["posterior_prediction"])

# Save the posterior predictions to a CSV file
posterior_df.to_csv("posterior_predictions.csv", index=False)

print(f"Posterior predictions saved to 'posterior_predictions.csv'")
