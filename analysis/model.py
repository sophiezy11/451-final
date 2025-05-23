import cmdstanpy
import pandas as pd
import numpy as np
import os

# Load your training dataset
df = pd.read_csv("../features/training_data.csv")

# Prepare the data dictionary for Stan
stan_data = {
    'N': len(df),  # Number of data points
    'motif_score': df["motif_score"].values,
    'chip_signal': df["chip_signal"].values,
    'y': df["label"].values
}

# Different models to try
model_versions = {
    'baseline': 'model_baseline.stan',
    'wide_prior': 'model_wide.stan',
    'tight_prior': 'model_tight.stan',
    'shifted_prior': 'model_shifted.stan'
}

# Folder to save everything
output_dir = "posterior_outputs"
os.makedirs(output_dir, exist_ok=True)

# Collect summary stats
summary_records = []

for version_name, stan_file in model_versions.items():
    print(f"Running {version_name}...")
    # Compile and fit the model
    model = cmdstanpy.CmdStanModel(stan_file=stan_file)
    fit = model.sample(data=stan_data)

    # Save the full posterior draws
    fit.save_csvfiles(dir=os.path.join(output_dir, f"{version_name}_cmdstan"))

    # Save posterior summaries
    summary = fit.summary()
    summary.to_csv(os.path.join(output_dir, f"summary_{version_name}.csv"))

    # Collect a few key parameter summaries for plotting later
    for param in ["alpha", "beta_motif", "mu_0", "mu_1"]:
        mean_val = summary.loc[param, "Mean"]
        sd_val = summary.loc[param, "StdDev"]
        summary_records.append({
            "version": version_name,
            "param": param,
            "mean": mean_val,
            "sd": sd_val
        })

# Save the summary across all experiments
summary_df = pd.DataFrame(summary_records)
summary_df.to_csv(os.path.join(output_dir, "all_posterior_summaries.csv"), index=False)

print("Finished all model fits!")

# For each data point, calculate the posterior prediction for binding probability
# Using the posterior samples of alpha, beta_motif, mu_0, mu_1, sigma_0, sigma_1


# # Calculate the posterior probabilities (probability of binding given motif and chip signal)
# # Logistic regression for the probability of binding
# logit_values = alpha_mean + beta_motif_mean * df["motif_score"].values
# binding_probs = 1 / (1 + np.exp(-logit_values))  # Sigmoid function to get probabilities

# # Now calculate the likelihood for each data point (using conditional Gaussian)
# chip_signal_log = np.log1p(df["chip_signal"].values)
# likelihood_0 = (1 / (sigma_0_mean * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((chip_signal_log - mu_0_mean) / sigma_0_mean) ** 2)
# likelihood_1 = (1 / (sigma_1_mean * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((chip_signal_log - mu_1_mean) / sigma_1_mean) ** 2)

# # Combine prior and likelihood to calculate posterior (using Bayes' rule)
# posterior_binding = (binding_probs * likelihood_1) / (binding_probs * likelihood_1 + (1 - binding_probs) * likelihood_0)

# # Convert posterior predictions to a DataFrame
# posterior_df = pd.DataFrame(posterior_binding, columns=["posterior_prediction"])

# # Save the posterior predictions to a CSV file
# posterior_df.to_csv("posterior_predictions.csv", index=False)

# print(f"Posterior predictions saved to 'posterior_predictions.csv'")
