import pandas as pd
import os

# Load the full dataset
df = pd.read_csv("../features/train_dataset.csv")

# Sample 5000 entries randomly
sampled_df = df.sample(n=5000, random_state=42)

# Define the path to the original dataset directory
original_file_path = "../features/train_dataset.csv"
directory = os.path.dirname(original_file_path)

# Save the sampled data to the same directory as train_dataset.csv
sampled_df.to_csv(os.path.join(directory, "training_data.csv"), index=False)

print("Sampled data saved to 'training_data.csv'.")
