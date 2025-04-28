import cmdstanpy

# Install CmdStan
cmdstanpy.install_cmdstan()

# Verify installation
print("CmdStan path:", cmdstanpy.cmdstan_path())