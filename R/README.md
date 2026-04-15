library(SynPALM)
library(Matrix)

# 1. Initialize Simulation Parameters
n_obs <- 1000        # Number of samples with observed primary phenotype
missing_rate <- 0.5  # 50% missingness
maf <- 0.3           # Minor Allele Frequency
rho <- 0.5           # Correlation for siblings

# 2. Generate Genotypes
# Creates a standardized genotype matrix for the full cohort
G_matrix <- generate_genotypes_rho(
  n_obs = n_obs, 
  miss = missing_rate, 
  n_snps = 1, 
  maf = maf, 
  rho = rho
)

# 3. Create a Genetic Relatedness Matrix (GRM)
n_all <- length(G_matrix)
GRM <- as(diag(n_all), "dgCMatrix") 

# 4. Define Step 1 Parameters (Mocked for this example)
pre_vars <- list(
  chol_Sigma22 = diag(n_all), 
  L_cond = diag(n_obs),
  Sigma12_Sigma22inv = matrix(0.5, n_obs, n_all),
  Sigma12_oracle_Sigma22inv = matrix(0.5, n_all, n_all),
  L_cond_oracle = diag(n_all)
)

# 5. Simulate Phenotypes (DGP)
sim_data <- DGP_GRMr_gmodel_step3_rho(
  n_obs = n_obs, 
  miss = missing_rate, 
  pve_g = 0.005, 
  Genotype_each = G_matrix, 
  pre_vars_step2 = pre_vars, 
  GRM = GRM
)

# 6. Run the Syn-PALM Score Test
results <- score_test_SynSurrG_power_sim_one_g_to_one_Y(sim_data)

# 7. Print Results
print(results)
