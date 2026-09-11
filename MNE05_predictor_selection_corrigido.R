# Ecological niche modelling
# by Fernando Paulino and Carolina Loss
# JAN - 2025
# Fixed - SEP/2025 (layer selection by name, overwrite, more robust writing)

# Load the packages required for the analysis
library(usdm)    # Para análise de multicolinearidade (VIF)
library(dismo)   # Para modelagem de distribuição de espécies
library(terra)   # Para manipulação de dados espaciais

# 1. LOAD OCCURRENCE DATA
# Read the CSV file with occurrence records
pts <- read.csv("./data/occour/Caimans.csv")
# Standardise species names (replace spaces with underscores)
pts$sp <- gsub(" ", "_", pts$sp)
# Get the unique list of species/population names
# (renamed from "names" to "species_names" so it does not shadow the names() function)
species_names <- unique(pts$sp)

# 2. DEFINE MAIN SPECIES FOR VARIABLE SELECTION
# Caiman_latirostris = species as a whole (aggregates GAC+GPI / the 5 units
# units), so it serves as the reference for calibrating VIF over the
# the species' full environmental range.
sp_principal <- "Caiman_latirostris"

# 3. MAIN SPECIES PROCESSING (VARIABLE SELECTION)
# Set the path to the main species' environmental predictors
pred_folder <- paste0("./output/predictors/predmod/predmodcrop/", sp_principal, "/")
# List ASC (predictor) files, already with full paths (avoids manual paste)
pred_files <- list.files(path = pred_folder, pattern = '*.asc$', full.names = TRUE)

# Load the predictor rasters
r_principal_full <- rast(pred_files)

# Select the initial layer subset (the 16 candidate variables)
# -----------------------------------------------------------------------
# IMPORTANT FIX: numeric indices (1, 5, 6, 7...) only point to the correct variable WITHIN the Caiman_latirostris folder, because list.files() sorts files alphabetically and the order can differ between species/population folders (e.g., if a file is missing, or a name is slightly different)
initial_indices <- c(1, 5, 6, 7, 8, 14, 15, 16, 17, 18, 21, 22, 24, 25, 26, 27)
initial_var_names <- names(r_principal_full)[initial_indices]

cat("Variáveis candidatas iniciais (capturadas por nome a partir de",
    sp_principal, "):\n")
print(initial_var_names)

r_principal <- subset(r_principal_full, initial_var_names)

# Sample background points for analysis
set.seed(2023) # Garante reprodutibilidade
bg_principal <- spatSample(r_principal, 1000, "random", na.rm = TRUE, as.df = TRUE)

# MULTICOLLINEARITY ANALYSIS (VIF)
# Run VIF analysis with a threshold of 10 (removes variables with VIF > 10)
v <- vifstep(bg_principal, th = 10)
# Get the names of the selected variables
variaveis_selecionadas <- v@results$Variables

cat("Variáveis selecionadas após VIF:\n")
print(variaveis_selecionadas)

# SAVE ANALYSIS RESULTS
# Create the output directory if it does not exist
# Fix: added recursive = TRUE (other calls in the script already used it)
dir.create("./output/vif/", recursive = TRUE, showWarnings = FALSE)
# Save the correlation matrix
write.csv(v@corMatrix, paste0("./output/vif/vif_cormatrix_", sp_principal, ".csv"), fileEncoding = "UTF-8")
# Save the VIF results
write.csv(v@results, paste0("./output/vif/vif_results_", sp_principal, ".csv"), fileEncoding = "UTF-8")

# 4. APPLY THE SAME SELECTED VARIABLES TO ALL SPECIES/POPULATIONS
for (i in 1:length(species_names)) {
  especie_atual <- species_names[i]

  # Load predictors for the current species/population
  pred_folder <- paste0("./output/predictors/predmod/predmodcrop/", especie_atual, "/")
  pred_files <- list.files(path = pred_folder, pattern = '*.asc$', full.names = TRUE)

  r_full <- rast(pred_files)

  # Fix: selects the same initial subset by NAME (not by position), ensuring the correct variable is picked even if the alphabetical file order in this folder differs from the main species folder
  faltando <- setdiff(initial_var_names, names(r_full))
  if (length(faltando) > 0) {
    warning("Espécie/população ", especie_atual, " não tem as camadas: ",
            paste(faltando, collapse = ", "), " - pulando esta espécie.")
    next
  }
  r <- subset(r_full, initial_var_names)

  # Select only the variables that passed VIF for the main species
  p <- subset(r, variaveis_selecionadas)

  # Fix: ensures there are no duplicate layer names before writing
  # (two layers with the same name would overwrite each other inconsistently
  # in the per-layer writeRaster below)
  if (any(duplicated(names(p)))) {
    warning("Camadas com nomes duplicados em ", especie_atual,
            " - verifique initial_var_names / variaveis_selecionadas. Pulando esta espécie.")
    next
  }

  # PREPARE OUTPUT
  # Create a dedicated directory for the species/population
  out.dir <- paste0("./output/vif/pred_vif/", especie_atual, "/")
  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

  # Save the rasters with the selected variables
  # Fix: writes layer by layer inside a tryCatch, so a problematic layer (e.g., an all-NA/NaN raster, which can happen with "asp" on very flat terrain) produces a clear warning instead of halting the whole loop before the remaining species/populations are processed
  for (j in 1:nlyr(p)) {
    layer_name <- names(p)[j]
    n_validas <- global(p[[j]], "notNA")[1, 1]
    if (n_validas == 0) {
      warning("Camada '", layer_name, "' de ", especie_atual,
              " está totalmente NA - não é possível salvar como .asc. Pulando esta camada.")
      next
    }
    tryCatch({
      writeRaster(p[[j]], paste0(out.dir, layer_name, ".asc"), overwrite = TRUE, NAflag = -9999)
    }, error = function(e) {
      message("Falhou ao salvar '", layer_name, "' para ", especie_atual, ": ",
              conditionMessage(e))
    })
  }

  cat("Concluído:", especie_atual, "\n")
}

cat("Processamento de seleção de variáveis (VIF) concluído para todas as espécies/populações.\n")
