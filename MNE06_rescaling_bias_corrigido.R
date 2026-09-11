# Ecological niche modelling
# by Fernando Paulino and Carolina Loss
# JAN - 2025
# Fixed - SEP/2025 (dir.create with a check, overwrite/NAflag in writeRaster,
# per-species robustness via tryCatch, spelling fixed in the plotting block)

# Load the 'terra' package for spatial data handling
library(terra)

# Define a function to rescale values to a new range
rescale <- function(x, x.min = NULL, x.max = NULL, new.min = 0, new.max = 1) {
  if (is.null(x.min)) x.min = min(x)  # Se mínimo não for fornecido, usa o mínimo dos dados
  if (is.null(x.max)) x.max = max(x)  # Se máximo não for fornecido, usa o máximo dos dados
  new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))  # Fórmula de reescalonamento
}

## Load the sampling-bias file generated in QGIS
# The bias file indicates areas with higher sampling probability
b <- rast("./data/bias.tif")

## Set the NA-mask path for all species
na_folder <- "./output/predNA/"
na_path <- list.files(path = na_folder, pattern = '*.tiff$')  # Lista arquivos TIFF
na_files <- paste(na_folder, na_path, sep = "")  # Cria caminhos completos

# Fix: the "[*] CRS do not match" warning seen at runtime indicates that bias.tif and the predNA masks are in DIFFERENT reference systems
maskNA_ref <- rast(na_files[1])
cat("CRS do bias.tif:\n"); print(crs(b, describe = TRUE))
cat("CRS das máscaras predNA:\n"); print(crs(maskNA_ref, describe = TRUE))

if (crs(b) != crs(maskNA_ref)) {
  cat("CRS diferentes detectados - reprojetando bias.tif para o CRS das máscaras...\n")
  b <- project(b, maskNA_ref)
}

# Extract species names from the file names
# (renamed from "names" to "species_names" so it does not shadow the names() function)
species_names <- gsub(na_folder, "", na_files)  # Remove caminho
species_names <- gsub("_predNA.tiff", "", species_names)  # Remove sufixo

# Fix: ensures the output folder exists before the loop, and stops immediately with a clear message if it cannot be created
out_dir <- "./output/bias/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) {
  stop("Não foi possível criar ", out_dir,
       " - verifique permissões ou se já existe um arquivo com esse nome no lugar da pasta.")
}

## Loop to process each species
for (i in 1:length(species_names)) {

  sp.name <- species_names[i]

  tryCatch({

    # Load NA mask for the current species
    maskNA <- rast(na_files[i])  # Máscara com 1 (dentro da área) e NA (fora)

    # Prepare the bias file:
    # 1. Resample the bias layer to match the mask's resolution/extent
    bNA <- resample(b, maskNA)
    # 2. Apply mask to restrict to the study area
    bNA <- bNA * maskNA

    # Create the inverse mask (0 inside the area, NA outside)
    mask0 <- maskNA - 1  # Transforma 1 em 0, mantém NA

    # Combine the inverse mask with bias:
    # - Where there is bias, keep the bias value
    # - Where there is no bias (NA), use 0 from the inverse mask
    bias <- sum(mask0, bNA, na.rm = TRUE)

    # Ensure the final bias layer is restricted to the study area
    bias <- bias * maskNA

    # Fix: if for any reason the bias layer ends up entirely NA
    # for that species/population, minmax() returns Inf/-Inf and rescale()
    # would silently produce invalid values; better to warn and skip.
    mm <- minmax(bias)
    if (!is.finite(mm[1]) || !is.finite(mm[2])) {
      warning("Camada de viés totalmente NA para ", sp.name, " - pulando esta espécie.")
      next  # next funciona normalmente aqui, mesmo dentro do tryCatch,
            # because it is a control-flow signal, not a condition
            # caught by tryCatch handlers (unlike a bare return()
            # outside a function, which would error since there is no function to return from)
    }

    # Rescale bias to the 1-100 range (for easier interpretation)
    b.scale <- rescale(x = bias,
                        x.min = mm[1],   # Mínimo original
                        x.max = mm[2],   # Máximo original
                        new.min = 1,     # Novo mínimo
                        new.max = 100)   # Novo máximo

    # Save the processed bias file
    # overwrite = TRUE avoids the "file already exists" error; NAflag = -9999 avoids invalid NODATA values
    writeRaster(b.scale, paste0(out_dir, sp.name, "_bias.tiff"),
                overwrite = TRUE, NAflag = -9999)

    cat("Concluído:", sp.name, "\n")

  }, error = function(e) {
    message("Falhou ao processar o viés de '", sp.name, "': ", conditionMessage(e))
  })
}

cat("Processamento de viés amostral concluído para todas as espécies/populações.\n")

## Visualisation of results (example for 3 species)
# Fix: "Caimamn_gpi" was the original typo in the CSV, already
# fixed at the source (see the MNE04/MNE05 fix); the correct name, used
# consistently in the predictor folders, is "Caiman_gpi".
r  <- rast(paste0(out_dir, "Caiman_gpi_bias.tiff"))          # Mapa de viés para Caiman_gpi
r2 <- rast(paste0(out_dir, "Caiman_flumi_bias.tiff"))        # Mapa de viés para Caiman_flumi
r3 <- rast(paste0(out_dir, "Caiman_latirostris_bias.tiff"))  # Mapa de viés para Caiman_latirostris

# Plot the bias maps
plot(r)   # Mapa de viés para Caiman_gpi
plot(r2)  # Mapa de viés para Caiman_flumi
plot(r3)  # Mapa de viés para Caiman_latirostris
