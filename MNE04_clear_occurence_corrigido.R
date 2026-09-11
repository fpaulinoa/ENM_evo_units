# Ecological niche modelling
# by Fernando Paulino and Carolina Loss
# JAN - 2025
# Fixed - SEP/2025 (functional warning for species with few records,
# dir.create on output folders, diagnostics for unmatched names)

# Load the 'terra' package for spatial data handling
library(terra)

## Processing of NA-predictor files (created earlier)

# Set the path to the NA-mask folder
pred_folder <- paste0("./output/predNA/")

# List all NA-mask TIFF files
pred_path <- list.files(path = pred_folder, pattern = '*.tiff$')

# Build full file paths
pred_files <- paste(pred_folder, pred_path, sep = "")

# Extract species names from the file names:
# 1. Strip the file path
# 2. Strip the '_predNA.tiff' suffix
# (renamed from "names" to "species_names" so it does not shadow the names() function)
species_names <- gsub("./output/predNA/", "", pred_files)
species_names <- gsub("_predNA.tiff", "", species_names)

## Load the species' occurrence data
xy_all <- read.csv("./data/occour/Caimans.csv")
# Standardise scientific names, replacing spaces with underscores
xy_all$sp <- gsub(" ", "_", xy_all$sp)

## Fix: ensures the output folders exist before the loop
# (previously, the script assumed they already existed; if not, write.csv
# would fail with "cannot open the connection" with no clue why)
dir.create("./output/ptsNA/", recursive = TRUE, showWarnings = FALSE)
dir.create("./output/pts_thin/", recursive = TRUE, showWarnings = FALSE)

# Fix: checks that the folders actually exist after dir.create
if (!dir.exists("./output/ptsNA/")) {
  stop("Não foi possível criar ./output/ptsNA/ - verifique permissões ou se já existe um arquivo com esse nome no lugar da pasta.")
}
if (!dir.exists("./output/pts_thin/")) {
  stop("Não foi possível criar ./output/pts_thin/ - verifique permissões ou se já existe um arquivo com esse nome no lugar da pasta.")
}

## Fix: diagnostics for names that do not match between predNA and the CSV,
# run BEFORE the main loop, to warn upfront about any
# species being silently skipped due to a missing match
nomes_sem_predNA <- setdiff(unique(xy_all$sp), species_names)
nomes_sem_csv    <- setdiff(species_names, unique(xy_all$sp))

if (length(nomes_sem_predNA) > 0) {
  warning("Espécie(s) no CSV sem máscara predNA correspondente (não serão processadas aqui): ",
          paste(nomes_sem_predNA, collapse = ", "))
}
if (length(nomes_sem_csv) > 0) {
  warning("Máscara(s) predNA sem nenhuma correspondência no CSV de ocorrências: ",
          paste(nomes_sem_csv, collapse = ", "))
}

## Main loop to process each species
for (i in 1:length(pred_files)) {

  sp.name <- species_names[i]

  # Load the NA mask for the current species
  predNA <- rast(pred_files[i])

  # Filter occurrence points for the current species
  xy <- xy_all[xy_all$sp == sp.name, ]

  # Fix: explicitly warns if there are no occurrence points
  # matching that mask, instead of letting the rest of the block run
  # silently over an empty data frame
  if (nrow(xy) == 0) {
    warning("Nenhum registro de ocorrência casou com '", sp.name,
            "' - verifique grafia/espaços entre o nome do arquivo predNA e xy_all$sp. Pulando.")
    next
  }

  # Convert points to a spatial (vector) object with WGS84 CRS
  pts <- vect(xy, crs = "epsg:4326", geom = c("lon", "lat"))

  # Spatial processing:
  # 1. Extract NA-mask values for each point
  pts$mask <- extract(predNA, pts)[, 2]
  # 2. Identify the raster cell corresponding to each point
  pts$cell <- cells(predNA, pts)[, 2]
  # 3. Remove duplicate records within the same cell (thinning)
  pts <- pts[!duplicated(pts$cell), ]

  # Identify and save points outside the valid-predictor area (NA)
  pts_out <- pts[is.na(pts$mask), ]  # Pontos com NA na máscara
  ptsNA <- as.data.frame(pts_out)    # Converte para dataframe
  file_name <- paste0("./output/ptsNA/ptsNA_", sp.name, ".csv")
  write.csv(ptsNA, file_name, fileEncoding = "UTF-8")

  # Keep only points with valid predictor values
  pts <- pts[!is.na(pts$mask), ]

  # Prepare the final data for modelling:
  # 1. Extract coordinates
  coords <- crds(pts)
  # 2. Convert to data frame and join with coordinates
  pts_xy <- as.data.frame(pts)
  pts_xy <- cbind(pts_xy, coords)

  # Check the minimum number of records for modelling (>= 3)
  # Fix: the else branch now actually emits the message (previously it was just a
  # "dead" paste0() call, which built the string and did nothing with it; so
  # species with few records used to be skipped in complete silence)
  if (nrow(pts_xy) >= 3) {
    file_name <- paste0("./output/pts_thin/xy_", sp.name, ".csv")
    write.csv(pts_xy, file_name, fileEncoding = "UTF-8")
    cat("Concluído:", sp.name, "-", nrow(pts_xy), "registros após thinning\n")
  } else {
    warning(sp.name, " tem menos de 3 registros (", nrow(pts_xy),
            ") após o thinning espacial. Não será modelada.")
  }
}

cat("Processamento de limpeza de ocorrências concluído.\n")
