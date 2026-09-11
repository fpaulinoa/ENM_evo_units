# Ecological niche modelling
# by Fernando Paulino and Carolina Loss
# JAN - 2025
# Fixed - SEP/2025 (critical cell/xy/bias-value misalignment bug,
# .asc file pattern instead of .tiff, dir.create on the output folders,
# per-species tryCatch)

# Load the terra package, used for spatial data handling
library(terra)

base_dir <- "./"  # ajuste aqui se os dados estiverem em outra pasta

# List the directory names inside the VIF-filtered predictors folder
# (renamed from "names" to "species_names" so it does not shadow the names() function)
species_names <- list.files(paste0(base_dir, "output/vif/pred_vif/"))

# Fix: ensures the 3 output folders exist before the loop, with
# explicit check (same pattern used in MNE04/MNE06)
out_dirs <- c(paste0(base_dir, "output/xySWD/"),
              paste0(base_dir, "output/envSWD/"),
              paste0(base_dir, "output/biasSWD/"))
for (d in out_dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(d)) {
    stop("Não foi possível criar ", d,
         " - verifique permissões ou se já existe um arquivo com esse nome no lugar da pasta.")
  }
}

# Loop over each set of variables
for (i in 1:length(species_names)) {

  sp.name <- species_names[i]

  tryCatch({

    # Set the folder holding the current species' predictors
    pred.folder <- paste0(base_dir, "output/vif/pred_vif/", sp.name)

    # Fix: pattern was misspelled "patter" (only worked due to R's partial matching) and the pattern was '*.tiff$', but MNE05 writes VIF-selected predictors as .asc, not .tiff
    pred.files <- list.files(pred.folder, pattern = '*.asc$', full.names = TRUE)

    if (length(pred.files) == 0) {
      warning("Nenhum preditor .asc encontrado em ", pred.folder, " - pulando ", sp.name)
      next  # next funciona corretamente aqui, mesmo dentro do tryCatch (ver nota no MNE06)
    }

    # Load the predictor rasters
    pred <- rast(pred.files)

    # Load the corresponding bias raster
    bias.path <- paste0(base_dir, "output/bias/", sp.name, "_bias.tiff")
    if (!file.exists(bias.path)) {
      warning("Raster de viés não encontrado para ", sp.name, ": ", bias.path, " - pulando.")
      next  # next funciona corretamente aqui, mesmo dentro do tryCatch (ver nota no MNE06)
    }
    bias <- rast(bias.path)

    # ---------------------------------------------------------------------
    # CRITICAL FIX: the original code extracted "cells(bias)" (ALL cells, including NA cells outside the study area) separately from "values(bias, na.rm = TRUE)" (only valid cells, with NA already removed)
    vals_all <- values(bias)[, 1]
    valid_idx <- which(!is.na(vals_all))

    cell_ids <- valid_idx                          # números de célula válidos
    xy <- xyFromCell(bias, valid_idx)              # coordenadas correspondentes
    prob <- vals_all[valid_idx]                    # valores de viés correspondentes

    KDEpts <- data.frame(cell = cell_ids, x = xy[, 1], y = xy[, 2], prob = prob)

    ## Load the occurrence records (already spatially thinned)
    occ.path <- paste0(base_dir, "output/pts_thin/xy_", sp.name, ".csv")
    if (!file.exists(occ.path)) {
      warning("Arquivo de ocorrências thinned não encontrado para ", sp.name,
              ": ", occ.path, " - pulando.")
      next  # next funciona corretamente aqui, mesmo dentro do tryCatch (ver nota no MNE06)
    }
    occ <- read.csv(occ.path)

    ## Extract the cell number and coordinates matching the occurrence points
    o <- extract(bias, occ[, 5:6], cells = TRUE, xy = TRUE)

    ## Extract predictor values for the occurrence points
    xySWD <- as.data.frame(extract(pred, o[, 3]))

    ## Add the species name to the data
    species <- rep(sp.name, nrow(o))

    ## Combine species name, cells, and predictor values into a single data frame
    xySWD <- cbind(species, o[, 3:5], xySWD)

    ## Save the file with the extracted variables for the occurrence points
    write.csv(xySWD, paste0(base_dir, "output/xySWD/xySWD_", sp.name, ".csv"))

    ## Remove occurrence points from the background-point set
    KDEpts <- subset(KDEpts, !is.element(KDEpts$cell, o$cell))

    ## Extract predictor values for all background points
    biasKDE_all <- as.data.frame(extract(pred, KDEpts[, 1]))

    ## Flag these points as "background"
    bg <- rep("background", nrow(KDEpts))

    ## Combine point type (background), coordinates and predictors
    biasKDE_all <- cbind(bg, KDEpts, biasKDE_all)

    ## Drop the probability column from the final environmental set
    envSWD <- biasKDE_all[, -5]

    ## Save the environmental data for the background points
    write.csv(envSWD, paste0(base_dir, "output/envSWD/envSWD_", sp.name, ".csv"))

    # If the number of background points exceeds 10,000:
    if (nrow(biasKDE_all) > 10000) {

      ## Sample 10,000 background points, weighted by bias probability
      biasSWD <- biasKDE_all[sample(seq(1:nrow(biasKDE_all)),
                                     size = 10000,
                                     replace = FALSE,
                                     prob = biasKDE_all[, "prob"]),
                              c(1:4, 6:ncol(biasKDE_all))]  # remove a coluna "prob"

      write.csv(biasSWD, paste0(base_dir, "output/biasSWD/biasSWD_", sp.name, ".csv"))

    } else {
      # If fewer than 10,000 background points are available
      message(sp.name, " tem menos de 10000 pontos de background (", nrow(biasKDE_all), ")")

      n <- round(nrow(biasKDE_all) * .8, digits = 0)

      biasSWD <- biasKDE_all[sample(seq(1:nrow(biasKDE_all)),
                                     size = n,
                                     replace = FALSE,
                                     prob = biasKDE_all[, "prob"]),
                              c(1:4, 6:ncol(biasKDE_all))]

      write.csv(biasSWD, paste0(base_dir, "output/biasSWD/biasSWD_", sp.name, ".csv"))
    }

    cat("Concluído:", sp.name, "\n")

  }, error = function(e) {
    message("Falhou ao processar '", sp.name, "': ", conditionMessage(e))
  })
}

cat("Processamento de background/xySWD/envSWD/biasSWD concluído.\n")
