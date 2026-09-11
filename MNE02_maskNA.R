# Ecological niche modelling
# by Fernando Paulino and Carolina Loss
# JAN - 2025
# Fixed - SEP/2025 (species_names no longer shadows base R's names() function;
# added dir.create, overwrite/NAflag, and per-species tryCatch)

library(terra)

base_dir <- "./"

pts <- read.csv(paste0(base_dir, "data/occour/Caimans.csv"))
pts$sp <- gsub(" ", "_", pts$sp)
species_names <- unique(pts$sp)

out_dir <- paste0(base_dir, "output/predNA/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (i in 1:length(species_names)) {

  tryCatch({

    sp.name <- species_names[i]
    pred_folder <- paste0(base_dir, "output/predictors/predmod/predmodcrop/", sp.name, "/")
    pred_files <- list.files(pred_folder, pattern = "\\.asc$", full.names = TRUE)

    if (length(pred_files) == 0) {
      warning("No cropped predictors found for ", sp.name, " in ", pred_folder, " - skipping.")
      next
    }

    pred0 <- rast(pred_files)

    # Build a single NA mask: replace 0 with 100 to avoid division artefacts,
    # normalise each layer to 1 (so all valid cells become 1), then multiply
    # all layers together, so any NA cell in any layer propagates to the
    # final mask (na.rm = FALSE is intentional here).
    pred <- subst(pred0, 0, 100)
    pred <- pred / pred
    predNA <- prod(pred, na.rm = FALSE)

    plot(predNA, main = paste0("NA mask: ", sp.name))

    writeRaster(predNA, paste0(out_dir, sp.name, "_predNA.tiff"), overwrite = TRUE, NAflag = -9999)

    cat("NA mask written for:", sp.name, "\n")

  }, error = function(e) {
    message("Failed to process ", species_names[i], ": ", conditionMessage(e))
  })
}
