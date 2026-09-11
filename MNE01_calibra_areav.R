# Ecological niche modelling
# by Fernando Paulino and Carolina Loss
# JAN - 2025
# Fixed - SEP/2025

library(terra)

base_dir <- "./"

bio_folder <- paste0(base_dir, "output/predictors/predmod/")
bio_files <- list.files(bio_folder, pattern = "\\.asc$", full.names = TRUE)

if (length(bio_files) == 0) {
  stop("No .asc predictor files found in ", bio_folder,
       " - run the earlier variable-assembly steps first.")
}

bio <- rast(bio_files)

# Calibration-area polygons, one per species/evolutionary unit, identified
# by the CLADO field in the shapefile attribute table.
b <- vect(paste0(base_dir, "data/shapes/calibrarea.shp"))
if (is.na(crs(b)) || crs(b) == "") {
  crs(b) <- "EPSG:4326"
}
plot(b)

for (i in 1:length(b)) {

  tryCatch({

    v <- b[i]
    id <- v$CLADO
    if (is.null(id) || is.na(id) || id == "") {
      warning("Polygon ", i, " has no valid CLADO identifier - skipping.")
      next
    }

    c <- crop(bio, v)
    m <- mask(c, v)

    out.dir <- paste0(base_dir, "output/predictors/predmod/predmodcrop/", id, "/")
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

    writeRaster(m, paste0(out.dir, names(m), ".asc"), overwrite = TRUE, NAflag = -9999)

    cat("Cropped predictors written for:", id, "\n")

  }, error = function(e) {
    message("Failed to process polygon ", i, ": ", conditionMessage(e))
  })
}

# Quick visual check on one of the cropped outputs (adjust the species name
# if Caiman_latirostris is not part of your calibration set).
check_path <- paste0(base_dir, "output/predictors/predmod/predmodcrop/Caiman_latirostris/bio1.asc")
if (file.exists(check_path)) {
  plot(rast(check_path))
} else {
  message("Visual check skipped: ", check_path, " not found.")
}
