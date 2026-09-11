# Ecological niche modelling
# by Fernando Paulino and Carolina Loss
# JAN - 2025
# Fixed - SEP/2025 (critical bug: the loop variable was named "rast", which
# shadows the terra::rast() function itself after the first iteration -
# rast(raster) would then try to call a SpatRaster object as a function and
# fail on the second basin/raster combination; renamed to r_layer throughout)

library(terra)

base_dir <- "./"

path_rasters <- paste0(base_dir, "output/predictors/Vars/")
path_shapes <- paste0(base_dir, "data/shapes/")
path_output <- paste0(base_dir, "output/predictors/")

bacias <- list.dirs(path_shapes, recursive = FALSE)
raster_files <- list.files(path_rasters, pattern = "\\.tif$", full.names = TRUE)

if (length(bacias) == 0) stop("No basin subfolders found in ", path_shapes)
if (length(raster_files) == 0) stop("No .tif files found in ", path_rasters)

for (bacia in bacias) {

  tryCatch({

    nome_bacia <- basename(bacia)
    output_bacia <- file.path(path_output, paste0(nome_bacia, "_rast"))
    dir.create(output_bacia, recursive = TRUE, showWarnings = FALSE)

    shp_files <- list.files(bacia, pattern = "\\.shp$", full.names = TRUE)
    if (length(shp_files) == 0) {
      warning("No .shp file found in ", bacia, " - skipping this basin.")
      next
    }
    shape <- vect(shp_files[1])

    for (raster_path in raster_files) {

      tryCatch({

        nome_raster <- basename(raster_path)
        r_layer <- rast(raster_path)

        if (!same.crs(r_layer, shape)) {
          shape <- project(shape, crs(r_layer))
        }

        r_crop <- crop(r_layer, shape)
        r_mask <- mask(r_crop, shape)

        writeRaster(r_mask, file.path(output_bacia, nome_raster), overwrite = TRUE, NAflag = -9999)

      }, error = function(e) {
        message("Failed on raster ", basename(raster_path), " for basin ", nome_bacia,
                ": ", conditionMessage(e))
      })
    }

    cat("Finished basin:", nome_bacia, "\n")

  }, error = function(e) {
    message("Failed to process basin ", basename(bacia), ": ", conditionMessage(e))
  })
}

cat("Done. Rasters cropped and organised by basin.\n")
