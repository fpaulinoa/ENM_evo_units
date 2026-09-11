# Ecological niche modelling
# by Fernando Paulino and Carolina Loss
# JAN - 2025
# Fixed - SEP/2025

library(dismo)  # provides biovars(); expects RasterStack input, not SpatRaster
library(raster)  # kept only for compatibility with dismo::biovars()
library(terra)

base_dir <- "./"  # adjust if data are stored elsewhere

tmin_dir <- paste0(base_dir, "data/climate/tmin/")
tmax_dir <- paste0(base_dir, "data/climate/tmax/")
precip_dir <- paste0(base_dir, "data/climate/pr/")

out_dir <- paste0(base_dir, "output/predictors/bioclim/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tmin_files <- sort(list.files(tmin_dir, pattern = "\\.tif$", full.names = TRUE))
tmax_files <- sort(list.files(tmax_dir, pattern = "\\.tif$", full.names = TRUE))
precip_files <- sort(list.files(precip_dir, pattern = "\\.tif$", full.names = TRUE))

if (length(tmin_files) == 0 || length(tmax_files) == 0 || length(precip_files) == 0) {
  stop("No .tif files found in one of the input folders. Check tmin_dir/tmax_dir/precip_dir.")
}

tmin_stack <- stack(tmin_files)
tmax_stack <- stack(tmax_files)
precip_stack <- stack(precip_files)

cat("Layers found - tmin:", nlayers(tmin_stack),
    "tmax:", nlayers(tmax_stack),
    "precip:", nlayers(precip_stack), "\n")

# If more than 12 layers are present (i.e., multiple years of monthly data),
# average across years to obtain a single climatological monthly series.
# Guards against a layer count that is not an exact multiple of 12, which
# would otherwise silently misalign months across years.
average_to_12_months <- function(s, label) {
  n <- nlayers(s)
  if (n == 12) return(s)
  if (n %% 12 != 0) {
    stop("Layer count for ", label, " (", n, ") is not a multiple of 12; ",
         "cannot safely average into a monthly climatology.")
  }
  stackApply(s, indices = rep(1:12, times = n / 12), fun = mean)
}

tmin_monthly <- average_to_12_months(tmin_stack, "tmin")
tmax_monthly <- average_to_12_months(tmax_stack, "tmax")
precip_monthly <- average_to_12_months(precip_stack, "precip")

bioclim_vars <- biovars(prec = precip_monthly, tmin = tmin_monthly, tmax = tmax_monthly)

for (i in 1:19) {
  writeRaster(bioclim_vars[[i]], paste0(out_dir, "bio", i, ".tif"), overwrite = TRUE)
}

cat("Bioclimatic variables (bio1-bio19) written to", out_dir, "\n")
