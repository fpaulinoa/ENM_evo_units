# Ecological niche modelling
# by Fernando Paulino and Carolina Loss
# JAN - 2025
# Fixed - SEP/2025 (Java memory sized for the available RAM, dir.create
# silent, per-species/run tryCatch, resume capability that skips what has already
# was processed, predictor selection by name, JAVA_HOME auto-detected)

# ---------------------------------------------------------------------------
# NOTE: the two lines below (java.parameters and JAVA_HOME) only take effect
# if this is a NEW R session, with no Java package loaded yet.
# If MNE08 has already run in this same session, close R completely and reopen
# again before running this script, or the memory setting below
# will be silently ignored.
# ---------------------------------------------------------------------------

# Fix: -Xmx8g reserved 8GB of heap for Java alone, but the machine has 8,192 MiB (8GB) of TOTAL RAM
options(java.parameters = "-Xmx3g")

# Fix: same JAVA_HOME auto-detection used in MNE08, so as not to
# depending on prior manual configuration
if (Sys.getenv("JAVA_HOME") == "") {
  java_home_detectado <- tryCatch(
    system("/usr/libexec/java_home", intern = TRUE),
    error = function(e) ""
  )
  if (length(java_home_detectado) > 0 && java_home_detectado != "") {
    Sys.setenv(JAVA_HOME = java_home_detectado)
    cat("JAVA_HOME definido automaticamente para:", java_home_detectado, "\n")
  }
}

# Load the required packages
library(dismo)
library(raster)

base_dir <- "./"  # ajuste aqui se os dados estiverem em outra pasta

# Read the CSV file with occurrences for all species
spp.xy <- read.csv(paste0(base_dir, "data/occour/Caimans.csv"))
spp.xy$sp <- gsub(" ", "_", spp.xy$sp)

# Build a list of unique species names
spp <- unique(spp.xy$sp)

# Set the output folder names for the 20 modelling replicates
run.dir <- paste0("run_", 1:20)

# Metadata columns that are NOT environmental variables (same principle as
# MNE08, selects by name, not by position)
metadata_cols_occ <- c("X", "species", "cell", "x", "y")
metadata_cols_bg  <- c("X", "bg", "cell", "x", "y")

# Start a loop over each species
for (i in 1:length(spp)) {

  sp.name <- spp[i]

  # Fix: checks that the 3 required files exist BEFORE entering
  # in the 20-run loop; avoids starting to process a species and only failing
  # on the first run due to a missing file
  xy.path <- paste0(base_dir, "output/xySWD/xySWD_", sp.name, ".csv")
  env.path <- paste0(base_dir, "output/biasSWD/biasSWD_", sp.name, ".csv")
  f.path <- paste0(base_dir, "output/Fclass_best/ENMeval_best_", sp.name, ".csv")

  if (!file.exists(xy.path) || !file.exists(env.path) || !file.exists(f.path)) {
    message("Arquivos necessários não encontrados para ", sp.name,
            " (xySWD/biasSWD/Fclass_best) - pulando esta espécie.")
    next
  }

  # Load the CSV with the best parameters selected via ENMeval
  f <- read.csv(f.path)
  # Fix: always uses the first row, protecting against the rare case of
  # an exact AICc tie in MNE08, which would leave f with more than one row and
  # would break the "if (fc == 'L')" check below with "condition has length > 1"
  fc <- f$fc[1]
  rm_value <- f$rm[1]

  # Set the MaxEnt feature-class arguments based on the ENMeval result
  if (fc == 'L') {
    fclass <- c("linear=true", "quadratic=false", "hinge=false", "product=false", "threshold=false")
  } else if (fc == 'LQ') {
    fclass <- c("linear=true", "quadratic=true", "hinge=false", "product=false", "threshold=false")
  } else {
    fclass <- c("linear=true", "quadratic=true", "hinge=true", "product=false", "threshold=false")
  }

  rm_arg <- paste0("betamultiplier=", rm_value)

  # Set the path to the species' cropped-predictor folder
  species_folder <- paste0(base_dir, "output/predictors/predmod/predmodcrop/", sp.name, "/")
  pred.files <- list.files(species_folder, pattern = '*.asc$', full.names = TRUE)

  if (length(pred.files) == 0) {
    message("Nenhum preditor encontrado em ", species_folder, " - pulando ", sp.name)
    next
  }
  pred.all <- stack(pred.files)

  # Load the presence-point SWD file
  pts <- read.csv(xy.path)
  pred_cols_occ <- setdiff(colnames(pts), metadata_cols_occ)
  pts <- pts[, pred_cols_occ, drop = FALSE]

  # Load the pseudo-absence/background SWD file
  env <- read.csv(env.path)
  pred_cols_bg <- setdiff(colnames(env), metadata_cols_bg)
  env <- env[, pred_cols_bg, drop = FALSE]

  # Combine presence and pseudo-absence data into a single data frame
  swd <- rbind(pts, env)
  swd.v <- c(rep(1, nrow(pts)), rep(0, nrow(env)))

  # Remove the word "current" from column names, for standardisation
  colnames(swd) <- gsub("current", "", colnames(swd))

  # Select the predictors matching the SWD columns
  p.names <- colnames(swd)
  pred <- subset(pred.all, p.names)

  # Start a loop over each replicate (run_1 to run_20)
  for (j in 1:length(run.dir)) {

    dir.name <- paste0(base_dir, "output/Models Maxent/", sp.name, "/", run.dir[j], "/")
    predict_filename <- paste0(dir.name, sp.name, "_prediction.tif")

    # Fix: resume capability
    if (file.exists(predict_filename)) {
      cat("Já existe, pulando:", sp.name, run.dir[j], "\n")
      next
    }

    tryCatch({

      dir.create(dir.name, recursive = TRUE, showWarnings = FALSE)

      # Set all arguments for MaxEnt modelling
      args <- c('randomseed=true',
                'outputformat=logistic',
                'replicatetype=crossvalidate',
                'replicates=5',
                fclass,
                rm_arg,
                'writebackgroundpredictions=true')

      # Fit the MaxEnt model using presences and pseudo-absences
      model <- maxent(x = swd, p = swd.v, path = dir.name, args = args)

      # Spatially project the fitted model
      predict(object = model, x = pred, filename = predict_filename,
              na.rm = TRUE, format = 'GTiff', overwrite = TRUE, progress = 'text')

      cat("Concluído:", sp.name, run.dir[j], "\n")

    }, error = function(e) {
      message("Falhou em ", sp.name, " ", run.dir[j], ": ", conditionMessage(e))
    })
  }

  print(paste0("Modelagem concluída para a espécie: ", sp.name))
}

print("Processamento completo para todas as espécies.")
