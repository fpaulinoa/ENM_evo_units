# Ecological niche modelling
# by Fernando Paulino and Carolina Loss
# JAN - 2025
# Fixed - SEP/2025 (critical bug: x/y being treated as an environmental variable in
# ENMeval; conditional package installation; dir.create; per-species tryCatch)

base_dir <- "./"  # ajuste aqui se os dados estiverem em outra pasta

# Fix: explicitly sets JAVA_HOME BEFORE anything else, instead of relying on it already being set in the session (this was the cause of the "object 'Library' not found" error that used to hang ENMevaluate)
if (Sys.getenv("JAVA_HOME") == "") {
  java_home_detectado <- tryCatch(
    system("/usr/libexec/java_home", intern = TRUE),
    error = function(e) ""
  )
  if (length(java_home_detectado) > 0 && java_home_detectado != "") {
    Sys.setenv(JAVA_HOME = java_home_detectado)
    cat("JAVA_HOME definido automaticamente para:", java_home_detectado, "\n")
  } else {
    warning("Não foi possível detectar JAVA_HOME automaticamente - ",
            "se o ENMevaluate falhar com 'object not found', defina manualmente com Sys.setenv(JAVA_HOME = '...').")
  }
}

# Install packages only if not already available, avoiding reinstallation on every run
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
if (!requireNamespace("ENMeval", quietly = TRUE)) {
  remotes::install_github("danlwarren/ENMeval")
}
if (!requireNamespace("rJava", quietly = TRUE)) {
  install.packages("rJava")
}

# Load the required packages
library(ENMeval)
library(terra)
library(dismo)
library(rJava)

# Locate the maxent.jar path, required to run MaxEnt from R
system.file("java", "maxent.jar", package = "dismo")

# Check the installed Java version
system("java -version")
.jinit()
Sys.getenv("JAVA_HOME")  # Mostra o caminho do Java configurado no sistema

# Note on parameter combinations that can be tuned in MaxEnt
# ('linear', 'quadratic', 'hinge', 'product', 'threshold', 'betamultiplier')

# List background file names, stripping prefixes and extensions
# to work with clean names
# (renamed from "names" to "species_names" so it does not shadow the names() function)
species_names <- list.files(paste0(base_dir, "output/biasSWD/"))
species_names <- gsub("biasSWD_", "", species_names)
species_names <- gsub(".csv", "", species_names)

# Fix: ensures the output folders exist before the loop
out_dirs <- c(paste0(base_dir, "output/Fclass_best/"),
              paste0(base_dir, "output/Fclass_results/"))
for (d in out_dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(d)) {
    stop("Não foi possível criar ", d,
         " - verifique permissões ou se já existe um arquivo com esse nome no lugar da pasta.")
  }
}

# Metadata columns that are NOT environmental variables; any column that
# not in this list is treated as a predictor. This replaces the earlier cut by
# position (xySWD[,4:ncol(xySWD)]), which was including "x" and "y" (the
# geographic coordinates) as if they were environmental variables.
metadata_cols_occ <- c("X", "species", "cell", "x", "y")
metadata_cols_bg  <- c("X", "bg", "cell", "x", "y")

# Loop to run ENMevaluate for each species or group
for (i in 1:length(species_names)) {

  sp.name <- species_names[i]

  tryCatch({

    # Read the occurrence SWD file
    pts.file <- paste0(base_dir, "output/xySWD/xySWD_", sp.name, ".csv")
    xySWD <- read.csv(pts.file)

    # Fix: selects predictors by NAME, excluding the
    # known metadata columns, instead of cutting by numeric position
    pred_cols_occ <- setdiff(colnames(xySWD), metadata_cols_occ)
    xySWD <- xySWD[, pred_cols_occ, drop = FALSE]

    # Read the background SWD file (pseudo-absences)
    bg.file <- paste0(base_dir, "output/biasSWD/biasSWD_", sp.name, ".csv")
    bgSWD <- read.csv(bg.file)

    pred_cols_bg <- setdiff(colnames(bgSWD), metadata_cols_bg)
    bgSWD <- bgSWD[, pred_cols_bg, drop = FALSE]

    # Fix: warns if, for any reason, the occurrence and background column sets do not match (e.g., a variable present in one file but not the other)
    if (!setequal(colnames(xySWD), colnames(bgSWD))) {
      warning("Colunas de preditores não batem entre ocorrência e background para ",
              sp.name, ". Ocorrência: ", paste(colnames(xySWD), collapse = ", "),
              " | Background: ", paste(colnames(bgSWD), collapse = ", "),
              " - pulando esta espécie.")
      next
    }

    # If there are at least 5 occurrence records
    if (nrow(xySWD) >= 5) {
      # Calibrate models using ENMevaluate with:
      # - different combinations of "feature classes" (L, LQ, LQH)
      # - regularisation values (rm = 1, 2, 3)
      # - using MaxEnt (maxent.jar)
      # - random partitioning into folds ("randomkfold")
      # Fix: parallel = FALSE (instead of TRUE/numCores=7)
      eval <- ENMevaluate(occs = xySWD, bg = bgSWD, parallel = FALSE,
                          tune.args = list(fc = c("L", "LQ", "LQH"), rm = 1:3),
                          algorithm = "maxent.jar", partitions = "randomkfold")

      # Additional notes:
      # - fc = types of variable transformations (linear, quadratic, hinge)
      # - rm = regularisation values (the higher, the smoother the model)
      # - Common, recommended parameters following Elith et al. 2011

      # Select the best model based on the lowest AICc (corrected information criterion)
      bestmod <- which(eval@results$AICc == min(eval@results$AICc, na.rm = TRUE))
      f <- eval@results[bestmod, ]

      # Save the results:
      # - best model's parameters
      write.csv(f, paste0(base_dir, "output/Fclass_best/ENMeval_best_", sp.name, ".csv"))

      # - all evaluation results
      write.csv(eval@results, paste0(base_dir, "output/Fclass_results/ENMeval_evalResults_", sp.name, ".csv"))

      cat("Concluído:", sp.name, "-", f$fc[1], "/ rm =", f$rm[1], "\n")

    } else {
      # If the species has fewer than 5 records
      message(sp.name, " tem menos de 5 registros de ocorrência - não será calibrada.")
    }

  }, error = function(e) {
    message("Falhou ao calibrar ENMeval para '", sp.name, "': ", conditionMessage(e))
  })
}

cat("Calibração ENMeval concluída para todas as espécies/populações.\n")
