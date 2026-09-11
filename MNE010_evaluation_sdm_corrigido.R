# Ecological niche modelling
# by Fernando Paulino and Carolina Loss
# JAN - 2025
# Fixed - SEP/2025 (removed the xlsx/Java dependency, filter for only
# real directories, per-replicate tryCatch, full paths, incremental saving
# incremental saving so progress is not lost if something fails midway)

library(terra)
library(enmSdmX)

base_dir <- "./"  # ajuste aqui se os dados estiverem em outra pasta

### Set the path to the generated models
output_path <- paste0(base_dir, "output/Models Maxent/")

# Fix: explicitly filters only directories (avoids breaking if there is a
# algum arquivo solto na pasta, como .DS_Store no Mac)
all_items <- list.files(output_path, full.names = TRUE)
species_list <- basename(all_items[file.info(all_items)$isdir])

cat("Espécies/populações encontradas:", length(species_list), "\n")
print(species_list)

# Create the directory to store evaluations, if it does not exist
evaluation_path <- paste0(base_dir, "Results New/evaluation/maxent")
dir.create(evaluation_path, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(evaluation_path)) {
  stop("Não foi possível criar ", evaluation_path,
       " - verifique permissões ou se já existe um arquivo com esse nome no lugar da pasta.")
}

# Initialise an empty data frame to store evaluation metrics
eval_df <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(eval_df) <- c("species", "Replica", "Bin.Prob", "AUC", "wAUC", "TSS", "CBI")

# Counters for the final report
n_ok <- 0
n_falhas <- 0

# Main loop: iterates over each species
for (species in species_list) {
  cat("Processando espécie:", species, "\n")
  species_path <- paste0(output_path, species, "/")

  # Loop over each run (run_1 to run_20)
  for (run_num in 1:20) {
    run_folder <- paste0(species_path, "run_", run_num)

    # Verifica se a pasta da rodada existe
    if (!dir.exists(run_folder)) {
      cat("Aviso: Subpasta", run_folder, "não encontrada. Pulando...\n")
      next
    }

    # Carrega o arquivo de resultados do Maxent
    maxent_results_file <- paste0(run_folder, "/maxentResults.csv")
    if (!file.exists(maxent_results_file)) {
      cat("Aviso: maxentResults.csv não encontrado em", run_folder, "pulando...\n")
      next
    }
    maxres <- read.csv(maxent_results_file)

    # List the presence-prediction files for each replicate
    replicates <- list.files(run_folder, pattern = "_samplePredictions.csv")

    for (replicate in replicates) {
      replicate_name <- sub("_samplePredictions.csv", "", replicate)

      # Fix: each replicate is processed inside a tryCatch; a failure
      # isolated (incomplete file, missing column, corrupted replicate)
      # does not bring down the other ~800 iterations still to be processed
      tryCatch({

        bg_file <- paste0(run_folder, "/", replicate_name, "_backgroundPredictions.csv")
        pres_file <- paste0(run_folder, "/", replicate_name, "_samplePredictions.csv")

        if (!file.exists(bg_file)) {
          warning("Arquivo de background não encontrado para ", species, "/", run_folder,
                  "/", replicate_name, " - pulando esta réplica.")
          n_falhas <<- n_falhas + 1
          next  # next funciona corretamente aqui, mesmo dentro do tryCatch (ver nota no MNE06/07)
        }

        # Load the presence and background prediction data
        presence <- read.csv(pres_file)
        background <- read.csv(bg_file)

        if (!"Logistic.prediction" %in% colnames(presence) || !"Logistic" %in% colnames(background)) {
          warning("Colunas esperadas não encontradas para ", species, "/", run_folder,
                  "/", replicate_name, " - pulando esta réplica.")
          n_falhas <<- n_falhas + 1
          next  # next funciona corretamente aqui, mesmo dentro do tryCatch (ver nota no MNE06/07)
        }

        # Extract the predictions
        predPres <- presence$Logistic.prediction
        predBg <- background$Logistic

        # Extrai o valor de AUC e o limiar (threshold) dos resultados do Maxent
        auc <- maxres[maxres[, 1] == replicate_name, "Test.AUC"]
        threshold <- maxres[maxres[, 1] == replicate_name, "Maximum.test.sensitivity.plus.specificity.Logistic.threshold"]

        if (length(auc) == 0 || length(threshold) == 0) {
          warning("Réplica '", replicate_name, "' não encontrada em maxentResults.csv de ",
                  run_folder, " - pulando.")
          n_falhas <<- n_falhas + 1
          next  # next funciona corretamente aqui, mesmo dentro do tryCatch (ver nota no MNE06/07)
        }

        # Compute evaluation metrics: TSS, weighted AUC and CBI
        tss <- enmSdmX::evalTSS(pres = predPres, contrast = predBg, thresholds = threshold)
        auc_eval <- enmSdmX::evalAUC(pres = predPres, contrast = predBg)
        cbi <- enmSdmX::evalContBoyce(pres = predPres, contrast = predBg)

        # Build a row with the metrics and append it to the overall data frame
        eval_row <- data.frame(species, replicate_name, threshold, auc, auc_eval, tss, cbi)
        colnames(eval_row) <- c("species", "Replica", "Bin.Prob", "AUC", "wAUC", "TSS", "CBI")
        eval_df <<- rbind(eval_df, eval_row)
        n_ok <<- n_ok + 1

      }, error = function(e) {
        message("Falhou ao avaliar ", species, "/", run_folder, "/", replicate_name,
                ": ", conditionMessage(e))
        n_falhas <<- n_falhas + 1
      })
    }

    # Fix: saves progress after every processed run (not only at the end),
    # so work already done across ~800 iterations is not lost if something
    # interrompa o script no meio do caminho
    write.csv(eval_df, paste0(evaluation_path, "/eval_df_progresso.csv"), row.names = FALSE)
  }
}

cat("Avaliação concluída:", n_ok, "réplicas processadas com sucesso,", n_falhas, "falharam/foram puladas.\n")

# Print all evaluations performed
print("Avaliações completas:")
print(eval_df)

# Filter valid replicates according to the established criteria
eval_df$decision <- ifelse(eval_df$wAUC >= 0.5 & eval_df$TSS >= 0.0 & eval_df$CBI >= 0.4, "include", "remove")
valid_replicates <- subset(eval_df, decision == "include")

# Print the replicates considered valid
print("Replicados válidos:")
print(valid_replicates)

# If there are valid replicates, compute the per-species metric averages
if (nrow(valid_replicates) > 0) {
  metric_average <- aggregate(valid_replicates[, 3:7], by = list(valid_replicates$species), FUN = mean)
  colnames(metric_average) <- c("species", "Bin.Prob", "AUC", "wAUC", "TSS", "CBI")

  # Fix: removed the dependency on the xlsx package (which uses rJava); CSV only,
  # which opens normally in Excel/Numbers/Google Sheets without needing Java
  write.csv(metric_average, file = paste0(evaluation_path, "/validation_average.csv"), row.names = FALSE)
  write.csv(eval_df, file = paste0(evaluation_path, "/eval_df_completo.csv"), row.names = FALSE)

  # Print the mean metrics
  print("Métricas médias por espécie:")
  print(metric_average)
} else {
  message("Nenhum replicado válido encontrado. Verifique os critérios ou as métricas calculadas.")
}

cat("Processamento do MNE010 completo.\n")
