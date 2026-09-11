# Ecological niche modelling - GAM (algorithm complementary to MaxEnt and RF)
# by Fernando Paulino and Carolina Loss
# SEP - 2025
# Updated - SEP/2025 (base_dir paths, permutation importance and
# per-variable response curves, saved to Results New/Caiman_latirostris/)
#
# FOCUS: Article 1, Caiman_latirostris only (whole species). The 7 units
# units/aggregates are left for Article 2, in a separate script
# (MNE09c_gam_unidades.R), kept for when this project begins.
#
# Reuses the same SWD files used for MaxEnt (MNE09) and for Random
# Forest (MNE09b). Follows the same validation design: 20 runs x 5-fold
# cross-validation, to allow direct comparison across the three algorithms.

library(terra)
library(mgcv)
library(dismo)
library(enmSdmX)
library(ggplot2)

set.seed(2025)

species <- "Caiman_latirostris"
base_dir <- "./"  # ajuste aqui se os dados estiverem em outra pasta

# ---------------------------------------------------------------------------
# 1. LOAD ALREADY-PREPARED SWD DATA (same files as MaxEnt/RF)
# ---------------------------------------------------------------------------
xy.path <- paste0(base_dir, "output/xySWD/xySWD_", species, ".csv")
pts <- read.csv(xy.path)
pts <- pts[, 6:ncol(pts)]

env.path <- paste0(base_dir, "output/biasSWD/biasSWD_", species, ".csv")
env <- read.csv(env.path)
env <- env[, 6:ncol(env)]

colnames(pts) <- gsub("current", "", colnames(pts))
colnames(env) <- gsub("current", "", colnames(env))

p.names <- colnames(pts)

# ---------------------------------------------------------------------------
# 2. LOAD PREDICTOR RASTERS
# ---------------------------------------------------------------------------
species_folder <- paste0(base_dir, "output/predictors/predmod/predmodcrop/", species, "/")
pred.files <- list.files(species_folder, pattern = '*.asc$', full.names = TRUE)
pred.all <- rast(pred.files)
pred <- subset(pred.all, p.names)

# ---------------------------------------------------------------------------
# 3. OUTPUT FOLDERS
# ---------------------------------------------------------------------------
out.base <- paste0(base_dir, "output/models_gam/", species, "/")
dir.create(out.base, recursive = TRUE, showWarnings = FALSE)

eval.dir <- paste0(base_dir, "Results New/evaluation/gam/")
dir.create(eval.dir, recursive = TRUE, showWarnings = FALSE)

results_dir <- paste0(base_dir, "Results New/", species, "/")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 4. BUILD THE GAM FORMULA DYNAMICALLY FROM THE VIF-SELECTED PREDICTORS
# ---------------------------------------------------------------------------
smooth.terms <- paste0("s(", p.names, ", k = 3, bs = 'cr')", collapse = " + ")
gam.formula <- as.formula(paste("presence ~", smooth.terms))
print(gam.formula)

n.runs <- 20
k.folds <- 5
grid.n <- 25

eval_df <- data.frame(species = character(), run = integer(), fold = integer(),
                       Bin.Prob = numeric(), wAUC = numeric(),
                       TSS = numeric(), CBI = numeric())

# Accumulate the mean rasters from each valid run, to build the final map of
# GAM (mean + SD across runs) directly in this script
run_rasters_validos <- list()

importance_df <- data.frame(run = integer(), fold = integer(),
                             variavel = character(), importancia = numeric())
curves_df <- data.frame(run = integer(), fold = integer(),
                         variavel = character(), x = numeric(), y = numeric())

# ---------------------------------------------------------------------------
# Function: response curve for ONE variable, holding the others fixed at the
# the training-data median
# ---------------------------------------------------------------------------
response_curve_gam <- function(model, train.data, var.name, grid.n) {
  medianas <- sapply(train.data[, setdiff(colnames(train.data), c("presence", "weight"))],
                      median, na.rm = TRUE)
  grid.vals <- seq(min(train.data[[var.name]], na.rm = TRUE),
                    max(train.data[[var.name]], na.rm = TRUE),
                    length.out = grid.n)

  newdata <- as.data.frame(matrix(rep(medianas, grid.n), nrow = grid.n, byrow = TRUE))
  colnames(newdata) <- names(medianas)
  newdata[[var.name]] <- grid.vals

  preds <- predict(model, newdata = newdata, type = "response")
  data.frame(x = grid.vals, y = preds)
}

# ---------------------------------------------------------------------------
# Function: permutation importance, same logic as RF/MaxEnt
# ---------------------------------------------------------------------------
permutation_importance_gam <- function(model, test.pres, test.bg, baseline.wauc, var.names) {
  resultado <- data.frame(variavel = character(), importancia = numeric())
  for (v in var.names) {
    test.pres.perm <- test.pres
    test.bg.perm <- test.bg
    test.pres.perm[[v]] <- sample(test.pres.perm[[v]])
    test.bg.perm[[v]] <- sample(test.bg.perm[[v]])

    pred.pres.perm <- predict(model, newdata = test.pres.perm, type = "response")
    pred.bg.perm <- predict(model, newdata = test.bg.perm, type = "response")
    wauc.perm <- enmSdmX::evalAUC(pres = pred.pres.perm, contrast = pred.bg.perm)

    queda <- max(0, baseline.wauc - wauc.perm)
    resultado <- rbind(resultado, data.frame(variavel = v, importancia = queda))
  }
  if (sum(resultado$importancia) > 0) {
    resultado$importancia <- 100 * resultado$importancia / sum(resultado$importancia)
  }
  resultado
}

# ---------------------------------------------------------------------------
# 5. MAIN LOOP: 20 runs x 5-fold cross-validation
# ---------------------------------------------------------------------------
for (run in 1:n.runs) {

  run.dir <- paste0(out.base, "run_", run, "/")
  dir.create(run.dir, recursive = TRUE, showWarnings = FALSE)

  fold.pres <- dismo::kfold(pts, k = k.folds)
  fold.bg   <- dismo::kfold(env, k = k.folds)

  fold.preds <- vector("list", k.folds)

  for (fold in 1:k.folds) {

    train.idx.p <- fold.pres != fold
    train.idx.b <- fold.bg   != fold

    train.data <- rbind(pts[train.idx.p, , drop = FALSE],
                         env[train.idx.b, , drop = FALSE])
    train.data$presence <- c(rep(1, sum(train.idx.p)), rep(0, sum(train.idx.b)))

    test.pres <- pts[!train.idx.p, , drop = FALSE]
    test.bg   <- env[!train.idx.b, , drop = FALSE]

    n.pres <- sum(train.idx.p)
    n.bg   <- sum(train.idx.b)
    w.pres <- 1
    w.bg   <- n.pres / n.bg
    train.data$weight <- ifelse(train.data$presence == 1, w.pres, w.bg)

    gam.model <- gam(gam.formula, data = train.data, family = binomial(link = "logit"),
                      weights = train.data$weight, method = "REML")

    pred.pres <- predict(gam.model, newdata = test.pres, type = "response")
    pred.bg   <- predict(gam.model, newdata = test.bg,   type = "response")

    th.grid <- seq(0, 1, by = 0.01)
    sens.spec <- sapply(th.grid, function(th) {
      tpr <- mean(pred.pres >= th)
      tnr <- mean(pred.bg   <  th)
      tpr + tnr
    })
    best.th <- th.grid[which.max(sens.spec)]

    wauc <- enmSdmX::evalAUC(pres = pred.pres, contrast = pred.bg)
    tss  <- enmSdmX::evalTSS(pres = pred.pres, contrast = pred.bg, thresholds = best.th)
    cbi  <- enmSdmX::evalContBoyce(pres = pred.pres, contrast = pred.bg)

    eval_df <- rbind(eval_df, data.frame(species = species, run = run, fold = fold,
                                          Bin.Prob = best.th, wAUC = wauc,
                                          TSS = tss, CBI = cbi))

    if (wauc >= 0.5 && tss >= 0.0 && cbi >= 0.4) {

      imp <- permutation_importance_gam(gam.model, test.pres, test.bg, wauc, p.names)
      imp$run <- run
      imp$fold <- fold
      importance_df <- rbind(importance_df, imp[, c("run", "fold", "variavel", "importancia")])

      for (v in p.names) {
        curva <- response_curve_gam(gam.model, train.data, v, grid.n)
        curva$run <- run
        curva$fold <- fold
        curva$variavel <- v
        curves_df <- rbind(curves_df, curva[, c("run", "fold", "variavel", "x", "y")])
      }
    }

    fold.preds[[fold]] <- predict(pred, gam.model, type = "response", na.rm = TRUE)

    print(list(run = run, fold = fold, wAUC = wauc, TSS = tss, CBI = cbi))
  }

  run.stack <- rast(fold.preds)
  run.mean <- mean(run.stack)
  writeRaster(run.mean, paste0(run.dir, species, "_prediction.tif"), overwrite = TRUE)

  metricas_run <- eval_df[eval_df$run == run, ]
  run_valida <- with(metricas_run, mean(wAUC) >= 0.5 & mean(TSS) >= 0.0 & mean(CBI) >= 0.4)
  if (isTRUE(run_valida)) {
    run_rasters_validos[[length(run_rasters_validos) + 1]] <- run.mean
  }

  metricas_run <- eval_df[eval_df$run == run, ]
  run_valida <- with(metricas_run, mean(wAUC) >= 0.5 & mean(TSS) >= 0.0 & mean(CBI) >= 0.4)
  if (isTRUE(run_valida)) {
    run_rasters_validos[[length(run_rasters_validos) + 1]] <- run.mean
  }

  # Also save this run's specific subset inside its own folder
  # run_X/ folder, same logic used for RF
  write.csv(metricas_run, paste0(run.dir, "gam_evaluation_run_", run, ".csv"), row.names = FALSE)
  importance_run <- importance_df[importance_df$run == run, ]
  curves_run <- curves_df[curves_df$run == run, ]
  write.csv(importance_run, paste0(run.dir, "gam_importance_run_", run, ".csv"), row.names = FALSE)
  write.csv(curves_run, paste0(run.dir, "gam_curves_run_", run, ".csv"), row.names = FALSE)

  write.csv(importance_df, paste0(eval.dir, "gam_importance_progresso.csv"), row.names = FALSE)
  write.csv(curves_df, paste0(eval.dir, "gam_curves_progresso.csv"), row.names = FALSE)

  print(paste0("GAM concluído - run ", run, " de ", n.runs))
}

# ---------------------------------------------------------------------------
# 6. SAVE EVALUATIONS AND FILTER VALID REPLICATES (same criteria as MNE010)
# ---------------------------------------------------------------------------
write.csv(eval_df, paste0(eval.dir, "gam_evaluation_", species, ".csv"), row.names = FALSE)

eval_df$decision <- ifelse(eval_df$wAUC >= 0.5 & eval_df$TSS >= 0.0 & eval_df$CBI >= 0.4,
                            "include", "remove")
valid_replicates <- subset(eval_df, decision == "include")

if (nrow(valid_replicates) > 0) {
  metric_average <- data.frame(
    species  = species,
    Bin.Prob = mean(valid_replicates$Bin.Prob),
    wAUC     = mean(valid_replicates$wAUC),
    TSS      = mean(valid_replicates$TSS),
    CBI      = mean(valid_replicates$CBI)
  )
  write.csv(metric_average, paste0(eval.dir, "gam_validation_average_", species, ".csv"),
            row.names = FALSE)
  print("Métricas médias (GAM):")
  print(metric_average)
} else {
  message("Nenhuma réplica válida (GAM) encontrada para ", species, ". Verifique os critérios.")
}

# ---------------------------------------------------------------------------
# 7. IMPORTANCE SUMMARY AND RESPONSE CURVES
# ---------------------------------------------------------------------------
if (nrow(importance_df) > 0) {

  resumo_importancia <- aggregate(importancia ~ variavel, data = importance_df,
                                   FUN = function(x) c(media = mean(x), sd = sd(x)))
  resumo_importancia <- do.call(data.frame, resumo_importancia)
  colnames(resumo_importancia) <- c("variavel", "media", "sd")
  write.csv(resumo_importancia, paste0(results_dir, species, "_gam_importancia.csv"), row.names = FALSE)

  resumo_curvas <- aggregate(y ~ variavel + x, data = curves_df,
                              FUN = function(v) c(media = mean(v), sd = sd(v)))
  resumo_curvas <- do.call(data.frame, resumo_curvas)
  colnames(resumo_curvas) <- c("variavel", "x", "media", "sd")
  write.csv(resumo_curvas, paste0(results_dir, species, "_gam_curvas_resposta.csv"), row.names = FALSE)

  g_imp <- ggplot(resumo_importancia, aes(x = reorder(variavel, -media), y = media)) +
    geom_bar(stat = "identity", fill = "darkorange") +
    geom_errorbar(aes(ymin = pmax(0, media - sd), ymax = media + sd), width = 0.25) +
    labs(title = paste("Importância das Variáveis (GAM) -", species),
         x = "Variável", y = "Importância (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(results_dir, species, "_gam_importancia.png"), g_imp, width = 10, height = 6, dpi = 300)

  g_curvas <- ggplot(resumo_curvas, aes(x = x, y = media)) +
    geom_ribbon(aes(ymin = pmax(0, media - sd), ymax = pmin(1, media + sd)), fill = "grey70", alpha = 0.5) +
    geom_line(color = "black") +
    facet_wrap(~ variavel, scales = "free_x") +
    labs(title = paste("Curvas de Resposta (GAM) -", species),
         x = "Valor da variável", y = "Adequabilidade") +
    theme_minimal()
  ggsave(paste0(results_dir, species, "_gam_curvas_resposta.png"), g_curvas, width = 12, height = 8, dpi = 300)

  cat("Importância e curvas de resposta (GAM) salvas em:", results_dir, "\n")
} else {
  message("Nenhuma réplica válida suficiente para calcular importância/curvas de resposta (GAM).")
}

print("Processamento GAM completo para Caiman_latirostris.")

# ---------------------------------------------------------------------------
# 8. GAM MEAN AND SD MAP (across valid runs), with elevation cut-off
# ---------------------------------------------------------------------------
if (length(run_rasters_validos) > 0) {

  gam_all <- rast(run_rasters_validos)
  gam_mean <- mean(gam_all)
  gam_sd <- app(gam_all, sd)

  writeRaster(gam_mean, paste0(results_dir, species, "_gam_mean.tif"), overwrite = TRUE)
  writeRaster(gam_sd, paste0(results_dir, species, "_gam_sd.tif"), overwrite = TRUE)

  cat("Mapa médio/SD do GAM salvo em:", results_dir, "(", length(run_rasters_validos),
      "de", n.runs, "runs válidas)\n")

  elev_path <- paste0(base_dir, "output/vif/pred_vif/", species, "/elev.asc")
  if (file.exists(elev_path)) {
    elev_raster <- rast(elev_path)
    elev_mask <- elev_raster <= 800

    if (!compareGeom(elev_mask, gam_mean, stopOnError = FALSE)) {
      elev_mask <- resample(elev_mask, gam_mean, method = "near")
    }

    gam_mean_800m <- mask(gam_mean, elev_mask, maskvalues = c(0, NA))
    gam_sd_800m   <- mask(gam_sd, elev_mask, maskvalues = c(0, NA))

    writeRaster(gam_mean_800m, paste0(results_dir, species, "_gam_mean_800m.tif"), overwrite = TRUE)
    writeRaster(gam_sd_800m, paste0(results_dir, species, "_gam_sd_800m.tif"), overwrite = TRUE)

    cat("Corte por elevação (>800m) do GAM salvo em:", results_dir, "\n")
  } else {
    warning("Camada 'elev' não encontrada para ", species,
            " - mapas do GAM salvos SEM o corte de elevação. Caminho esperado: ", elev_path)
  }

} else {
  message("Nenhuma run válida (GAM) para gerar mapa médio/SD de ", species, ".")
}
