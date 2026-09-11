# Ecological niche modelling - GAM (algorithm complementary to MaxEnt and RF)
# by Fernando Paulino and Carolina Loss
# SEP - 2025
# Updated - SEP/2025 (generalised to run on the 7 evolutionary units/
# aggregates separately, Caiman_flumi, Caiman_gac, Caiman_gpi, Caiman_nocaa,
# Caiman_norma, Caiman_paran, Caiman_saof, each with its own set
# maps, permutation importance and response curves, saved to
# Results New/<unit>/. Caiman_latirostris is excluded from this script.)

library(terra)
library(mgcv)
library(dismo)
library(enmSdmX)
library(ggplot2)

set.seed(2025)

base_dir <- "./"  # ajuste aqui se os dados estiverem em outra pasta
models_dir <- paste0(base_dir, "output/models/")

all_items <- list.files(models_dir, full.names = TRUE)
species_list <- basename(all_items[file.info(all_items)$isdir])
species_list <- setdiff(species_list, "Caiman_latirostris")

cat("Unidades a processar no GAM:", length(species_list), "\n")
print(species_list)

grid.n <- 25
n.runs <- 20
k.folds <- 5

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
# OUTER LOOP: one evolutionary unit/aggregate at a time
# ---------------------------------------------------------------------------
for (species in species_list) {

  tryCatch({

    cat("\n========== Iniciando GAM para", species, "==========\n")

    xy.path <- paste0(base_dir, "output/xySWD/xySWD_", species, ".csv")
    env.path <- paste0(base_dir, "output/biasSWD/biasSWD_", species, ".csv")

    if (!file.exists(xy.path) || !file.exists(env.path)) {
      warning("Arquivos SWD não encontrados para ", species, " - pulando.")
      next
    }

    pts <- read.csv(xy.path)
    pts <- pts[, 6:ncol(pts)]
    env <- read.csv(env.path)
    env <- env[, 6:ncol(env)]

    colnames(pts) <- gsub("current", "", colnames(pts))
    colnames(env) <- gsub("current", "", colnames(env))
    p.names <- colnames(pts)

    species_folder <- paste0(base_dir, "output/predictors/predmod/predmodcrop/", species, "/")
    pred.files <- list.files(species_folder, pattern = '*.asc$', full.names = TRUE)
    if (length(pred.files) == 0) {
      warning("Nenhum preditor encontrado para ", species, " - pulando.")
      next
    }
    pred.all <- rast(pred.files)
    pred <- subset(pred.all, p.names)

    out.base <- paste0(base_dir, "output/models_gam/", species, "/")
    dir.create(out.base, recursive = TRUE, showWarnings = FALSE)

    eval.dir <- paste0(base_dir, "output/evaluation_gam/")
    dir.create(eval.dir, recursive = TRUE, showWarnings = FALSE)

    results_dir <- paste0(base_dir, "Results New/", species, "/")
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

    # GAM formula built dynamically for this unit's predictors
    smooth.terms <- paste0("s(", p.names, ", k = 3, bs = 'cr')", collapse = " + ")
    gam.formula <- as.formula(paste("presence ~", smooth.terms))

    eval_df <- data.frame(species = character(), run = integer(), fold = integer(),
                           Bin.Prob = numeric(), wAUC = numeric(),
                           TSS = numeric(), CBI = numeric())
    importance_df <- data.frame(run = integer(), fold = integer(),
                                 variavel = character(), importancia = numeric())
    curves_df <- data.frame(run = integer(), fold = integer(),
                             variavel = character(), x = numeric(), y = numeric())

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
        train.data$weight <- ifelse(train.data$presence == 1, 1, n.pres / n.bg)

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

        print(list(species = species, run = run, fold = fold, wAUC = wauc, TSS = tss, CBI = cbi))
      }

      run.stack <- rast(fold.preds)
      run.mean <- mean(run.stack)
      writeRaster(run.mean, paste0(run.dir, species, "_prediction.tif"), overwrite = TRUE)

      write.csv(importance_df, paste0(eval.dir, "gam_importance_progresso_", species, ".csv"), row.names = FALSE)
      write.csv(curves_df, paste0(eval.dir, "gam_curves_progresso_", species, ".csv"), row.names = FALSE)

      cat("GAM concluído -", species, "- run", run, "de", n.runs, "\n")
    }

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
      print(metric_average)
    } else {
      message("Nenhuma réplica válida (GAM) encontrada para ", species, ".")
    }

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

      cat("Importância e curvas de resposta (GAM) salvas para", species, "\n")
    } else {
      message("Nenhuma réplica válida suficiente para importância/curvas de resposta (GAM) - ", species)
    }

  }, error = function(e) {
    message("Falhou o processamento GAM de ", species, ": ", conditionMessage(e))
  })
}

cat("\nProcessamento GAM completo para as unidades evolutivas/agregados (excluindo Caiman_latirostris).\n")
