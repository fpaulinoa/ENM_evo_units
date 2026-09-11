# MaxEnt response curves - Caiman_latirostris
# by Fernando Paulino and Carolina Loss
# SEP - 2025
#
# Dedicated script, added later, only to extract response curves
# MaxEnt partial-dependence curves; the original MNE09 was not run with
# 'responsecurves=true', and this information matters for comparing
# directly with the curves already generated for RF (MNE09b) and GAM (MNE09c).
#
# Instead of relying on MaxEnt's internal response-curve output format
# curves, this script uses the SAME partial-dependence approach already
# applied to RF and GAM: varies one variable at a time (holding the
# others at the training-data median) and uses dismo::predict() with a
# a data.frame (not a raster) to obtain the prediction; this ensures the
# curves from the three algorithms are directly comparable, computed
# in exactly the same way.
#
# Restricted to Caiman_latirostris; reuses the fc/rm already chosen by
# ENMeval (MNE08), does not repeat the tuning. Does not rewrite the prediction rasters
# per run (they already exist under output/Models Maxent/); it only computes and saves the
# response curves and permutation importance (for symmetry
# methodological symmetry with RF/GAM).

options(java.parameters = "-Xmx3g")
if (Sys.getenv("JAVA_HOME") == "") {
  java_home_detectado <- tryCatch(system("/usr/libexec/java_home", intern = TRUE), error = function(e) "")
  if (length(java_home_detectado) > 0 && java_home_detectado != "") {
    Sys.setenv(JAVA_HOME = java_home_detectado)
    cat("JAVA_HOME definido automaticamente para:", java_home_detectado, "\n")
  }
}

library(dismo)
library(rJava)
library(enmSdmX)
library(ggplot2)

species <- "Caiman_latirostris"
base_dir <- "./"  # ajuste aqui se os dados estiverem em outra pasta

results_dir <- paste0(base_dir, "Results New/", species, "/")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Load the already-prepared data (same as MNE09/MNE09b/MNE09c)
# ---------------------------------------------------------------------------
xy.path <- paste0(base_dir, "output/xySWD/xySWD_", species, ".csv")
env.path <- paste0(base_dir, "output/biasSWD/biasSWD_", species, ".csv")
f.path <- paste0(base_dir, "output/Fclass_best/ENMeval_best_", species, ".csv")

pts <- read.csv(xy.path)
metadata_cols_occ <- c("X", "species", "cell", "x", "y")
pts <- pts[, setdiff(colnames(pts), metadata_cols_occ), drop = FALSE]

env <- read.csv(env.path)
metadata_cols_bg <- c("X", "bg", "cell", "x", "y")
env <- env[, setdiff(colnames(env), metadata_cols_bg), drop = FALSE]

colnames(pts) <- gsub("current", "", colnames(pts))
colnames(env) <- gsub("current", "", colnames(env))
p.names <- colnames(pts)

# Reuses the fc/rm already chosen by ENMeval; does not repeat the tuning
f <- read.csv(f.path)
fc <- f$fc[1]
rm_value <- f$rm[1]

if (fc == 'L') {
  fclass <- c("linear=true", "quadratic=false", "hinge=false", "product=false", "threshold=false")
} else if (fc == 'LQ') {
  fclass <- c("linear=true", "quadratic=true", "hinge=false", "product=false", "threshold=false")
} else {
  fclass <- c("linear=true", "quadratic=true", "hinge=true", "product=false", "threshold=false")
}
rm_arg <- paste0("betamultiplier=", rm_value)

cat("Usando fc =", fc, "/ rm =", rm_value, "(já escolhido pelo ENMeval no MNE08)\n")

# ---------------------------------------------------------------------------
# Function: response curve for ONE variable for MaxEnt, holding the
# others fixed at the median, same logic used for RF (MNE09b) and GAM
# (MNE09c), but using dismo::predict() with a data.frame instead of a raster
# ---------------------------------------------------------------------------
response_curve_maxent <- function(model, train.data, var.name, grid.n = 25) {
  medianas <- sapply(train.data, median, na.rm = TRUE)
  grid.vals <- seq(min(train.data[[var.name]], na.rm = TRUE),
                    max(train.data[[var.name]], na.rm = TRUE),
                    length.out = grid.n)

  newdata <- as.data.frame(matrix(rep(medianas, grid.n), nrow = grid.n, byrow = TRUE))
  colnames(newdata) <- names(medianas)
  newdata[[var.name]] <- grid.vals

  preds <- predict(model, x = newdata, args = "outputformat=logistic")
  data.frame(x = grid.vals, y = preds)
}

# Permutation importance, same logic used for RF/GAM, to keep the
# methodologically symmetric comparison across the three algorithms
permutation_importance_maxent <- function(model, test.pres, test.bg, baseline.wauc, var.names) {
  resultado <- data.frame(variavel = character(), importancia = numeric())
  for (v in var.names) {
    test.pres.perm <- test.pres
    test.bg.perm <- test.bg
    test.pres.perm[[v]] <- sample(test.pres.perm[[v]])
    test.bg.perm[[v]] <- sample(test.bg.perm[[v]])

    pred.pres.perm <- predict(model, x = test.pres.perm, args = "outputformat=logistic")
    pred.bg.perm <- predict(model, x = test.bg.perm, args = "outputformat=logistic")
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
# 2. Loop: 20 runs x 5-fold, only to compute curves/importance, without
#    rewrite the prediction rasters (they already exist)
# ---------------------------------------------------------------------------
n.runs <- 20
k.folds <- 5

curves_df <- data.frame(run = integer(), fold = integer(),
                         variavel = character(), x = numeric(), y = numeric())
importance_df <- data.frame(run = integer(), fold = integer(),
                             variavel = character(), importancia = numeric())

for (run in 1:n.runs) {

  fold.pres <- dismo::kfold(pts, k = k.folds)
  fold.bg   <- dismo::kfold(env, k = k.folds)

  for (fold in 1:k.folds) {

    tryCatch({

      train.idx.p <- fold.pres != fold
      train.idx.b <- fold.bg   != fold

      train.data <- rbind(pts[train.idx.p, , drop = FALSE],
                           env[train.idx.b, , drop = FALSE])
      train.resp <- c(rep(1, sum(train.idx.p)), rep(0, sum(train.idx.b)))

      test.pres <- pts[!train.idx.p, , drop = FALSE]
      test.bg   <- env[!train.idx.b, , drop = FALSE]

      args <- c('randomseed=true', 'outputformat=logistic', fclass, rm_arg)

      model <- maxent(x = train.data, p = train.resp, args = args)

      pred.pres <- predict(model, x = test.pres, args = "outputformat=logistic")
      pred.bg   <- predict(model, x = test.bg,   args = "outputformat=logistic")

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

      if (wauc >= 0.5 && tss >= 0.0 && cbi >= 0.4) {

        imp <- permutation_importance_maxent(model, test.pres, test.bg, wauc, p.names)
        imp$run <- run
        imp$fold <- fold
        importance_df <- rbind(importance_df, imp[, c("run", "fold", "variavel", "importancia")])

        for (v in p.names) {
          curva <- response_curve_maxent(model, train.data, v)
          curva$run <- run
          curva$fold <- fold
          curva$variavel <- v
          curves_df <- rbind(curves_df, curva[, c("run", "fold", "variavel", "x", "y")])
        }
      }

      cat("run", run, "fold", fold, "- wAUC:", round(wauc, 3), "TSS:", round(tss, 3), "CBI:", round(cbi, 3), "\n")

    }, error = function(e) {
      message("Falhou em run ", run, " fold ", fold, ": ", conditionMessage(e))
    })
  }

  # saves progress after every run
  write.csv(curves_df, paste0(results_dir, species, "_maxent_curves_progresso.csv"), row.names = FALSE)
  write.csv(importance_df, paste0(results_dir, species, "_maxent_importance_progresso.csv"), row.names = FALSE)
  cat("Curvas MaxEnt - run", run, "de", n.runs, "concluída\n")
}

# ---------------------------------------------------------------------------
# 3. Final summary (mean + SD) and plots, same pattern as RF/GAM
# ---------------------------------------------------------------------------
if (nrow(curves_df) > 0) {

  resumo_curvas <- aggregate(y ~ variavel + x, data = curves_df,
                              FUN = function(v) c(media = mean(v), sd = sd(v)))
  resumo_curvas <- do.call(data.frame, resumo_curvas)
  colnames(resumo_curvas) <- c("variavel", "x", "media", "sd")
  write.csv(resumo_curvas, paste0(results_dir, species, "_maxent_curvas_resposta.csv"), row.names = FALSE)

  g_curvas <- ggplot(resumo_curvas, aes(x = x, y = media)) +
    geom_ribbon(aes(ymin = pmax(0, media - sd), ymax = pmin(1, media + sd)), fill = "grey70", alpha = 0.5) +
    geom_line(color = "black") +
    facet_wrap(~ variavel, scales = "free_x") +
    labs(title = paste("Curvas de Resposta (MaxEnt) -", species),
         x = "Valor da variável", y = "Adequabilidade") +
    theme_minimal()
  ggsave(paste0(results_dir, species, "_maxent_curvas_resposta.png"), g_curvas, width = 12, height = 8, dpi = 300)

  cat("Curvas de resposta do MaxEnt salvas em:", results_dir, "\n")
} else {
  message("Nenhuma réplica válida suficiente para calcular curvas de resposta (MaxEnt).")
}

if (nrow(importance_df) > 0) {
  resumo_importancia <- aggregate(importancia ~ variavel, data = importance_df,
                                   FUN = function(x) c(media = mean(x), sd = sd(x)))
  resumo_importancia <- do.call(data.frame, resumo_importancia)
  colnames(resumo_importancia) <- c("variavel", "media", "sd")
  write.csv(resumo_importancia, paste0(results_dir, species, "_maxent_importancia_permutacao.csv"), row.names = FALSE)
  cat("Importância por permutação (MaxEnt) salva - para comparação com a contribuição nativa já usada.\n")
}

print("Processamento de curvas de resposta do MaxEnt completo have been run for Caiman_latirostris.")
