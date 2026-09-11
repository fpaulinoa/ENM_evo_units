# Variable contribution by evolutionary unit
# by Fernando Paulino and Carolina Loss
# JAN - 2025
# Fixed - SEP/2025 (all 8 populations, without the "Caimamn_gpi" typo,
# filtering by valid runs using MNE010 metrics, per-file tryCatch,
# results saved to disk (plots and tables) instead of just print())

library(rvest)
library(dplyr)
library(ggplot2)
library(broom)

base_dir <- "./"  # ajuste aqui se os dados estiverem em outra pasta
models_dir <- paste0(base_dir, "output/models/")
out_dir <- paste0(base_dir, "Results New/comparacao_variaveis/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Fix: dynamically lists all populations present, instead of a
# fixed list of only 2 (with the "Caimamn_gpi" typo)
all_items <- list.files(models_dir, full.names = TRUE)
populacoes <- basename(all_items[file.info(all_items)$isdir])
cat("Populações a processar:", length(populacoes), "\n")
print(populacoes)

# ---------------------------------------------------------------------------
# Load the per-run quality metrics (generated in MNE010/reused
# in MNE011), to include here only the variable contribution from runs that
# passed the quality criteria (wAUC>=0.5, TSS>=0, CBI>=0.4)
# ---------------------------------------------------------------------------
eval_file <- paste0(base_dir, "output/evaluation/eval_df_completo.csv")
run_valid_lookup <- NULL
if (file.exists(eval_file)) {
  eval_df <- read.csv(eval_file)

  run_lookup <- data.frame(species = character(), Replica = character(), run_num = integer())
  for (populacao in populacoes) {
    species_dir <- paste0(models_dir, populacao, "/")
    run_dirs <- list.dirs(species_dir, recursive = FALSE)
    for (run_folder in run_dirs) {
      run_num <- as.integer(gsub(".*run_", "", run_folder))
      replicate_files <- list.files(run_folder, pattern = "_samplePredictions.csv")
      replicate_names <- sub("_samplePredictions.csv", "", replicate_files)
      if (length(replicate_names) > 0) {
        run_lookup <- rbind(run_lookup, data.frame(
          species = populacao, Replica = replicate_names, run_num = run_num
        ))
      }
    }
  }
  eval_df <- merge(eval_df, run_lookup, by = c("species", "Replica"), all.x = TRUE)
  run_metrics <- aggregate(cbind(wAUC, TSS, CBI) ~ species + run_num, data = eval_df, FUN = mean)
  run_metrics$run_valid <- with(run_metrics, wAUC >= 0.5 & TSS >= 0.0 & CBI >= 0.4)
  run_valid_lookup <- run_metrics
} else {
  warning("Não encontrei ", eval_file, " - seguindo SEM filtro de qualidade por run ",
          "(todas as runs serão incluídas). Rode o MNE010 antes para um resultado mais rigoroso.")
}

# Function to extract the variable-contribution table from the
# maxent.html (summary already averaged across that run's 5 replicates)
extrair_tabela_contribuicao <- function(caminho_arquivo) {
  html <- read_html(caminho_arquivo)
  tabela <- html %>%
    html_nodes("table") %>%
    html_table(fill = TRUE) %>%
    .[[1]]
  if (ncol(tabela) != 3) {
    stop("Tabela com formato inesperado (", ncol(tabela), " colunas) em ", caminho_arquivo)
  }
  colnames(tabela) <- c("Variavel", "Contribuicao_Percentual", "Importancia_Permutacao")
  return(tabela)
}

# Function to process all runs for a population
processar_populacao <- function(populacao, models_dir, run_valid_lookup) {

  sp_run_valid <- NULL
  if (!is.null(run_valid_lookup)) {
    sp_run_valid <- subset(run_valid_lookup, species == populacao)
  }

  linhas <- list()

  for (run_num in 1:20) {

    # Fix: skips runs that did not pass the quality criteria
    if (!is.null(sp_run_valid)) {
      ok <- sp_run_valid$run_valid[sp_run_valid$run_num == run_num]
      if (length(ok) == 0 || !isTRUE(ok)) next
    }

    arquivo <- paste0(models_dir, populacao, "/run_", run_num, "/maxent.html")

    # Fix: per-file tryCatch; a run with a corrupted/missing HTML
    # does not bring down processing of the other runs/populations
    tabela <- tryCatch(extrair_tabela_contribuicao(arquivo), error = function(e) {
      message("Falhou ao ler ", arquivo, ": ", conditionMessage(e))
      NULL
    })

    if (!is.null(tabela)) {
      tabela$Repeticao <- run_num
      linhas[[length(linhas) + 1]] <- tabela
    }
  }

  if (length(linhas) == 0) {
    warning("Nenhuma run válida/legível encontrada para ", populacao)
    return(NULL)
  }

  dados_populacao <- bind_rows(linhas)
  dados_populacao$Populacao <- populacao
  return(dados_populacao)
}

# Process all populations and consolidate
dados_finais <- lapply(populacoes, function(pop) {
  processar_populacao(pop, models_dir, run_valid_lookup)
}) %>% bind_rows()

if (nrow(dados_finais) == 0) {
  stop("Nenhum dado de contribuição de variável foi lido - confira se as runs do MNE09 geraram os maxent.html corretamente.")
}

# Save the raw data (one row per run x variable x population)
write.csv(dados_finais, paste0(out_dir, "contribuicao_variaveis_bruto.csv"), row.names = FALSE)

# Statistical summary: mean and standard deviation per variable and population
resumo_contribuicoes <- dados_finais %>%
  group_by(Populacao, Variavel) %>%
  summarise(
    Media_Contribuicao = mean(Contribuicao_Percentual, na.rm = TRUE),
    Desvio_Padrao = sd(Contribuicao_Percentual, na.rm = TRUE),
    n_runs = n(),
    .groups = "drop"
  )

write.csv(resumo_contribuicoes, paste0(out_dir, "contribuicao_variaveis_resumo.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Comparative plot (all populations together), saved to disk
# ---------------------------------------------------------------------------
grafico <- ggplot(resumo_contribuicoes, aes(x = Variavel, y = Media_Contribuicao, fill = Populacao)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = pmax(0, Media_Contribuicao - Desvio_Padrao),
                     ymax = Media_Contribuicao + Desvio_Padrao),
                position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Contribuição das Variáveis por População",
       x = "Variável", y = "Contribuição Média (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(out_dir, "contribuicao_variaveis_todas_populacoes.png"), grafico,
       width = 12, height = 7, dpi = 300)

# Fix: also saves one plot per population individually, since a
# single plot with all 8 populations gets cluttered/hard to read
for (pop in unique(resumo_contribuicoes$Populacao)) {
  sub_dados <- subset(resumo_contribuicoes, Populacao == pop)
  g_ind <- ggplot(sub_dados, aes(x = Variavel, y = Media_Contribuicao)) +
    geom_bar(stat = "identity", fill = "salmon") +
    geom_errorbar(aes(ymin = pmax(0, Media_Contribuicao - Desvio_Padrao),
                       ymax = Media_Contribuicao + Desvio_Padrao), width = 0.25) +
    labs(title = paste("Contribuição das Variáveis -", pop),
         x = "Variável", y = "Contribuição Média (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # creates the same per-species subfolder used in MNE011, to keep everything together
  pop_dir <- paste0(base_dir, "Results New/", pop, "/")
  dir.create(pop_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(pop_dir, pop, "_contribuicao_variaveis.png"), g_ind, width = 8, height = 6, dpi = 300)
}

# ---------------------------------------------------------------------------
# Statistical tests comparing populations, for each variable
# ---------------------------------------------------------------------------
resultados_anova <- tryCatch({
  dados_finais %>%
    group_by(Variavel) %>%
    do(tidy(aov(Contribuicao_Percentual ~ Populacao, data = .))) %>%
    ungroup()
}, error = function(e) {
  message("Falhou ao calcular ANOVA: ", conditionMessage(e))
  NULL
})

resultados_kruskal <- tryCatch({
  dados_finais %>%
    group_by(Variavel) %>%
    do(tidy(kruskal.test(Contribuicao_Percentual ~ Populacao, data = .))) %>%
    ungroup()
}, error = function(e) {
  message("Falhou ao calcular Kruskal-Wallis: ", conditionMessage(e))
  NULL
})

if (!is.null(resultados_anova)) {
  write.csv(resultados_anova, paste0(out_dir, "anova_contribuicao_variaveis.csv"), row.names = FALSE)
}
if (!is.null(resultados_kruskal)) {
  write.csv(resultados_kruskal, paste0(out_dir, "kruskal_contribuicao_variaveis.csv"), row.names = FALSE)
}

# ---------------------------------------------------------------------------
# Boxplot of one specific variable ("prox") comparing populations
# ---------------------------------------------------------------------------
variavel_exemplo <- "prox"
boxplot <- ggplot(dados_finais %>% filter(Variavel == variavel_exemplo),
                   aes(x = Populacao, y = Contribuicao_Percentual, fill = Populacao)) +
  geom_boxplot() +
  labs(title = paste("Contribuição da Variável", variavel_exemplo, "por População"),
       x = "População", y = "Contribuição Percentual") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(out_dir, "boxplot_", variavel_exemplo, "_por_populacao.png"), boxplot,
       width = 10, height = 6, dpi = 300)

cat("MNE012 concluído. Resultados salvos em:", out_dir, "\n")
cat("Gráficos individuais por população salvos dentro de cada subpasta em Results New/.\n")
