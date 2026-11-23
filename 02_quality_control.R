# 1. Подключение библиотек

library(edgeR)
library(ggplot2)


# 2. Создание директории результатов

dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)


# 3. Загрузка подготовленных данных

prepared_data <- readRDS(
  "data/processed/hcv_data.rds"
)

count_matrix <- prepared_data$counts
gene_annotation <- prepared_data$gene_annotation
sample_table <- prepared_data$sample_table


# 4. Проверка загруженных данных

stopifnot(
  identical(
    colnames(count_matrix),
    sample_table$Run
  )
)

stopifnot(
  all(count_matrix >= 0),
  all(count_matrix == round(count_matrix)),
  !anyNA(count_matrix)
)

cat("Размер матрицы счётчиков:\n")
print(dim(count_matrix))


# 5. Подготовка групп

sample_table$fibrosis_group <- factor(
  sample_table$fibrosis_group,
  levels = c("F0-F2", "F3-F4")
)

sample_table$sex <- factor(
  sample_table$`Sample Characteristic[sex]`
)


# 6. Создание объекта DGEList

y <- DGEList(
  counts = count_matrix,
  group = sample_table$fibrosis_group
)


# 7. Первичные показатели качества

library_size <- colSums(y$counts)

zero_fraction <- colMeans(
  y$counts == 0
)

qc_metrics <- data.frame(
  Run = colnames(y),
  Library_Size = library_size,
  Library_Size_Millions = library_size / 1e6,
  Zero_Fraction = zero_fraction
)

cat("\nРазмер библиотек:\n")
print(summary(qc_metrics$Library_Size))

cat("\nДоля нулевых значений:\n")
print(summary(qc_metrics$Zero_Fraction))


# 8. Фильтрация низкоэкспрессируемых генов

keep_genes <- filterByExpr(
  y,
  group = sample_table$fibrosis_group
)

cat("\nКоличество генов до фильтрации:\n")
print(nrow(y))

cat("Количество генов после фильтрации:\n")
print(sum(keep_genes))

cat("Удалено генов:\n")
print(sum(!keep_genes))

y <- y[keep_genes, , keep.lib.sizes = FALSE]


# 9. TMM-нормализация

y <- normLibSizes(
  y,
  method = "TMM"
)

cat("\nTMM normalization factors:\n")
print(y$samples$norm.factors)


# 10. Расчёт log-CPM

log_cpm <- cpm(
  y,
  log = TRUE,
  prior.count = 2
)


# 11. График распределения экспрессии

png(
  "results/log_cpm_distribution.png",
  width = 1200,
  height = 800,
  res = 120
)

boxplot(
  log_cpm,
  las = 2,
  outline = FALSE,
  main = "Log-CPM distribution after TMM normalization",
  xlab = "Samples",
  ylab = "Log-CPM"
)

dev.off()


# 12. График плотности распределения

log_cpm_df <- as.data.frame(log_cpm)

log_cpm_long <- stack(log_cpm_df)

colnames(log_cpm_long) <- c(
  "Log_CPM",
  "Sample"
)

density_plot <- ggplot(
  log_cpm_long,
  aes(
    x = Log_CPM,
    group = Sample
  )
) +
  geom_density(
    alpha = 0.15
  ) +
  labs(
    title = "Expression density after TMM normalization",
    x = "Log-CPM",
    y = "Density"
  ) +
  theme_minimal()

ggsave(
  "results/log_cpm_density.png",
  density_plot,
  width = 10,
  height = 7,
  dpi = 300
)


# 13. PCA

pca <- prcomp(
  t(log_cpm),
  center = TRUE,
  scale. = FALSE
)

pca_variance <- pca$sdev^2 /
  sum(pca$sdev^2) * 100

cat("\nДоля объяснённой дисперсии:\n")
print(round(pca_variance[1:5], 2))


# 14. Подготовка данных PCA

pca_data <- data.frame(
  Run = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

pca_data <- merge(
  pca_data,
  sample_table,
  by.x = "Run",
  by.y = "Run",
  sort = FALSE
)


# 15. PCA по стадиям фиброза

pca_fibrosis <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = fibrosis_group
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  labs(
    title = "PCA by fibrosis stage",
    x = paste0(
      "PC1 (",
      round(pca_variance[1], 2),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(pca_variance[2], 2),
      "%)"
    ),
    color = "Fibrosis"
  ) +
  theme_minimal()

ggsave(
  "results/pca_fibrosis.png",
  pca_fibrosis,
  width = 9,
  height = 7,
  dpi = 300
)


# 16. PCA по полу

pca_sex <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = sex
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  labs(
    title = "PCA by sex",
    x = paste0(
      "PC1 (",
      round(pca_variance[1], 2),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(pca_variance[2], 2),
      "%)"
    ),
    color = "Sex"
  ) +
  theme_minimal()

ggsave(
  "results/pca_sex.png",
  pca_sex,
  width = 9,
  height = 7,
  dpi = 300
)



# 17. Поиск потенциальных выбросов

# 17.1. Корреляция между образцами

sample_cor <- cor(
  log_cpm,
  method = "pearson"
)

median_correlation <- apply(
  sample_cor,
  2,
  function(x) median(x[x < 1])
)

correlation_threshold <- median(
  median_correlation
) - 3 * mad(median_correlation)

correlation_outlier <- (
  median_correlation < correlation_threshold
)


# 17.2. Отклонения по PCA

pc1_deviation <- abs(
  pca$x[, 1] - median(pca$x[, 1])
)

pc2_deviation <- abs(
  pca$x[, 2] - median(pca$x[, 2])
)

pc1_threshold <- 3 * mad(pca$x[, 1])
pc2_threshold <- 3 * mad(pca$x[, 2])

pc1_outlier <- pc1_deviation > pc1_threshold
pc2_outlier <- pc2_deviation > pc2_threshold


# 17.3. Объединение результатов

outlier_details <- data.frame(
  Run = colnames(log_cpm),
  Median_Correlation = median_correlation,
  Correlation_Outlier = correlation_outlier,
  PC1_Outlier = pc1_outlier,
  PC2_Outlier = pc2_outlier
)

outlier_details$Potential_Outlier <- (
  outlier_details$Correlation_Outlier |
    outlier_details$PC1_Outlier |
    outlier_details$PC2_Outlier
)

# Добавление характеристик образцов
outlier_details <- merge(
  outlier_details,
  qc_metrics,
  by = "Run",
  sort = FALSE
)

outlier_details <- merge(
  outlier_details,
  sample_table[
    ,
    c("Run", "fibrosis_group", "Sample Characteristic[sex]")
  ],
  by = "Run",
  sort = FALSE
)

# Итоговая таблица
outlier_details <- outlier_details[
  order(
    outlier_details$Potential_Outlier,
    decreasing = TRUE
  ),
]

cat("\nПричины обнаружения потенциальных выбросов:\n")

print(
  outlier_details[
    outlier_details$Potential_Outlier,
    c(
      "Run",
      "Library_Size",
      "Zero_Fraction",
      "Median_Correlation",
      "Correlation_Outlier",
      "PC1_Outlier",
      "PC2_Outlier",
      "fibrosis_group",
      "Sample Characteristic[sex]"
    )
  ]
)

# Сохранение подробной таблицы
write.csv(
  outlier_details,
  "results/outlier_details.csv",
  row.names = FALSE
)

# Список потенциальных выбросов
excluded_samples <- outlier_details$Run[
  outlier_details$Potential_Outlier
]



# 18. PCA без потенциальных выбросов

# Получаем список образцов, отмеченных алгоритмом
excluded_samples <- outlier_details$Run[
  outlier_details$Potential_Outlier
]

cat("\nКоличество потенциальных выбросов:\n")
print(length(excluded_samples))

cat("\nСписок потенциальных выбросов:\n")
print(excluded_samples)

# Выбираем образцы, которые не попали в список
keep_samples <- !colnames(log_cpm) %in% excluded_samples

if (sum(keep_samples) >= 3) {
  
  log_cpm_filtered <- log_cpm[
    ,
    keep_samples,
    drop = FALSE
  ]
  
  pca_no_outliers <- prcomp(
    t(log_cpm_filtered),
    center = TRUE,
    scale. = FALSE
  )
  
  pca_no_outliers_variance <- (
    pca_no_outliers$sdev^2 /
      sum(pca_no_outliers$sdev^2)
  ) * 100
  
  pca_no_outliers_data <- data.frame(
    Run = rownames(pca_no_outliers$x),
    PC1 = pca_no_outliers$x[, 1],
    PC2 = pca_no_outliers$x[, 2]
  )
  
  pca_no_outliers_data <- merge(
    pca_no_outliers_data,
    sample_table,
    by = "Run",
    sort = FALSE
  )
  
  pca_no_outliers_plot <- ggplot(
    pca_no_outliers_data,
    aes(
      x = PC1,
      y = PC2,
      color = fibrosis_group
    )
  ) +
    geom_point(
      size = 3,
      alpha = 0.8
    ) +
    labs(
      title = "PCA without potential outliers",
      x = paste0(
        "PC1 (",
        round(pca_no_outliers_variance[1], 2),
        "%)"
      ),
      y = paste0(
        "PC2 (",
        round(pca_no_outliers_variance[2], 2),
        "%)"
      ),
      color = "Fibrosis"
    ) +
    theme_minimal()
  
  ggsave(
    "results/pca_without_outliers.png",
    pca_no_outliers_plot,
    width = 9,
    height = 7,
    dpi = 300
  )
  
} else {
  
  warning(
    "Недостаточно образцов для PCA без потенциальных выбросов."
  )
  
}


# 19. Сохранение QC-метрик

write.csv(
  qc_metrics,
  "results/qc_metrics.csv",
  row.names = FALSE
)


# 20. Сохранение результатов нормализации

saveRDS(
  list(
    dge = y,
    log_cpm = log_cpm,
    pca = pca,
    qc_metrics = qc_metrics,
    potential_outliers = excluded_samples
  ),
  "data/processed/qc_results.rds"
)


# 21. Итог

cat("\nКонтроль качества завершён.\n")

cat("Генов после фильтрации:\n")
print(nrow(log_cpm))

cat("Образцов:\n")
print(ncol(log_cpm))

cat("Потенциальных выбросов:\n")
print(length(excluded_samples))

cat("\nСозданные файлы:\n")
print(
  list.files("results")
)
