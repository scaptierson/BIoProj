# 1. Создание необходимых директорий

dir.create(
  "data/processed",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results",
  recursive = TRUE,
  showWarnings = FALSE
)


# 2. Загрузка исходных данных

counts <- read.delim(
  "data/raw/raw-counts.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

design <- read.delim(
  "data/metadata/experiment-design.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)


# 3. Проверка размеров таблиц

cat("Размер таблицы counts:\n")
print(dim(counts))

cat("Размер таблицы design:\n")
print(dim(design))


# 4. Проверка структуры исходных данных

stopifnot(
  all(c("Gene ID", "Gene Name") %in% colnames(counts)),
  "Run" %in% colnames(design)
)

# Проверка отсутствующих значений
stopifnot(
  !anyNA(counts),
  !anyNA(design$Run)
)

# Проверка уникальности идентификаторов образцов
stopifnot(
  !anyDuplicated(colnames(counts)[-(1:2)]),
  !anyDuplicated(design$Run)
)


# 5. Сопоставление образцов

count_ids <- colnames(counts)[-(1:2)]
design_ids <- design$Run

# Проверка соответствия ID
missing_in_design <- setdiff(count_ids, design_ids)
missing_in_counts <- setdiff(design_ids, count_ids)

cat("Образцы, отсутствующие в design:\n")
print(missing_in_design)

cat("Образцы, отсутствующие в counts:\n")
print(missing_in_counts)

stopifnot(
  length(missing_in_design) == 0,
  length(missing_in_counts) == 0
)

# Выравнивание порядка образцов
design <- design[
  match(count_ids, design$Run),
]

stopifnot(
  identical(count_ids, design$Run)
)



# 6. Проверка исходных счётчиков

count_matrix_all <- as.matrix(
  counts[, -(1:2)]
)

storage.mode(count_matrix_all) <- "numeric"

stopifnot(
  is.numeric(count_matrix_all),
  all(is.finite(count_matrix_all)),
  all(count_matrix_all >= 0),
  all(count_matrix_all == round(count_matrix_all))
)

cat("Проверка исходных счётчиков пройдена.\n")


# 7. Отбор HCV генотипов 1 и 3

hcv_design <- design[
  design$`Sample Characteristic[infect]` %in% c(
    "Hepatitis C virus genotype 1",
    "Hepatitis C virus genotype 3"
  ),
]

cat("Количество образцов HCV:\n")
print(nrow(hcv_design))

cat("Распределение по генотипам:\n")
print(
  table(hcv_design$`Sample Characteristic[infect]`)
)


# 8. Объединение стадий фиброза

fibrosis_stage <- hcv_design$`Sample Characteristic[fibrosis stage]`

hcv_design$fibrosis_group <- ifelse(
  fibrosis_stage %in% c(
    "METAVIR fibrosis score 0",
    "METAVIR fibrosis score 1",
    "METAVIR fibrosis score 2"
  ),
  "F0-F2",
  ifelse(
    fibrosis_stage %in% c(
      "METAVIR fibrosis score 3",
      "METAVIR fibrosis score 4"
    ),
    "F3-F4",
    NA
  )
)

# Проверка отсутствующих групп
stopifnot(
  !anyNA(hcv_design$fibrosis_group)
)

cat("Распределение по стадиям фиброза:\n")
print(table(hcv_design$fibrosis_group))

cat("Распределение по исходным стадиям METAVIR:\n")
print(table(fibrosis_stage))


# 9. Проверка распределения по полу

cat("Распределение по полу:\n")
print(
  table(hcv_design$`Sample Characteristic[sex]`)
)

cat("Распределение пола внутри групп фиброза:\n")
print(
  table(
    hcv_design$fibrosis_group,
    hcv_design$`Sample Characteristic[sex]`
  )
)



# 10. Формирование матрицы счётчиков HCV

hcv_counts <- counts[
  ,
  c("Gene ID", "Gene Name", hcv_design$Run)
]

# Сохранение аннотации генов

gene_annotation <- hcv_counts[
  ,
  c("Gene ID", "Gene Name")
]

# Создание матрицы счётчиков

count_matrix <- as.matrix(
  hcv_counts[, -(1:2)]
)

storage.mode(count_matrix) <- "numeric"

# Назначение настоящих идентификаторов генов

rownames(count_matrix) <- as.character(
  gene_annotation$`Gene ID`
)

# Проверки идентификаторов

stopifnot(
  !anyNA(rownames(count_matrix)),
  !any(rownames(count_matrix) == ""),
  !anyDuplicated(rownames(count_matrix))
)

# Проверка соответствия образцов

stopifnot(
  identical(colnames(count_matrix), hcv_design$Run)
)

# Проверка соответствия аннотации

stopifnot(
  identical(
    rownames(count_matrix),
    as.character(gene_annotation$`Gene ID`)
  )
)

cat("Размер матрицы HCV:\n")
print(dim(count_matrix))

cat("\nПервые 10 идентификаторов генов:\n")
print(head(rownames(count_matrix), 10))

cat("\nПервые 10 названий генов:\n")
print(head(gene_annotation$`Gene Name`, 10))


# 11. Сохранение подготовленных данных

prepared_data <- list(
  counts = count_matrix,
  gene_annotation = gene_annotation,
  sample_table = hcv_design
)

saveRDS(
  prepared_data,
  "data/processed/hcv_data.rds"
)

write.csv(
  hcv_design,
  "results/hcv_sample_table.csv",
  row.names = FALSE
)


# 12. Итоговая проверка

cat("\nПодготовка данных завершена.\n")

cat("Количество образцов:\n")
print(ncol(count_matrix))

cat("Количество генов:\n")
print(nrow(count_matrix))

cat("Количество образцов по группам:\n")
print(table(hcv_design$fibrosis_group))

cat("\nСозданные файлы:\n")
print(
  c(
    "data/processed/hcv_data.rds",
    "results/hcv_sample_table.csv"
  )
)
