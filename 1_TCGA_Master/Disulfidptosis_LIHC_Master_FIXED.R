DATA_DIR  <- "."          # Directorio donde están los datos TCGA descargados
OUT_DIR   <- "resultados" # Directorio de salida (se crea si no existe)


# Crear carpetas de salida
dirs_needed <- c(
  file.path(OUT_DIR, "01_datos_procesados"),
  file.path(OUT_DIR, "02_scores_clusters"),
  file.path(OUT_DIR, "03_supervivencia"),
  file.path(OUT_DIR, "04_microambiente"),
  file.path(OUT_DIR, "05_expresion_diferencial"),
  file.path(OUT_DIR, "06_figuras"),
  file.path(OUT_DIR, "07_tablas_suplementarias")
)
invisible(lapply(dirs_needed, dir.create, recursive = TRUE, showWarnings = FALSE))

# Rutas cortas para comodidad
DOUT  <- file.path(OUT_DIR, "01_datos_procesados")
SOUT  <- file.path(OUT_DIR, "02_scores_clusters")
KOUT  <- file.path(OUT_DIR, "03_supervivencia")
MOUT  <- file.path(OUT_DIR, "04_microambiente")
EOUT  <- file.path(OUT_DIR, "05_expresion_diferencial")
FOUT  <- file.path(OUT_DIR, "06_figuras")
TOUT  <- file.path(OUT_DIR, "07_tablas_suplementarias")

## Paleta
BLUE_DARKEST <- "#03254C"   # Marino profundo — títulos premium
BLUE_DARK    <- "#1167B1"   # Azul real — color principal
BLUE_MID     <- "#2A9D8F"   # Teal vibrante — acento complementario
BLUE_LIGHT   <- "#74B9D4"   # Celeste — C2 / Low Score
BLUE_PALE    <- "#D0E8F5"   # Azul hielo — fondos suaves
ACCENT_CORAL <- "#E76F51"   # Coral cálido — resaltar significativo
GREY_LINE    <- "#DEE2E6"   # Gris líneas
GREY_TEXT    <- "#495057"   # Gris texto secundario

COL_HIGH <- BLUE_DARK
COL_LOW  <- BLUE_LIGHT
COL_C1   <- BLUE_DARK     # C1 = alto score → azul oscuro
COL_C2   <- BLUE_LIGHT    # C2 = bajo score → azul claro

## Tema ggplot2
theme_blue <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) %+replace%
    ggplot2::theme(
      # Tipografía
      text              = ggplot2::element_text(family = "sans", color = GREY_TEXT),
      plot.title        = ggplot2::element_text(
                            face = "bold", size = base_size + 3,
                            color = BLUE_DARKEST, hjust = 0,
                            margin = ggplot2::margin(b = 4)),
      plot.subtitle     = ggplot2::element_text(
                            size = base_size - 0.5, color = "grey50",
                            hjust = 0, margin = ggplot2::margin(b = 10)),
      plot.caption      = ggplot2::element_text(
                            size = base_size - 2, color = "grey65",
                            hjust = 1, margin = ggplot2::margin(t = 8)),
      # Ejes
      axis.title        = ggplot2::element_text(
                            face = "bold", size = base_size,
                            color = BLUE_DARKEST),
      axis.text         = ggplot2::element_text(
                            size = base_size - 1, color = GREY_TEXT),
      axis.line         = ggplot2::element_line(
                            color = GREY_LINE, linewidth = 0.6),
      axis.ticks        = ggplot2::element_line(
                            color = GREY_LINE, linewidth = 0.4),
      # Leyenda
      legend.title      = ggplot2::element_text(
                            face = "bold", size = base_size - 1,
                            color = BLUE_DARKEST),
      legend.text       = ggplot2::element_text(
                            size = base_size - 2, color = GREY_TEXT),
      legend.key.size   = ggplot2::unit(0.45, "cm"),
      legend.background = ggplot2::element_rect(
                            fill = "white", color = GREY_LINE,
                            linewidth = 0.3),
      legend.margin     = ggplot2::margin(4, 6, 4, 6),
      # Grid
      panel.grid.major  = ggplot2::element_line(
                            color = "#EEF2F7", linewidth = 0.35),
      panel.grid.minor  = ggplot2::element_blank(),
      # Fondos
      plot.background   = ggplot2::element_rect(fill = "white", color = NA),
      panel.background  = ggplot2::element_rect(fill = "white", color = NA),
      # Facetas
      strip.background  = ggplot2::element_rect(fill = BLUE_PALE, color = NA),
      strip.text        = ggplot2::element_text(
                            face = "bold", color = BLUE_DARKEST,
                            size = base_size - 1),
      # Márgenes
      plot.margin       = ggplot2::margin(14, 18, 12, 14)
    )
}

## figuras (PNG alta resolución + PDF)
guardar_figura <- function(plot, nombre, ancho = 12, alto = 9,
                           dpi = 300, dir = FOUT) {
  ruta_base <- file.path(dir, nombre)
  ggplot2::ggsave(paste0(ruta_base, ".png"),
                  plot  = plot,
                  width = ancho, height = alto,
                  dpi   = dpi, bg = "white")
  ggplot2::ggsave(paste0(ruta_base, ".pdf"),
                  plot  = plot,
                  width = ancho, height = alto,
                  device = grDevices::cairo_pdf)
  message("Guardé: ", basename(ruta_base))
}

message("Ya tengo la configuración global lista.")


message("Empiezo con la carga de los datos de TCGA-LIHC.")

## Instalar paquetes necesarios 
pkgs_bioc <- c("TCGAbiolinks", "SummarizedExperiment")
pkgs_cran <- c("data.table", "dplyr", "tidyr")

for (p in pkgs_bioc)
  if (!requireNamespace(p, quietly = TRUE))
    BiocManager::install(p, ask = FALSE, update = FALSE)

for (p in pkgs_cran)
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(data.table)
library(dplyr)

# Opción A: reutilizar datos ya descargados
se_path  <- file.path(DOUT, "TCGA_LIHC_SE_raw.rds")
clin_path <- file.path(DOUT, "TCGA_LIHC_clinical.csv")

if (!file.exists(se_path)) {
  message("Voy a descargar los datos de TCGA-LIHC, esto puede tardar un rato...")

  query_rna <- GDCquery(
    project       = "TCGA-LIHC",
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type   = c("Primary Tumor", "Solid Tissue Normal")
  )
  GDCdownload(query_rna, method = "api", files.per.chunk = 10,
              directory = DATA_DIR)
  se <- GDCprepare(query_rna, directory = DATA_DIR)
  saveRDS(se, se_path)
  message("Guardé el SummarizedExperiment en: ", se_path)
} else {
  message("Cargo el SE desde la caché: ", se_path)
  se <- readRDS(se_path)
}

# Datos clínicos
if (!file.exists(clin_path)) {
  clinical <- GDCquery_clinic(project = "TCGA-LIHC", type = "clinical")
  fwrite(clinical, clin_path)
} else {
  clinical <- fread(clin_path, data.table = FALSE)
}

# Normalización log2(FPKM+1)
fpkm_sym_path <- file.path(DOUT, "TCGA_LIHC_log2fpkm_SYMBOL.rds")
coldata_path  <- file.path(DOUT, "TCGA_LIHC_colData.rds")
rowdata_path  <- file.path(DOUT, "TCGA_LIHC_rowData.rds")

col_data <- as.data.frame(colData(se))
row_data <- as.data.frame(rowData(se))
saveRDS(col_data, coldata_path)
saveRDS(row_data, rowdata_path)

if (!file.exists(fpkm_sym_path)) {
  # Intentar fpkm_unstrand primero, luego calcular manualmente
  if ("fpkm_unstrand" %in% assayNames(se)) {
    fpkm_mat <- assay(se, "fpkm_unstrand")
  } else {
    message("No encontré 'fpkm_unstrand', así que voy a calcular el FPKM a partir de los counts...")
    counts_mat  <- assay(se, "unstranded")
    if (!"gene_length" %in% colnames(row_data) ||
        all(is.na(row_data$gene_length))) {
      stop(paste0(
        "No se encontró 'gene_length' en rowData. ",
        "Proporciona la longitud génica manualmente o usa un SE con ",
        "FPKM precalculado (assay 'fpkm_unstrand')."))
    }
    gene_lengths <- pmax(row_data$gene_length, 1, na.rm = TRUE)
    total_counts <- colSums(counts_mat)
    fpkm_mat <- sweep(counts_mat, 2, total_counts / 1e6, "/")
    fpkm_mat <- sweep(fpkm_mat, 1, gene_lengths / 1e3, "/")
    rm(counts_mat); gc()
  }

  log2_fpkm <- log2(fpkm_mat + 1)
  rm(fpkm_mat); gc()

  # Convertir ENSEMBL → SYMBOL
  ens_ids    <- gsub("\\..*", "", rownames(log2_fpkm))
  ens_lookup <- gsub("\\..*", "", row_data$gene_id)
  sym_vec    <- row_data$gene_name[match(ens_ids, ens_lookup)]
  rownames(log2_fpkm) <- sym_vec

  # Limpiar: quitar NA, duplicados y genes con poca expresión
  log2_fpkm <- log2_fpkm[!is.na(rownames(log2_fpkm)), ]
  log2_fpkm <- log2_fpkm[!duplicated(rownames(log2_fpkm)), ]
  keep       <- rowSums(log2_fpkm > 0.1) >= (0.10 * ncol(log2_fpkm))
  log2_fpkm  <- log2_fpkm[keep, ]

  saveRDS(log2_fpkm, fpkm_sym_path)
  rm(se, log2_fpkm); gc()
} else {
  message("El log2-FPKM con símbolos ya existe, así que omito la normalización.")
  rm(se); gc()
}

message("Genes que quedaron en la matriz filtrada: ",
        nrow(readRDS(fpkm_sym_path)))
message("Terminé el bloque 1.")


message("Ahora calculo el score de disulfidptosis.")

library(dplyr)

log2_fpkm <- readRDS(fpkm_sym_path)

# Genes del panel de disulfidptosis
genes_canonical <- c("SLC7A11", "SLC3A2", "SLC2A1", "LRPPRC",
                     "RNF213",  "NUBPL",  "OXSM",   "SHMT2",
                     "GYS1",    "FLII")
genes_extended  <- c(genes_canonical,
                     "G6PD", "TKT", "TALDO1", "PGLS", "PGD", "H6PD")

genes_can_ok <- genes_canonical[genes_canonical %in% rownames(log2_fpkm)]
genes_ext_ok <- genes_extended[genes_extended   %in% rownames(log2_fpkm)]

message("Genes canónicos encontrados: ",
        length(genes_can_ok), "/", length(genes_canonical))
message("Genes extendidos encontrados: ",
        length(genes_ext_ok), "/", length(genes_extended))
message("Genes canónicos que faltan: ",
        paste(setdiff(genes_canonical, genes_can_ok), collapse = ", "))

# Score por z-scores (método principal)
expr_z          <- t(scale(t(log2_fpkm)))
score_canonical <- colMeans(expr_z[genes_can_ok, , drop = FALSE], na.rm = TRUE)
score_extended  <- colMeans(expr_z[genes_ext_ok, , drop = FALSE], na.rm = TRUE)

# PCA score (método alternativo robusto)
sub_mat  <- t(log2_fpkm[genes_can_ok, , drop = FALSE])
keep_var <- apply(sub_mat, 2, var, na.rm = TRUE) > 0
sub_mat  <- sub_mat[, keep_var, drop = FALSE]
pca_res  <- prcomp(sub_mat, scale. = TRUE, center = TRUE)
pc1_score <- pca_res$x[, 1]
# Orientar PC1 para que correlacione positivamente con expresión media
if (cor(pc1_score, rowMeans(sub_mat), use = "complete.obs") < 0)
  pc1_score <- -pc1_score

# Varianza explicada por PC1
pve_pc1 <- round(summary(pca_res)$importance[2, 1] * 100, 1)
message("PC1 explica ", pve_pc1, "% de la varianza")

scores_df <- data.frame(
  sample_id                = colnames(log2_fpkm),
  Disulfidptosis_Canonical = as.numeric(score_canonical),
  Disulfidptosis_Extended  = as.numeric(score_extended),
  PCA_score                = as.numeric(pc1_score),
  stringsAsFactors         = FALSE
)

write.csv(scores_df,
          file.path(SOUT, "Scores_Disulfidptosis_todos.csv"),
          row.names = FALSE)
saveRDS(scores_df,
        file.path(SOUT, "Scores_Disulfidptosis_todos.rds"))

rm(expr_z); gc()
message("Terminé el bloque 2.")


message("Sigo con el clustering de pacientes.")

for (p in c("cluster", "factoextra", "ggplot2"))
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")

library(cluster)
library(ggplot2)

scores_df <- readRDS(file.path(SOUT, "Scores_Disulfidptosis_todos.rds"))
col_data  <- readRDS(coldata_path)

# Filtrar solo tumores primarios
tumor_ids    <- col_data$barcode[
  col_data$sample_type %in% c("Primary Tumor", "01")]
scores_tumor <- scores_df[scores_df$sample_id %in% tumor_ids, ]
message("Tumores primarios que voy a usar para el clustering: ", nrow(scores_tumor))

clust_mat <- scale(scores_tumor[,
                   c("Disulfidptosis_Canonical",
                     "Disulfidptosis_Extended",
                     "PCA_score")])
rownames(clust_mat) <- scores_tumor$sample_id

# k óptimo por silhouette
set.seed(42)
sil_width <- sapply(2:6, function(k) {
  km  <- kmeans(clust_mat, centers = k, nstart = 25, iter.max = 100)
  sil <- silhouette(km$cluster, dist(clust_mat))
  mean(sil[, 3])
})
k_optimal <- which.max(sil_width) + 1
message("El k óptimo según silhouette es: ", k_optimal)

# k-means final — siempre k = 2 (biológicamente interpretable en HCC:
# alto score / bajo score; k variable rompe figuras, Cox y validación
# que asumen exactamente C1 y C2)
k_use <- 2L
if (k_optimal != 2L)
  message("Silhouette sugirió k = ", k_optimal,
          ", pero uso k = 2 porque así está diseñado el análisis.")
set.seed(42)
km_final <- kmeans(clust_mat, centers = k_use,
                   nstart = 50, iter.max = 200)

scores_tumor$Cluster <- paste0("C", km_final$cluster)

# Garantizar que C1 = alto score (BLUE_DARK) y C2 = bajo score (BLUE_LIGHT)
med_by_cluster <- tapply(scores_tumor$Disulfidptosis_Canonical,
                         scores_tumor$Cluster, median, na.rm = TRUE)
highest_cluster <- names(which.max(med_by_cluster))
if (highest_cluster != "C1") {
  scores_tumor$Cluster <- ifelse(scores_tumor$Cluster == "C1",
                                 "C2", "C1")
  message("Reetiqueté los clusters: C1 queda como alto score y C2 como bajo score.")
}

write.csv(scores_tumor,
          file.path(SOUT, "Clusters_Disulfidptosis_tumores.csv"),
          row.names = FALSE)
saveRDS(scores_tumor,
        file.path(SOUT, "Clusters_Disulfidptosis_tumores.rds"))

# Guardar ancho silhouette para figura
sil_df <- data.frame(k = 2:6, sil = sil_width)
saveRDS(sil_df, file.path(SOUT, "Silhouette_widths.rds"))

message("Terminé el bloque 3.")


message("Ahora hago el análisis de supervivencia.")

for (p in c("survival", "survminer"))
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")

library(survival)
library(survminer)
library(dplyr)

clusters_df <- readRDS(file.path(SOUT, "Clusters_Disulfidptosis_tumores.rds"))
clinical    <- read.csv(clin_path, stringsAsFactors = FALSE)

# Estandarizar patient_id
clusters_df$patient_id <- toupper(substr(clusters_df$sample_id, 1, 12))
clinical$patient_id    <- toupper(clinical$submitter_id)

surv_df <- inner_join(clusters_df, clinical, by = "patient_id") %>%
  mutate(
    # Tiempo de seguimiento: días hasta muerte o último contacto
    OS_time = dplyr::case_when(
      !is.na(days_to_death) & days_to_death > 0 ~
        as.numeric(days_to_death),
      !is.na(days_to_last_follow_up) & days_to_last_follow_up > 0 ~
        as.numeric(days_to_last_follow_up),
      TRUE ~ NA_real_
    ),
    OS_status = dplyr::case_when(
      tolower(vital_status) %in% c("dead","deceased") ~ 1L,
      tolower(vital_status) %in% c("alive","living")  ~ 0L,
      TRUE ~ NA_integer_
    ),
    # Grupos por mediana del score canónico
    Score_group = ifelse(
      Disulfidptosis_Canonical >=
        median(Disulfidptosis_Canonical, na.rm = TRUE),
      "High", "Low"),
    # Variables clínicas simplificadas
    Age_group = factor(
      ifelse(age_at_index >= 60, ">=60", "<60"),
      levels = c("<60", ">=60")),
    Stage_simple = dplyr::case_when(
      grepl("Stage I$|Stage II$|Stage IIA$|Stage IIB$",
            ajcc_pathologic_stage, ignore.case = TRUE) ~ "Early",
      grepl("Stage III|Stage IV",
            ajcc_pathologic_stage, ignore.case = TRUE) ~ "Advanced",
      TRUE ~ NA_character_
    ),
    Stage_simple = factor(Stage_simple,
                          levels = c("Advanced", "Early"))
  ) %>%
  filter(OS_time > 0, !is.na(OS_time), !is.na(OS_status))

message("Pacientes que entran en el análisis de supervivencia: ", nrow(surv_df))
message("Eventos (muertes): ", sum(surv_df$OS_status, na.rm = TRUE))
message("Alto score: ",
        sum(surv_df$Score_group == "High", na.rm = TRUE),
        " | Bajo score: ",
        sum(surv_df$Score_group == "Low",  na.rm = TRUE))

saveRDS(surv_df, file.path(KOUT, "Datos_supervivencia.rds"))
write.csv(surv_df, file.path(KOUT, "Datos_supervivencia.csv"),
          row.names = FALSE)
message("Terminé el bloque 4.")


message("Ahora reviso el microambiente tumoral.")

if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes", repos = "https://cloud.r-project.org")
if (!requireNamespace("IOBR", quietly = TRUE)) {
  # Versión fijada para reproducibilidad — actualizar deliberadamente si se requiere
  remotes::install_github("IOBR/IOBR", ref = "v0.99.9", upgrade = "never")
}


library(IOBR)
library(dplyr)
library(tidyr)

log2_fpkm <- readRDS(fpkm_sym_path)
surv_df   <- readRDS(file.path(KOUT, "Datos_supervivencia.rds"))

tumor_ids  <- surv_df$sample_id
expr_tumor <- log2_fpkm[, colnames(log2_fpkm) %in% tumor_ids, drop = FALSE]
message("Muestras de tumor que voy a usar para la deconvolución: ", ncol(expr_tumor))

epic_path <- file.path(MOUT, "EPIC_normalizado.csv")

if (!file.exists(epic_path)) {
  epic_res <- deconvo_tme(
    eset   = expr_tumor,
    method = "epic",
    arrays = FALSE
  )  # perm no aplica a EPIC (parámetro exclusivo de CIBERSORT)
  write.csv(cibersort_res,
            file.path(MOUT, "EPIC_bruto.csv"),
            row.names = FALSE)

  # Identificar columnas de fracciones celulares EPIC (_EPIC suffix)
  frac_cols <- grep("_EPIC$",
                    colnames(epic_res), value = TRUE, ignore.case = FALSE)
  # EPIC no genera columnas de p-valor ni RMSE; frac_cols ya son solo fracciones

  # Renombrar columnas: mantener nombre limpio sin sufijo _EPIC
  cib_renamed <- epic_res
  new_names   <- gsub("_EPIC$", "", colnames(epic_res))
  new_names   <- trimws(new_names)
  colnames(cib_renamed) <- new_names

  # Recalcular frac_cols con nuevos nombres (sin sufijo _EPIC)
  frac_cols_clean <- gsub("_EPIC$", "", frac_cols)

  # Normalizar fracciones a suma = 1 por muestra
  fc_idx <- match(frac_cols_clean, colnames(cib_renamed))
  fc_idx <- fc_idx[!is.na(fc_idx)]
  row_sums <- rowSums(cib_renamed[, fc_idx], na.rm = TRUE)
  row_sums[row_sums == 0] <- 1
  cib_renamed[, fc_idx] <- cib_renamed[, fc_idx] / row_sums

  write.csv(cib_renamed, epic_path, row.names = FALSE)
  saveRDS(list(cib = cib_renamed,
               frac_cols = frac_cols_clean[!is.na(fc_idx)]),  # nombres limpios sin sufijo _EPIC
          file.path(MOUT, "EPIC_objeto.rds"))
} else {
  message("EPIC ya está calculado, lo cargo desde la caché.")
}

rm(log2_fpkm, expr_tumor); gc()
message("Terminé el bloque 5.")


message("Sigo con DEG, GSEA y ORA.")

pkgs_bioc6 <- c("limma", "clusterProfiler", "org.Hs.eg.db", "enrichplot")
for (p in pkgs_bioc6)
  if (!requireNamespace(p, quietly = TRUE))
    BiocManager::install(p, ask = FALSE, update = FALSE)

library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)

surv_df    <- readRDS(file.path(KOUT, "Datos_supervivencia.rds"))
log2_fpkm  <- readRDS(fpkm_sym_path)

# Alinear muestras disponibles
samples_use  <- surv_df$sample_id[!is.na(surv_df$Score_group)]
samples_use  <- samples_use[samples_use %in% colnames(log2_fpkm)]
expr_aligned <- log2_fpkm[, samples_use, drop = FALSE]
rm(log2_fpkm); gc()

group <- factor(
  ifelse(samples_use %in%
           surv_df$sample_id[surv_df$Score_group == "High"],
         "High", "Low"),
  levels = c("Low", "High")
)
message("Muestras: High = ",
        sum(group == "High"), ", Low = ", sum(group == "Low"))

# DEG con limma
design <- model.matrix(~ group)
fit    <- lmFit(expr_aligned, design)
fit    <- eBayes(fit)
deg_res <- topTable(fit, coef = 2, number = Inf,
                    adjust.method = "BH", sort.by = "logFC")
deg_res$gene <- rownames(deg_res)
rm(expr_aligned); gc()

n_up   <- sum(deg_res$adj.P.Val < 0.05 & deg_res$logFC >  1, na.rm = TRUE)
n_down <- sum(deg_res$adj.P.Val < 0.05 & deg_res$logFC < -1, na.rm = TRUE)
message("DEG significativos: Up = ", n_up, ", Down = ", n_down)

write.csv(deg_res,
          file.path(EOUT, "DEG_Alto_vs_Bajo_score.csv"),
          row.names = FALSE)

# GSEA GO Biological Process
gene_list <- sort(setNames(deg_res$logFC, deg_res$gene),
                  decreasing = TRUE)
gene_list <- gene_list[!is.na(names(gene_list))]

gsea_go_path <- file.path(EOUT, "GSEA_GO_BP.rds")
if (!file.exists(gsea_go_path)) {
  set.seed(42)
  gsea_go <- gseGO(
    geneList     = gene_list,
    ont          = "BP",
    OrgDb        = org.Hs.eg.db,
    keyType      = "SYMBOL",
    minGSSize    = 15,
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    eps          = 0,
    seed         = TRUE,
    verbose      = TRUE
  )
  saveRDS(gsea_go, gsea_go_path)
  write.csv(as.data.frame(gsea_go),
            file.path(EOUT, "GSEA_GO_BP_tabla.csv"),
            row.names = FALSE)
} else {
  gsea_go <- readRDS(gsea_go_path)
}
message("GSEA GO BP: ",
        nrow(as.data.frame(gsea_go)), " términos significativos")

# GSEA KEGG
gene_entrez <- bitr(names(gene_list),
                    fromType = "SYMBOL",
                    toType   = "ENTREZID",
                    OrgDb    = org.Hs.eg.db)
gl_entrez <- gene_list[gene_entrez$SYMBOL]
names(gl_entrez) <- gene_entrez$ENTREZID
gl_entrez <- sort(gl_entrez[!duplicated(names(gl_entrez))],
                  decreasing = TRUE)

gsea_kegg_path <- file.path(EOUT, "GSEA_KEGG.rds")
if (!file.exists(gsea_kegg_path)) {
  set.seed(42)
  gsea_kegg <- gseKEGG(
    geneList      = gl_entrez,
    organism      = "hsa",
    minGSSize     = 15,
    maxGSSize     = 500,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    eps           = 0,
    seed          = TRUE,
    verbose       = TRUE
  )
  saveRDS(gsea_kegg, gsea_kegg_path)
  write.csv(as.data.frame(gsea_kegg),
            file.path(EOUT, "GSEA_KEGG_tabla.csv"),
            row.names = FALSE)
} else {
  gsea_kegg <- readRDS(gsea_kegg_path)
}
message("GSEA KEGG: ",
        nrow(as.data.frame(gsea_kegg)), " rutas significativas")

# ORA GO BP
sig_genes <- deg_res$gene[
  !is.na(deg_res$adj.P.Val) &
  deg_res$adj.P.Val < 0.05 &
  abs(deg_res$logFC) > 1]
sig_genes <- sig_genes[!is.na(sig_genes)]

ora_go_path <- file.path(EOUT, "ORA_GO_BP.rds")
if (length(sig_genes) >= 10 && !file.exists(ora_go_path)) {
  ora_go <- enrichGO(
    gene          = sig_genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  saveRDS(ora_go, ora_go_path)
  write.csv(as.data.frame(ora_go),
            file.path(EOUT, "ORA_GO_BP_tabla.csv"),
            row.names = FALSE)
} else if (length(sig_genes) < 10) {
  warning("Tengo menos de 10 genes significativos para el ORA: ",
          length(sig_genes))
}

message("Terminé el bloque 6.")


message("Ahora armo las tablas suplementarias.")

if (!requireNamespace("openxlsx", quietly = TRUE))
  install.packages("openxlsx", repos = "https://cloud.r-project.org")
if (!requireNamespace("tibble",  quietly = TRUE))
  install.packages("tibble", repos = "https://cloud.r-project.org")

library(openxlsx)
library(dplyr)
library(stringr)
library(tibble)

surv_df   <- readRDS(file.path(KOUT, "Datos_supervivencia.rds"))
deg_res   <- read.csv(file.path(EOUT, "DEG_Alto_vs_Bajo_score.csv"),
                      stringsAsFactors = FALSE)
gsea_go   <- readRDS(file.path(EOUT, "GSEA_GO_BP.rds"))
gsea_kegg <- readRDS(file.path(EOUT, "GSEA_KEGG.rds"))

go_df   <- as.data.frame(gsea_go)
kegg_df <- as.data.frame(gsea_kegg)

# Cargar resultados EPIC
cib_obj   <- readRDS(file.path(MOUT, "EPIC_objeto.rds"))
cib_res   <- cib_obj$cib
frac_cols <- cib_obj$frac_cols

# Cox univariado y multivariado
library(survival)

cox_uni     <- coxph(Surv(OS_time, OS_status) ~ Disulfidptosis_Canonical,
                     data = surv_df)
cox_uni_sum <- summary(cox_uni)

surv_df2 <- surv_df %>%
  mutate(
    Stage_simple = factor(
      ifelse(grepl("Stage I$|Stage II$|Stage IIA$|Stage IIB$",
                   ajcc_pathologic_stage, ignore.case = TRUE),
             "Early", "Advanced"),
      levels = c("Advanced", "Early")),
    Age_group = factor(
      ifelse(age_at_index >= 60, ">=60", "<60"),
      levels = c("<60", ">=60"))
  )

cox_multi     <- coxph(
  Surv(OS_time, OS_status) ~
    Disulfidptosis_Canonical + Stage_simple + Age_group + gender,
  data = surv_df2
)
cox_multi_sum <- summary(cox_multi)

hr_uni  <- round(cox_uni_sum$conf.int[1, "exp(coef)"],   2)
ci_low  <- round(cox_uni_sum$conf.int[1, "lower .95"],    2)
ci_high <- round(cox_uni_sum$conf.int[1, "upper .95"],    2)
p_cox   <- cox_uni_sum$coefficients[1, "Pr(>|z|)"]
lr_test <- survdiff(Surv(OS_time, OS_status) ~ Score_group, data = surv_df)
p_lr    <- 1 - pchisq(lr_test$chisq, df = 1)

cat(sprintf(
  "\n[RESULTADO PRINCIPAL]\n  HR = %.2f (IC95%%: %.2f–%.2f)\n  p Cox = %.2e | p Log-rank = %.2e\n\n",
  hr_uni, ci_low, ci_high, p_cox, p_lr))

# Tabla Cox multivariado
cox_table <- as.data.frame(cox_multi_sum$conf.int) %>%
  rownames_to_column("Variable") %>%
  rename(HR = `exp(coef)`,
         HR_IC_inferior = `lower .95`,
         HR_IC_superior = `upper .95`) %>%
  mutate(p_valor = cox_multi_sum$coefficients[, "Pr(>|z|)"],
         across(where(is.numeric), ~ round(., 4))) %>%
  mutate(Variable = recode(Variable,
    "Disulfidptosis_Canonical" = "Score Disulfidptosis",
    "Stage_simpleEarly"        = "Estadio: Temprano vs Avanzado",
    "Age_group>=60"            = "Edad: ≥60 vs <60",
    "gendermale"               = "Sexo: Masculino vs Femenino"))

# Crear libro Excel
wb <- createWorkbook()

estilo_header <- createStyle(
  fontColour = "#FFFFFF",
  fgFill     = "#1167B1",
  halign     = "LEFT",
  fontName   = "Calibri",
  fontSize   = 11,
  textDecoration = "bold",
  border     = "Bottom",
  borderColour = "#03254C"
)

estilo_alt <- createStyle(fgFill = "#D0E8F5")

agregar_hoja <- function(wb, nombre, datos,
                         estilo_h = estilo_header,
                         estilo_a = estilo_alt) {
  addWorksheet(wb, nombre)
  writeData(wb, nombre, datos, headerStyle = estilo_h)
  # Filas alternas
  if (nrow(datos) > 1) {
    for (i in seq(2, nrow(datos), by = 2))
      addStyle(wb, nombre, style = estilo_a,
               rows = i + 1, cols = seq_len(ncol(datos)),
               gridExpand = TRUE)
  }
  setColWidths(wb, nombre, cols = seq_len(ncol(datos)),
               widths = "auto")
}

# S1 — Cohorte clínica
agregar_hoja(wb, "S1_Cohorte_Clinica",
  surv_df %>%
    dplyr::select(
      sample_id, patient_id, age_at_index, gender,
      ajcc_pathologic_stage, OS_time, OS_status,
      Disulfidptosis_Canonical, Disulfidptosis_Extended,
      PCA_score, Score_group, Cluster) %>%
    rename(
      Muestra = sample_id,
      Paciente = patient_id,
      Edad = age_at_index,
      Sexo = gender,
      Estadio_AJCC = ajcc_pathologic_stage,
      Tiempo_OS_dias = OS_time,
      Estado_OS = OS_status,
      Score_Canonico = Disulfidptosis_Canonical,
      Score_Extendido = Disulfidptosis_Extended,
      Score_PCA = PCA_score,
      Grupo_Score = Score_group
    ) %>%
    arrange(Cluster, Grupo_Score)
)

# S2 — Scores de disulfidptosis
agregar_hoja(wb, "S2_Scores_Disulfidptosis",
  surv_df %>%
    dplyr::select(sample_id, Disulfidptosis_Canonical,
                  Disulfidptosis_Extended, PCA_score,
                  Score_group, Cluster) %>%
    rename(
      Muestra = sample_id,
      Score_Canonico = Disulfidptosis_Canonical,
      Score_Extendido = Disulfidptosis_Extended,
      Score_PCA = PCA_score,
      Grupo_Score = Score_group
    )
)

# S3 — DEG completo
agregar_hoja(wb, "S3_DEG_Completo",
  deg_res %>%
    dplyr::select(gene, logFC, AveExpr, t, P.Value, adj.P.Val, B) %>%
    rename(Gen = gene,
           log2FC = logFC,
           Expresion_media = AveExpr,
           Estadistico_t = t,
           p_valor = P.Value,
           p_valor_ajustado_BH = adj.P.Val,
           Estadistico_B = B) %>%
    arrange(p_valor_ajustado_BH)
)

# S4 — DEG significativos
agregar_hoja(wb, "S4_DEG_Significativos",
  deg_res %>%
    filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
    dplyr::select(gene, logFC, AveExpr, t, P.Value, adj.P.Val, B) %>%
    rename(Gen = gene,
           log2FC = logFC,
           Expresion_media = AveExpr,
           Estadistico_t = t,
           p_valor = P.Value,
           p_valor_ajustado_BH = adj.P.Val,
           Estadistico_B = B) %>%
    arrange(p_valor_ajustado_BH)
)

# S5 — GSEA GO BP
agregar_hoja(wb, "S5_GSEA_GO_BiologicalProcess",
  go_df %>%
    dplyr::select(ID, Description, setSize, enrichmentScore,
                  NES, pvalue, p.adjust, qvalue, leading_edge) %>%
    rename(
      Termino_GO = ID,
      Descripcion = Description,
      Tamano_set = setSize,
      Score_enriquecimiento = enrichmentScore,
      NES_normalizado = NES,
      p_valor = pvalue,
      p_ajustado_BH = p.adjust,
      q_valor = qvalue,
      Genes_principales = leading_edge) %>%
    arrange(p_ajustado_BH)
)

# S6 — GSEA KEGG
agregar_hoja(wb, "S6_GSEA_KEGG_Rutas",
  kegg_df %>%
    dplyr::select(ID, Description, setSize, enrichmentScore,
                  NES, pvalue, p.adjust, qvalue, leading_edge) %>%
    rename(
      Ruta_KEGG = ID,
      Descripcion = Description,
      Tamano_set = setSize,
      Score_enriquecimiento = enrichmentScore,
      NES_normalizado = NES,
      p_valor = pvalue,
      p_ajustado_BH = p.adjust,
      q_valor = qvalue,
      Genes_principales = leading_edge) %>%
    arrange(p_ajustado_BH)
)

# S7 — ORA GO BP
if (file.exists(file.path(EOUT, "ORA_GO_BP.rds"))) {
  ora_go <- readRDS(file.path(EOUT, "ORA_GO_BP.rds"))
  agregar_hoja(wb, "S7_ORA_GO_BiologicalProcess",
    as.data.frame(ora_go) %>%
      dplyr::select(ID, Description, GeneRatio, BgRatio,
                    pvalue, p.adjust, qvalue, geneID, Count) %>%
      rename(
        Termino_GO = ID,
        Descripcion = Description,
        Razon_genes = GeneRatio,
        Razon_fondo = BgRatio,
        p_valor = pvalue,
        p_ajustado_BH = p.adjust,
        q_valor = qvalue,
        Genes = geneID,
        N_genes = Count) %>%
      arrange(p_ajustado_BH)
  )
}

# S8 — EPIC TME fracciones
agregar_hoja(wb, "S8_EPIC_TME_Fracciones",
  cib_res %>%
    left_join(
      surv_df[, c("sample_id", "Cluster", "Score_group")],
      by = setNames("sample_id",
                    names(cib_res)[1])   # primera columna = ID
    ) %>%
    dplyr::select(1, Cluster, Score_group,
                  all_of(frac_cols[frac_cols %in% colnames(cib_res)]))
)

# S9 — Cox multivariado
agregar_hoja(wb, "S9_Cox_Multivariado", cox_table)

xl_path <- file.path(TOUT,
  "Tablas_Suplementarias_Disulfidptosis_LIHC.xlsx")
saveWorkbook(wb, xl_path, overwrite = TRUE)
message("Guardé el Excel en: ", xl_path)
message("Terminé el bloque 7.")


message("Ahora genero las figuras para la publicación.")

for (p in c("ggplot2","ggrepel","patchwork","survminer",
            "stringr","scales","cluster"))
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")

library(ggplot2); library(ggrepel); library(patchwork)
library(survival); library(survminer)
library(dplyr);    library(tidyr);   library(stringr)
library(scales);   library(cluster)
library(clusterProfiler); library(enrichplot)
library(tibble)

# Cargar datos
surv_df   <- readRDS(file.path(KOUT, "Datos_supervivencia.rds"))
deg_res   <- read.csv(file.path(EOUT, "DEG_Alto_vs_Bajo_score.csv"),
                      stringsAsFactors = FALSE)
gsea_go   <- readRDS(file.path(EOUT, "GSEA_GO_BP.rds"))
gsea_kegg <- readRDS(file.path(EOUT, "GSEA_KEGG.rds"))
cib_obj   <- readRDS(file.path(MOUT, "EPIC_objeto.rds"))
cib_res   <- cib_obj$cib
frac_cols <- cib_obj$frac_cols
go_df     <- as.data.frame(gsea_go)
kegg_df   <- as.data.frame(gsea_kegg)

surv_df2 <- surv_df %>%
  mutate(
    Stage_simple = factor(
      ifelse(grepl("Stage I$|Stage II$|Stage IIA$|Stage IIB$",
                   ajcc_pathologic_stage, ignore.case = TRUE),
             "Early", "Advanced"),
      levels = c("Advanced", "Early")),
    Age_group = factor(
      ifelse(age_at_index >= 60, ">=60", "<60"),
      levels = c("<60", ">=60"))
  )

cox_multi     <- coxph(
  Surv(OS_time, OS_status) ~
    Disulfidptosis_Canonical + Stage_simple + Age_group + gender,
  data = surv_df2
)
cox_ms <- summary(cox_multi)


if (!requireNamespace("cowplot", quietly = TRUE))
  install.packages("cowplot")
library(cowplot)


# FIG 1 — KM por Score de Disulfidptosis
message("Generando la Fig 01, KM por Score...")
km_score   <- survfit(Surv(OS_time, OS_status) ~ Score_group,
                      data = surv_df)
p_km_score <- ggsurvplot(
  km_score, data = surv_df,
  palette        = c(COL_HIGH, COL_LOW),
  conf.int       = TRUE,
  conf.int.alpha = 0.12,
  pval           = TRUE,
  pval.size      = 4.5,
  pval.coord     = c(150, 0.10),
  risk.table     = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.26,
  risk.table.fontsize = 3.5,
  tables.theme   = theme_cleantable(),
  legend.labs    = c("High Score", "Low Score"),
  legend.title   = "Disulfidptosis Score",
  xlab           = "Time (days)",
  ylab           = "Overall Survival Probability",
  title          = "Overall Survival by Disulfidptosis Score",
  subtitle       = "TCGA-LIHC",  # se actualiza en la siguiente línea
  ggtheme        = theme_blue(),
  surv.median.line = "hv",
  size           = 1,
  censor.shape   = 124,
  censor.size    = 3,
  fontsize       = 3.8
)
p_km_score$plot <- p_km_score$plot +
  ggplot2::labs(subtitle = sprintf(
    "TCGA-LIHC (n = %d, events = %d)",
    nrow(surv_df), sum(surv_df$OS_status)))
p_km_score$table <- p_km_score$table +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = 10, face = "bold", color = BLUE_DARKEST))

km_score_grob <- cowplot::plot_grid(
  p_km_score$plot,
  p_km_score$table,
  ncol = 1,
  rel_heights = c(0.74, 0.26)
)
# fallback: guardar con print estándar si la conversión falla
tryCatch({
  guardar_figura(km_score_grob, "Fig01_KM_Score_Disulfidptosis",
                 ancho = 10, alto = 9)
}, error = function(e) {
  png(file.path(FOUT, "Fig01_KM_Score_Disulfidptosis.png"),
      width = 3000, height = 2700, res = 300)
  print(p_km_score)
  dev.off()
  pdf(file.path(FOUT, "Fig01_KM_Score_Disulfidptosis.pdf"),
      width = 10, height = 9)
  print(p_km_score)
  dev.off()
})

# FIG 2 — KM por Cluster
message("Generando la Fig 02, KM por Cluster...")
km_clust   <- survfit(Surv(OS_time, OS_status) ~ Cluster,
                      data = surv_df)
n_clusters <- nlevels(factor(surv_df$Cluster))
pal_cluster <- setNames(
  c(BLUE_DARKEST, BLUE_DARK, BLUE_LIGHT, BLUE_PALE)[seq_len(n_clusters)],
  paste0("C", seq_len(n_clusters))
)
labs_cluster <- paste0("Cluster C", seq_len(n_clusters))

p_km_clust <- ggsurvplot(
  km_clust, data = surv_df,
  palette        = pal_cluster,
  conf.int       = TRUE,
  conf.int.alpha = 0.12,
  pval           = TRUE,
  pval.size      = 4.5,
  pval.coord     = c(150, 0.10),
  risk.table     = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.26,
  risk.table.fontsize = 3.5,
  tables.theme   = theme_cleantable(),
  legend.labs    = labs_cluster,
  legend.title   = "Disulfidptosis Cluster",
  xlab           = "Time (days)",
  ylab           = "Overall Survival Probability",
  title          = "Overall Survival by Disulfidptosis Cluster",
  subtitle       = sprintf("TCGA-LIHC (n = %d)", nrow(surv_df)),
  ggtheme        = theme_blue(),
  surv.median.line = "hv",
  size           = 1,
  censor.shape   = 124,
  censor.size    = 3,
  fontsize       = 3.8
)
p_km_clust$table <- p_km_clust$table +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = 10, face = "bold", color = BLUE_DARKEST))

png(file.path(FOUT, "Fig02_KM_Cluster.png"),
    width = 3000, height = 2700, res = 300)
print(p_km_clust); dev.off()
pdf(file.path(FOUT, "Fig02_KM_Cluster.pdf"), width = 10, height = 9)
print(p_km_clust); dev.off()

# FIG 3 — Violin + Boxplot: Scores por Cluster
message("Generando la Fig 03, violín de scores por cluster...")
scores_long <- surv_df %>%
  dplyr::select(sample_id, Cluster,
                Disulfidptosis_Canonical,
                Disulfidptosis_Extended,
                PCA_score) %>%
  pivot_longer(
    cols      = c(Disulfidptosis_Canonical,
                  Disulfidptosis_Extended, PCA_score),
    names_to  = "Score_type",
    values_to = "Value") %>%
  mutate(Score_type = recode(Score_type,
    "Disulfidptosis_Canonical" = "Canonical\nScore",
    "Disulfidptosis_Extended"  = "Extended\nScore",
    "PCA_score"                = "PCA\nScore"))

pval_df <- scores_long %>%
  group_by(Score_type) %>%
  summarise(
    p     = tryCatch(wilcox.test(Value ~ Cluster)$p.value,
                     error = function(e) NA_real_),
    y_pos = max(Value, na.rm = TRUE) * 1.08,
    .groups = "drop") %>%
  mutate(label = dplyr::case_when(
    is.na(p)    ~ "n.s.",
    p < 0.001   ~ "p < 0.001",
    p < 0.01    ~ sprintf("p = %.3f", p),
    p < 0.05    ~ sprintf("p = %.3f", p),
    TRUE        ~ "n.s."))

p_boxplot <- ggplot(scores_long,
                    aes(x = Cluster, y = Value, fill = Cluster)) +
  geom_violin(alpha = 0.65, trim = TRUE,
              linewidth = 0.25, color = "white") +
  geom_boxplot(width = 0.16, outlier.shape = 21,
               outlier.size = 1.0, outlier.alpha = 0.4,
               linewidth = 0.5, color = "grey25", fill = "white") +
  geom_text(data = pval_df,
            aes(x = 1.5, y = y_pos, label = label),
            inherit.aes = FALSE, size = 3.2,
            color = "grey30", fontface = "italic") +
  scale_fill_manual(
    values = setNames(c(BLUE_DARK, BLUE_LIGHT, BLUE_MID, BLUE_PALE),
                      paste0("C", 1:4)),
    guide = "none") +
  facet_wrap(~ Score_type, scales = "free_y", nrow = 1) +
  labs(
    title   = "Disulfidptosis Scores by Cluster",
    subtitle = "Wilcoxon test | TCGA-LIHC Primary Tumors",
    x       = "Cluster",
    y       = "Score (z-normalized)",
    caption = "Violin = distribution density; box = IQR; whiskers = 1.5×IQR"
  ) +
  theme_blue(base_size = 12)

guardar_figura(p_boxplot, "Fig03_Violin_Scores_por_Cluster",
               ancho = 12, alto = 7)

# FIG 4 — PCA Clusters
message("Generando la Fig 04, PCA de clusters...")
pca_data <- surv_df %>%
  dplyr::select(sample_id, Cluster,
                Disulfidptosis_Canonical,
                Disulfidptosis_Extended, PCA_score) %>%
  filter(complete.cases(.))

pca_mat  <- scale(pca_data[, c("Disulfidptosis_Canonical",
                                "Disulfidptosis_Extended",
                                "PCA_score")])
pca_res  <- prcomp(pca_mat, center = FALSE, scale. = FALSE)
pca_df   <- as.data.frame(pca_res$x[, 1:2])
pca_df$Cluster <- pca_data$Cluster
var_exp  <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster)) +
  stat_ellipse(aes(fill = Cluster), geom = "polygon",
               alpha = 0.07, level = 0.95,
               linewidth = 0.9, linetype = "solid",
               show.legend = FALSE) +
  geom_point(size = 2.2, alpha = 0.78, stroke = 0) +
  scale_color_manual(
    values = setNames(c(BLUE_DARKEST, BLUE_LIGHT, BLUE_MID, BLUE_PALE),
                      paste0("C", 1:4)),
    labels = paste0("C", 1:n_clusters,
                    c(" (High Score)", " (Low Score)",
                      "", "")[1:n_clusters])) +
  scale_fill_manual(
    values = setNames(c(BLUE_DARKEST, BLUE_LIGHT, BLUE_MID, BLUE_PALE),
                      paste0("C", 1:4)),
    guide = "none") +
  labs(
    title   = "PCA — Disulfidptosis Clusters",
    subtitle = "95% confidence ellipses | TCGA-LIHC Primary Tumors",
    x       = sprintf("PC1 (%s%%)", var_exp[1]),
    y       = sprintf("PC2 (%s%%)", var_exp[2]),
    color   = "Cluster"
  ) +
  theme_blue(base_size = 12) +
  theme(legend.position = c(0.88, 0.12))

guardar_figura(p_pca, "Fig04_PCA_Clusters", ancho = 9, alto = 8)

# FIG 5 — Silhouette Width
message("Generando la Fig 05, silhouette...")
sil_df <- readRDS(file.path(SOUT, "Silhouette_widths.rds"))

p_sil <- ggplot(sil_df, aes(x = k, y = sil)) +
  geom_area(fill = BLUE_PALE, alpha = 0.5) +
  geom_line(color = BLUE_DARK, linewidth = 1.1) +
  geom_point(aes(color = k == sil_df$k[which.max(sil_df$sil)]),
             size = 4.5, shape = 21,
             fill = "white", stroke = 1.8) +
  scale_color_manual(values = c("FALSE" = BLUE_LIGHT,
                                "TRUE"  = BLUE_DARKEST),
                     guide = "none") +
  geom_vline(xintercept = sil_df$k[which.max(sil_df$sil)],
             linetype = "dashed", color = BLUE_MID,
             linewidth = 0.7) +
  geom_label(
    data = sil_df[which.max(sil_df$sil), ],
    aes(label = sprintf("k = %d\n(optimal)", k)),
    hjust = -0.15, vjust = 0.5, size = 3.5,
    color = BLUE_DARKEST, fill = "white",
    label.size = 0.35, label.padding = unit(0.25, "lines")) +
  scale_x_continuous(breaks = 2:6) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  labs(
    title   = "Optimal Number of Clusters — Silhouette Analysis",
    subtitle = "k-means clustering | Disulfidptosis Score | TCGA-LIHC",
    x       = "Number of Clusters (k)",
    y       = "Average Silhouette Width",
    caption = "Higher silhouette width indicates better-defined clusters"
  ) +
  theme_blue(base_size = 12)

guardar_figura(p_sil, "Fig05_Silhouette", ancho = 9, alto = 7)

# FIG 6 — Forest Plot Cox Multivariado
message("Generando la Fig 06, forest plot de Cox...")
n_tab <- list(
  Disulfidptosis_Canonical = nrow(surv_df2),
  Stage_simpleEarly        = sum(surv_df2$Stage_simple == "Early",  na.rm = TRUE),
  `Age_group>=60`          = sum(surv_df2$Age_group    == ">=60",   na.rm = TRUE),
  gendermale               = sum(surv_df2$gender        == "male",  na.rm = TRUE)
)

forest_df <- as.data.frame(cox_ms$conf.int) %>%
  rownames_to_column("term") %>%
  rename(HR = `exp(coef)`, lower = `lower .95`, upper = `upper .95`) %>%
  mutate(
    p_val    = cox_ms$coefficients[, "Pr(>|z|)"],
    sig      = dplyr::case_when(
                 p_val < 0.001 ~ "***",
                 p_val < 0.01  ~ "**",
                 p_val < 0.05  ~ "*",
                 TRUE          ~ ""),
    n_label  = unlist(n_tab[term]),
    Variable = c("Disulfidptosis Score",
                 "Stage: Early vs Advanced",
                 "Age: ≥60 vs <60",
                 "Sex: Male vs Female"),
    ci_label = sprintf("%.2f (%.2f–%.2f)%s",
                       HR, lower, upper, sig),
    highlight = term == "Disulfidptosis_Canonical"
  )

ref_rows <- data.frame(
  term      = c("Stage_ref", "Age_ref", "gender_ref"),
  HR        = 1, lower = NA, upper = NA,
  p_val     = NA, sig = "",
  n_label   = c(sum(surv_df2$Stage_simple == "Advanced", na.rm = TRUE),
                sum(surv_df2$Age_group    == "<60",       na.rm = TRUE),
                sum(surv_df2$gender        == "female",   na.rm = TRUE)),
  Variable  = c("Stage: Advanced (ref)",
                "Age: <60 (ref)",
                "Sex: Female (ref)"),
  ci_label  = "reference",
  highlight = FALSE,
  stringsAsFactors = FALSE
)

forest_full <- bind_rows(
  forest_df[1, ], ref_rows[1, ], forest_df[2, ],
  ref_rows[2, ],  forest_df[3, ], ref_rows[3, ],
  forest_df[4, ]
) %>%
  mutate(
    row_id = rev(seq_len(n())),
    is_ref = ci_label == "reference"
  )

x_max <- max(forest_full$upper, na.rm = TRUE) * 2.6

p_forest <- ggplot(forest_full, aes(y = row_id)) +
  # Filas de fondo alternas
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin  = row_id - 0.5,
                ymax  = row_id + 0.5,
                fill  = as.factor(row_id %% 2)),
            alpha = 0.15, show.legend = FALSE) +
  scale_fill_manual(values = c("0" = "white", "1" = BLUE_PALE)) +
  # Línea de referencia HR = 1
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey55", linewidth = 0.65) +
  # IC 95%
  geom_errorbarh(
    data = dplyr::filter(forest_full, !is_ref, !is.na(lower)),
    aes(xmin = lower, xmax = upper, color = highlight),
    height = 0.28, linewidth = 1.0) +
  # Puntos HR
  geom_point(
    data = dplyr::filter(forest_full, !is_ref),
    aes(x = HR, color = highlight,
        size  = highlight, shape = highlight)) +
  # Marcador de referencia
  geom_point(
    data = dplyr::filter(forest_full, is_ref),
    aes(x = 1), color = "grey55", size = 2.8, shape = 18) +
  # Etiquetas de variables (izquierda)
  geom_text(
    aes(x     = -0.06 * x_max,
        label = Variable,
        fontface = ifelse(highlight, "bold", "plain"),
        color   = highlight),
    hjust = 1, size = 3.5) +
  # N por fila
  geom_text(
    aes(x     = x_max * 0.64,
        label = ifelse(!is.na(n_label), paste0("N=", n_label), "")),
    hjust = 0.5, size = 3.0, color = "grey40") +
  # HR (IC) texto
  geom_text(
    aes(x    = x_max * 0.79,
        label = ci_label,
        fontface = ifelse(highlight, "bold", "plain")),
    hjust = 0, size = 3.1, color = "grey20") +
  scale_color_manual(
    values = c("TRUE" = BLUE_DARKEST, "FALSE" = BLUE_DARK),
    guide  = "none") +
  scale_size_manual(
    values = c("TRUE" = 5.5, "FALSE" = 3.5),
    guide  = "none") +
  scale_shape_manual(
    values = c("TRUE" = 18, "FALSE" = 15),
    guide  = "none") +
  scale_x_continuous(
    limits = c(-x_max * 0.42, x_max),
    breaks = c(0.25, 0.5, 1, 2, 3, 5),
    expand = c(0, 0)) +
  # Encabezados de columnas
  annotate("text", x = x_max * 0.64,
           y = max(forest_full$row_id) + 0.8,
           label = "N", hjust = 0.5, size = 3.7,
           fontface = "bold", color = BLUE_DARKEST) +
  annotate("text", x = x_max * 0.79,
           y = max(forest_full$row_id) + 0.8,
           label = "HR (95% CI)", hjust = 0, size = 3.7,
           fontface = "bold", color = BLUE_DARKEST) +
  labs(
    title   = "Multivariable Cox Regression — TCGA-LIHC",
    subtitle = sprintf("Events: %d  |  C-index: %.3f  |  %s",
                       sum(surv_df$OS_status, na.rm = TRUE),
                       cox_ms$concordance["C"],
                       ifelse(any(forest_df$p_val < 0.05, na.rm = TRUE),
                              "★ Significant variables", "")),
    x       = "Hazard Ratio (95% CI)",
    y       = NULL
  ) +
  theme_blue(base_size = 11) +
  theme(
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.line.y      = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin      = margin(14, 18, 12, 65)
  )

guardar_figura(p_forest, "Fig06_ForestPlot_Cox",
               ancho = 14, alto = 7)

# FIG 7 — TME Violin por Score_group
message("Generando la Fig 07, violín de TME por score...")

excluir <- c("sample", "Tumor", "Uncharacterized", "RMSE", "Correlation",
             "tumor.purity", "otherCells")  # otherCells = ~68% → excluir de figuras principales
frac_ok <- colnames(cib_res)[!colnames(cib_res) %in% excluir &
                               sapply(cib_res, is.numeric)]
top8_cells <- cib_res %>%
  dplyr::select(all_of(frac_ok)) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(),
               names_to  = "cell",
               values_to = "mean") %>%
  slice_max(mean, n = 8) %>%
  pull(cell)

id_col <- colnames(cib_res)[1]   # primera columna es el ID

cib_score <- cib_res %>%
  dplyr::select(all_of(c(id_col, top8_cells))) %>%
  left_join(
    surv_df[, c("sample_id", "Score_group")],
    by = setNames("sample_id", id_col)) %>%
  filter(!is.na(Score_group)) %>%
  pivot_longer(-c(all_of(id_col), Score_group),
               names_to  = "CellType",
               values_to = "Fraction") %>%
  mutate(CellType = str_to_title(
    str_replace_all(CellType, "_", " ")))

p_tme_score <- ggplot(cib_score,
                      aes(x = Score_group, y = Fraction,
                          fill = Score_group)) +
  geom_violin(alpha = 0.70, trim = TRUE,
              linewidth = 0.25, color = "white") +
  geom_boxplot(width = 0.14, outlier.shape = NA,
               fill = "white", linewidth = 0.45,
               color = "grey25") +
  scale_fill_manual(
    values = c("High" = COL_HIGH, "Low" = COL_LOW),
    name   = "Score Group") +
  facet_wrap(~ CellType, scales = "free_y", nrow = 2) +
  labs(
    title    = "Tumor Microenvironment by Disulfidptosis Score",
    subtitle = "High vs Low — EPIC Deconvolution (IOBR) | TCGA-LIHC",
    x        = "Score Group",
    y        = "Cell Fraction (normalized)",
    caption  = "EPIC cell-fraction estimates; fractions normalized per sample (sum = 1)"
  ) +
  theme_blue(base_size = 10) +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 9))

guardar_figura(p_tme_score, "Fig07_TME_Violin_Score",
               ancho = 16, alto = 8)

# FIG 8 — TME Violin por Cluster
message("Generando la Fig 08, violín de TME por cluster...")
cib_clust <- cib_res %>%
  dplyr::select(all_of(c(id_col, top8_cells))) %>%
  left_join(
    surv_df[, c("sample_id", "Cluster")],
    by = setNames("sample_id", id_col)) %>%
  filter(!is.na(Cluster)) %>%
  pivot_longer(-c(all_of(id_col), Cluster),
               names_to  = "CellType",
               values_to = "Fraction") %>%
  mutate(CellType = str_to_title(
    str_replace_all(CellType, "_", " ")))

p_tme_clust <- ggplot(cib_clust,
                      aes(x = Cluster, y = Fraction, fill = Cluster)) +
  geom_violin(alpha = 0.70, trim = TRUE,
              linewidth = 0.25, color = "white") +
  geom_boxplot(width = 0.14, outlier.shape = NA,
               fill = "white", linewidth = 0.45, color = "grey25") +
  scale_fill_manual(
    values = setNames(c(BLUE_DARKEST, BLUE_LIGHT, BLUE_MID, BLUE_PALE),
                      paste0("C", 1:4)),
    name = "Cluster") +
  facet_wrap(~ CellType, scales = "free_y", nrow = 2) +
  labs(
    title    = "Tumor Microenvironment by Disulfidptosis Cluster",
    subtitle = "EPIC Deconvolution (IOBR) | TCGA-LIHC",
    x        = "Cluster",
    y        = "Cell Fraction (normalized)"
  ) +
  theme_blue(base_size = 10) +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 9))

guardar_figura(p_tme_clust, "Fig08_TME_Violin_Cluster",
               ancho = 16, alto = 8)

# FIG 9 — Correlación TME × Score (barplot Spearman)
message("Generando la Fig 09, correlación de TME...")
cib_wide <- cib_res %>%
  dplyr::select(all_of(c(id_col, frac_ok))) %>%
  left_join(
    surv_df[, c("sample_id", "Disulfidptosis_Canonical")],
    by = setNames("sample_id", id_col)) %>%
  filter(!is.na(Disulfidptosis_Canonical))

cor_df <- lapply(frac_ok, function(cell) {
  x  <- cib_wide[[cell]]
  y  <- cib_wide$Disulfidptosis_Canonical
  ok <- complete.cases(x, y)
  if (sum(ok) < 10) return(NULL)
  ct <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  clean_name <- str_to_title(str_replace_all(cell, "_", " "))
  data.frame(
    CellType  = clean_name,
    rho       = as.numeric(ct$estimate),
    p_value   = ct$p.value,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows() %>%
  mutate(
    sig       = ifelse(p_value < 0.05, "p < 0.05", "n.s."),
    rho_label = ifelse(p_value < 0.05,
                       sprintf("%.2f", rho), "")
  ) %>%
  arrange(rho)

p_corr <- ggplot(cor_df,
                 aes(x    = rho,
                     y    = reorder(CellType, rho),
                     fill = rho,
                     alpha = sig)) +
  geom_col(width = 0.70, color = NA) +
  geom_vline(xintercept = 0,
             color = "grey40", linewidth = 0.55) +
  geom_text(aes(label = rho_label,
                hjust = ifelse(rho >= 0, -0.2, 1.2)),
            size = 2.9, color = "grey20") +
  scale_fill_gradient2(
    low      = BLUE_LIGHT,
    mid      = "white",
    high     = BLUE_DARKEST,
    midpoint = 0,
    name     = "Spearman ρ") +
  scale_alpha_manual(
    values = c("p < 0.05" = 1, "n.s." = 0.35),
    name   = "Significance") +
  labs(
    title    = "Spearman Correlation: Disulfidptosis Score vs Immune Cells",
    subtitle = "EPIC Deconvolution (IOBR) | TCGA-LIHC  |  ρ shown for p < 0.05",
    x        = "Spearman ρ",
    y        = NULL
  ) +
  xlim(min(cor_df$rho, na.rm = TRUE) * 1.35,
       max(cor_df$rho, na.rm = TRUE) * 1.35) +
  theme_blue(base_size = 11) +
  theme(legend.position = "right")

guardar_figura(p_corr, "Fig09_TME_Correlacion_Spearman",
               ancho = 13, alto = 10)

# FIG 10 — Volcano DEG
message("Generando la Fig 10, volcano de DEG...")
genes_canonical <- c("SLC7A11","SLC3A2","SLC2A1","LRPPRC","RNF213",
                     "NUBPL","OXSM","SHMT2","GYS1","FLII")

deg_plot <- deg_res %>%
  filter(!is.na(adj.P.Val), !is.na(logFC)) %>%
  mutate(
    neg_log10_p = -log10(adj.P.Val + 1e-300),
    status = dplyr::case_when(
      adj.P.Val < 0.05 & logFC >  1 ~ "Up",
      adj.P.Val < 0.05 & logFC < -1 ~ "Down",
      TRUE                           ~ "NS"),
    label = ifelse(
      status != "NS" &
        (neg_log10_p > quantile(neg_log10_p[status != "NS"],
                                0.88, na.rm = TRUE) |
           abs(logFC) > quantile(abs(logFC[status != "NS"]),
                                 0.88, na.rm = TRUE)),
      gene, NA_character_)
  )
# Forzar etiqueta de genes canónicos
deg_plot$label[deg_plot$gene %in% genes_canonical] <-
  deg_plot$gene[deg_plot$gene %in% genes_canonical]

n_up_v   <- sum(deg_plot$status == "Up",   na.rm = TRUE)
n_down_v <- sum(deg_plot$status == "Down", na.rm = TRUE)

p_volcano <- ggplot(
  deg_plot[order(deg_plot$status == "NS"), ],  # NS al fondo
  aes(x = logFC, y = neg_log10_p,
      color = status, size = status)) +
  geom_point(alpha = 0.55, stroke = 0) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey50",
             linewidth = 0.55) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed", color = "grey50",
             linewidth = 0.55) +
  geom_label_repel(
    aes(label = label),
    size           = 2.8,
    max.overlaps   = 25,
    box.padding    = 0.45,
    label.padding  = 0.2,
    label.size     = 0.2,
    fill           = "white",
    color          = "grey15",
    segment.color  = "grey60",
    segment.size   = 0.3,
    na.rm          = TRUE,
    seed           = 42) +
  scale_color_manual(
    values = c("Up" = BLUE_DARKEST,
               "Down" = BLUE_LIGHT,
               "NS"   = "grey82"),
    labels = c(
      "Up"   = sprintf("Up-regulated (n = %d)", n_up_v),
      "Down" = sprintf("Down-regulated (n = %d)", n_down_v),
      "NS"   = "Not significant"),
    name   = NULL) +
  scale_size_manual(
    values = c("Up" = 1.8, "Down" = 1.8, "NS" = 0.7),
    guide  = "none") +
  labs(
    title    = "Differential Gene Expression",
    subtitle = "High vs Low Disulfidptosis Score  |  FDR < 0.05, |log₂FC| > 1  |  TCGA-LIHC",
    x        = "log₂ Fold Change",
    y        = "-log₁₀ (adjusted p-value)",
    caption  = "Labels: top 12% by significance or fold change + canonical disulfidptosis genes"
  ) +
  theme_blue(base_size = 12) +
  theme(legend.position = c(0.80, 0.92),
        legend.background = element_rect(
          fill = "white", color = GREY_LINE, linewidth = 0.3))

guardar_figura(p_volcano, "Fig10_Volcano_DEG",
               ancho = 12, alto = 10)

# FIG 11 — GSEA GO BP Dotplot
message("Generando la Fig 11, dotplot de GSEA GO...")
top_go <- bind_rows(
  go_df %>% arrange(desc(NES)) %>% slice_head(n = 15),
  go_df %>% arrange(NES)       %>% slice_head(n = 15)
) %>%
  distinct(ID, .keep_all = TRUE) %>%
  mutate(
    Direction   = ifelse(NES > 0, "Activated", "Suppressed"),
    Description = str_wrap(Description, 48)
  )

p_go_dot <- ggplot(top_go,
                   aes(x = NES, y = reorder(Description, NES),
                       size = setSize, color = p.adjust)) +
  geom_point(alpha = 0.90) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.55) +
  scale_color_gradient(
    low   = BLUE_DARKEST,
    high  = BLUE_PALE,
    name  = "FDR",
    guide = guide_colorbar(reverse = TRUE)) +
  scale_size_continuous(range = c(2, 9),
                        name  = "Gene Set\nSize") +
  facet_grid(Direction ~ ., scales = "free_y", space = "free") +
  labs(
    title    = "GSEA — GO Biological Process",
    subtitle = "Top 15 Activated & Suppressed  |  High vs Low Score  |  TCGA-LIHC",
    x        = "Normalized Enrichment Score (NES)",
    y        = NULL
  ) +
  theme_blue(base_size = 10) +
  theme(axis.text.y   = element_text(size = 8),
        legend.position = "right")

guardar_figura(p_go_dot, "Fig11_GSEA_GO_BP_Dotplot",
               ancho = 16, alto = 14)

# FIG 12 — GSEA KEGG Dotplot
message("Generando la Fig 12, dotplot de GSEA KEGG...")
top_kegg <- bind_rows(
  kegg_df %>% arrange(desc(NES)) %>% slice_head(n = 15),
  kegg_df %>% arrange(NES)       %>% slice_head(n = 15)
) %>%
  distinct(ID, .keep_all = TRUE) %>%
  mutate(
    Direction   = ifelse(NES > 0, "Activated", "Suppressed"),
    Description = str_wrap(Description, 42)
  )

p_kegg_dot <- ggplot(top_kegg,
                     aes(x = NES,
                         y = reorder(Description, NES),
                         size = setSize,
                         color = p.adjust)) +
  geom_point(alpha = 0.90) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.55) +
  scale_color_gradient(
    low   = BLUE_DARKEST,
    high  = BLUE_PALE,
    name  = "FDR",
    guide = guide_colorbar(reverse = TRUE)) +
  scale_size_continuous(range = c(2, 9),
                        name  = "Gene Set\nSize") +
  facet_grid(Direction ~ ., scales = "free_y", space = "free") +
  labs(
    title    = "GSEA — KEGG Pathway Enrichment",
    subtitle = "Top 15 Activated & Suppressed  |  High vs Low Score  |  TCGA-LIHC",
    x        = "Normalized Enrichment Score (NES)",
    y        = NULL
  ) +
  theme_blue(base_size = 10) +
  theme(axis.text.y   = element_text(size = 8),
        legend.position = "right")

guardar_figura(p_kegg_dot, "Fig12_GSEA_KEGG_Dotplot",
               ancho = 16, alto = 14)

# FIG 13 — GSEA GO Enrichment Curves
message("Generando la Fig 13, curvas de enriquecimiento GSEA...")
top_go_ids <- c(
  go_df %>% arrange(desc(NES)) %>% slice_head(n = 3) %>% pull(ID),
  go_df %>% arrange(NES)       %>% slice_head(n = 3) %>% pull(ID)
)
# Filtrar IDs que existen en el objeto
top_go_ids <- top_go_ids[top_go_ids %in% gsea_go@result$ID]

blue_pal <- c(BLUE_DARKEST, BLUE_DARK, BLUE_LIGHT,
              "#74B9D4", "#A8D8EA", BLUE_MID)

tryCatch({
  p_gsea_curves <- gseaplot2(
    gsea_go,
    geneSetID   = top_go_ids,
    title       = "GSEA Enrichment — Top GO Biological Process Terms",
    pvalue_table = TRUE,
    color        = blue_pal[seq_along(top_go_ids)],
    base_size    = 11
  )
  guardar_figura(p_gsea_curves, "Fig13_GSEA_GO_Enrichment_Curves",
                 ancho = 14, alto = 10)
}, error = function(e) {
  warning("La Fig13 falló: ", conditionMessage(e))
})

# FIG 14 — ORA GO BP
if (file.exists(file.path(EOUT, "ORA_GO_BP.rds"))) {
  message("Generando la Fig 14, ORA GO BP...")
  ora_go <- readRDS(file.path(EOUT, "ORA_GO_BP.rds"))
  ora_df <- as.data.frame(ora_go) %>%
    slice_max(Count, n = 20) %>%
    mutate(Description = str_wrap(Description, 42))

  p_ora <- ggplot(ora_df,
                  aes(x    = Count,
                      y    = reorder(Description, Count),
                      fill = p.adjust)) +
    geom_col(width = 0.72, alpha = 0.92,
             color = "white", linewidth = 0.2) +
    scale_fill_gradient(
      low   = BLUE_DARKEST,
      high  = BLUE_PALE,
      name  = "FDR",
      guide = guide_colorbar(reverse = TRUE)) +
    labs(
      title    = "Over-Representation Analysis — GO Biological Process",
      subtitle = "Significant DEGs (FDR < 0.05, |log₂FC| > 1)  |  TCGA-LIHC",
      x        = "Gene Count",
      y        = NULL
    ) +
    theme_blue(base_size = 11) +
    theme(axis.text.y    = element_text(size = 9),
          legend.position = "right")

  guardar_figura(p_ora, "Fig14_ORA_GO_BP",
                 ancho = 14, alto = 11)
}

message("Terminé el bloque 8. Guardé las 14 figuras.")


message("Por último, corro la validación automática.")

PASS    <- function(msg) cat(sprintf("  OK: %s\n", msg))
FAIL    <- function(msg) cat(sprintf("  Falló: %s\n", msg))
SECTION <- function(t)   cat(sprintf("\n%s\n", t))

errores <- list()

# KM Score
SECTION("Kaplan-Meier por Score")
tryCatch({
  n_high <- sum(surv_df$Score_group == "High", na.rm = TRUE)
  n_low  <- sum(surv_df$Score_group == "Low",  na.rm = TRUE)
  lrt    <- survdiff(Surv(OS_time, OS_status) ~ Score_group,
                     data = surv_df)
  p_val  <- 1 - pchisq(lrt$chisq, df = 1)
  PASS(sprintf("High = %d, Low = %d", n_high, n_low))
  # Verificar tamaño mínimo de cohorte
  if (nrow(surv_df) >= 100)
    PASS(sprintf("Cohorte adecuada: %d pacientes en supervivencia", nrow(surv_df)))
  else
    FAIL(sprintf("Cohorte pequeña: solo %d pacientes (se esperan ≥ 100 en LIHC)", nrow(surv_df)))
  if (p_val < 0.05)
    PASS(sprintf("Log-rank p = %.2e (significativo)", p_val))
  else
    FAIL(sprintf("Log-rank p = %.4f (no significativo)", p_val))
}, error = function(e) {
  FAIL(e$message)
  errores[["KM_Score"]] <<- e$message
})

# Cox multivariado
SECTION("Cox Forest (multivariado)")
tryCatch({
  hr    <- cox_ms$conf.int["Disulfidptosis_Canonical", "exp(coef)"]
  p_cox <- cox_ms$coefficients["Disulfidptosis_Canonical", "Pr(>|z|)"]
  ci_l  <- cox_ms$conf.int["Disulfidptosis_Canonical", "lower .95"]
  ci_u  <- cox_ms$conf.int["Disulfidptosis_Canonical", "upper .95"]
  PASS(sprintf("HR = %.2f (%.2f–%.2f), p = %.2e",
               hr, ci_l, ci_u, p_cox))
  c_idx <- cox_ms$concordance["C"]
  PASS(sprintf("C-index = %.3f", c_idx))
  if (c_idx >= 0.55)
    PASS(sprintf("C-index ≥ 0.55 — discriminación adecuada (%.3f)", c_idx))
  else
    FAIL(sprintf("C-index bajo: %.3f (se esperan ≥ 0.55)", c_idx))
  if (hr > 1 & p_cox < 0.05)
    PASS("HR > 1 y p < 0.05 — score es factor de riesgo independiente")
  else
    FAIL(sprintf("HR = %.2f no significativo (p = %.3f)", hr, p_cox))
}, error = function(e) {
  FAIL(e$message)
  errores[["Cox"]] <<- e$message
})

# TME — nombres limpios
SECTION("TME: Limpieza de nombres de células")
tryCatch({
  dirty <- cor_df$CellType[grepl("_EPIC|EPIC",
                                  cor_df$CellType,
                                  ignore.case = FALSE)]
  if (length(dirty) == 0)
    PASS("Todos los nombres de tipos celulares están limpios")
  else
    FAIL(paste("Nombres sucios:", paste(dirty, collapse = ", ")))
}, error = function(e) {
  FAIL(e$message)
  errores[["TME_Nombres"]] <<- e$message
})

# PCA — dirección de clusters
SECTION("PCA: Orientación de clusters")
tryCatch({
  med_c1 <- median(surv_df$Disulfidptosis_Canonical[
                     surv_df$Cluster == "C1"], na.rm = TRUE)
  med_c2 <- median(surv_df$Disulfidptosis_Canonical[
                     surv_df$Cluster == "C2"], na.rm = TRUE)
  if (med_c1 > med_c2)
    PASS(sprintf("C1 (%.3f) > C2 (%.3f) — orientación correcta",
                 med_c1, med_c2))
  else
    FAIL(sprintf("C1 (%.3f) no mayor que C2 (%.3f)", med_c1, med_c2))
}, error = function(e) {
  FAIL(e$message)
  errores[["PCA_Dir"]] <<- e$message
})

# GSEA — resultados no vacíos
SECTION("GSEA: Resultados no vacíos")
tryCatch({
  n_go   <- nrow(as.data.frame(gsea_go))
  n_kegg <- nrow(as.data.frame(gsea_kegg))
  if (n_go > 0)
    PASS(sprintf("GSEA GO BP: %d términos significativos", n_go))
  else
    FAIL("GSEA GO BP devolvió 0 términos — revisar gene_list y parámetros")
  if (n_kegg > 0)
    PASS(sprintf("GSEA KEGG: %d rutas significativas", n_kegg))
  else
    FAIL("GSEA KEGG devolvió 0 rutas — revisar conversión SYMBOL→ENTREZID")
}, error = function(e) {
  FAIL(e$message)
  errores[["GSEA_Resultados"]] <<- e$message
})

# Verificar archivos de salida
SECTION("Archivos de salida")
archivos_esperados <- c(
  file.path(KOUT, "Datos_supervivencia.rds"),
  file.path(EOUT, "DEG_Alto_vs_Bajo_score.csv"),
  file.path(EOUT, "GSEA_GO_BP.rds"),
  file.path(EOUT, "GSEA_KEGG.rds"),
  file.path(TOUT, "Tablas_Suplementarias_Disulfidptosis_LIHC.xlsx")
)
for (f in archivos_esperados) {
  if (file.exists(f))
    PASS(basename(f))
  else
    FAIL(paste("Falta:", f))
}

# Figuras
n_figs <- length(list.files(FOUT, pattern = "\\.png$"))
PASS(sprintf("%d figuras PNG generadas en %s", n_figs, FOUT))

cat("\n")
if (length(errores) == 0) {
  message("Validación completa, no detecté problemas.")
} else {
  message("Detecté ", length(errores), " problema(s): ",
          paste(names(errores), collapse = ", "))
}

# Resumen final
cat("\n")
message("Terminé el análisis de disulfidptosis en TCGA-LIHC.")
message("Las figuras (PNG + PDF) quedaron en: ", FOUT)
message("  Fig01: KM Score Disulfidptosis")
message("  Fig02: KM Cluster")
message("  Fig03: Violin Scores por Cluster")
message("  Fig04: PCA Clusters")
message("  Fig05: Silhouette Width")
message("  Fig06: Forest Plot Cox Multivariado")
message("  Fig07: TME Violin por Score")
message("  Fig08: TME Violin por Cluster")
message("  Fig09: Correlación Spearman TME vs Score")
message("  Fig10: Volcano DEG")
message("  Fig11: GSEA GO BP Dotplot")
message("  Fig12: GSEA KEGG Dotplot")
message("  Fig13: GSEA GO Enrichment Curves")
message("  Fig14: ORA GO BP Barplot")
message("")
message("Las tablas Excel (9 hojas) quedaron en: ", TOUT)
message("  Tablas_Suplementarias_Disulfidptosis_LIHC.xlsx")
message("")
message("Datos intermedios:")
message("  ", DOUT, ": SE raw, log2FPKM, colData, rowData")
message("  ", SOUT, ": Scores, Clusters, Silhouette")
message("  ", KOUT, ": Datos de supervivencia")
message("  ", MOUT, ": resultados de la deconvolución EPIC")
message("  ", EOUT, ": DEG, GSEA, ORA")
