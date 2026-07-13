
OUT_DIR <- "resultados_fase4"


dirs_fase4 <- c(
  file.path(OUT_DIR, "01_geo_datos"),
  file.path(OUT_DIR, "02_scores_geo"),
  file.path(OUT_DIR, "03_supervivencia_geo"),
  file.path(OUT_DIR, "04_nomograma"),
  file.path(OUT_DIR, "05_checkpoints"),
  file.path(OUT_DIR, "06_figuras_fase4")
)
invisible(lapply(dirs_fase4, dir.create, recursive = TRUE, showWarnings = FALSE))

GEO_DIR  <- file.path(OUT_DIR, "01_geo_datos")
SGEO_DIR <- file.path(OUT_DIR, "02_scores_geo")
KGEO_DIR <- file.path(OUT_DIR, "03_supervivencia_geo")
NOMO_DIR <- file.path(OUT_DIR, "04_nomograma")
CHKP_DIR <- file.path(OUT_DIR, "05_checkpoints")
FIG4_DIR <- file.path(OUT_DIR, "06_figuras_fase4")

BLUE_DARKEST <- "#03254C"
BLUE_DARK    <- "#1167B1"
BLUE_MID     <- "#2A9D8F"
BLUE_LIGHT   <- "#74B9D4"
BLUE_PALE    <- "#D0E8F5"
ACCENT_CORAL <- "#E76F51"
GREY_LINE    <- "#DEE2E6"
GREY_TEXT    <- "#495057"

theme_blue <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) %+replace%
    ggplot2::theme(
      text              = ggplot2::element_text(family = "sans", color = GREY_TEXT),
      plot.title        = ggplot2::element_text(face = "bold", size = base_size + 3,
                            color = BLUE_DARKEST, hjust = 0,
                            margin = ggplot2::margin(b = 4)),
      plot.subtitle     = ggplot2::element_text(size = base_size - 0.5,
                            color = "grey50", hjust = 0,
                            margin = ggplot2::margin(b = 10)),
      plot.caption      = ggplot2::element_text(size = base_size - 2,
                            color = "grey65", hjust = 1,
                            margin = ggplot2::margin(t = 8)),
      axis.title        = ggplot2::element_text(face = "bold", size = base_size,
                            color = BLUE_DARKEST),
      axis.text         = ggplot2::element_text(size = base_size - 1,
                            color = GREY_TEXT),
      axis.line         = ggplot2::element_line(color = GREY_LINE, linewidth = 0.6),
      axis.ticks        = ggplot2::element_line(color = GREY_LINE, linewidth = 0.4),
      legend.title      = ggplot2::element_text(face = "bold", size = base_size - 1,
                            color = BLUE_DARKEST),
      legend.text       = ggplot2::element_text(size = base_size - 2,
                            color = GREY_TEXT),
      legend.key.size   = ggplot2::unit(0.45, "cm"),
      legend.background = ggplot2::element_rect(fill = "white", color = GREY_LINE,
                            linewidth = 0.3),
      panel.grid.major  = ggplot2::element_line(color = "#EEF2F7", linewidth = 0.35),
      panel.grid.minor  = ggplot2::element_blank(),
      plot.background   = ggplot2::element_rect(fill = "white", color = NA),
      panel.background  = ggplot2::element_rect(fill = "white", color = NA),
      strip.background  = ggplot2::element_rect(fill = BLUE_PALE, color = NA),
      strip.text        = ggplot2::element_text(face = "bold", color = BLUE_DARKEST,
                            size = base_size - 1),
      plot.margin       = ggplot2::margin(14, 18, 12, 14)
    )
}

guardar_figura <- function(plot, nombre, ancho = 12, alto = 9,
                           dpi = 300, dir = FIG4_DIR) {
  ruta_base <- file.path(dir, nombre)
  ggplot2::ggsave(paste0(ruta_base, ".png"), plot = plot,
                  width = ancho, height = alto, dpi = dpi, bg = "white")
  ggplot2::ggsave(paste0(ruta_base, ".pdf"), plot = plot,
                  width = ancho, height = alto,
                  device = grDevices::cairo_pdf)
  message("Guardé la figura: ", basename(ruta_base))
}

message("Listo, configuración de la Fase 4 cargada.")



message("Empiezo con la extracción de datos de TCGA.")

TCGA_KOUT <- "resultados/03_supervivencia"
TCGA_SOUT <- "resultados/02_scores_clusters"
TCGA_EOUT <- "resultados/05_expresion_diferencial"
TCGA_DOUT <- "resultados/01_datos_procesados"

surv_tcga <- readRDS(file.path(TCGA_KOUT, "Datos_supervivencia.rds"))

datos_tcga_nomo <- surv_tcga[, c(
  "sample_id", "patient_id",
  "Disulfidptosis_Canonical",
  "Score_group",
  "Cluster",
  "OS_time", "OS_status",
  "age_at_index",
  "ajcc_pathologic_stage",
  "Stage_simple"
)]

if ("alpha_fetoprotein_outcome_value" %in% colnames(surv_tcga)) {
  datos_tcga_nomo$AFP_log <- log2(
    as.numeric(surv_tcga$alpha_fetoprotein_outcome_value) + 1
  )
} else {
  datos_tcga_nomo$AFP_log <- NA_real_
  message("No encontré la variable de AFP en los datos clínicos de TCGA, así que la omito del nomograma.")
}

write.csv(datos_tcga_nomo,
          file.path(NOMO_DIR, "TCGA_datos_nomograma.csv"),
          row.names = FALSE)
message("Extraje los datos de TCGA para ", nrow(datos_tcga_nomo), " pacientes.")

expr_tcga <- readRDS(file.path(TCGA_DOUT, "TCGA_LIHC_log2fpkm_SYMBOL.rds"))
message("Cargué la matriz de TCGA: ", nrow(expr_tcga), " genes por ",
        ncol(expr_tcga), " muestras.")

message("Sigo ahora con la descarga de GEO.")

pkgs_geo <- c("GEOquery", "Biobase", "limma")
for (p in pkgs_geo)
  if (!requireNamespace(p, quietly = TRUE))
    BiocManager::install(p, ask = FALSE, update = FALSE)

for (p in c("dplyr", "data.table"))
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")

library(GEOquery)
library(Biobase)
library(limma)
library(dplyr)
library(data.table)

cargar_gse <- function(gse_id, dir = GEO_DIR) {
  rds_path <- file.path(dir, paste0(gse_id, ".rds"))
  if (file.exists(rds_path)) {
    message("Cargando ", gse_id, " desde la caché.")
    return(readRDS(rds_path))
  }
  message("Descargando ", gse_id, " desde GEO, esto puede tardar unos minutos.")
  gse <- getGEO(gse_id, destdir = dir, GSEMatrix = TRUE, AnnotGPL = TRUE)
  saveRDS(gse, rds_path)
  gse
}

gse14520_list <- cargar_gse("GSE14520")
gse14520 <- gse14520_list[[1]]

expr_14520 <- exprs(gse14520)
pheno_14520 <- pData(gse14520)
feat_14520  <- fData(gse14520)

sym_col_14520 <- grep("gene.symbol|gene symbol|symbol",
                       colnames(feat_14520),
                       ignore.case = TRUE, value = TRUE)[1]
message("En GSE14520 la columna de símbolo génico que detecté es '", sym_col_14520, "'.")

syms_14520 <- toupper(trimws(feat_14520[[sym_col_14520]]))
keep_14520  <- syms_14520 != "" & !is.na(syms_14520) & !grepl("^---", syms_14520)
expr_14520  <- expr_14520[keep_14520, ]
syms_14520  <- syms_14520[keep_14520]

expr_14520_sym <- limma::avereps(expr_14520, ID = syms_14520)
message("GSE14520 queda con ", nrow(expr_14520_sym), " genes únicos y ",
        ncol(expr_14520_sym), " muestras.")

saveRDS(expr_14520_sym, file.path(GEO_DIR, "GSE14520_expr_SYMBOL.rds"))
saveRDS(pheno_14520,    file.path(GEO_DIR, "GSE14520_phenoData.rds"))

gse76427_list <- cargar_gse("GSE76427")
gse76427 <- gse76427_list[[1]]

expr_76427  <- exprs(gse76427)
pheno_76427 <- pData(gse76427)
feat_76427  <- fData(gse76427)

sym_col_76427 <- grep("gene.symbol|gene symbol|symbol",
                       colnames(feat_76427),
                       ignore.case = TRUE, value = TRUE)[1]
message("  GSE76427: columna símbolo detectada = '", sym_col_76427, "'")

syms_76427 <- toupper(trimws(feat_76427[[sym_col_76427]]))
keep_76427  <- syms_76427 != "" & !is.na(syms_76427) & !grepl("^---", syms_76427)
expr_76427  <- expr_76427[keep_76427, ]
syms_76427  <- syms_76427[keep_76427]

expr_76427_sym <- limma::avereps(expr_76427, ID = syms_76427)
message("  GSE76427: ", nrow(expr_76427_sym), " genes únicos, ",
        ncol(expr_76427_sym), " muestras")

saveRDS(expr_76427_sym, file.path(GEO_DIR, "GSE76427_expr_SYMBOL.rds"))
saveRDS(pheno_76427,    file.path(GEO_DIR, "GSE76427_phenoData.rds"))

message("Termino con la descarga de GEO, ahora calculo el score de disulfidptosis en estos datasets.")








genes_canonical <- c("SLC7A11", "SLC3A2", "SLC2A1", "LRPPRC",
                     "RNF213",  "NUBPL",  "OXSM",   "SHMT2",
                     "GYS1",    "FLII")

calcular_score_geo <- function(expr_mat, dataset_name) {
  genes_ok <- genes_canonical[genes_canonical %in% rownames(expr_mat)]
  message("En ", dataset_name, " encontré ",
          length(genes_ok), " de ", length(genes_canonical), " genes del panel.")
  if (length(genes_ok) < 5)
    stop("Menos de 5 genes del panel en ", dataset_name,
         ". Revisar normalización o plataforma.")

  expr_z <- t(scale(t(expr_mat)))

  score <- colMeans(expr_z[genes_ok, , drop = FALSE], na.rm = TRUE)
  df <- data.frame(
    sample_id   = colnames(expr_mat),
    Score_Disulf = as.numeric(score),
    Score_group  = ifelse(score >= median(score, na.rm = TRUE), "High", "Low"),
    stringsAsFactors = FALSE
  )
  df
}

scores_14520 <- calcular_score_geo(expr_14520_sym, "GSE14520")
scores_76427 <- calcular_score_geo(expr_76427_sym, "GSE76427")

write.csv(scores_14520, file.path(SGEO_DIR, "GSE14520_scores.csv"), row.names = FALSE)
write.csv(scores_76427, file.path(SGEO_DIR, "GSE76427_scores.csv"), row.names = FALSE)

message("Guardé los scores de GEO, sigo con el análisis de supervivencia.")

for (p in c("survival", "survminer", "ggplot2"))
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")

library(survival)
library(survminer)
library(ggplot2)

extraer_os_geo <- function(pheno_df, dataset_name) {
  cn <- tolower(colnames(pheno_df))

  t_col <- grep("os.time|overall.survival.time|surv.*time|time.*survival|
                  days.*death|survival.*months|time.to.death|
                  os.*month|month.*survival",
                 cn, value = TRUE)[1]

  # Estado: buscar columnas de evento
  s_col <- grep("os.event|os.status|vital.*status|status.*os|
                  death.event|survival.*status|event.*os",
                 cn, value = TRUE)[1]

  if (is.na(t_col) || is.na(s_col)) {
    message("No pude detectar automáticamente las columnas de OS en ", dataset_name, ".")
    message("Estas son las columnas disponibles en phenoData:")
    print(colnames(pheno_df))
    message("Voy a ajustar manualmente las líneas correspondientes.")
    return(NULL)
  }

  orig_tcol <- colnames(pheno_df)[match(t_col, cn)]
  orig_scol <- colnames(pheno_df)[match(s_col, cn)]

  df <- data.frame(
    sample_id = rownames(pheno_df),
    OS_time   = suppressWarnings(as.numeric(pheno_df[[orig_tcol]])),
    OS_status = suppressWarnings(as.numeric(pheno_df[[orig_scol]])),
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$OS_time) & !is.na(df$OS_status) & df$OS_time > 0, ]
  message(dataset_name, " tiene ", nrow(df), " muestras con OS disponible.")
  df
}

os_14520 <- extraer_os_geo(pheno_14520, "GSE14520")
os_76427 <- extraer_os_geo(pheno_76427, "GSE76427")

os_76427 <- data.frame(
  sample_id = rownames(pheno_76427),
  OS_time   = as.numeric(pheno_76427[["duryears_os:ch1"]]) * 12,
  OS_status = as.numeric(pheno_76427[["event_os:ch1"]])
) %>% filter(!is.na(OS_time), !is.na(OS_status), OS_time > 0)
message("Con las columnas correctas identificadas, GSE76427 queda con ", nrow(os_76427), " muestras con OS.")

construir_surv_geo <- function(scores_df, os_df, pheno_df, dataset_name) {
  if (is.null(os_df)) {
    message("No tengo datos de OS para ", dataset_name, ", hay que revisarlo manualmente.")
    return(NULL)
  }

  df <- inner_join(scores_df, os_df, by = "sample_id")

  cn <- tolower(colnames(pheno_df))
  age_col <- grep("age|edad", cn, value = TRUE)[1]
  if (!is.na(age_col)) {
    orig <- colnames(pheno_df)[match(age_col, cn)]
    df$age <- suppressWarnings(as.numeric(pheno_df[df$sample_id, orig]))
  } else {
    df$age <- NA_real_
  }

  stage_col <- grep("stage|estadio|tnm", cn, value = TRUE)[1]
  if (!is.na(stage_col)) {
    orig <- colnames(pheno_df)[match(stage_col, cn)]
    df$stage <- pheno_df[df$sample_id, orig]
  } else {
    df$stage <- NA_character_
  }

  df$Dataset <- dataset_name
  df
}

surv_14520 <- construir_surv_geo(scores_14520, os_14520, pheno_14520, "GSE14520")
surv_76427 <- construir_surv_geo(scores_76427, os_76427, pheno_76427, "GSE76427")

if (!is.null(surv_14520))
  write.csv(surv_14520, file.path(KGEO_DIR, "GSE14520_supervivencia.csv"), row.names = FALSE)
if (!is.null(surv_76427))
  write.csv(surv_76427, file.path(KGEO_DIR, "GSE76427_supervivencia.csv"), row.names = FALSE)

km_geo <- function(surv_df, dataset_name) {
  if (is.null(surv_df) || nrow(surv_df) < 20) {
    message(dataset_name, " no tiene suficientes muestras para hacer el KM, lo salto.")
    return(invisible(NULL))
  }

  fit <- survfit(Surv(OS_time, OS_status) ~ Score_group, data = surv_df)
  lrt <- survdiff(Surv(OS_time, OS_status) ~ Score_group, data = surv_df)
  p_val <- round(1 - pchisq(lrt$chisq, df = 1), 4)
  p_lab  <- ifelse(p_val < 0.001, "p < 0.001",
                   paste0("p = ", format(p_val, digits = 3)))

  p_km <- ggsurvplot(
    fit,
    data         = surv_df,
    palette      = c(BLUE_DARK, BLUE_LIGHT),
    risk.table   = TRUE,
    pval         = p_lab,
    pval.coord   = c(0, 0.08),
    conf.int     = FALSE,
    xlab         = "Tiempo (meses)",
    ylab         = "Probabilidad de supervivencia global",
    legend.labs  = c("High Score", "Low Score"),
    legend.title = "Score Disulfidptosis",
    title        = paste0("Curvas KM — ", dataset_name),
    ggtheme      = theme_blue(base_size = 12),
    fontsize     = 4,
    risk.table.col = "strata",
    surv.median.line = "hv",
    tables.theme = theme_cleantable()
  )

png(file.path(FIG4_DIR, paste0("KM_", dataset_name, ".png")),
    width = 3000, height = 2600, res = 300)
print(p_km)
dev.off()

pdf_path <- file.path(FIG4_DIR, paste0("KM_", dataset_name, ".pdf"))
grDevices::cairo_pdf(pdf_path, width = 10, height = 8.7)
print(p_km)
dev.off()

message("Guardé la curva KM de ", dataset_name, " en PNG y PDF.")

  cox1 <- coxph(Surv(OS_time, OS_status) ~ Score_Disulf, data = surv_df)
  message("Cox univariado de ", dataset_name, ":")
  print(summary(cox1)$conf.int)
}

km_geo(surv_14520, "GSE14520")
km_geo(surv_76427, "GSE76427")

message("Termino con las curvas de supervivencia en GEO, sigo con el nomograma clínico.")

if (!requireNamespace("rms", quietly = TRUE))
  install.packages("rms", repos = "https://cloud.r-project.org")

library(rms)
library(ggplot2)

nomo_df <- read.csv(file.path(NOMO_DIR, "TCGA_datos_nomograma.csv"),
                    stringsAsFactors = FALSE) %>%
  filter(OS_time > 0, !is.na(OS_time), !is.na(OS_status),
         !is.na(Disulfidptosis_Canonical),
         !is.na(age_at_index)) %>%
  mutate(
    Stage_bin = case_when(
      Stage_simple == "Advanced" ~ 1,
      Stage_simple == "Early"    ~ 0,
      TRUE ~ NA_real_
    ),
    Age_cont = as.numeric(age_at_index),
    Score    = Disulfidptosis_Canonical
  )

if (all(is.na(nomo_df$AFP_log))) {
  message("AFP no está disponible en TCGA, así que construyo el nomograma sin ella.")
  vars_nomo <- c("Score", "Stage_bin", "Age_cont")
} else {
  nomo_df$AFP_log[is.na(nomo_df$AFP_log)] <- median(nomo_df$AFP_log, na.rm = TRUE)
  vars_nomo <- c("Score", "Stage_bin", "AFP_log", "Age_cont")
}

nomo_df <- nomo_df %>% filter(!is.na(Stage_bin))
message("Me quedan ", nrow(nomo_df), " pacientes con estadio disponible.")
dd <- datadist(nomo_df[, vars_nomo])
options(datadist = "dd")

formula_nomo <- as.formula(
  paste0("Surv(OS_time, OS_status) ~ ",
         paste(vars_nomo, collapse = " + "))
)
cox_nomo <- cph(formula_nomo, data = nomo_df,
                x = TRUE, y = TRUE, surv = TRUE)

message("El C-index de concordancia del modelo es ",
        round(cox_nomo$stats["C"], 3), ".")

cal_1yr <- calibrate(cox_nomo, cmethod = "KM", method = "boot",
                     u = 365, B = 100)
cal_3yr <- calibrate(cox_nomo, cmethod = "KM", method = "boot",
                     u = 365 * 3, B = 100)

png(file.path(FIG4_DIR, "Calibracion_Nomograma.png"),
    width = 2800, height = 2400, res = 300)
par(mfrow = c(1, 2),
    col.axis = BLUE_DARKEST, col.lab = BLUE_DARK,
    font.lab = 2, las = 1, bty = "l")
plot(cal_1yr, main = "Calibración 1 año",
     col  = BLUE_DARK,
     lwd  = 2,
     subtitles = FALSE)
plot(cal_3yr, main = "Calibración 3 años",
     col  = BLUE_DARK,
     lwd  = 2,
     subtitles = FALSE)
dev.off()

grDevices::cairo_pdf(file.path(FIG4_DIR, "Calibracion_Nomograma.pdf"),
                     width = 9.3, height = 8)
par(mfrow = c(1, 2),
    col.axis = BLUE_DARKEST, col.lab = BLUE_DARK,
    font.lab = 2, las = 1, bty = "l")
plot(cal_1yr, main = "Calibración 1 año",  col = BLUE_DARK, lwd = 2, subtitles = FALSE)
plot(cal_3yr, main = "Calibración 3 años", col = BLUE_DARK, lwd = 2, subtitles = FALSE)
dev.off()

message("Guardé las curvas de calibración.")

surv_base <- survest(cox_nomo, times = c(365, 365*3), what = "survival")
S1 <- surv_base$surv[1]
S3 <- surv_base$surv[2]
message("S0 a 1 año = ", round(S1, 4), ", S0 a 3 años = ", round(S3, 4), ".")

fun_1yr <- function(lp) 1 - S1^exp(lp)
fun_3yr <- function(lp) 1 - S3^exp(lp)

nomo_obj <- nomogram(
  cox_nomo,
  fun      = list(fun_1yr, fun_3yr),
  fun.at   = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
  funlabel = c("Mort. 1 año", "Mort. 3 años"),
  lp       = FALSE
)

png(file.path(FIG4_DIR, "Nomograma_TCGA.png"),
    width = 3400, height = 2000, res = 300)
plot(nomo_obj, lmgp = 0.2, cex.axis = 0.75, col.grid = gray(0.9))
title(main = "Nomograma — Score Disulfidptosis · TCGA-LIHC",
      col.main = BLUE_DARKEST, font.main = 2, cex.main = 1.1)
dev.off()

grDevices::cairo_pdf(file.path(FIG4_DIR, "Nomograma_TCGA.pdf"),
                     width = 11.3, height = 6.7)
plot(nomo_obj, lmgp = 0.2, cex.axis = 0.75, col.grid = gray(0.9))
title(main = "Nomograma — Score Disulfidptosis · TCGA-LIHC",
      col.main = BLUE_DARKEST, font.main = 2, cex.main = 1.1)
dev.off()
message("Guardé el nomograma en PDF.")

message("Guardé el nomograma.")

cox_score_only <- cph(Surv(OS_time, OS_status) ~ Score,
                      data = nomo_df, x = TRUE, y = TRUE)
c_score <- round(cox_score_only$stats["C"], 3)
c_full  <- round(cox_nomo$stats["C"], 3)

cat(sprintf("C-index solo con el score: %.3f\n", c_score))
cat(sprintf("C-index con el modelo completo: %.3f\n", c_full))

message("Termino con el nomograma, sigo con los checkpoints inmunes y el TMB.")

library(ggplot2)
library(dplyr)

checkpoints <- c(
  "CD274",    # PD-L1
  "CTLA4",    # CTLA-4
  "PDCD1",    # PD-1
  "HAVCR2",   # TIM-3
  "LAG3",     # LAG-3
  "TIGIT",    # TIGIT
  "CD276",    # B7-H3
  "VTCN1"     # B7-H4
)

surv_tcga <- readRDS("resultados/03_supervivencia/Datos_supervivencia.rds")

correlacion_checkpoints <- function(expr_mat, surv_df_in,
                                    dataset_label, genes_cp = checkpoints) {

  genes_ok <- genes_cp[genes_cp %in% rownames(expr_mat)]

  muestras_comunes <- intersect(colnames(expr_mat), surv_df_in$sample_id)
  if (length(muestras_comunes) < 20) {
    message("Muy pocas muestras en común para ", dataset_label, ", lo salto.")
    return(NULL)
  }

  score_aligned <- surv_df_in$Disulfidptosis_Canonical[
    match(muestras_comunes, surv_df_in$sample_id)
  ]

  res <- lapply(genes_ok, function(g) {
    expr_g <- as.numeric(expr_mat[g, muestras_comunes])
    ct <- cor.test(score_aligned, expr_g, method = "spearman", exact = FALSE)
    data.frame(Gene     = g,
               rho      = round(ct$estimate, 3),
               p_value  = ct$p.value,
               Dataset  = dataset_label,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
}

cor_tcga <- correlacion_checkpoints(expr_tcga, surv_tcga, "TCGA-LIHC")

if (!is.null(surv_14520)) {
  expr_14520_mat <- readRDS(file.path(GEO_DIR, "GSE14520_expr_SYMBOL.rds"))
  surv_14520_cp <- surv_14520 %>%
    rename(Disulfidptosis_Canonical = Score_Disulf)
  cor_14520 <- correlacion_checkpoints(expr_14520_mat, surv_14520_cp, "GSE14520")
} else {
  cor_14520 <- NULL
}

if (!is.null(surv_76427)) {
  expr_76427_mat <- readRDS(file.path(GEO_DIR, "GSE76427_expr_SYMBOL.rds"))
  surv_76427_cp <- surv_76427 %>%
    rename(Disulfidptosis_Canonical = Score_Disulf)
  cor_76427 <- correlacion_checkpoints(expr_76427_mat, surv_76427_cp, "GSE76427")
} else {
  cor_76427 <- NULL
}

cor_all <- bind_rows(cor_tcga, cor_14520, cor_76427) %>%
  mutate(
    FDR        = p.adjust(p_value, method = "BH"),
    Significativo = ifelse(FDR < 0.05, "FDR < 0.05", "ns"),
    Etiqueta   = case_when(
      Gene == "CD274"  ~ "PD-L1 (CD274)",
      Gene == "CTLA4"  ~ "CTLA-4",
      Gene == "PDCD1"  ~ "PD-1",
      Gene == "HAVCR2" ~ "TIM-3",
      Gene == "LAG3"   ~ "LAG-3",
      Gene == "TIGIT"  ~ "TIGIT",
      Gene == "CD276"  ~ "B7-H3",
      Gene == "VTCN1"  ~ "B7-H4",
      TRUE ~ Gene
    )
  )

write.csv(cor_all,
          file.path(CHKP_DIR, "Checkpoints_correlacion_todos_datasets.csv"),
          row.names = FALSE)

library(tidyr)

heatmap_df <- cor_all %>%
  select(Etiqueta, Dataset, rho) %>%
  pivot_wider(names_from = Dataset, values_from = rho, values_fill = NA) %>%
  as.data.frame()

rownames(heatmap_df) <- heatmap_df$Etiqueta
heatmap_mat <- as.matrix(heatmap_df[, -1])

hm_long <- cor_all %>%
  select(Etiqueta, Dataset, rho, Significativo) %>%
  mutate(Etiqueta = factor(Etiqueta, levels = rev(sort(unique(Etiqueta)))))

p_heatmap <- ggplot(hm_long, aes(x = Dataset, y = Etiqueta, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = ifelse(Significativo == "FDR < 0.05", "*", "")),
    color = "white", size = 5, fontface = "bold"
  ) +
  scale_fill_gradient2(
    low      = BLUE_PALE,
    mid      = "white",
    high     = BLUE_DARKEST,
    midpoint = 0,
    limits   = c(-0.7, 0.7),
    name     = "Spearman ρ",
    guide    = guide_colorbar(reverse = FALSE,
                              barwidth = 1.2, barheight = 8)
  ) +
  labs(
    title    = "Correlación: Score Disulfidptosis vs Checkpoints Inmunes",
    subtitle = "TCGA-LIHC · GSE14520 · GSE76427  |  * FDR < 0.05",
    x        = NULL,
    y        = NULL
  ) +
  theme_blue(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 30, hjust = 1, face = "bold"),
    panel.grid   = element_blank(),
    legend.position = "right"
  )

guardar_figura(p_heatmap, "Fig_Checkpoints_Heatmap", ancho = 10, alto = 8)

checkpoints_key <- c("CD274", "CTLA4")
muestras_com <- intersect(colnames(expr_tcga), surv_tcga$sample_id)

dot_df <- data.frame(
  sample_id = muestras_com,
  Score     = surv_tcga$Disulfidptosis_Canonical[
                match(muestras_com, surv_tcga$sample_id)],
  PD_L1     = as.numeric(expr_tcga["CD274", muestras_com]),
  CTLA4     = as.numeric(expr_tcga["CTLA4", muestras_com]),
  Score_group = surv_tcga$Score_group[
                  match(muestras_com, surv_tcga$sample_id)]
)

p_pdl1 <- ggplot(dot_df, aes(x = Score, y = PD_L1, color = Score_group)) +
  geom_point(alpha = 0.55, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, color = BLUE_DARK,
              fill = BLUE_PALE, linewidth = 1.1) +
  scale_color_manual(values = c("High" = BLUE_DARK, "Low" = BLUE_LIGHT),
                     name = "Score") +
  labs(
    title    = "PD-L1 (CD274) vs Score Disulfidptosis",
    subtitle = "TCGA-LIHC  |  Correlación de Spearman",
    x        = "Score de Disulfidptosis",
    y        = "PD-L1 (log2 FPKM+1)"
  ) +
  annotate("text",
           x = quantile(dot_df$Score, 0.05, na.rm = TRUE),
           y = max(dot_df$PD_L1, na.rm = TRUE) * 0.95,
           label = paste0("ρ = ",
                          round(cor(dot_df$Score, dot_df$PD_L1,
                                    method = "spearman",
                                    use = "complete.obs"), 3)),
           color = BLUE_DARKEST, fontface = "bold", size = 4.5, hjust = 0) +
  theme_blue(base_size = 12)

guardar_figura(p_pdl1, "Fig_PD-L1_vs_Score_TCGA")

p_ctla4 <- ggplot(dot_df, aes(x = Score, y = CTLA4, color = Score_group)) +
  geom_point(alpha = 0.55, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, color = BLUE_DARK,
              fill = BLUE_PALE, linewidth = 1.1) +
  scale_color_manual(values = c("High" = BLUE_DARK, "Low" = BLUE_LIGHT),
                     name = "Score") +
  labs(
    title    = "CTLA-4 vs Score Disulfidptosis",
    subtitle = "TCGA-LIHC  |  Correlación de Spearman",
    x        = "Score de Disulfidptosis",
    y        = "CTLA-4 (log2 FPKM+1)"
  ) +
  annotate("text",
           x = quantile(dot_df$Score, 0.05, na.rm = TRUE),
           y = max(dot_df$CTLA4, na.rm = TRUE) * 0.95,
           label = paste0("ρ = ",
                          round(cor(dot_df$Score, dot_df$CTLA4,
                                    method = "spearman",
                                    use = "complete.obs"), 3)),
           color = BLUE_DARKEST, fontface = "bold", size = 4.5, hjust = 0) +
  theme_blue(base_size = 12)

guardar_figura(p_ctla4, "Fig_CTLA4_vs_Score_TCGA")

maf_path <- "resultados/01_datos_procesados/TCGA_LIHC.maf"

if (file.exists(maf_path)) {
  if (!requireNamespace("maftools", quietly = TRUE))
    BiocManager::install("maftools", ask = FALSE, update = FALSE)

  library(maftools)

  laml <- read.maf(maf   = maf_path,
                   clinicalData = write.csv(surv_tcga,
                                            tempfile(), row.names = FALSE))

  tmb_df <- tmb(laml, captureSize = 35.8) %>%
    mutate(Tumor_Sample_Barcode_12 = toupper(substr(Tumor_Sample_Barcode, 1, 12)))

  tmb_surv <- inner_join(
    tmb_df,
    surv_tcga %>% mutate(patient_id = toupper(substr(sample_id, 1, 12))),
    by = c("Tumor_Sample_Barcode_12" = "patient_id")
  ) %>% mutate(TMB_log2 = log2(total_perMB + 1))

  p_tmb <- ggplot(tmb_surv, aes(x = Score_group, y = TMB_log2,
                                  fill = Score_group)) +
    geom_boxplot(alpha = 0.85, outlier.shape = 21,
                 outlier.alpha = 0.5, width = 0.5) +
    geom_jitter(width = 0.12, alpha = 0.3, size = 1,
                color = GREY_TEXT) +
    scale_fill_manual(values = c("High" = BLUE_DARK, "Low" = BLUE_LIGHT),
                      guide = "none") +
    ggpubr::stat_compare_means(
      comparisons = list(c("High", "Low")),
      method      = "wilcox.test",
      label       = "p.signif",
      size        = 5,
      color       = BLUE_DARKEST
    ) +
    labs(
      title    = "Tumor Mutational Burden vs Score Disulfidptosis",
      subtitle = "TCGA-LIHC  |  Test Wilcoxon",
      x        = "Grupo Score",
      y        = "TMB (log₂ mutaciones/Mb)"
    ) +
    theme_blue(base_size = 12)

  guardar_figura(p_tmb, "Fig_TMB_vs_Score_TCGA")

} else {
  message("No encontré el MAF en: ", maf_path)
  message("Para calcular el TMB con maftools hace falta descargar el MAF de TCGA-LIHC con:")
  message("  query_maf <- GDCquery(project='TCGA-LIHC',")
  message("    data.category='Simple Nucleotide Variation',")
  message("    data.type='Masked Somatic Mutation')")
  message("  GDCdownload(query_maf, method='api')")
  message("  maf_data <- GDCprepare(query_maf)")
  message("Por ahora omito el bloque de TMB.")
}

message("Termino con los checkpoints inmunes y el TMB.")

n_figs4 <- length(list.files(FIG4_DIR, pattern = "\\.png$"))

message("Ya terminé la Fase 4, la validación externa en HCC.")
message("Generé ", n_figs4, " figuras en PNG y PDF, guardadas en ", FIG4_DIR, ".")
message("Los datos de GEO quedaron en ", GEO_DIR, ".")
message("Los scores de GEO quedaron en ", SGEO_DIR, ".")
message("La supervivencia de GEO quedó en ", KGEO_DIR, ".")
message("El nomograma quedó en ", NOMO_DIR, ".")
message("Los checkpoints quedaron en ", CHKP_DIR, ".")
