


OUT_DIR  <- "resultados"



DOUT <- file.path(OUT_DIR, "01_datos_procesados")
SOUT <- file.path(OUT_DIR, "02_scores_clusters")
KOUT <- file.path(OUT_DIR, "03_supervivencia")
MOUT <- file.path(OUT_DIR, "04_microambiente")
EOUT <- file.path(OUT_DIR, "05_expresion_diferencial")
FOUT <- file.path(OUT_DIR, "06_figuras")
TOUT <- file.path(OUT_DIR, "07_tablas_suplementarias")

fpkm_sym_path <- file.path(DOUT, "TCGA_LIHC_log2fpkm_SYMBOL.rds")

# Paleta de colores
BLUE_DARKEST <- "#03254C"
BLUE_DARK    <- "#1167B1"
BLUE_MID     <- "#2A9D8F"
BLUE_LIGHT   <- "#74B9D4"
BLUE_PALE    <- "#D0E8F5"
ACCENT_CORAL <- "#E76F51"
GREY_LINE    <- "#DEE2E6"
GREY_TEXT    <- "#495057"
COL_HIGH     <- BLUE_DARK
COL_LOW      <- BLUE_LIGHT

# Tema ggplot2
theme_blue <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) %+replace%
    ggplot2::theme(
      text              = ggplot2::element_text(family = "sans", color = GREY_TEXT),
      plot.title        = ggplot2::element_text(face = "bold", size = base_size + 3,
                            color = BLUE_DARKEST, hjust = 0, margin = ggplot2::margin(b = 4)),
      plot.subtitle     = ggplot2::element_text(size = base_size - 0.5, color = "grey50",
                            hjust = 0, margin = ggplot2::margin(b = 10)),
      plot.caption      = ggplot2::element_text(size = base_size - 2, color = "grey65",
                            hjust = 1, margin = ggplot2::margin(t = 8)),
      axis.title        = ggplot2::element_text(face = "bold", size = base_size, color = BLUE_DARKEST),
      axis.text         = ggplot2::element_text(size = base_size - 1, color = GREY_TEXT),
      axis.line         = ggplot2::element_line(color = GREY_LINE, linewidth = 0.6),
      axis.ticks        = ggplot2::element_line(color = GREY_LINE, linewidth = 0.4),
      legend.title      = ggplot2::element_text(face = "bold", size = base_size - 1, color = BLUE_DARKEST),
      legend.text       = ggplot2::element_text(size = base_size - 2, color = GREY_TEXT),
      legend.key.size   = ggplot2::unit(0.45, "cm"),
      legend.background = ggplot2::element_rect(fill = "white", color = GREY_LINE, linewidth = 0.3),
      legend.margin     = ggplot2::margin(4, 6, 4, 6),
      panel.grid.major  = ggplot2::element_line(color = "#EEF2F7", linewidth = 0.35),
      panel.grid.minor  = ggplot2::element_blank(),
      plot.background   = ggplot2::element_rect(fill = "white", color = NA),
      panel.background  = ggplot2::element_rect(fill = "white", color = NA),
      strip.background  = ggplot2::element_rect(fill = BLUE_PALE, color = NA),
      strip.text        = ggplot2::element_text(face = "bold", color = BLUE_DARKEST, size = base_size - 1),
      plot.margin       = ggplot2::margin(14, 18, 12, 14)
    )
}

guardar_figura <- function(plot, nombre, ancho = 12, alto = 9, dpi = 300, dir = FOUT) {
  ruta_base <- file.path(dir, nombre)
  ggplot2::ggsave(paste0(ruta_base, ".png"), plot = plot,
                  width = ancho, height = alto, dpi = dpi, bg = "white")
  ggplot2::ggsave(paste0(ruta_base, ".pdf"), plot = plot,
                  width = ancho, height = alto, device = grDevices::cairo_pdf)
  message("Guardado: ", basename(ruta_base))
}

message("Configuración cargada.")

message("Deconvolución con EPIC")

if (!requireNamespace("remotes",  quietly = TRUE))
  install.packages("remotes", repos = "https://cloud.r-project.org")
if (!requireNamespace("IOBR", quietly = TRUE))
  remotes::install_github("IOBR/IOBR", ref = "v0.99.9", upgrade = "never")

library(IOBR)
library(dplyr)
library(tidyr)

log2_fpkm <- readRDS(fpkm_sym_path)
surv_df   <- readRDS(file.path(KOUT, "Datos_supervivencia.rds"))

tumor_ids  <- surv_df$sample_id
expr_tumor <- log2_fpkm[, colnames(log2_fpkm) %in% tumor_ids, drop = FALSE]
message("Muestras de tumor para deconvolución: ", ncol(expr_tumor))

epic_path <- file.path(MOUT, "EPIC_normalizado.csv")

old_cache <- file.path(MOUT, "CIBERSORT_normalizado.csv")
if (file.exists(old_cache) && !file.exists(epic_path)) {
  message("Encuentro un caché anterior con nombre CIBERSORT_normalizado.csv.")
  message("Reejecuto con EPIC, tarda unos minutos.")
}

epic_res <- deconvo_tme(
  eset   = expr_tumor,
  method = "epic",
  arrays = FALSE
)

write.csv(epic_res, file.path(MOUT, "EPIC_bruto.csv"), row.names = FALSE)
message("EPIC bruto guardado.")

frac_cols <- grep("_EPIC$", colnames(epic_res), value = TRUE, ignore.case = FALSE)
message("Columnas EPIC detectadas: ", paste(frac_cols, collapse = ", "))

cib_renamed <- epic_res
new_names   <- gsub("_EPIC$", "", colnames(epic_res))
new_names   <- trimws(new_names)
colnames(cib_renamed) <- new_names

frac_cols_clean <- gsub("_EPIC$", "", frac_cols)

fc_idx   <- match(frac_cols_clean, colnames(cib_renamed))
fc_idx   <- fc_idx[!is.na(fc_idx)]
row_sums <- rowSums(cib_renamed[, fc_idx], na.rm = TRUE)
row_sums[row_sums == 0] <- 1
cib_renamed[, fc_idx] <- cib_renamed[, fc_idx] / row_sums

write.csv(cib_renamed, epic_path, row.names = FALSE)
saveRDS(
  list(cib = cib_renamed, frac_cols = frac_cols_clean[!is.na(fc_idx)]),
  file.path(MOUT, "EPIC_objeto.rds")
)

message("Columnas de fracciones limpias: ", paste(frac_cols_clean, collapse = ", "))
message("Verificación suma = 1 (primeras 3 muestras): ",
        paste(round(rowSums(cib_renamed[1:3, fc_idx]), 4), collapse = " | "))

rm(log2_fpkm, expr_tumor); gc()
message("Listo con la deconvolución.")

message("Armando las tablas suplementarias")

for (p in c("openxlsx", "tibble", "stringr"))
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")

library(openxlsx); library(dplyr); library(stringr); library(tibble)
library(survival)

surv_df   <- readRDS(file.path(KOUT, "Datos_supervivencia.rds"))
deg_res   <- read.csv(file.path(EOUT, "DEG_Alto_vs_Bajo_score.csv"), stringsAsFactors = FALSE)
gsea_go   <- readRDS(file.path(EOUT, "GSEA_GO_BP.rds"))
gsea_kegg <- readRDS(file.path(EOUT, "GSEA_KEGG.rds"))
go_df     <- as.data.frame(gsea_go)
kegg_df   <- as.data.frame(gsea_kegg)

cib_obj   <- readRDS(file.path(MOUT, "EPIC_objeto.rds"))   # <-- EPIC
cib_res   <- cib_obj$cib
frac_cols <- cib_obj$frac_cols

surv_df2 <- surv_df %>%
  mutate(
    Stage_simple = factor(
      ifelse(grepl("Stage I$|Stage II$|Stage IIA$|Stage IIB$",
                   ajcc_pathologic_stage, ignore.case = TRUE), "Early", "Advanced"),
      levels = c("Advanced", "Early")),
    Age_group = factor(ifelse(age_at_index >= 60, ">=60", "<60"), levels = c("<60", ">=60"))
  )

cox_multi     <- coxph(Surv(OS_time, OS_status) ~
                         Disulfidptosis_Canonical + Stage_simple + Age_group + gender,
                       data = surv_df2)
cox_multi_sum <- summary(cox_multi)

cox_uni     <- coxph(Surv(OS_time, OS_status) ~ Disulfidptosis_Canonical, data = surv_df)
cox_uni_sum <- summary(cox_uni)
hr_uni  <- round(cox_uni_sum$conf.int[1, "exp(coef)"], 2)
ci_low  <- round(cox_uni_sum$conf.int[1, "lower .95"],  2)
ci_high <- round(cox_uni_sum$conf.int[1, "upper .95"],  2)
p_cox   <- cox_uni_sum$coefficients[1, "Pr(>|z|)"]
lr_test <- survdiff(Surv(OS_time, OS_status) ~ Score_group, data = surv_df)
p_lr    <- 1 - pchisq(lr_test$chisq, df = 1)
cat(sprintf("HR = %.2f (IC95%%: %.2f-%.2f), p Cox = %.2e, p Log-rank = %.2e\n",
            hr_uni, ci_low, ci_high, p_cox, p_lr))

cox_table <- as.data.frame(cox_multi_sum$conf.int) %>%
  rownames_to_column("Variable") %>%
  rename(HR = `exp(coef)`, HR_IC_inferior = `lower .95`, HR_IC_superior = `upper .95`) %>%
  mutate(p_valor = cox_multi_sum$coefficients[, "Pr(>|z|)"],
         across(where(is.numeric), ~ round(., 4))) %>%
  mutate(Variable = recode(Variable,
    "Disulfidptosis_Canonical" = "Score Disulfidptosis",
    "Stage_simpleEarly"        = "Estadio: Temprano vs Avanzado",
    "Age_group>=60"            = "Edad: ≥60 vs <60",
    "gendermale"               = "Sexo: Masculino vs Femenino"))

estilo_header <- createStyle(fontColour = "#FFFFFF", fgFill = "#1167B1", halign = "LEFT",
                              fontName = "Calibri", fontSize = 11, textDecoration = "bold",
                              border = "Bottom", borderColour = "#03254C")
estilo_alt    <- createStyle(fgFill = "#D0E8F5")

agregar_hoja <- function(wb, nombre, datos) {
  addWorksheet(wb, nombre)
  writeData(wb, nombre, datos, headerStyle = estilo_header)
  if (nrow(datos) > 1)
    for (i in seq(2, nrow(datos), by = 2))
      addStyle(wb, nombre, style = estilo_alt, rows = i + 1,
               cols = seq_len(ncol(datos)), gridExpand = TRUE)
  setColWidths(wb, nombre, cols = seq_len(ncol(datos)), widths = "auto")
}

wb <- createWorkbook()

agregar_hoja(wb, "S1_Cohorte_Clinica",
  surv_df %>%
    dplyr::select(sample_id, patient_id, age_at_index, gender,
                  ajcc_pathologic_stage, OS_time, OS_status,
                  Disulfidptosis_Canonical, Disulfidptosis_Extended,
                  PCA_score, Score_group, Cluster) %>%
    rename(Muestra = sample_id, Paciente = patient_id, Edad = age_at_index,
           Sexo = gender, Estadio_AJCC = ajcc_pathologic_stage,
           Tiempo_OS_dias = OS_time, Estado_OS = OS_status,
           Score_Canonico = Disulfidptosis_Canonical,
           Score_Extendido = Disulfidptosis_Extended,
           Score_PCA = PCA_score, Grupo_Score = Score_group) %>%
    arrange(Cluster, Grupo_Score))

agregar_hoja(wb, "S2_Scores_Disulfidptosis",
  surv_df %>%
    dplyr::select(sample_id, Disulfidptosis_Canonical, Disulfidptosis_Extended,
                  PCA_score, Score_group, Cluster) %>%
    rename(Muestra = sample_id, Score_Canonico = Disulfidptosis_Canonical,
           Score_Extendido = Disulfidptosis_Extended, Score_PCA = PCA_score,
           Grupo_Score = Score_group))

agregar_hoja(wb, "S3_DEG_Completo",
  deg_res %>%
    dplyr::select(gene, logFC, AveExpr, t, P.Value, adj.P.Val, B) %>%
    rename(Gen = gene, log2FC = logFC, Expresion_media = AveExpr,
           Estadistico_t = t, p_valor = P.Value,
           p_valor_ajustado_BH = adj.P.Val, Estadistico_B = B) %>%
    arrange(p_valor_ajustado_BH))

agregar_hoja(wb, "S4_DEG_Significativos",
  deg_res %>% filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
    dplyr::select(gene, logFC, AveExpr, t, P.Value, adj.P.Val, B) %>%
    rename(Gen = gene, log2FC = logFC, Expresion_media = AveExpr,
           Estadistico_t = t, p_valor = P.Value,
           p_valor_ajustado_BH = adj.P.Val, Estadistico_B = B) %>%
    arrange(p_valor_ajustado_BH))

agregar_hoja(wb, "S5_GSEA_GO_BiologicalProcess",
  go_df %>%
    dplyr::select(ID, Description, setSize, enrichmentScore, NES,
                  pvalue, p.adjust, qvalue, leading_edge) %>%
    rename(Termino_GO = ID, Descripcion = Description, Tamano_set = setSize,
           Score_enriquecimiento = enrichmentScore, NES_normalizado = NES,
           p_valor = pvalue, p_ajustado_BH = p.adjust, q_valor = qvalue,
           Genes_principales = leading_edge) %>%
    arrange(p_ajustado_BH))

agregar_hoja(wb, "S6_GSEA_KEGG_Rutas",
  kegg_df %>%
    dplyr::select(ID, Description, setSize, enrichmentScore, NES,
                  pvalue, p.adjust, qvalue, leading_edge) %>%
    rename(Ruta_KEGG = ID, Descripcion = Description, Tamano_set = setSize,
           Score_enriquecimiento = enrichmentScore, NES_normalizado = NES,
           p_valor = pvalue, p_ajustado_BH = p.adjust, q_valor = qvalue,
           Genes_principales = leading_edge) %>%
    arrange(p_ajustado_BH))

if (file.exists(file.path(EOUT, "ORA_GO_BP.rds"))) {
  ora_go <- readRDS(file.path(EOUT, "ORA_GO_BP.rds"))
  agregar_hoja(wb, "S7_ORA_GO_BiologicalProcess",
    as.data.frame(ora_go) %>%
      dplyr::select(ID, Description, GeneRatio, BgRatio,
                    pvalue, p.adjust, qvalue, geneID, Count) %>%
      rename(Termino_GO = ID, Descripcion = Description,
             Razon_genes = GeneRatio, Razon_fondo = BgRatio,
             p_valor = pvalue, p_ajustado_BH = p.adjust,
             q_valor = qvalue, Genes = geneID, N_genes = Count) %>%
      arrange(p_ajustado_BH))
}

agregar_hoja(wb, "S8_EPIC_TME_Fracciones",
  cib_res %>%
    left_join(surv_df[, c("sample_id", "Cluster", "Score_group")],
              by = setNames("sample_id", names(cib_res)[1])) %>%
    dplyr::select(1, Cluster, Score_group,
                  all_of(frac_cols[frac_cols %in% colnames(cib_res)])))

agregar_hoja(wb, "S9_Cox_Multivariado", cox_table)

xl_path <- file.path(TOUT, "Tablas_Suplementarias_Disulfidptosis_LIHC.xlsx")
saveWorkbook(wb, xl_path, overwrite = TRUE)
message("Excel guardado, con la hoja S8_EPIC_TME_Fracciones")
message("Listo con las tablas.")

message("Ahora las figuras de TME (07, 08 y 09)")

for (p in c("ggplot2", "ggrepel", "patchwork", "survminer",
            "stringr", "scales", "cluster", "tibble"))
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")

library(ggplot2); library(dplyr); library(tidyr); library(stringr)
library(tibble)

surv_df  <- readRDS(file.path(KOUT, "Datos_supervivencia.rds"))
cib_obj  <- readRDS(file.path(MOUT, "EPIC_objeto.rds"))
cib_res  <- cib_obj$cib
frac_cols <- cib_obj$frac_cols

# Excluyo ID, Tumor, Uncharacterized y otherCells (esta última ronda el 68%
# y aplasta la escala de las demás células en las figuras)
excluir <- c("sample", "Tumor", "Uncharacterized", "RMSE", "Correlation",
             "tumor.purity", "otherCells", "NKcells")
frac_ok <- colnames(cib_res)[!colnames(cib_res) %in% excluir & sapply(cib_res, is.numeric)]

top_cells <- cib_res %>%
  dplyr::select(all_of(frac_ok)) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "cell", values_to = "mean") %>%
  slice_max(mean, n = 7) %>%
  pull(cell)

id_col <- colnames(cib_res)[1]
message("Tipos celulares en las figuras: ", paste(top_cells, collapse = ", "))

message("Generando la Fig 07: violín de TME por score")
cib_score <- cib_res %>%
  dplyr::select(all_of(c(id_col, top_cells))) %>%
  left_join(surv_df[, c("sample_id", "Score_group")],
            by = setNames("sample_id", id_col)) %>%
  filter(!is.na(Score_group)) %>%
  pivot_longer(-c(all_of(id_col), Score_group),
               names_to = "CellType", values_to = "Fraction") %>%
  mutate(CellType = str_to_title(str_replace_all(CellType, "_", " ")))

cell_labels <- c(
  "Bcells"      = "B cells",
  "Cafs"        = "CAFs",
  "Cd4 Tcells"  = "CD4\u207a T cells",
  "Cd8 Tcells"  = "CD8\u207a T cells",
  "Endothelial" = "Endothelial",
  "Macrophages" = "Macrophages",
  "Nkcells"     = "NK cells"
)

mutate(CellType = dplyr::recode(CellType, !!!cell_labels))


p_tme_score <- ggplot(cib_score,
                      aes(x = Score_group, y = Fraction, fill = Score_group)) +
  geom_violin(alpha = 0.70, trim = TRUE, linewidth = 0.25, color = "white") +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white",
               linewidth = 0.45, color = "grey25") +
  scale_fill_manual(values = c("High" = COL_HIGH, "Low" = COL_LOW),
                    name = "Score Group") +
  facet_wrap(~ CellType, scales = "free_y", nrow = 2) +
  labs(
    title    = "Tumor Microenvironment by Disulfidptosis Score",
    subtitle = "High vs Low — EPIC Deconvolution (IOBR) | TCGA-LIHC",
    x        = "Score Group",
    y        = "Cell Fraction (normalized)",
    caption  = "EPIC cell-fraction estimates; fractions normalized per sample (sum = 1)"
  ) +
  theme_blue(base_size = 10) +
  theme(legend.position = "top", axis.text.x = element_text(size = 9))

guardar_figura(p_tme_score, "Fig07_TME_Violin_Score", ancho = 16, alto = 8)

message("Generando la Fig 08: violín de TME por cluster")
cib_clust <- cib_res %>%
  dplyr::select(all_of(c(id_col, top_cells))) %>%
  left_join(surv_df[, c("sample_id", "Cluster")],
            by = setNames("sample_id", id_col)) %>%
  filter(!is.na(Cluster)) %>%
  pivot_longer(-c(all_of(id_col), Cluster),
               names_to = "CellType", values_to = "Fraction") %>%
  mutate(CellType = str_to_title(str_replace_all(CellType, "_", " ")))

p_tme_clust <- ggplot(cib_clust,
                      aes(x = Cluster, y = Fraction, fill = Cluster)) +
  geom_violin(alpha = 0.70, trim = TRUE, linewidth = 0.25, color = "white") +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white",
               linewidth = 0.45, color = "grey25") +
  scale_fill_manual(
    values = setNames(c(BLUE_DARKEST, BLUE_LIGHT, BLUE_MID, BLUE_PALE), paste0("C", 1:4)),
    name = "Cluster") +
  facet_wrap(~ CellType, scales = "free_y", nrow = 2) +
  labs(
    title    = "Tumor Microenvironment by Disulfidptosis Cluster",
    subtitle = "EPIC Deconvolution (IOBR) | TCGA-LIHC",
    x        = "Cluster",
    y        = "Cell Fraction (normalized)"
  ) +
  theme_blue(base_size = 10) +
  theme(legend.position = "top", axis.text.x = element_text(size = 9))

guardar_figura(p_tme_clust, "Fig08_TME_Violin_Cluster", ancho = 16, alto = 8)

message("Generando la Fig 09: correlación de Spearman entre TME y score")
cib_wide <- cib_res %>%
  dplyr::select(all_of(c(id_col, frac_ok))) %>%
  left_join(surv_df[, c("sample_id", "Disulfidptosis_Canonical")],
            by = setNames("sample_id", id_col)) %>%
  filter(!is.na(Disulfidptosis_Canonical))

cor_df <- lapply(frac_ok, function(cell) {
  x  <- cib_wide[[cell]]
  y  <- cib_wide$Disulfidptosis_Canonical
  ok <- complete.cases(x, y)
  if (sum(ok) < 10) return(NULL)
  ct <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  data.frame(CellType = str_to_title(str_replace_all(cell, "_", " ")),
             rho = as.numeric(ct$estimate), p_value = ct$p.value,
             stringsAsFactors = FALSE)
}) %>% bind_rows() %>%
  mutate(sig       = ifelse(p_value < 0.05, "p < 0.05", "n.s."),
         rho_label = ifelse(p_value < 0.05, sprintf("%.2f", rho), "")) %>%
  arrange(rho)

p_corr <- ggplot(cor_df, aes(x = rho, y = reorder(CellType, rho),
                              fill = rho, alpha = sig)) +
  geom_col(width = 0.70, color = NA) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.55) +
  geom_text(aes(label = rho_label,
                hjust = ifelse(rho >= 0, -0.2, 1.2)), size = 2.9, color = "grey20") +
  scale_fill_gradient2(low = BLUE_LIGHT, mid = "white", high = BLUE_DARKEST,
                       midpoint = 0, name = "Spearman ρ") +
  scale_alpha_manual(values = c("p < 0.05" = 1, "n.s." = 0.35), name = "Significance") +
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

guardar_figura(p_corr, "Fig09_TME_Correlacion_Spearman", ancho = 13, alto = 10)

message("Figuras 07, 08 y 09 regeneradas con las etiquetas de EPIC.")
message("Excel actualizado con la hoja S8_EPIC_TME_Fracciones.")
message("Archivos generados en:")
message("Figuras: ", FOUT)
message("Excel: ", file.path(TOUT, "Tablas_Suplementarias_Disulfidptosis_LIHC.xlsx"))
message("Datos: ", file.path(MOUT, "EPIC_objeto.rds"))
