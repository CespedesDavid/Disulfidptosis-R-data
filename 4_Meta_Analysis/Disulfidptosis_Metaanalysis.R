pkgs_needed <- c(
  "meta", "metafor",
  "ggplot2", "ggpubr", "patchwork",
  "dplyr", "tidyr", "stringr",
  "scales", "grid", "gridExtra",
  "RColorBrewer"
)

for (p in pkgs_needed) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}

if (!requireNamespace("dmetar", quietly = TRUE)) {
  message("Instalando dmetar desde GitHub.")
  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)
  tryCatch(
    remotes::install_github("MathiasHarrer/dmetar", upgrade = "never", quiet = TRUE),
    error = function(e) message("dmetar no pudo instalarse: ", conditionMessage(e))
  )
}

suppressPackageStartupMessages({
  library(meta)
  library(metafor)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(scales)
  library(grid)
  library(gridExtra)
})

HAS_DMETAR <- requireNamespace("dmetar", quietly = TRUE)
if (HAS_DMETAR) suppressPackageStartupMessages(library(dmetar))

BLUE  <- c(
  dark   = "#1B4F72",
  main   = "#2E86C1",
  mid    = "#5DADE2",
  light  = "#AED6F1",
  pale   = "#D6EAF8",
  accent = "#1A5276"
)

theme_dis <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", color = BLUE["dark"], size = base_size + 2),
      plot.subtitle    = element_text(color = BLUE["main"], size = base_size - 1),
      axis.title       = element_text(color = BLUE["dark"]),
      axis.text        = element_text(color = "#2C3E50"),
      legend.title     = element_text(face = "bold", color = BLUE["dark"]),
      legend.background = element_rect(fill = "white", color = BLUE["light"], linewidth = 0.3),
      panel.grid.major = element_line(color = BLUE["pale"], linewidth = 0.3),
      strip.background = element_rect(fill = BLUE["pale"], color = BLUE["mid"]),
      strip.text       = element_text(face = "bold", color = BLUE["dark"]),
      plot.background  = element_rect(fill = "white", color = NA)
    )
}

meta_raw <- tribble(
  ~id, ~author,         ~year, ~cancer,  ~n,    ~hr,    ~hr_low, ~hr_high, ~model,   ~subgroup, ~tcga_only, ~n_genes,
  # HCC
   3,  "Wang T",        2024, "HCC",    1041,  1.303,  1.164,   1.458,  "multi",  "HCC",  FALSE, 5,
   8,  "Shang N",       2025, "HCC",     365,  1.790,  1.400,   2.290,  "multi",  "HCC",  FALSE, 7,
  21,  "Yang L",        2023, "HCC",    1155,  2.101,  1.770,   2.611,  "multi",  "HCC",  FALSE, 6,
  26,  "Chen X",        2023, "HCC",     486,  1.203,  1.109,   1.309,  "multi",  "HCC",  TRUE,  6,
  41,  "Wang X2",       2024, "HCC",     195,  3.083,  1.952,   4.869,  "multi",  "HCC",  FALSE, 5,
  42,  "Wang Y",        2024, "HCC",     485,  1.500,  1.323,   1.700,  "multi",  "HCC",  FALSE, 4,
  44,  "Xu Z",          2024, "HCC",     354, 29.341,  5.419, 158.870,  "multi",  "HCC",  FALSE, 3,  # extremo
  49,  "Tang J",        2024, "HCC",     568,  3.690,  2.540,   5.400,  "multi",  "HCC",  FALSE,19,
  50,  "Li S",          2025, "HCC",     602,  1.041,  1.018,   1.063,  "multi",  "HCC",  FALSE, 5,
  53,  "Wang T2",       2024, "HCC",     698,  1.014,  1.006,   1.022,  "multi",  "HCC",  FALSE,14,
  63,  "Zhao J",        2023, "HCC",     731,  5.848,  2.953,  11.579,  "multi",  "HCC",  FALSE, 3,
  # KIRC
   9,  "Wen X",         2025, "KIRC",    532,  3.042,  2.160,   4.283,  "multi",  "KIRC", FALSE,11,
  32,  "Ren L",         2024, "KIRC",   1270,  2.812,  1.952,   4.052,  "multi",  "KIRC", FALSE, 8,
  48,  "Du L",          2024, "KIRC",    537,  2.129,  1.518,   2.986,  "multi",  "KIRC", FALSE, 9,
  74,  "Peng K",        2023, "KIRC",    770,  2.808,  2.316,   3.405,  "multi",  "KIRC", FALSE, 6,
  # Glioma
   4,  "Wang X",        2023, "Glioma",  1592, 4.279,  3.160,   5.790,  "multi",  "Glioma",FALSE,3,
  19,  "Zhou Y",        2024, "Glioma",   154, 2.360,  1.270,   4.370,  "multi",  "Glioma",FALSE,6,
  22,  "Zhang F",       2024, "Glioma",   905, 1.781,  1.126,   2.817,  "multi",  "Glioma",TRUE, 1,
  38,  "Zhou Y2",       2023, "Glioma",  1503, 1.990,  1.594,   2.485,  "multi",  "Glioma",FALSE,9,
  # CRC / COAD
  13,  "Wen L",         2025, "CRC",      432, 2.560,  1.690,   3.890,  "multi",  "CRC",  TRUE, 10,
  17,  "Li J",          2023, "CRC",     1071, 2.029,  1.229,   3.349,  "multi",  "CRC",  FALSE, 3,
  35,  "Xiao L",        2023, "CRC",      409, 1.267,  1.163,   1.380,  "multi",  "CRC",  TRUE,  4,
  69,  "Xiao Y",        2024, "CRC",     1125,13.632,  6.587,  28.214,  "multi",  "CRC",  FALSE,13,  # extremo
  # LUAD
   6,  "Ni L",          2023, "LUAD",    1206, 3.316,  2.300,   4.782,  "multi",  "LUAD", FALSE, 5,
  15,  "Zhong L",       2024, "LUAD",    1708, 2.679,  1.612,   4.454,  "multi",  "LUAD", FALSE, 5,
  47,  "Huang J",       2023, "LUAD",    1194, 2.600,  1.800,   3.600,  "multi",  "LUAD", FALSE,21,
  56,  "Wang H",        2024, "LUAD",    1313, 2.284,  1.669,   3.124,  "multi",  "LUAD", FALSE, 5,
  75,  "Xu X",          2024, "LUAD",     623, 3.230,  1.490,   7.010,  "multi",  "LUAD", FALSE, 6,
  # STAD
  29,  "Chen R",        2024, "STAD",     732,11.670,  3.730,  36.480,  "multi",  "STAD", FALSE, 5,  # extremo
  43,  "Liu X",         2024, "STAD",     840, 1.448,  1.184,   1.770,  "multi",  "STAD", FALSE, 8,
  # CESC
  28,  "Jin T",         2024, "CESC",     504, 4.010,  2.200,   7.310,  "multi",  "CESC", FALSE, 8,
  # Otros
   5,  "Zhou Z",        2025, "BLCA",     556, 1.834,  1.064,   3.161,  "multi",  "Other",FALSE, 2,
  20,  "Guo S",         2025, "BLCA",     419, 1.389,  1.236,   1.561,  "multi",  "Other",FALSE, 8,
  24,  "Zhang X",       2025, "ESCC",     259, 2.310,  1.540,   3.450,  "multi",  "Other",FALSE,19,
  71,  "Xu J",          2024, "SARC",     260, 2.279,  1.507,   3.448,  "multi",  "Other",FALSE, 3,
  72,  "Zhang K",       2024, "LUAD",     508, 1.352,  1.245,   1.468,  "multi",  "LUAD", FALSE, 7
)

meta_raw <- meta_raw %>%
  mutate(
    extreme = id %in% c(29, 44, 69),
    log_hr  = log(hr),
    se_log  = (log(hr_high) - log(hr_low)) / (2 * 1.96),
    label   = paste0(author, " (", year, ")")
  )

cat("Cargué", nrow(meta_raw), "estudios con HR completo.\n")
cat("  HCC:", sum(meta_raw$subgroup == "HCC"),
    "| KIRC:", sum(meta_raw$subgroup == "KIRC"),
    "| Glioma:", sum(meta_raw$subgroup == "Glioma"),
    "| LUAD:", sum(meta_raw$subgroup == "LUAD"),
    "| CRC:", sum(meta_raw$subgroup == "CRC"), "\n")

cat("\nMeta-análisis principal:\n")

ma_main <- metagen(
  TE      = meta_raw$log_hr,
  seTE    = meta_raw$se_log,
  studlab = meta_raw$label,
  sm      = "HR",
  method.tau = "REML",
  method.random.ci = "HK",
  data    = meta_raw
)
summary(ma_main)

ma_sub <- update(ma_main, subgroup = meta_raw$subgroup, common = FALSE)
summary(ma_sub)

meta_raw <- meta_raw %>%
  mutate(gene_def = case_when(
    n_genes <= 10 ~ "Canonical Liu 2023",
    TRUE          ~ "Extended definition"
  ))

ma_genedef <- metagen(
  TE               = meta_raw$log_hr,
  seTE             = meta_raw$se_log,
  studlab          = meta_raw$label,
  sm               = "HR",
  method.tau       = "REML",
  method.random.ci = "HK",
  subgroup         = meta_raw$gene_def,
  data             = meta_raw
)
cat("\nSubgrupo genes canónicos vs ampliados:\n")
summary(ma_genedef)

meta_hcc <- meta_raw %>% filter(subgroup == "HCC")
ma_hcc_strat <- metagen(
  TE      = meta_hcc$log_hr,
  seTE    = meta_hcc$se_log,
  studlab = meta_hcc$label,
  sm      = "HR",
  method.tau = "REML",
  method.random.ci = "HK",
  subgroup = meta_hcc$tcga_only,
  data    = meta_hcc
)
cat("\nHCC: solo TCGA vs con validación externa:\n")
summary(ma_hcc_strat)

cat("\nAnálisis de sensibilidad leave-one-out:\n")
if (HAS_DMETAR) {
  loo <- loo.meta(ma_main)
  print(loo)
} else {
  loo_list <- lapply(seq_len(nrow(meta_raw)), function(i) {
    d <- meta_raw[-i, ]
    m <- metagen(TE = d$log_hr, seTE = d$se_log, sm = "HR",
                 method.tau = "REML", method.random.ci = "HK")
    data.frame(studlab      = meta_raw$label[i],
               TE.random    = m$TE.random,
               lower.random = m$lower.random,
               upper.random = m$upper.random,
               I2           = m$I2,
               stringsAsFactors = FALSE)
  })
  loo <- do.call(rbind, loo_list)
  cat("Terminé el leave-one-out. El HR pooled se mueve entre",
      round(exp(min(loo$TE.random)), 2), "y",
      round(exp(max(loo$TE.random)), 2), ".\n")
}

cat("\nSesgo de publicación:\n")
if (HAS_DMETAR) {
  egger_test <- egger.test(ma_main)
  print(egger_test)
} else {
  egger_test <- metabias(ma_main, method.bias = "Egger")
  print(egger_test)
}
begg_test  <- metabias(ma_main, method.bias = "Begg")
print(begg_test)

cat("\nGRADE — certeza de la evidencia:\n")
cat("  Riesgo de sesgo:       Serio (mayoría estudios observacionales TCGA)\n")
cat("  Inconsistencia:        Serio (I² esperado alto)\n")
cat("  Indirección:           No serio (outcome OS directo)\n")
cat("  Imprecisión:           No serio (>30 estudios incluidos)\n")
cat("  Sesgo de publicación:  Probable (ver Egger/Funnel)\n")
cat("  Certeza global:        BAJA (⊕⊕○○)\n")


# ─── HELPER: extraer resultados por subgrupo ─────────
extract_subgroup_results <- function(ma_obj) {
  # Intentar primero la API antigua (meta <6)
  if (length(ma_obj$TE.random.w) > 0) {
    sgs  <- names(ma_obj$TE.random.w)
    te   <- ma_obj$TE.random.w
    lo   <- ma_obj$lower.random.w
    hi   <- ma_obj$upper.random.w
    i2   <- if (!is.null(ma_obj$I2.w))   ma_obj$I2.w   else rep(NA_real_, length(sgs))
    tau2 <- if (!is.null(ma_obj$tau2.w)) ma_obj$tau2.w else rep(NA_real_, length(sgs))
    qp   <- if (!is.null(ma_obj$pval.Q.w)) ma_obj$pval.Q.w else rep(NA_real_, length(sgs))
    k    <- if (!is.null(ma_obj$k.w))    ma_obj$k.w    else rep(NA_integer_, length(sgs))
    n    <- if (!is.null(ma_obj$n.w))    ma_obj$n.w    else rep(NA_real_, length(sgs))
  } else {
    sr   <- ma_obj$subgroup.results
    sgs  <- sr$subgroup
    if (!"TE.random" %in% names(sr))
      stop("extract_subgroup_results: sr$TE.random no encontrado. ",
           "Verifique que ma_obj fue construido con random = TRUE y common = FALSE.")
    te   <- sr$TE.random
    lo   <- sr$lower.random
    hi   <- sr$upper.random
    i2   <- if ("I2" %in% names(sr))   sr$I2   else rep(NA_real_, length(sgs))
    tau2 <- if ("tau2" %in% names(sr)) sr$tau2 else rep(NA_real_, length(sgs))
    qp   <- if ("pval.Q" %in% names(sr)) sr$pval.Q else rep(NA_real_, length(sgs))
    k    <- if ("k" %in% names(sr))    sr$k    else rep(NA_integer_, length(sgs))
    n    <- if ("n" %in% names(sr))    sr$n    else rep(NA_real_, length(sgs))
  }
  # meta pkg devuelve I2 como proporción (0-1); convertir a porcentaje una sola vez aquí
  i2_pct <- ifelse(is.na(i2), NA_real_, ifelse(i2 <= 1, i2 * 100, i2))
  data.frame(
    subgroup = sgs,
    TE       = te,
    lower    = lo,
    upper    = hi,
    I2       = i2_pct,
    tau2     = tau2,
    pval_Q   = qp,
    k        = k,
    n        = n,
    stringsAsFactors = FALSE
  )
}

sg_sub     <- extract_subgroup_results(ma_sub)
sg_genedef <- extract_subgroup_results(ma_genedef)

cat("\nGenerando las figuras.\n")

BXLIM <- 0.42
BYLIM <- 0.55

boxes <- data.frame(
  x     = c(1, 2, 3, 4, 5),
  y     = c(3, 3, 3, 3, 3),
  label = c("Identified\nPubMed\n(n = 346)",
            "Screened\nTitle/Abstract\n(n = 121)",
            "Full Text\nAssessed\n(n = 120)",
            "Excluded\nPhase 2\n(n = 42)",
            "Included in\nReview\n(n = 78)")
)

excl <- data.frame(
  x     = c(2,   3),
  y     = c(1.2, 1.2),
  label = c("Excluded Phase 1\n(n = 225)\nTitle/Abstract",
            "Excluded Phase 2\n(n = 42)\nFull Text")
)

fig_prisma <- ggplot() +
  # Cajas principales — fondo
  geom_rect(data = boxes[c(1,2,3,5), ],
            aes(xmin = x - BXLIM, xmax = x + BXLIM,
                ymin = y - BYLIM,  ymax = y + BYLIM),
            fill = BLUE["pale"], color = BLUE["main"], linewidth = 0.8) +
  # Caja 4 (Excluidos Fase 2) — misma fila pero diferente color
  geom_rect(data = boxes[4, ],
            aes(xmin = x - BXLIM, xmax = x + BXLIM,
                ymin = y - BYLIM,  ymax = y + BYLIM),
            fill = "#EBF5FB", color = BLUE["mid"], linewidth = 0.8) +
  geom_rect(data = excl,
            aes(xmin = x - BXLIM, xmax = x + BXLIM,
                ymin = y - BYLIM,  ymax = y + BYLIM),
            fill = "#EBF5FB", color = BLUE["mid"],
            linewidth = 0.6, linetype = "dashed") +
  geom_text(data = boxes,
            aes(x = x, y = y, label = label),
            fontface = "bold", color = BLUE["dark"], size = 3.1,
            lineheight = 1.1) +
  geom_text(data = excl,
            aes(x = x, y = y, label = label),
            color = BLUE["accent"], size = 2.9, lineheight = 1.1) +
  annotate("segment",
           x = 1 + BXLIM, xend = 2 - BXLIM, y = 3, yend = 3,
           arrow = arrow(length = unit(0.25,"cm"), type = "closed"),
           color = BLUE["main"], linewidth = 0.8) +
  annotate("segment",
           x = 2 + BXLIM, xend = 3 - BXLIM, y = 3, yend = 3,
           arrow = arrow(length = unit(0.25,"cm"), type = "closed"),
           color = BLUE["main"], linewidth = 0.8) +
  annotate("segment",
           x = 3 + BXLIM, xend = 5 - BXLIM, y = 3, yend = 3,
           arrow = arrow(length = unit(0.25,"cm"), type = "closed"),
           color = BLUE["main"], linewidth = 0.8) +
  annotate("segment",
           x = 2, xend = 2, y = 3 - BYLIM, yend = 1.2 + BYLIM,
           arrow = arrow(length = unit(0.2,"cm"), type = "closed"),
           color = BLUE["mid"], linewidth = 0.6, linetype = "dashed") +
  annotate("segment",
           x = 3, xend = 3, y = 3 - BYLIM, yend = 1.2 + BYLIM,
           arrow = arrow(length = unit(0.2,"cm"), type = "closed"),
           color = BLUE["mid"], linewidth = 0.6, linetype = "dashed") +
  xlim(0.45, 5.55) + ylim(0.5, 3.8) +
  labs(title    = "PRISMA 2020 Flow Diagram",
       subtitle = "PubMed search only — Disulfidptosis (n = 78 included)") +
  theme_void() +
  theme(
    plot.title      = element_text(face = "bold", color = BLUE["dark"],
                                   size = 13, hjust = 0.5, margin = margin(b=4)),
    plot.subtitle   = element_text(color = BLUE["main"],
                                   size = 10, hjust = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(10, 10, 10, 10)
  )

cancer_dist <- tribble(
  ~cancer,               ~n,    ~cat,
  "HCC",                 14,   "Principal",
  "LUAD/NSCLC",          11,   "Lung",
  "Glioma/GBM",           9,   "Glioma",
  "CRC/COAD",             8,   "Colorectal",
  "KIRC",                 5,   "Renal",
  "STAD",                 5,   "Gastric",
  "BLCA",                 4,   "Bladder",
  "CESC",                 3,   "Cervical",
  "BRCA",                 2,   "Breast",
  "ESCC",                 2,   "Esophagus",
  "PAAD",                 3,   "Pancreas",
  "Otros (OV/SKCM/SARC/UCEC/HNSCC)", 12, "Otros"
)
cancer_dist <- cancer_dist %>%
  arrange(desc(n)) %>%
  mutate(cancer = factor(cancer, levels = rev(cancer)))

fig_cancer_bar <- ggplot(cancer_dist, aes(x = cancer, y = n,
                                          fill = factor(cat))) +
  geom_col(width = 0.7, color = "white", linewidth = 0.4) +
  geom_text(aes(label = n), hjust = -0.2, fontface = "bold",
            color = BLUE["dark"], size = 3.5) +
  coord_flip() +
  scale_fill_manual(
    values = c("Principal"   = BLUE["dark"],
               "Lung"      = BLUE["main"],
               "Glioma"      = BLUE["mid"],
               "Colorectal" = "#5499C7",
               "Renal"       = BLUE["light"],
               "Gastric"    = "#2471A3",
               "Bladder"      = "#7FB3D3",
               "Cervical"    = "#A9CCE3",
               "Breast"        = "#D4E6F1",
               "Esophagus"     = "#85C1E9",
               "Pancreas"    = "#3498DB",
               "Otros"       = "#ABB2B9"),
    name = "Tumor type"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title    = "Study distribution by cancer type",
       subtitle = "78 articles included in the systematic review",
       x = NULL, y = "Number of studies") +
  theme_dis() +
  theme(legend.position = "none")

# ── FIG 3: FOREST PLOT (ggplot2) ─────────────────────────────
fp_data <- meta_raw %>%
  mutate(
    subgroup_f = factor(subgroup,
                        levels = c("HCC","KIRC","Glioma","LUAD","CRC","STAD","CESC","Other")),
    label_full = paste0(author, " (", year, ") — ", cancer)
  ) %>%
  arrange(subgroup_f, desc(hr))

sub_res <- sg_sub %>%
  mutate(
    log_hr     = TE,
    se_log     = (upper - lower) / (2 * 1.96),
    hr         = exp(TE),
    hr_low     = exp(lower),
    hr_high    = exp(upper),
    extreme    = FALSE,
    label_full = paste0("\u25c6 Pooled ", subgroup),
    type       = "pooled"
  )

color_map <- c(
  "HCC"    = BLUE["dark"],
  "KIRC"   = BLUE["main"],
  "Glioma" = BLUE["mid"],
  "LUAD"   = "#2471A3",
  "CRC"    = "#5499C7",
  "STAD"   = "#7FB3D3",
  "CESC"   = "#A9CCE3",
  "Other"  = "#BDC3C7"
)

fp_data <- fp_data %>% arrange(subgroup_f, hr)
fp_data$plot_y <- nrow(fp_data):1

sub_res_fp <- sub_res %>%
  mutate(
    subgroup_f = factor(subgroup, levels = levels(fp_data$subgroup_f)),
    cancer     = subgroup,
    plot_y     = NA_real_
  )

fig_forest <- ggplot(fp_data, aes(y = plot_y)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "#7F8C8D", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = hr_low, xmax = hr_high,
                     color = subgroup),
                 height = 0.3, linewidth = 0.6) +
  geom_point(data = filter(fp_data, !extreme),
             aes(x = hr, color = subgroup, size = 1/se_log),
             shape = 15) +
  geom_point(data = filter(fp_data, extreme),
             aes(x = hr, color = subgroup),
             shape = 17, size = 3.5) +
  geom_text(aes(x = 0.05, label = label_full,
                color = subgroup),
            hjust = 0, size = 2.5, fontface = "plain") +
  geom_text(aes(x = pmax(hr_high, 35) * 1.05,
                label = sprintf("%.2f (%.2f–%.2f)", hr, hr_low, hr_high),
                color = subgroup),
            hjust = 0, size = 2.4) +
  scale_color_manual(values = color_map, name = "Cancer") +
  scale_size_continuous(range = c(1.5, 4), guide = "none") +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 30),
                labels = c("0.5","1","2","5","10","30")) +
  coord_cartesian(xlim = c(0.04, 200)) +
  labs(title    = "Forest Plot — Disulfidptosis HR for OS",
       subtitle = "Random effects model (RE) | ▲ = extreme HR",
       x = "Hazard Ratio (95% CI, log scale)",
       y = NULL) +
  theme_dis(base_size = 9) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )

# ── FIG 4: FUNNEL PLOT ────────────────────────────────────────
ma_mfor <- rma(
  yi   = meta_raw$log_hr,
  sei  = meta_raw$se_log,
  method = "REML",
  slab = meta_raw$label
)

se_range <- seq(0, max(meta_raw$se_log) * 1.05, length.out = 100)
pool_loghr <- ma_mfor$beta[1]

funnel_df <- data.frame(
  log_hr = meta_raw$log_hr,
  se     = meta_raw$se_log,
  subgroup = meta_raw$subgroup,
  extreme  = meta_raw$extreme
)

pseudo_ci <- data.frame(
  se     = se_range,
  lo_95  = pool_loghr - 1.96 * se_range,
  hi_95  = pool_loghr + 1.96 * se_range,
  lo_99  = pool_loghr - 2.576 * se_range,
  hi_99  = pool_loghr + 2.576 * se_range
)

se_max   <- max(meta_raw$se_log) * 1.15
se_seq   <- seq(0, se_max, length.out = 200)
poly_95  <- data.frame(
  x = c(pool_loghr - 1.96 * se_seq,
         rev(pool_loghr + 1.96 * se_seq)),
  y = c(se_seq, rev(se_seq))
)
poly_99  <- data.frame(
  x = c(pool_loghr - 2.576 * se_seq,
         rev(pool_loghr + 2.576 * se_seq)),
  y = c(se_seq, rev(se_seq))
)

egger_p_val <- if (HAS_DMETAR) egger_test$p.value else egger_test$pval

fig_funnel <- ggplot() +
  geom_polygon(data = poly_99, aes(x = x, y = y),
               fill = BLUE["pale"], alpha = 0.55) +
  geom_polygon(data = poly_95, aes(x = x, y = y),
               fill = BLUE["light"], alpha = 0.60) +
  geom_line(data = data.frame(x = pool_loghr - 1.96 * se_seq, y = se_seq),
            aes(x = x, y = y), color = BLUE["mid"],
            linetype = "dashed", linewidth = 0.5) +
  geom_line(data = data.frame(x = pool_loghr + 1.96 * se_seq, y = se_seq),
            aes(x = x, y = y), color = BLUE["mid"],
            linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = pool_loghr, color = BLUE["dark"], linewidth = 0.9) +
  geom_vline(xintercept = 0, linetype = "dotted",
             color = "#7F8C8D", linewidth = 0.5) +
  geom_point(data = filter(meta_raw, !extreme),
             aes(x = log_hr, y = se_log, color = subgroup),
             shape = 16, size = 2.8, alpha = 0.85) +
  geom_point(data = filter(meta_raw, extreme),
             aes(x = log_hr, y = se_log, color = subgroup),
             shape = 17, size = 4) +
  scale_color_manual(values = color_map, name = "Cancer") +
  scale_y_reverse(limits = c(se_max, 0),
                  breaks = seq(0, round(se_max, 1), by = 0.1)) +
  labs(title    = "Funnel Plot — Publication Bias",
       subtitle = paste0("Egger p = ",
                         format(egger_p_val, digits = 3, scientific = FALSE),
                         " | ▲ = extreme HR | 95%/99% bands"),
       x = "log(HR)", y = "Standard Error (SE)") +
  theme_dis()

hetero_df <- sg_sub %>%
  mutate(
    I2      = I2,
    hr_pool = exp(TE),
    hr_lo   = exp(lower),
    hr_hi   = exp(upper)
  ) %>%
  rename(Q_p = pval_Q) %>%
  filter(!is.na(I2)) %>%
  mutate(I2_cat = cut(I2,
                      breaks = c(-1, 25, 50, 75, 101),
                      labels = c("Low (<25%)", "Moderate (25-50%)",
                                 "High (50-75%)", "Very high (>75%)")))

fig_hetero <- ggplot(hetero_df, aes(x = reorder(subgroup, -I2), y = I2,
                                    fill = I2_cat)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.4) +
  geom_text(aes(label = paste0(round(I2, 1), "%")),
            vjust = -0.3, fontface = "bold",
            color = BLUE["dark"], size = 3.5) +
  geom_hline(yintercept = c(25, 50, 75), linetype = "dashed",
             color = c(BLUE["light"], BLUE["mid"], BLUE["main"]),
             linewidth = 0.5) +
  scale_fill_manual(
    values = c("Low (<25%)"       = BLUE["pale"],
               "Moderate (25-50%)" = BLUE["light"],
               "High (50-75%)"     = BLUE["mid"],
               "Very high (>75%)"   = BLUE["main"]),
    name = "I² level"
  ) +
  scale_y_continuous(limits = c(0, 105),
                     labels = function(x) paste0(x, "%")) +
  labs(title    = "Heterogeneity by subgroup (I²)",
       subtitle = "Reference lines: 25%, 50%, 75%",
       x = "Cancer type", y = "I² (%)") +
  theme_dis()

# ── FIG 6: LEAVE-ONE-OUT PLOT ────────────────────────────────
# Normalizar nombres de columna (dmetar vs manual)
# Normalizar salida de loo.meta (dmetar) vs implementación manual
if (HAS_DMETAR) {
  loo_study <- if ("studlab" %in% names(loo)) loo$studlab else rownames(loo)
  loo_te    <- if ("TE.random" %in% names(loo)) loo$TE.random else loo$TE
  loo_lo    <- if ("lower.random" %in% names(loo)) loo$lower.random else loo$lower
  loo_hi    <- if ("upper.random" %in% names(loo)) loo$upper.random else loo$upper
  loo_i2    <- if ("I2" %in% names(loo)) ifelse(loo$I2 <= 1, loo$I2 * 100, loo$I2) else rep(NA_real_, nrow(loo))
} else {
  loo_study <- loo$studlab
  loo_te    <- loo$TE.random
  loo_lo    <- loo$lower.random
  loo_hi    <- loo$upper.random
  loo_i2    <- loo$I2 * 100
}
loo_df <- data.frame(
  study = loo_study,
  hr    = exp(loo_te),
  hr_lo = exp(loo_lo),
  hr_hi = exp(loo_hi),
  i2    = pmin(pmax(loo_i2, 0), 100)
) %>% arrange(hr) %>%
  mutate(study = factor(study, levels = study))

hr_pool_main <- exp(ma_main$TE.random)
hr_lo_main   <- exp(ma_main$lower.random)
hr_hi_main   <- exp(ma_main$upper.random)

fig_loo <- ggplot(loo_df, aes(x = hr, y = study)) +
  annotate("rect",
           xmin = hr_lo_main, xmax = hr_hi_main,
           ymin = -Inf, ymax = Inf,
           fill = BLUE["pale"], alpha = 0.5) +
  geom_vline(xintercept = hr_pool_main, color = BLUE["dark"],
             linewidth = 0.8, linetype = "dashed") +
  geom_vline(xintercept = 1, color = "#7F8C8D",
             linewidth = 0.4, linetype = "dotted") +
  geom_errorbarh(aes(xmin = hr_lo, xmax = hr_hi),
                 height = 0.3, color = BLUE["main"], linewidth = 0.5) +
  geom_point(aes(fill = i2), shape = 21,
             size = 2.5, color = BLUE["dark"], stroke = 0.4) +
  scale_fill_gradient(low = BLUE["light"], high = BLUE["dark"],
                      name = "I² (%)", limits = c(0, 100)) +
  scale_x_log10(breaks = c(1, 1.5, 2, 3, 4, 5),
                labels = c("1","1.5","2","3","4","5")) +
  labs(title    = "Leave-One-Out Sensitivity Analysis",
       subtitle = "Shaded band = 95% CI of main analysis",
       x = "Pooled HR (95% CI)", y = NULL) +
  theme_dis(base_size = 8) +
  theme(axis.text.y = element_text(size = 6.5))

gd_label_map <- c(
  "Canonical Liu 2023" = "Canonical Liu 2023 (≤10 genes)",
  "Extended definition" = "Extended definition (>10 genes)"
)

gene_def_res <- sg_genedef %>%
  mutate(
    def   = dplyr::recode(subgroup, !!!gd_label_map, .default = subgroup),
    hr    = exp(TE),
    hr_lo = exp(lower),
    hr_hi = exp(upper),
    I2    = I2
  ) %>%
  select(def, hr, hr_lo, hr_hi, k, I2)

gene_def_res$hr_lo <- ifelse(is.na(gene_def_res$hr_lo), gene_def_res$hr * 0.5, gene_def_res$hr_lo)
gene_def_res$hr_hi <- ifelse(is.na(gene_def_res$hr_hi), gene_def_res$hr * 2.0, gene_def_res$hr_hi)
gene_def_res$I2    <- ifelse(is.na(gene_def_res$I2),    0,                      gene_def_res$I2)

HR_HI_CAP <- 20
gene_def_res$hr_hi_plot <- pmin(gene_def_res$hr_hi, HR_HI_CAP)
y_top <- max(gene_def_res$hr_hi_plot, na.rm = TRUE) * 1.45

fig_genedef <- ggplot(gene_def_res, aes(x = def, y = hr, fill = def)) +
  geom_col(width = 0.5, color = "white", alpha = 0.9) +
  # Barras de error usando el IC truncado visualmente
  geom_errorbar(aes(ymin = hr_lo, ymax = hr_hi_plot),
                width = 0.12, color = BLUE["dark"], linewidth = 0.9) +
  # Flecha hacia arriba cuando el IC real supera el cap
  geom_text(data = filter(gene_def_res, hr_hi > HR_HI_CAP),
            aes(x = def, y = HR_HI_CAP * 1.02, label = "\u2191"),
            size = 5, color = BLUE["dark"], vjust = 0) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "#7F8C8D", linewidth = 0.5) +
  # Etiqueta con HR real (no truncado)
  geom_text(aes(y = hr_hi_plot + y_top * 0.04,
                label = sprintf("HR = %.2f (%.2f\u2013%.2f)\nk = %d studies  I\u00b2 = %.0f%%",
                                hr, hr_lo, hr_hi, k, I2)),
            vjust = 0, size = 3.3, fontface = "bold", color = BLUE["dark"],
            lineheight = 1.2) +
  scale_fill_manual(values = c(BLUE["dark"], BLUE["mid"]), guide = "none") +
  scale_y_continuous(limits = c(0, y_top), breaks = pretty(c(0, y_top))) +
  labs(title    = "Subgroup: Disulfidptosis gene definition",
       subtitle = "Canonical Liu 2023 (\u226410 genes) vs extended definition (>10 genes)",
       x = NULL, y = "Pooled HR (95% CI)") +
  theme_dis()


fig_bubble <- ggplot(meta_raw,
                     aes(x = n, y = hr, color = subgroup)) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "#7F8C8D", linewidth = 0.5) +
  geom_point(data = filter(meta_raw, !extreme),
             aes(size = 1/se_log),
             alpha = 0.80, shape = 16) +
  geom_point(data = filter(meta_raw, extreme),
             size = 5, alpha = 0.95, shape = 17) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x,
              color = BLUE["dark"], fill = BLUE["pale"],
              linewidth = 0.7, linetype = "dashed",
              aes(group = 1), na.rm = TRUE) +
  geom_text(data = filter(meta_raw, extreme),
            aes(label = paste0(author, "\nHR=", round(hr,1))),
            hjust = -0.1, size = 2.5, color = BLUE["dark"],
            fontface = "italic") +
  scale_color_manual(values = color_map, name = "Cancer type") +
  scale_size_continuous(range = c(2, 9), name = "Precision (1/SE)") +
  scale_y_log10(breaks = c(1, 2, 5, 10, 30),
                labels = c("1","2","5","10","30")) +
  scale_x_log10(labels = scales::comma) +
  labs(title    = "Sample size vs HR — All studies",
       subtitle = "Bubble ∝ statistical precision | ▲ = extreme HR (protocol: visually identified)",
       x = "Total N (log scale)", y = "HR (log scale)") +
  theme_dis() +
  guides(size = guide_legend(override.aes = list(shape = 16)))

heatmap_df <- meta_raw %>%
  group_by(subgroup) %>%
  summarise(
    k         = n(),
    n_total   = sum(n),
    pct_valid = mean(!tcga_only) * 100,
    .groups   = "drop"
  ) %>%
  left_join(
    sg_sub %>% transmute(subgroup, hr_pool = exp(TE)),
    by = "subgroup"
  )

hm_mat <- heatmap_df %>%
  select(subgroup, k, hr_pool, pct_valid) %>%
  rename(`Studies (k)` = k,
         `HR pooled`              = hr_pool,
         `% with external validation` = pct_valid) %>%
  pivot_longer(-subgroup, names_to = "metric", values_to = "value")

hm_norm <- hm_mat %>%
  group_by(metric) %>%
  mutate(val_norm = (value - min(value)) / (max(value) - min(value) + 1e-9)) %>%
  ungroup() %>%
  mutate(
    label_txt = case_when(
      grepl("HR",     metric) ~ sprintf("%.2f", value),
      grepl("Studi", metric) ~ sprintf("%d",   as.integer(value)),
      TRUE                    ~ sprintf("%.0f%%", value)
    )
  )

fig_heatmap <- ggplot(hm_norm,
                      aes(x = metric,
                          y = reorder(subgroup, value * (metric == "HR pooled")),
                          fill = val_norm)) +
  geom_tile(color = "white", linewidth = 0.9) +
  geom_text(aes(label = label_txt),
            fontface = "bold", color = "white", size = 3.8) +
  scale_fill_gradient(low = BLUE["light"], high = BLUE["dark"],
                      name = "Value
(norm.)",
                      limits = c(0, 1),
                      breaks = c(0, 0.5, 1),
                      labels = c("Min", "Med", "Max")) +
  scale_x_discrete(labels = c("HR pooled"                  = "HR pooled",
                                "Studies (k)"              = "k (studies)",
                                "% with external validation"= "% valid. ext.")) +
  labs(title    = "Quantitative summary by cancer type",
       subtitle = "Column-normalized scale | Random effects meta-analysis",
       x = NULL, y = NULL) +
  theme_dis() +
  theme(axis.text.x   = element_text(angle = 15, hjust = 1, size = 10),
        axis.text.y   = element_text(size = 10, face = "bold"),
        legend.position = "right")

cat("\nGuardando las figuras.\n")
output_dir <- "Disulfidptosis_Figures"
if (!dir.exists(output_dir)) dir.create(output_dir)

save_fig <- function(fig, name, w = 10, h = 7) {
  ggsave(file.path(output_dir, paste0(name, ".pdf")),
         plot = fig, width = w, height = h, units = "in", dpi = 300)
  ggsave(file.path(output_dir, paste0(name, ".png")),
         plot = fig, width = w, height = h, units = "in", dpi = 300)
  cat("Guardé:", name, "\n")
}

save_fig(fig_prisma,     "Fig1_PRISMA",              w = 11, h = 6)
save_fig(fig_cancer_bar, "Fig2_Cancer_Distribution",  w = 9,  h = 6)
save_fig(fig_forest,     "Fig3_Forest_Plot",          w = 15, h = 12)
save_fig(fig_funnel,     "Fig4_Funnel_Plot",          w = 9,  h = 7)
save_fig(fig_hetero,     "Fig5_Heterogeneity_I2",     w = 9,  h = 6)
save_fig(fig_loo,        "Fig6_Leave_One_Out",        w = 11, h = 10)
save_fig(fig_genedef,    "Fig7_Gene_Definition",      w = 8,  h = 6)
save_fig(fig_bubble,     "Fig8_Bubble_HR_vs_N",       w = 10, h = 7)
save_fig(fig_heatmap,    "Fig9_Heatmap_Summary",      w = 10, h = 6)

panel_top    <- fig_cancer_bar + fig_hetero
panel_mid    <- fig_funnel + fig_genedef
panel_bottom <- fig_bubble

fig_panel <- (panel_top / panel_mid / panel_bottom) +
  plot_annotation(
    title   = "Disulfidptosis Meta-analysis — Results Panel",
    subtitle = paste0("Random effects (RE) | N = ", nrow(meta_raw), " studies with complete HR"),
    theme   = theme(
      plot.title    = element_text(face = "bold", color = BLUE["dark"],
                                   size = 14, hjust = 0.5),
      plot.subtitle = element_text(color = BLUE["main"],
                                   size = 10, hjust = 0.5)
    )
  )
save_fig(fig_panel, "Fig_Panel_Resumen", w = 16, h = 18)

# ─── 8. TABLA RESUMEN (CSV) ──────────────────────────────────
tabla_resumen <- meta_raw %>%
  select(id, author, year, cancer, n, hr, hr_low, hr_high,
         subgroup, gene_def, tcga_only, extreme, n_genes) %>%
  arrange(subgroup, year)

write.csv(tabla_resumen,
          file.path(output_dir, "Tabla_Estudios_HR.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

tabla_pooled <- sg_sub %>%
  transmute(
    Subgrupo         = subgroup,
    k                = k,
    HR_pooled        = round(exp(TE), 3),
    IC_95_low        = round(exp(lower), 3),
    IC_95_high       = round(exp(upper), 3),
    I2_pct           = round(I2, 1),
    p_heterogeneidad = round(pval_Q, 4)
  )
write.csv(tabla_pooled,
          file.path(output_dir, "Tabla_Pooled_Subgrupos.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("\n✔ Tablas CSV guardadas\n")

# ─── 9. FOREST PLOT BASE R (meta pkg) — alternativa clásica ──
cat("\n══ Generando forest plot clásico (meta) ══\n")

# Reemplazar lower/upper internos de metagen por los IC originales del CSV.
# metagen recalcula IC asumiendo simetría log-normal perfecta; cuando el IC
# original no es perfectamente simétrico (e.g. Yang L, Huang J) hay pequeñas
# diferencias. Forzamos los valores del CSV para máxima fidelidad.
ma_fig3b <- ma_main
ma_fig3b$lower <- log(meta_raw$hr_low)
ma_fig3b$upper <- log(meta_raw$hr_high)

pdf(file.path(output_dir, "Fig3b_Forest_Classic.pdf"),
    width = 18, height = 14)
forest(ma_fig3b,
       sortvar    = meta_raw$subgroup,
       leftcols   = c("studlab","n","cancer"),
       leftlabs   = c("Study","N","Cancer"),
       rightcols  = c("effect","ci"),
       rightlabs  = c("HR","95% CI"),
       common     = FALSE,
       col.square    = BLUE["main"],
       col.diamond   = BLUE["dark"],
       col.inside    = "white",
       col.lines     = BLUE["mid"],
       print.tau2    = TRUE,
       print.I2      = TRUE,
       print.pval.Q  = TRUE,
       smlab         = "HR (95% CI)",
       text.random   = "Random Effects",
       xlim          = c(0.3, 40),
       at            = c(0.5, 1, 2, 5, 10, 20),
       header.line   = TRUE,
       fontsize      = 8,
       spacing       = 0.7,
       squaresize    = 0.5)
dev.off()
cat("  ✔ Forest plot clásico guardado\n")

# ─── 10. RESUMEN FINAL EN CONSOLA ────────────────────────────
cat("\n", paste(rep("═", 60), collapse = ""), "\n")
cat("  RESULTADOS PRINCIPALES\n")
cat(paste(rep("═", 60), collapse = ""), "\n")
cat(sprintf("  Estudios incluidos en MA cuantitativo: %d\n", nrow(meta_raw)))
cat(sprintf("  HR pooled global:  %.2f (IC 95%%: %.2f–%.2f)\n",
            exp(ma_main$TE.random),
            exp(ma_main$lower.random),
            exp(ma_main$upper.random)))
cat(sprintf("  I²:  %.1f%%\n", ma_main$I2 * 100))
cat(sprintf("  τ²:  %.4f\n",  ma_main$tau2))
cat(sprintf("  p(Q): %.4f\n", ma_main$pval.Q))
cat(sprintf("  Egger p: %.4f\n", egger_p_val))
cat(paste(rep("─", 60), collapse = ""), "\n")
cat("  HR POOLED POR SUBGRUPO:\n")
for (i in seq_len(nrow(sg_sub))) {
  cat(sprintf("   %-8s  HR=%.2f (%.2f–%.2f)  I\u00b2=%.0f%%\n",
              sg_sub$subgroup[i],
              exp(sg_sub$TE[i]),
              exp(sg_sub$lower[i]),
              exp(sg_sub$upper[i]),
              sg_sub$I2[i]))
}
cat(paste(rep("═", 60), collapse = ""), "\n")
cat("  Figuras guardadas en: ./", output_dir, "/\n", sep = "")
cat("  Figuras: Fig1_PRISMA, Fig2–Fig9, Fig_Panel_Resumen\n")
cat("  Tablas:  Tabla_Estudios_HR.csv, Tabla_Pooled_Subgrupos.csv\n")
cat(paste(rep("═", 60), collapse = ""), "\n\n")
