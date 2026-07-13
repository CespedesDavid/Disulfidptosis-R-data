# Paquetes requeridos
required_packages <- c("ggplot2", "dplyr", "tidyr", "forcats", "scales",
                       "patchwork", "ggtext", "showtext")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, quiet = TRUE)
}
invisible(lapply(required_packages, install_if_missing))

library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(scales)
library(patchwork)
library(ggtext)

# Fuente de Google
if (requireNamespace("showtext", quietly = TRUE)) {
  library(showtext)
  font_add_google("Inter", "inter")
  showtext_auto()
  base_font <- "inter"
} else {
  base_font <- "sans"
}

# Datos

# Conteos por dominio
domain_data <- data.frame(
  Domain = c("D1: Participation", "D2: Attrition",
             "D3: Prognostic Factor", "D4: Outcome",
             "D5: Confounding", "D6: Analysis & Reporting",
             "Overall Risk"),
  Low      = c(78, 75, 53, 78, 70, 28, 52),
  Moderate = c( 0,  0, 18,  0,  6, 37, 18),
  High     = c( 0,  0,  7,  0,  2, 13,  8),
  Unclear  = c( 0,  3,  0,  0,  0,  0,  0)
)

# Año de publicación
year_data <- data.frame(
  Year  = c(2023, 2024, 2025, 2026),
  Count = c(19, 35, 19, 5)
)

# Tipo de cáncer (se agrupan los tipos poco frecuentes)
cancer_raw <- data.frame(
  CancerType = c("HCC", "LUAD", "Colon/Colorectal", "Glioma/GBM/LGG",
                 "Bladder", "Gastric", "ccRCC/KIRC", "NSCLC",
                 "Pan-cancer", "Ovarian", "Cervical", "Pancreatic",
                 "ESCC", "Melanoma", "HNSCC", "Other"),
  Count = c(15, 12, 8, 9, 4, 4, 4, 2, 4, 2, 2, 3, 2, 2, 2, 3)
)

# Riesgo general
overall_data <- data.frame(
  Risk  = c("Low", "Moderate", "High"),
  Count = c(52, 18, 8),
  Pct   = c(52/78, 18/78, 8/78)
)

# Paleta de azules
blue_low      <- "#1A6FA8"   # azul técnico profundo → Low
blue_moderate <- "#56B4D3"   # azul cian medio       → Moderate
blue_high     <- "#C8E6F5"   # azul cielo pálido     → High
blue_unclear  <- "#8ECAE6"   # azul suave            → Unclear

domain_colors <- c("Low" = blue_low, "Moderate" = blue_moderate,
                   "High" = blue_high, "Unclear"  = blue_unclear)

gradient_blues <- c("#0D4F8B", "#1A6FA8", "#2E8EC4", "#56B4D3",
                    "#8ECAE6", "#B9DFF2", "#D9EEF9", "#EAF5FC")

# Tema compartido
theme_tech <- function(base_size = 12) {
  theme_minimal(base_size = base_size, base_family = base_font) +
    theme(
      plot.background    = element_rect(fill = "#F0F6FB", color = NA),
      panel.background   = element_rect(fill = "#F8FCFE", color = NA),
      panel.grid.major   = element_line(color = "#D0E8F5", linewidth = 0.4),
      panel.grid.minor   = element_blank(),
      axis.text          = element_text(color = "#1B3A52", size = rel(0.85)),
      axis.title         = element_text(color = "#1B3A52", face = "bold",
                                        size = rel(0.9)),
      plot.title         = element_text(color = "#0D2B40", face = "bold",
                                        size = rel(1.15), margin = margin(b = 6)),
      plot.subtitle      = element_text(color = "#3B6B8C", size = rel(0.85),
                                        margin = margin(b = 10)),
      plot.caption       = element_text(color = "#6A93AB", size = rel(0.7),
                                        hjust = 0, margin = margin(t = 8)),
      legend.background  = element_rect(fill = "#EAF5FC", color = NA),
      legend.key         = element_rect(fill = NA),
      legend.title       = element_text(color = "#1B3A52", face = "bold",
                                        size = rel(0.85)),
      legend.text        = element_text(color = "#1B3A52", size = rel(0.8)),
      strip.text         = element_text(color = "#0D2B40", face = "bold"),
      plot.margin        = margin(14, 14, 10, 14)
    )
}

# Función auxiliar para las etiquetas de porcentaje
pct_label <- function(n, total = 78) paste0(n, "\n(", round(100 * n / total), "%)")

# Figura 1: barras apiladas de riesgo de sesgo por dominio

fig1_long <- domain_data %>%
  pivot_longer(cols = c(Low, Moderate, High, Unclear),
               names_to = "Risk", values_to = "Count") %>%
  filter(Count > 0) %>%
  mutate(
    Risk   = factor(Risk, levels = c("High", "Unclear", "Moderate", "Low")),
    Domain = fct_rev(factor(Domain, levels = domain_data$Domain)),
    Pct    = Count / 78
  )

fig1 <- ggplot(fig1_long, aes(x = Pct, y = Domain, fill = Risk)) +
  geom_col(width = 0.65, position = "stack") +
  geom_text(aes(label = ifelse(Count > 2, paste0(Count), "")),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 3.2,
            family = base_font) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     expand  = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = domain_colors,
                    breaks = c("Low", "Moderate", "High", "Unclear"),
                    guide  = guide_legend(reverse = TRUE)) +
  labs(
    title    = "Figure 1. Risk of Bias per QUIPS Domain",
    subtitle = "78 prognostic studies on disulfidptosis-related gene signatures in cancer",
    x        = "Proportion of studies (%)",
    y        = NULL,
    fill     = "Risk level",
    caption  = "QUIPS = Quality In Prognosis Studies (Hayden et al., Ann Intern Med 2006)."
  ) +
  theme_tech()

# Figura 2: riesgo de sesgo general (dona/polar)

overall_data2 <- overall_data %>%
  mutate(
    Risk   = factor(Risk, levels = c("Low", "Moderate", "High")),
    label  = paste0(Count, "\n(", round(100 * Pct), "%)"),
    ymax   = cumsum(Pct),
    ymin   = lag(ymax, default = 0),
    ymid   = (ymin + ymax) / 2
  )

fig2 <- ggplot(overall_data2, aes(ymax = ymax, ymin = ymin,
                                   xmax = 4, xmin = 2.4, fill = Risk)) +
  geom_rect(color = "white", linewidth = 0.6) +
  geom_text(aes(x = 3.65, y = ymid, label = label),
            color = "white", fontface = "bold", size = 3.8,
            family = base_font) +
  annotate("text", x = 0, y = 0, label = "Overall\nRisk",
           color = "#0D4F8B", fontface = "bold", size = 5,
           family = base_font, hjust = 0.5) +
  coord_polar(theta = "y") +
  xlim(0, 4) +
  scale_fill_manual(values = c("Low" = blue_low,
                                "Moderate" = blue_moderate,
                                "High" = blue_high)) +
  labs(
    title    = "Figure 2. Overall Risk of Bias Distribution",
    subtitle = "n = 78 studies",
    fill     = "Risk level",
    caption  = "Two-thirds of included studies were rated as Low overall risk."
  ) +
  theme_tech() +
  theme(
    axis.text  = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

# Figura 3: estudios por año de publicación

fig3 <- ggplot(year_data, aes(x = factor(Year), y = Count, fill = Count)) +
  geom_col(width = 0.55, show.legend = FALSE) +
  geom_text(aes(label = Count), vjust = -0.5, fontface = "bold",
            color = "#0D4F8B", size = 4.2, family = base_font) +
  scale_fill_gradient(low = "#8ECAE6", high = "#0D4F8B") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "Figure 3. Publication Year of Included Studies",
    subtitle = "Rapid growth in disulfidptosis prognostic research (2023–2026)",
    x        = "Publication Year",
    y        = "Number of Studies",
    caption  = "2026 includes studies published through the search date."
  ) +
  theme_tech() +
  theme(panel.grid.major.x = element_blank())

# Figura 4: tipos de cáncer

cancer_plot <- cancer_raw %>%
  arrange(Count) %>%
  mutate(CancerType = factor(CancerType, levels = CancerType))

fig4 <- ggplot(cancer_plot, aes(x = Count, y = CancerType)) +
  geom_segment(aes(xend = 0, yend = CancerType, color = Count),
               linewidth = 1.4) +
  geom_point(aes(color = Count), size = 4.5) +
  geom_text(aes(label = Count), nudge_x = 0.5, hjust = 0,
            color = "#0D2B40", fontface = "bold", size = 3.2,
            family = base_font) +
  scale_color_gradient(low = "#8ECAE6", high = "#0D4F8B",
                       guide = "none") +
  scale_x_continuous(breaks = seq(0, 16, 4),
                     expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = "Figure 4. Distribution of Cancer Types",
    subtitle = "HCC and LUAD are the most frequently studied tumor types",
    x        = "Number of Studies",
    y        = NULL,
    caption  = "HCC = Hepatocellular carcinoma; LUAD = Lung adenocarcinoma;\nGlioma/GBM/LGG = combined glioma subtypes; ccRCC/KIRC = clear-cell renal carcinoma."
  ) +
  theme_tech() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "#D0E8F5", linewidth = 0.35))

# Exportación de las figuras

output_dir <- "."   # cambiar a la ruta deseada si hace falta

ggsave(file.path(output_dir, "Fig1_RoB_per_domain.pdf"),
       plot = fig1, width = 10, height = 6, dpi = 300, bg = "#F0F6FB")

ggsave(file.path(output_dir, "Fig2_Overall_RoB.pdf"),
       plot = fig2, width = 7, height = 7, dpi = 300, bg = "#F0F6FB")

ggsave(file.path(output_dir, "Fig3_Publication_Year.pdf"),
       plot = fig3, width = 8, height = 5, dpi = 300, bg = "#F0F6FB")

ggsave(file.path(output_dir, "Fig4_Cancer_Types.pdf"),
       plot = fig4, width = 9, height = 7, dpi = 300, bg = "#F0F6FB")

# Panel combinado 2x2
combined <- (fig1 + fig2) / (fig3 + fig4) +
  plot_annotation(
    title   = "QUIPS Risk of Bias Assessment — Disulfidptosis Prognostic Studies",
    subtitle = "Systematic Review · n = 78 studies",
    caption  = "QUIPS = Quality In Prognosis Studies. HCC = Hepatocellular carcinoma; LUAD = Lung adenocarcinoma.",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 15,
                                   color = "#0D2B40", family = base_font),
      plot.subtitle = element_text(size = 10, color = "#3B6B8C",
                                   family = base_font),
      plot.caption  = element_text(size = 8, color = "#6A93AB",
                                   family = base_font),
      plot.background = element_rect(fill = "#E8F3FA", color = NA)
    )
  )

ggsave(file.path(output_dir, "QUIPS_Panel_Combined.pdf"),
       plot = combined, width = 18, height = 13, dpi = 300, bg = "#E8F3FA")

message("Listo, ya guardé las cinco figuras: Fig1_RoB_per_domain.pdf, Fig2_Overall_RoB.pdf, Fig3_Publication_Year.pdf, Fig4_Cancer_Types.pdf y QUIPS_Panel_Combined.pdf.")
