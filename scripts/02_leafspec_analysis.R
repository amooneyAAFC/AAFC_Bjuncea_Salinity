# scripts/02_leafspec_analysis.R
# LeafSpec analysis (selected traits): outliers, t-tests, summary stats, H2, PCA figure
# Inputs:
#   data/master_run_1.xlsx
#   data/master_run_2.xlsx

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(lme4)
  library(emmeans)
  
  library(ggplot2)
  library(ggrepel)
  library(ggnewscale)
  library(FactoMineR)
  library(missMDA)
  library(patchwork)
})

# -----------------------------
# Config
# -----------------------------
run1_file <- file.path("data", "master_run_1.xlsx")
run2_file <- file.path("data", "master_run_2.xlsx")
stopifnot(file.exists(run1_file), file.exists(run2_file))

leafspec_traits <- c("ARI1","ARI2","CCI","CPHLT","CRI1","CRI2","NDVI","WBI")

sheets_run1 <- list(
  list(sheet = 1, week = "Week 3", treatment = "Saline"),
  list(sheet = 2, week = "Week 3", treatment = "Control")
)

sheets_run2 <- list(
  list(sheet = 1, week = "Week 3", treatment = "Saline"),
  list(sheet = 2, week = "Week 3", treatment = "Control"),
  list(sheet = 5, week = "Week 6", treatment = "Saline"),
  list(sheet = 6, week = "Week 6", treatment = "Control")
)

# -----------------------------
# Helpers
# -----------------------------
load_data_leafspec <- function(file, sheet_mapping) {
  map_dfr(sheet_mapping, function(s) {
    raw <- read_excel(file, sheet = s$sheet)
    
    value_cols <- grep("^Value", names(raw), value = TRUE)
    title_cols <- grep("^Title", names(raw), value = TRUE)
    
    if (length(value_cols) != length(title_cols)) {
      stop("Mismatch between number of Title and Value columns in sheet ", s$sheet)
    }
    
    map_dfr(seq_along(value_cols), function(i) {
      raw %>%
        select(Genotype, !!title_cols[i], !!value_cols[i]) %>%
        rename(Title = !!title_cols[i], Value = !!value_cols[i]) %>%
        mutate(
          Value = as.numeric(Value),
          Week = s$week,
          Treatment = s$treatment
        )
    }) %>%
      filter(!is.na(Title), !is.na(Value))
  })
}

remove_outliers_iqr <- function(df, value_col = "Value") {
  x <- df[[value_col]]
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lb <- q1 - 1.5 * iqr
  ub <- q3 + 1.5 * iqr
  df %>% filter(.data[[value_col]] >= lb, .data[[value_col]] <= ub)
}

run_ttests_by_week <- function(df, traits) {
  weeks <- sort(unique(df$Week))
  expand_grid(Week = weeks, Trait = traits) %>%
    pmap_dfr(function(Week, Trait) {
      sub <- df %>% filter(Week == !!Week, Title == !!Trait)
      
      if (nrow(sub) < 2 || n_distinct(sub$Treatment) < 2) {
        return(tibble(Week = Week, Trait = Trait, p_value = NA_real_, t_value = NA_real_,
                      mean_saline = NA_real_, mean_control = NA_real_))
      }
      
      tt <- t.test(Value ~ Treatment, data = sub)
      tibble(
        Week = Week,
        Trait = Trait,
        p_value = tt$p.value,
        t_value = unname(tt$statistic),
        mean_saline  = mean(sub$Value[sub$Treatment == "Saline"], na.rm = TRUE),
        mean_control = mean(sub$Value[sub$Treatment == "Control"], na.rm = TRUE)
      )
    })
}

calculate_H2_leafspec <- function(df, trait) {
  sub <- df %>% filter(Title == trait) %>%
    mutate(
      Genotype = as.factor(Genotype),
      Treatment = as.factor(Treatment)
    )
  
  if (n_distinct(sub$Genotype) < 2 || n_distinct(sub$Treatment) < 2) {
    return(tibble(Trait = trait, H2 = NA_real_, Sigma_G = NA_real_, vdelta_BLUE = NA_real_))
  }
  
  model_re <- tryCatch(lmer(Value ~ Treatment + (1 | Genotype), data = sub),
                       error = function(e) NULL)
  if (is.null(model_re)) {
    return(tibble(Trait = trait, H2 = NA_real_, Sigma_G = NA_real_, vdelta_BLUE = NA_real_))
  }
  
  vc <- as.data.frame(VarCorr(model_re))
  sigma_g2 <- vc$vcov[vc$grp == "Genotype"][1]
  
  model_fe <- tryCatch(lm(Value ~ Treatment + Genotype, data = sub),
                       error = function(e) NULL)
  if (is.null(model_fe)) {
    return(tibble(Trait = trait, H2 = NA_real_, Sigma_G = sigma_g2, vdelta_BLUE = NA_real_))
  }
  
  emm <- tryCatch(as.data.frame(summary(emmeans(model_fe, ~ Genotype))),
                  error = function(e) NULL)
  if (is.null(emm) || !"SE" %in% names(emm)) {
    return(tibble(Trait = trait, H2 = NA_real_, Sigma_G = sigma_g2, vdelta_BLUE = NA_real_))
  }
  
  v_delta <- 2 * mean(emm$SE^2, na.rm = TRUE)
  H2 <- sigma_g2 / (sigma_g2 + v_delta / 2)
  
  tibble(Trait = trait, H2 = H2, Sigma_G = sigma_g2, vdelta_BLUE = v_delta)
}

# -----------------------------
# Load + clean: Run 1
# -----------------------------
leafspec_run1 <- load_data_leafspec(run1_file, sheets_run1) %>%
  mutate(Run = "Run1") %>%
  filter(Title %in% leafspec_traits) %>%
  select(Genotype, Run, Week, Treatment, Title, Value) %>%
  group_by(Title) %>%
  group_modify(~ remove_outliers_iqr(.x, "Value")) %>%
  ungroup()

# -----------------------------
# Load + clean: Run 2
# -----------------------------
leafspec_run2 <- load_data_leafspec(run2_file, sheets_run2) %>%
  mutate(Run = "Run2") %>%
  filter(Title %in% leafspec_traits) %>%
  select(Genotype, Run, Week, Treatment, Title, Value) %>%
  group_by(Week, Treatment, Title) %>%
  group_modify(~ remove_outliers_iqr(.x, "Value")) %>%
  ungroup()

# -----------------------------
# Combined dataset
# -----------------------------
leafspec_all <- bind_rows(leafspec_run1, leafspec_run2)

# -----------------------------
# Stats outputs
# -----------------------------
ttests <- leafspec_all %>%
  group_by(Run) %>%
  group_modify(~ run_ttests_by_week(.x, leafspec_traits)) %>%
  ungroup()

summary_stats <- leafspec_all %>%
  group_by(Run, Week, Title, Treatment) %>%
  summarise(
    n = n(),
    min = min(Value, na.rm = TRUE),
    max = max(Value, na.rm = TRUE),
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    .groups = "drop"
  )

heritability <- map_dfr(leafspec_traits, ~ calculate_H2_leafspec(leafspec_all, .x))

# -----------------------------
# Save outputs
# -----------------------------
dir.create("outputs", showWarnings = FALSE)

write.csv(ttests, file.path("outputs", "leafspec_ttests.csv"), row.names = FALSE)
write.csv(summary_stats, file.path("outputs", "leafspec_summary_stats.csv"), row.names = FALSE)
write.csv(heritability, file.path("outputs", "leafspec_heritability.csv"), row.names = FALSE)

message("LeafSpec analysis complete. Outputs saved to /outputs.")

# -----------------------------
# PCA figure (manuscript style)
# -----------------------------
plot_leafspec_pca_style <- function(data, week_name = NULL, arrow_scale = 8, label_size = 5) {
  d <- data
  if (!is.null(week_name)) {
    d <- d %>% filter(Week == week_name)
  }
  
  wide <- d %>%
    group_by(Genotype, Treatment, Title) %>%
    summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Title, values_from = Value)
  
  meta <- wide %>% select(Genotype, Treatment)
  traits <- wide %>% select(where(is.numeric)) %>% mutate(across(everything(), as.numeric))
  traits <- traits[, !sapply(traits, function(x) any(is.nan(x) | is.infinite(x)))]
  
  ncp <- min(estim_ncpPCA(traits, method = "Regularized")$ncp, ncol(traits))
  imp <- imputePCA(traits, ncp = ncp)$completeObs
  pca <- PCA(imp, graph = FALSE)
  
  ind_df <- as.data.frame(pca$ind$coord) %>% mutate(Treatment = meta$Treatment)
  var_df <- as.data.frame(pca$var$coord) %>%
    rownames_to_column("Variable") %>%
    mutate(
      Contribution = rowSums(pca$var$contrib[, 1:2]),
      xend = Dim.1 * arrow_scale,
      yend = Dim.2 * arrow_scale
    )
  
  x_limits <- range(c(ind_df$Dim.1, var_df$xend)) * 1.2
  y_limits <- range(c(ind_df$Dim.2, var_df$yend)) * 1.2
  
  ggplot() +
    stat_ellipse(data = ind_df, aes(Dim.1, Dim.2, color = Treatment),
                 level = 0.75, size = 1) +
    scale_color_manual(values = c(Control = "deeppink4", Saline = "orange2"),
                       name = "Treatment") +
    geom_point(data = ind_df, aes(Dim.1, Dim.2, shape = Treatment),
               color = "grey86", size = 2, alpha = 0.7) +
    scale_shape_manual(values = c(Control = 17, Saline = 16), name = "Treatment") +
    ggnewscale::new_scale_color() +
    geom_segment(
      data = var_df,
      aes(x = 0, y = 0, xend = xend, yend = yend, color = Contribution),
      arrow = arrow(length = unit(0.25, "cm")),
      size = 0.8
    ) +
    geom_text_repel(
      data = var_df,
      aes(xend, yend, label = Variable),
      color = "black",
      size = label_size,
      max.overlaps = Inf,
      force = 20,
      nudge_x = 0.2,
      nudge_y = 0.2,
      segment.size = 0.4,
      box.padding = 0.3,
      point.padding = 0.3
    ) +
    scale_color_gradientn(colors = c("steelblue4", "#008ECE", "#54BFB7"),
                          name = "Contribution") +
    labs(
      x = paste0("PC1 (", round(pca$eig[1, 2], 1), "%)"),
      y = paste0("PC2 (", round(pca$eig[2, 2], 1), "%)")
    ) +
    coord_cartesian(xlim = x_limits, ylim = y_limits, expand = FALSE) +
    theme_classic(base_size = 14) +
    theme(
      text = element_text(family = "Times New Roman"),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(size = 14, color = "black"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      plot.background = element_blank(),
      panel.border = element_blank()
    )
}

pca_week3 <- plot_leafspec_pca_style(leafspec_all, week_name = "Week 3") + theme(legend.position = "none")
pca_week6 <- plot_leafspec_pca_style(leafspec_all, week_name = "Week 6") + theme(legend.position = "none")
pca_combined <- plot_leafspec_pca_style(leafspec_all)

pca_leafspec_panel <- pca_week3 + pca_week6 + pca_combined +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("A", "B", "C"))) &
  theme(legend.position = "right")

pca_leafspec_panel <- pca_leafspec_panel +
  plot_annotation(theme = theme(plot.background = element_rect(color = "black", fill = NA, size = 2)))

ggsave(
  filename = file.path("outputs", "Fig7_LeafSpec_PCA.jpeg"),
  plot = pca_leafspec_panel,
  width = 16, height = 8, dpi = 300
)
