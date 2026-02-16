# scripts/01_photosynq_analysis.R
# PhotosynQ analysis: outliers, t-tests, summary stats, H2
# Inputs:
#   data/master_run1.xlsx
#   data/master_run2.xlsx

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
run1_file <- file.path("data", "master_run1.xlsx")
run2_file <- file.path("data", "master_run2.xlsx")
stopifnot(file.exists(run1_file), file.exists(run2_file))

photosynq_traits_raw <- c(
  "SPAD",
  "FmPrime", "FoPrime", "FvP_over_FmP",
  "qL",
  "ECSt mAU",
  "gH+",
  "vH+",
  "NPQt",
  "Phi2",
  "PhiNPQ",
  "PhiNO",
  "PS1 Active Centers",
  "PS1 Open Centers",
  "PS1 Oxidized Centers",
  "PS1 Over Reduced Centers"
)

photosynq_traits <- make.names(photosynq_traits_raw)

sheets_run1 <- list(
  list(sheet = 3, week = "Week 3", treatment = "Saline"),
  list(sheet = 4, week = "Week 3", treatment = "Control"),
  list(sheet = 5, week = "Week 6", treatment = "Saline"),
  list(sheet = 6, week = "Week 6", treatment = "Control")
)

sheets_run2 <- list(
  list(sheet = 3, week = "Week 2", treatment = "Saline"),
  list(sheet = 4, week = "Week 2", treatment = "Control"),
  list(sheet = 7, week = "Week 5", treatment = "Saline"),
  list(sheet = 8, week = "Week 5", treatment = "Control")
)

# -----------------------------
# Helpers
# -----------------------------
average_if_csv <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) return(NA_real_)
    val <- as.character(val)
    if (grepl(",", val)) {
      nums <- suppressWarnings(as.numeric(strsplit(val, ",")[[1]]))
      return(mean(nums, na.rm = TRUE))
    }
    suppressWarnings(as.numeric(val))
  })
}

load_photosynq <- function(file, sheet_mapping, trait_names_raw) {
  map_dfr(sheet_mapping, function(s) {
    df <- read_excel(file, sheet = s$sheet) %>%
      mutate(Week = s$week, Treatment = s$treatment)
    
    names(df) <- make.names(names(df))
    
    for (tr in make.names(trait_names_raw)) {
      if (tr %in% names(df)) {
        df[[tr]] <- average_if_csv(df[[tr]])
      }
    }
    
    df %>%
      select(any_of(c("Genotype", "Week", "Treatment", make.names(trait_names_raw))))
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

run_ttests_by_week <- function(df) {
  weeks <- sort(unique(df$Week))
  traits <- sort(unique(df$Trait))
  
  expand_grid(Week = weeks, Trait = traits) %>%
    pmap_dfr(function(Week, Trait) {
      sub <- df %>% filter(Week == !!Week, Trait == !!Trait)
      
      if (nrow(sub) < 2 || n_distinct(sub$Treatment) < 2) {
        return(tibble(Week = Week, Trait = Trait, p_value = NA_real_, t_value = NA_real_,
                      mean_saline = NA_real_, mean_control = NA_real_))
      }
      
      vals_saline  <- sub$Value[sub$Treatment == "Saline"]
      vals_control <- sub$Value[sub$Treatment == "Control"]
      
      if (length(unique(vals_saline)) < 2 || length(unique(vals_control)) < 2) {
        return(tibble(
          Week = Week, Trait = Trait,
          p_value = NA_real_, t_value = NA_real_,
          mean_saline = mean(vals_saline, na.rm = TRUE),
          mean_control = mean(vals_control, na.rm = TRUE)
        ))
      }
      
      tt <- tryCatch(t.test(Value ~ Treatment, data = sub, var.equal = FALSE),
      error = function(e) NULL)
      if (is.null(tt)) {
        return(tibble(
          Week = Week, Trait = Trait,
          p_value = NA_real_, t_value = NA_real_,
          mean_saline = mean(vals_saline, na.rm = TRUE),
          mean_control = mean(vals_control, na.rm = TRUE)
        ))
      }
      
      tibble(
        Week = Week,
        Trait = Trait,
        p_value = tt$p.value,
        t_value = unname(tt$statistic),
        mean_saline = mean(vals_saline, na.rm = TRUE),
        mean_control = mean(vals_control, na.rm = TRUE)
      )
    })
}

calculate_H2 <- function(df, trait) {
  sub <- df %>% filter(Trait == trait) %>%
    mutate(
      Genotype = as.factor(Genotype),
      Treatment = as.factor(Treatment),
      Week = as.factor(Week)
    )
  
  if (n_distinct(sub$Genotype) < 2 || n_distinct(sub$Treatment) < 2) {
    return(tibble(Trait = trait, H2 = NA_real_, Sigma_G = NA_real_, vdelta_BLUE = NA_real_))
  }
  
  model_re <- tryCatch(
    lmer(Value ~ Treatment + (1 | Genotype), data = sub),
    error = function(e) NULL
  )
  if (is.null(model_re)) {
    return(tibble(Trait = trait, H2 = NA_real_, Sigma_G = NA_real_, vdelta_BLUE = NA_real_))
  }
  
  vc <- as.data.frame(VarCorr(model_re))
  sigma_g2 <- vc$vcov[vc$grp == "Genotype"][1]
  
  model_fe <- tryCatch(
    lm(Value ~ Treatment + Genotype, data = sub),
    error = function(e) NULL
  )
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
photo_run1_wide <- load_photosynq(run1_file, sheets_run1, photosynq_traits_raw) %>%
  mutate(Run = "Run1")

photo_run1_long <- photo_run1_wide %>%
  pivot_longer(cols = any_of(photosynq_traits), names_to = "Trait", values_to = "Value") %>%
  filter(!is.na(Value)) %>%
  group_by(Trait) %>%
  group_modify(~ remove_outliers_iqr(.x, "Value")) %>%
  ungroup()

# -----------------------------
# Load + clean: Run 2
# -----------------------------
photo_run2_wide <- load_photosynq(run2_file, sheets_run2, photosynq_traits_raw) %>%
  mutate(Run = "Run2")

photo_run2_long <- photo_run2_wide %>%
  pivot_longer(cols = any_of(photosynq_traits), names_to = "Trait", values_to = "Value") %>%
  filter(!is.na(Value)) %>%
  group_by(Trait) %>%
  group_modify(~ remove_outliers_iqr(.x, "Value")) %>%
  ungroup()

# -----------------------------
# Standardize weeks for combined analysis
# -----------------------------
photo_run2_long_std <- photo_run2_long %>%
  mutate(Week = recode(Week, "Week 2" = "Week 3", "Week 5" = "Week 6"))

photo_all <- bind_rows(
  photo_run1_long,
  photo_run2_long_std
) %>% mutate(Week = factor(Week, levels = c("Week 3", "Week 6")))

# -----------------------------
# Outputs
# -----------------------------
ttests_run1 <- run_ttests_by_week(photo_run1_long)
ttests_run2 <- run_ttests_by_week(photo_run2_long)
ttests_all  <- run_ttests_by_week(photo_all)

summary_run1 <- photo_run1_long %>%
  group_by(Week, Trait, Treatment) %>%
  summarise(n = n(), min = min(Value, na.rm = TRUE), max = max(Value, na.rm = TRUE),
            mean = mean(Value, na.rm = TRUE), sd = sd(Value, na.rm = TRUE), .groups = "drop")

summary_run2 <- photo_run2_long %>%
  group_by(Week, Trait, Treatment) %>%
  summarise(n = n(), min = min(Value, na.rm = TRUE), max = max(Value, na.rm = TRUE),
            mean = mean(Value, na.rm = TRUE), sd = sd(Value, na.rm = TRUE), .groups = "drop")

summary_all <- photo_all %>%
  group_by(Week, Trait, Treatment) %>%
  summarise(n = n(), min = min(Value, na.rm = TRUE), max = max(Value, na.rm = TRUE),
            mean = mean(Value, na.rm = TRUE), sd = sd(Value, na.rm = TRUE), .groups = "drop")

heritability_all <- map_dfr(sort(unique(photo_all$Trait)), ~ calculate_H2(photo_all, .x))

# -----------------------------
# Save outputs
# -----------------------------
dir.create("outputs", showWarnings = FALSE)

write.csv(ttests_run1, file.path("outputs", "photosynq_ttests_run1.csv"), row.names = FALSE)
write.csv(ttests_run2, file.path("outputs", "photosynq_ttests_run2.csv"), row.names = FALSE)
write.csv(ttests_all,  file.path("outputs", "photosynq_ttests_combined.csv"), row.names = FALSE)

write.csv(summary_run1, file.path("outputs", "photosynq_summary_run1.csv"), row.names = FALSE)
write.csv(summary_run2, file.path("outputs", "photosynq_summary_run2.csv"), row.names = FALSE)
write.csv(summary_all,  file.path("outputs", "photosynq_summary_combined.csv"), row.names = FALSE)

write.csv(heritability_all, file.path("outputs", "photosynq_heritability_combined.csv"), row.names = FALSE)

message("PhotosynQ analysis complete. Outputs saved to /outputs.")

# -----------------------------
# PCA
# -----------------------------

photo_all_with_reps <- photo_all %>%
  group_by(Genotype, Week, Treatment, Run, Trait) %>%
  mutate(Rep = row_number()) %>%
  ungroup()

combined_photosynq_wide <- photo_all_with_reps %>%
  pivot_wider(
    id_cols = c(Genotype, Week, Treatment, Run, Rep),
    names_from = Trait,
    values_from = Value
  )

plot_photosynq_pca_style <- function(data, week_name = NULL, arrow_scale = 3, label_size = 5) {
  
  d <- data
  if (!is.null(week_name)) {
    d <- d %>% filter(Week == week_name)
  }
  
  d$Treatment <- factor(d$Treatment, levels = c("Control", "Saline"))
  
  meta <- d %>% select(Genotype, Treatment)
  
  traits <- d %>%
    select(all_of(photosynq_traits)) %>%
    mutate(across(everything(), as.numeric))
  
  traits <- traits[, !sapply(traits, function(x) any(is.nan(x) | is.infinite(x)))]
  
  ncp <- min(estim_ncpPCA(traits, method = "Regularized")$ncp, ncol(traits))
  imp <- imputePCA(traits, ncp = ncp)$completeObs
  
  pca <- PCA(imp, graph = FALSE)
  
  ind_df <- as.data.frame(pca$ind$coord) %>%
    mutate(Treatment = meta$Treatment)
  
  var_df <- as.data.frame(pca$var$coord) %>%
    tibble::rownames_to_column("Variable") %>%
    mutate(
      Contribution = rowSums(pca$var$contrib[, 1:2]),
      xend = Dim.1 * arrow_scale,
      yend = Dim.2 * arrow_scale
    )
  
  x_limits <- range(c(ind_df$Dim.1, var_df$xend)) * 1.2
  y_limits <- range(c(ind_df$Dim.2, var_df$yend)) * 1.2
  
  ggplot() +
    stat_ellipse(data = ind_df,
                 aes(Dim.1, Dim.2, color = Treatment),
                 level = 0.75,
                 size = 1) +
    
    geom_point(data = ind_df,
               aes(Dim.1, Dim.2, shape = Treatment),
               color = "grey86",
               size = 2.5,
               alpha = 0.8) +
    
    scale_color_manual(values = c(Control = "deeppink4",
                                  Saline = "orange2"),
                       name = "Treatment") +
    
    scale_shape_manual(values = c(Control = 17,
                                  Saline = 16),
                       name = "Treatment") +
    
    ggnewscale::new_scale_color() +
    
    geom_segment(data = var_df,
                 aes(x = 0, y = 0, xend = xend, yend = yend,
                     color = Contribution),
                 arrow = arrow(length = unit(0.25, "cm")),
                 size = 0.8) +
    
    geom_text_repel(data = var_df,
                    aes(xend, yend, label = Variable),
                    color = "black",
                    size = label_size,
                    max.overlaps = Inf,
                    force = 20,
                    nudge_x = 0.2,
                    nudge_y = 0.2,
                    segment.size = 0.4,
                    box.padding = 0.3,
                    point.padding = 0.3) +
    
    scale_color_gradientn(colors = c("steelblue4",
                                     "#008ECE",
                                     "#54BFB7"),
                          name = "Contribution") +
    
    labs(
      x = paste0("PC1 (", round(pca$eig[1, 2], 1), "%)"),
      y = paste0("PC2 (", round(pca$eig[2, 2], 1), "%)")
    ) +
    
    coord_cartesian(xlim = x_limits,
                    ylim = y_limits,
                    expand = FALSE) +
    
    theme_classic(base_size = 14) +
    theme(
      text = element_text(family = "Times New Roman"),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(size = 14, color = "black"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
}

pca_week3   <- plot_photosynq_pca_style(combined_photosynq_wide, "Week 3") + theme(legend.position = "none")
pca_week6   <- plot_photosynq_pca_style(combined_photosynq_wide, "Week 6") + theme(legend.position = "none")
pca_combined <- plot_photosynq_pca_style(combined_photosynq_wide)

pca_final <- pca_week3 + pca_week6 + pca_combined +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = list(c("A", "B", "C"))) &
  theme(legend.position = "right")

pca_final <- pca_final +
  patchwork::plot_annotation(
    theme = theme(plot.background = element_rect(color = "black",
                                                 fill = NA,
                                                 size = 2))
  )

ggsave(
  filename = file.path("outputs", "pca_PhotosynQ.jpeg"),
  plot = pca_final,
  width = 16,
  height = 8,
  dpi = 300
)

