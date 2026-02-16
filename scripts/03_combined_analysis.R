# scripts/03_composite_scores_and_sti.R
# Combined physiology + yield: run-adjusted z-scores + weighted composite ranking + yield STI
# Inputs:
#   data/agronomydata_copy.xlsx
# Requires objects from scripts 01 and 02 already in memory:
#   - photo_run1_wide, photo_run2_wide
#   - leafspec_run1, leafspec_run2
# Outputs:
#   outputs/composite_scores_weighted.csv
#   outputs/yield_sti.csv

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(lme4)
})

# -----------------------------
# Config
# -----------------------------
agron_file <- file.path("data", "agronomydata_copy.xlsx")
stopifnot(file.exists(agron_file))

# -----------------------------
# Helpers
# -----------------------------
extract_list_first_numeric <- function(x) {
  if (!is.list(x)) return(x)
  sapply(x, function(v) if (length(v) >= 1) v[1] else NA_real_)
}

# -----------------------------
# Yield data (Mean_Yield per Genotype x Treatment)
# -----------------------------
yield_sheets <- list(
  list(sheet = 16, treatment = "Saline"),
  list(sheet = 17, treatment = "Control")
)

yield_data <- bind_rows(lapply(yield_sheets, function(info) {
  read_excel(agron_file, sheet = info$sheet, skip = 2, col_names = FALSE) %>%
    setNames(c("Tank", "Section", "Plant", "Accession", "Name", "# of Plants",
               "Pod Count", "Seed Wt. (g)", "# of diseased plants", "Seed wt per pod (g)",
               "...11", "...12", "Samples Remaining", "257", "Treatment_col",
               "Rep", "Pot", "Approximate Shattering Loss (%)", "Shatter Adjusted Seed wt per pod (g)")) %>%
    mutate(
      Treatment = info$treatment,
      `# of Plants` = as.numeric(`# of Plants`),
      `Seed Wt. (g)` = as.numeric(`Seed Wt. (g)`),
      Yield_per_Plant = if_else(
        Treatment == "Saline" & !is.na(`# of Plants`) & `# of Plants` > 0,
        `Seed Wt. (g)` / `# of Plants`,
        `Seed Wt. (g)`
      )
    )
})) %>%
  group_by(Accession, Treatment) %>%
  summarise(Mean_Yield = mean(Yield_per_Plant, na.rm = TRUE), .groups = "drop") %>%
  rename(Genotype = Accession) %>%
  filter(!is.na(Mean_Yield) & !is.nan(Mean_Yield))

# -----------------------------
# Yield STI
# -----------------------------
yield_control <- yield_data %>%
  filter(Treatment == "Control") %>%
  select(Genotype, Mean_Yield) %>%
  rename(Yield_Control = Mean_Yield)

yield_saline <- yield_data %>%
  filter(Treatment == "Saline") %>%
  select(Genotype, Mean_Yield) %>%
  rename(Yield_Saline = Mean_Yield)

yield_sti <- inner_join(yield_control, yield_saline, by = "Genotype")

mean_yield_control <- mean(yield_control$Yield_Control, na.rm = TRUE)

yield_sti <- yield_sti %>%
  mutate(STI = (Yield_Saline * Yield_Control) / (mean_yield_control^2)) %>%
  arrange(desc(STI))

# -----------------------------
# Build combined physiology (LeafSpec + PhotosynQ)
# -----------------------------
leafspec_all <- bind_rows(leafspec_run1, leafspec_run2)
photosynq_all_wide <- bind_rows(photo_run1_wide, photo_run2_wide)

leafspec_wide <- leafspec_all %>%
  select(Genotype, Run, Week, Treatment, Title, Value) %>%
  pivot_wider(names_from = Title, values_from = Value)

photosynq_wide <- photosynq_all_wide %>%
  pivot_wider(
    id_cols = c(Genotype, Run, Week, Treatment),
    names_from = Trait,
    values_from = Value
  )

physiology <- full_join(
  leafspec_wide,
  photosynq_wide,
  by = c("Genotype", "Run", "Week", "Treatment")
)

master_data <- physiology %>%
  left_join(yield_data, by = c("Genotype", "Treatment"))

master <- master_data %>%
  mutate(across(where(is.list), extract_list_first_numeric))

# -----------------------------
# Run-adjusted z-scores + Composite score
# -----------------------------
trait_cols <- master %>%
  select(where(is.numeric)) %>%
  select(-any_of(c("Mean_Yield"))) %>%
  names()

adj_list <- list()

for (trait in trait_cols) {
  rows_with_data <- !is.na(master[[trait]])
  df_subset <- master[rows_with_data, ]
  
  if (nrow(df_subset) == 0) next
  
  model <- tryCatch(
    lmer(as.formula(paste(trait, "~ (1|Run)")), data = df_subset, REML = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(model)) next
  
  residuals_vec <- resid(model)
  if (length(residuals_vec) != nrow(df_subset)) next
  
  adj_values <- rep(NA_real_, nrow(master))
  adj_values[rows_with_data] <- residuals_vec
  adj_list[[paste0(trait, "_adj")]] <- adj_values
}

master <- bind_cols(master, as.data.frame(adj_list))

adj_trait_cols <- names(adj_list)

z_scores <- master %>%
  select(all_of(adj_trait_cols)) %>%
  mutate(across(everything(), ~ scale(.x)[, 1]))

colnames(z_scores) <- paste0(adj_trait_cols, "_z")
master <- bind_cols(master, z_scores)

master <- master %>%
  rowwise() %>%
  mutate(
    Composite_Score = mean(c_across(all_of(paste0(adj_trait_cols, "_z"))), na.rm = TRUE)
  ) %>%
  ungroup()

# -----------------------------
# Weighted composite (0.4 Week 3, 0.6 Week 6)
# -----------------------------
composite_scores <- master %>%
  select(Genotype, Run, Week, Treatment, Composite_Score)

week3_scores <- composite_scores %>%
  filter(Week == "Week 3") %>%
  select(Genotype, Run, Treatment, Composite_Score) %>%
  rename(Composite_Week3 = Composite_Score)

week6_scores <- composite_scores %>%
  filter(Week == "Week 6") %>%
  select(Genotype, Run, Treatment, Composite_Score) %>%
  rename(Composite_Week6 = Composite_Score)

combined_weighted <- full_join(
  week3_scores,
  week6_scores,
  by = c("Genotype", "Run", "Treatment")
) %>%
  mutate(
    Composite_Score_Weighted = case_when(
      !is.na(Composite_Week3) & !is.na(Composite_Week6) ~ 0.4 * Composite_Week3 + 0.6 * Composite_Week6,
      !is.na(Composite_Week3) &  is.na(Composite_Week6) ~ Composite_Week3,
      is.na(Composite_Week3) & !is.na(Composite_Week6) ~ Composite_Week6,
      TRUE ~ NA_real_
    )
  ) %>%
  arrange(Run, Treatment, desc(Composite_Score_Weighted))

# -----------------------------
# Save outputs
# -----------------------------
dir.create("outputs", showWarnings = FALSE)

write.csv(
  combined_weighted,
  file.path("outputs", "composite_scores_weighted.csv"),
  row.names = FALSE
)

write.csv(
  yield_sti,
  file.path("outputs", "yield_sti.csv"),
  row.names = FALSE
)

message("Composite scores + STI complete. Outputs saved to /outputs.")
