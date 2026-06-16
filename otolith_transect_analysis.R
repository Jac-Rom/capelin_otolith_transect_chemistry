# =============================================================================
# DOMINANT MIGRATION STRATEGIES WITHIN THE GULF OF ST. LAWRENCE CAPELIN STOCK
# INFERRED FROM OTOLITH CHEMISTRY TRANSECTS
#
# Authors : Romaric Jac, Olivier Le Pape, Elisabeth van Beveren,
#           Mathieu Boudreau, Lola Coussau, Pascal Sirois,
#           Dominique Robert, Pablo Brosset
#
# Journal : Canadian Journal of Fisheries and Aquatic Sciences
#           Small Pelagic Fish Symposium (SPF-2026) special issue
#
# ─────────────────────────────────────────────────────────────────────────────
# DATA AVAILABILITY
# ─────────────────────────────────────────────────────────────────────────────
# Otolith chemistry data are not publicly archived and are available upon
# reasonable request to the corresponding author. The following input files
# are required to run this script:
#
#   transect_capelin_clean.csv       – Cleaned otolith transect data (one row
#                                      per measurement point along each transect).
#                                      Key columns: Individual, Distance_to_core,
#                                      Length_at_dist, Length, Sex, Survey,
#                                      Region_caught, Region_pred, and the seven
#                                      elemental ratios below.
#
#   otolith_margin_clean_region.csv  – Otolith edge measurements used to train
#                                      the QDA model from Jac et al., 2026. Key columns: Region and
#                                      the seven elemental ratios below.
#
#   FL_nGSL.csv                      – Bottom trawl data of the northern Gulf of
#                                      St. Lawrence (n-GSL) capelin. Key columns:
#                                      year, latitude, longitude, length (mm),
#                                      number.caught.
#
#   FL_sGSL.csv                      – Bottom trawl data of the southern Gulf of
#                                      St. Lawrence (s-GSL) capelin. Key columns:
#                                      year, length (cm), number.caught.
#
# Elemental ratio columns used throughout (ICP-MS laser-ablation, Ca-normalised):
#   Li7_Ca, B11_Ca, Mg25_Ca, K39_Ca, Zn64_Ca, Sr88_Ca, Ba138_Ca
#
# ─────────────────────────────────────────────────────────────────────────────
# SCRIPT OVERVIEW
# ─────────────────────────────────────────────────────────────────────────────
# 1.  Libraries & global parameters
# 2.  Load & prepare data
# 3.  QDA model (region assignment from otolith edge chemistry)
# 4.  Correlation plot among elemental ratios (correlation.png)
# 5.  Transect grouping & region attribution across resolutions
# 6.  State-distribution plot  (State_distribution_plot.png)
# 7.  Permanent-transition detection & summary plots
#       – Transition diagram        (Transition_diagram.png)
#       – Barplot n° of transitions (Barplot_nombre_transitions.png)
# 8.  Sensitivity analysis on transition threshold
#       – Total transitions vs threshold   (Total_transitions_vs_threshold.png)
#       – Proportion per class vs threshold(Proportion_vs_threshold.png)
# 9.  Filtering individuals for clustering
#       – Undetermined-proportion histogram (undetermined.png)
# 10. Multiple Correspondence Analysis (MCA)
#       – Variance barplot (variances.png)
# 11. k-means clustering on MCA coordinates
#       – Cluster scatter plot       (Clusters.png)
#       – Region-per-cluster barplot (Barplot_clusters.png)
#       – Sex-per-cluster barplot    (Barplot_clusters_sex.png)
# 12. Synthetic cluster profiles – bubble plot (plots_cluster_majority_bis/)
# 13. Migration shift analysis (cluster 3)
#       – SBI vs South proportions  (shift_SBI_sud.png)
#       – Fork-length distributions (density_NGSL_SGSL.png)
# =============================================================================


# =============================================================================
# 1.  LIBRARIES & GLOBAL PARAMETERS
# =============================================================================

library(MASS)        # QDA
library(dplyr)       # data manipulation
library(tidyr)       # pivoting
library(ggplot2)     # figures
library(purrr)       # functional iteration
library(readr)       # CSV import
library(cluster)     # clustering utilities
library(factoextra)  # cluster diagnostics (fviz_*)
library(FactoMineR)  # MCA
library(fpc)         # cluster stats
library(ggsci)       # colour palettes
library(colorspace)  # lighten()
library(ggnewscale)  # new_scale_fill() in ggplot2
library(ggforce)     # geom_circle (loaded for compatibility)
library(nnet)        # multinom()
library(patchwork)   # panel figures

# ── Elemental ratio columns used in all QDA calls ────────────────────────────
CHEM_COLS <- c("Li7_Ca", "B11_Ca", "Mg25_Ca", "K39_Ca",
               "Zn64_Ca", "Sr88_Ca", "Ba138_Ca")

# ── Region colour palette (used in every figure) ─────────────────────────────
REGION_COLORS <- c(
  "estuary"      = "#FF9900",
  "strait"       = "#FE0000",
  "south"        = "#5D97EB",
  "unidentified" = "#CA92DF"
)
REGION_LEVELS <- c("unidentified", "strait", "south", "estuary")

# ── Resolution grid (grouping steps for the transect binning) ────────────────
RESOLUTIONS <- list(
  list(n_group = 1, label = "5",  dir = "transect_plots_5_clean"),
  list(n_group = 2, label = "10", dir = "transect_plots_10_clean"),
  list(n_group = 3, label = "15", dir = "transect_plots_15_clean"),
  list(n_group = 4, label = "20", dir = "transect_plots_20_clean"),
  list(n_group = 8, label = "40", dir = "transect_plots_40_clean")
)

# ── Analysis resolution (15 µm bins, chosen after sensitivity checks) ─────────
TARGET_RES <- "15"

# ── Transition-detection threshold (minimum consecutive observations) ─────────
MIN_RUN_THRESHOLD <- 6

# ── Length cut-off for clustering (same for all sexes, 70th-percentile filter
#    applied separately; hard cut at 12.9 cm) ─────────────────────────────────
LENGTH_CUTOFF <- 12.9

# ── Minimum number of individuals per length bin to keep a bin ───────────────
MIN_IND_PER_BIN <- 11

# ── Number of k-means clusters ────────────────────────────────────────────────
N_CLUSTERS <- 4

# ── Shift breakpoints (fork length, cm) used to annotate migration figures ───
SHIFT_VALUES <- data.frame(
  label   = c("Shift 1", "Shift 2"),
  cross_x = c(11.16, 12.9)
)

# ── Set seed for reproducibility ─────────────────────────────────────────────
set.seed(123)


# =============================================================================
# 2.  LOAD & PREPARE DATA
# =============================================================================

# ── Raw transect measurements ────────────────────────────────────────────────
transect_raw <- read.csv("transect_capelin_clean.csv",
                         na.strings = "NA", dec = ".")

# ── Otolith edge measurements (training data for QDA) ───────────────────────
#    Exclude the "north" region: it is not represented in the transect dataset
#    and was excluded from the QDA to avoid classification artefacts.
edge_data <- read.csv("otolith_margin_clean_region.csv",
                      na.strings = "NA") %>%
  filter(Region != "north")

# ── Fork-length survey data (for migration figure) ───────────────────────────
n_gsl_raw <- read.csv("FL_nGSL.csv", na.strings = "NA", dec = ".")
s_gsl_raw <- read.csv("FL_sGSL.csv", na.strings = "NA", dec = ".")


# =============================================================================
# 3.  QDA MODEL — REGION ASSIGNMENT FROM OTOLITH EDGE CHEMISTRY
# =============================================================================
# A Quadratic Discriminant Analysis (QDA) is fitted on the labelled edge
# measurements. The trained model is then applied to each measurement point
# along the transects (except the outermost 40 µm which are assigned directly
# from the capture region — ground truth).

# Keep only rows with a known region label
edge_train <- edge_data %>%
  filter(!is.na(Region)) %>%
  dplyr::select(all_of(CHEM_COLS), Region)

# Fit QDA
qda_model <- qda(Region ~ ., data = edge_train)


# =============================================================================
# 4.  CORRELATION PLOT AMONG ELEMENTAL RATIOS  →  correlation.png
# =============================================================================
# Visual check of collinearity between the seven elemental ratios used in the
# QDA. Helps justify the choice of predictors and flag redundancy.

corr_matrix <- cor(edge_train %>% dplyr::select(all_of(CHEM_COLS)),
                   use = "pairwise.complete.obs")

corr_long <- as.data.frame(as.table(corr_matrix)) %>%
  rename(Var1 = Var1, Var2 = Var2, Correlation = Freq)

correlation_plot <- ggplot(corr_long, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(Correlation, 2)), size = 3.5) +
  scale_fill_gradient2(low = "#5D97EB", mid = "white", high = "#FE0000",
                       midpoint = 0, limits = c(-1, 1)) +
  theme_bw(base_size = 14) +
  theme(
    axis.title       = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid       = element_blank(),
    legend.position  = "right",
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.75)
  )

ggsave("correlation.png", correlation_plot,
       width = 1800, height = 1500, dpi = 300, units = "px")


# =============================================================================
# 5.  TRANSECT GROUPING & REGION ATTRIBUTION
# =============================================================================
# Two helper functions are applied to each resolution level:
#   group_transect()          – bins consecutive measurements into n-µm windows
#   predict_transect_region() – assigns a region to each bin via QDA or
#                               ground-truth (outermost 40 µm)

# ── 5a. Grouping function ────────────────────────────────────────────────────
# The first 20 µm (core region) are kept at 1-bin resolution to preserve fine
# early-life chemistry. The remainder is grouped in windows of n_group bins.
# Non-chemical columns use the first value of each window (metadata); chemical
# columns (positions 18 and 22–28) are averaged.

group_transect <- function(transect_df, n_group) {

  # Core (≤ 20 µm): bin width = 1
  core_part <- transect_df %>%
    filter(Distance_to_core <= 20) %>%
    group_by(Individual) %>%
    mutate(bin = row_number() - 1) %>%
    group_by(Individual, bin) %>%
    summarise(
      across(-c(18, 22:28), ~ first(.x)),
      across( c(18, 22:28), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::select(-bin)

  # Rest (> 20 µm): bin width = n_group
  rest_part <- transect_df %>%
    filter(Distance_to_core > 20) %>%
    group_by(Individual) %>%
    mutate(bin = (row_number() - 1) %/% n_group) %>%
    group_by(Individual, bin) %>%
    summarise(
      across(-c(18, 22:28), ~ first(.x)),
      across( c(18, 22:28), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::select(-bin)

  bind_rows(core_part, rest_part) %>%
    arrange(Individual, Distance_to_core)
}

# ── 5b. Region-attribution function ─────────────────────────────────────────
# Outermost 40 µm  → ground truth (Region_caught, i.e. capture location)
# All other points → QDA prediction; points with max posterior < 0.50 are
#                    labelled "unidentified".

predict_transect_region <- function(df) {
  df <- df %>% mutate(Region_assigned = NA_character_)

  is_edge <- df$Distance_to_core >= (max(df$Distance_to_core) - 40)
  is_core <- !is_edge

  # Ground truth for edge
  df$Region_assigned[is_edge] <- df$Region_caught[is_edge]

  # QDA for core + mid-transect
  if (any(is_core)) {
    pred      <- predict(qda_model, df[is_core, CHEM_COLS])
    max_prob  <- apply(pred$posterior, 1, max)
    assigned  <- ifelse(
      max_prob >= 0.5,
      colnames(pred$posterior)[apply(pred$posterior, 1, which.max)],
      "unidentified"
    )
    df$Region_assigned[is_core] <- assigned
  }
  df
}

# ── 5c. Apply both functions across all resolutions ──────────────────────────
transect_by_res <- RESOLUTIONS %>%
  map(function(res) {
    grouped <- group_transect(transect_raw, n_group = res$n_group)
    grouped %>%
      group_by(Individual) %>%
      group_modify(~ predict_transect_region(.x)) %>%
      ungroup()
  }) %>%
  setNames(paste0("res_", c("5", "10", "15", "20", "40")))

# ── 5d. Named shortcuts ──────────────────────────────────────────────────────
transect_5  <- transect_by_res$res_5
transect_10 <- transect_by_res$res_10
transect_15 <- transect_by_res$res_15   # ← main analysis resolution
transect_20 <- transect_by_res$res_20
transect_40 <- transect_by_res$res_40

# ── 5e. Summary tables (assignment counts per resolution) ───────────────────
for (nm in names(transect_by_res)) {
  cat("\n─── Resolution:", nm, "───\n")
  df <- transect_by_res[[nm]]
  print(table(df$Region_assigned))
  print(round(100 * prop.table(table(df$Region_assigned)), 2))
  print(table(df$Region_assigned, df$Region_caught))
  print(df %>% summarise(n_individuals = n_distinct(Individual)))
}

# ── Working dataset (15 µm resolution, excluding the "north" region) ─────────
transect <- transect_15 %>% filter(Region_assigned != "north")


# =============================================================================
# 6.  STATE-DISTRIBUTION PLOT  →  State_distribution_plot.png
# =============================================================================
# For each 0.25-cm fork-length bin, the dominant (modal) region per individual
# is computed; the bar chart shows the proportion of individuals in each region,
# and the overlaid line shows the number of individuals contributing to that bin.
# Bins with fewer than MIN_IND_PER_BIN individuals are excluded.

# Modal region per individual × length bin
ind_region_per_bin <- transect %>%
  mutate(length_bin = floor(Length_at_dist / 0.25) * 0.25) %>%
  group_by(Individual, length_bin) %>%
  summarise(
    Region_assigned = names(which.max(table(Region_assigned))),
    .groups = "drop"
  )

# Proportions per bin
region_prop <- ind_region_per_bin %>%
  count(length_bin, Region_assigned) %>%
  group_by(length_bin) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Number of unique individuals per bin
ind_count_per_bin <- ind_region_per_bin %>%
  group_by(length_bin) %>%
  summarise(n_ind = n_distinct(Individual), .groups = "drop")

# Keep only bins above the minimum-individual threshold
valid_bins <- ind_count_per_bin %>%
  filter(n_ind >= MIN_IND_PER_BIN) %>%
  pull(length_bin)

region_prop_filt  <- region_prop      %>% filter(length_bin %in% valid_bins)
ind_count_filt    <- ind_count_per_bin %>% filter(length_bin %in% valid_bins)
scale_factor      <- 1 / max(ind_count_filt$n_ind)

state_distribution_plot <- ggplot() +
  geom_col(
    data  = region_prop_filt,
    aes(x = length_bin, y = prop, fill = Region_assigned),
    width = 0.25, colour = "black", linewidth = 0.2
  ) +
  geom_line(
    data = ind_count_filt,
    aes(x = length_bin, y = n_ind * scale_factor, group = 1),
    colour = "black", linewidth = 2, lineend = "round"
  ) +
  scale_x_continuous(
    breaks = c(0, 5, 10, 15),
    labels = c("0", "5", "10", "15"),
    limits = c(min(valid_bins) - 0.5, max(valid_bins) + 0.5)
  ) +
  scale_y_continuous(
    name     = "Rel. Freq.",
    limits   = c(0, 1),
    labels   = scales::percent_format(),
    sec.axis = sec_axis(~ . / scale_factor, name = "Number of unique individuals")
  ) +
  scale_fill_manual(values = REGION_COLORS) +
  theme_bw(base_size = 16) +
  theme(
    axis.title.x     = element_blank(),
    axis.text.x      = element_blank(),
    axis.text.y      = element_blank(),
    axis.title.y     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "none",
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.75)
  )

ggsave("State_distribution_plot.png", state_distribution_plot,
       width = 1500, height = 1500, dpi = 300, units = "px")


# =============================================================================
# 7.  PERMANENT-TRANSITION DETECTION
# =============================================================================
# A "permanent transition" is defined as a change from one region to another
# where both the preceding and following stretches contain ≥ MIN_RUN_THRESHOLD
# consecutive observations assigned to the same region. This conservative rule
# filters out brief, transient excursions.

# ── 7a. Detection function ───────────────────────────────────────────────────
detect_permanent_transitions <- function(df, min_run = MIN_RUN_THRESHOLD) {
  df <- df %>% arrange(Individual, Distance_to_core)
  results <- list()

  for (ind in unique(df$Individual)) {
    sub     <- df %>% filter(Individual == ind)
    rle_res <- rle(sub$Region_assigned)

    # Runs satisfying the minimum-length criterion
    valid_runs  <- rle_res$lengths >= min_run
    stable_regs <- rle_res$values[valid_runs]

    if (sum(valid_runs) < 2) next  # need at least two stable states

    # Record every change between consecutive stable regions
    for (i in seq_len(length(stable_regs) - 1)) {
      if (stable_regs[i] != stable_regs[i + 1]) {
        results[[length(results) + 1]] <- data.frame(
          Individual = ind,
          from       = stable_regs[i],
          to         = stable_regs[i + 1]
        )
      }
    }
  }

  if (length(results) == 0)
    return(data.frame(Individual = character(),
                      from = character(), to = character()))
  bind_rows(results)
}

transitions <- detect_permanent_transitions(transect)

# ── 7b. Transition summary ───────────────────────────────────────────────────
# Count and proportion of each directed transition (from → to)
trans_summary <- transitions %>%
  count(from, to, name = "n_trans") %>%
  mutate(prop = n_trans / sum(n_trans))

# Number of permanent transitions per individual
trans_per_ind <- transitions %>%
  count(Individual, name = "n_transitions")

# Count of individuals by number of transitions (including 0)
n_ind_total <- n_distinct(transect$Individual)

trans_count_summary <- trans_per_ind %>%
  count(n_transitions, name = "n_ind") %>%
  bind_rows(tibble(n_transitions = 0L,
                   n_ind = n_ind_total - sum(.$n_ind))) %>%
  arrange(n_transitions)

# ── 7c. Transition diagram  →  Transition_diagram.png ───────────────────────
# Directed graph where nodes are regions and edge width encodes the proportion
# of transitions following that path.

# Node positions (fixed layout: estuary top-left, strait top-right,
#                               south   bot-left, unidentified bot-right)
all_regions <- sort(unique(c(transitions$from, transitions$to)))
nodes <- data.frame(
  region = all_regions,
  x = case_when(all_regions == "estuary"      ~ -0.8,
                all_regions == "strait"        ~  0.8,
                all_regions == "south"         ~ -0.8,
                all_regions == "unidentified"  ~  0.8),
  y = case_when(all_regions == "estuary"      ~  0.8,
                all_regions == "strait"        ~  0.8,
                all_regions == "south"         ~ -0.8,
                all_regions == "unidentified"  ~ -0.8)
)

# Edges: arrow tip shortened to 85 % of the segment so it sits near the node
edges <- trans_summary %>%
  mutate(
    x    = nodes$x[match(from, nodes$region)],
    y    = nodes$y[match(from, nodes$region)],
    xend = nodes$x[match(to,   nodes$region)],
    yend = nodes$y[match(to,   nodes$region)],
    x2   = x + (xend - x) * 0.85,
    y2   = y + (yend - y) * 0.85
  )

transition_diagram <- ggplot() +
  geom_curve(
    data = edges,
    aes(x = x, y = y, xend = x2, yend = y2,
        linewidth = prop, color = from),
    curvature = 0.25,
    arrow     = arrow(type = "closed", length = unit(0.4, "cm")),
    lineend   = "round"
  ) +
  geom_point(
    data = nodes,
    aes(x = x, y = y, fill = region),
    shape = 21, size = 20, colour = "black", stroke = 1.5
  ) +
  scale_linewidth(range = c(0.5, 4)) +
  expand_limits(x = c(-1, 1), y = c(-1, 1)) +
  scale_fill_manual(values  = REGION_COLORS) +
  scale_color_manual(values = REGION_COLORS) +
  theme_bw(base_size = 16) +
  theme(
    axis.title       = element_blank(),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "none",
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.75)
  )

ggsave("Transition_diagram.png", transition_diagram,
       width = 1500, height = 1500, dpi = 300, units = "px")

# ── 7d. Barplot: number of transitions per individual  →  Barplot_nombre_transitions.png

barplot_transitions <- ggplot(trans_count_summary,
                              aes(x = factor(n_transitions), y = n_ind)) +
  geom_col(fill = "steelblue", colour = "black", width = 0.8) +
  scale_x_discrete(name = "Number of permanent transitions") +
  scale_y_continuous(name   = "Number of individuals",
                     expand = expansion(mult = c(0, 0.1))) +
  theme_bw(base_size = 16) +
  theme(
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "none",
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.75)
  )

ggsave("Barplot_nombre_transitions.png", barplot_transitions,
       width = 1500, height = 1500, dpi = 300, units = "px")


# =============================================================================
# 8.  SENSITIVITY ANALYSIS — TRANSITION THRESHOLD
# =============================================================================
# The detection threshold (minimum consecutive observations defining a "stable"
# stretch) is varied from 3 to 12 to assess how results depend on this choice.
# The chosen threshold (MIN_RUN_THRESHOLD = 6) is highlighted in the figures.

THRESHOLDS <- 3:12

# ── 8a. Run detection for each threshold ────────────────────────────────────
sensitivity_results <- lapply(THRESHOLDS, function(thr) {

  trans <- detect_permanent_transitions(transect, min_run = thr)

  # Per-individual transition count (including 0)
  all_inds <- data.frame(Individual = unique(transect$Individual))
  if (nrow(trans) > 0) {
    per_ind <- trans %>% count(Individual, name = "n_transitions")
  } else {
    per_ind <- data.frame(Individual = character(), n_transitions = integer())
  }
  per_ind_full <- all_inds %>%
    left_join(per_ind, by = "Individual") %>%
    mutate(n_transitions = coalesce(n_transitions, 0L))

  # Mobility classes
  mob <- per_ind_full %>%
    mutate(mobility_class = case_when(
      n_transitions == 0 ~ "0",
      n_transitions == 1 ~ "1",
      n_transitions == 2 ~ "2",
      TRUE               ~ "3+"
    )) %>%
    count(mobility_class) %>%
    mutate(prop = n / n_ind_total, threshold = thr)

  list(threshold = thr, n_trans = nrow(trans), mob_class = mob)
})

# ── 8b. Assemble plotting data frames ────────────────────────────────────────
df_total_trans <- data.frame(
  threshold = THRESHOLDS,
  n_trans   = sapply(sensitivity_results, `[[`, "n_trans")
)

df_mobility <- bind_rows(lapply(sensitivity_results, `[[`, "mob_class")) %>%
  mutate(mobility_class = factor(mobility_class, levels = c("0", "1", "2", "3+")))

palette_mob <- c("0" = "#76B7B2", "1" = "#E15759",
                 "2" = "#B07AA1", "3+" = "#59A14F")

# ── 8c. Figure: total transitions vs threshold  →  Total_transitions_vs_threshold.png

p_total_sens <- ggplot(df_total_trans, aes(x = threshold, y = n_trans)) +
  geom_line(colour = "grey40", linewidth = 1) +
  geom_point(shape = 21, size = 4, fill = "white", colour = "black", stroke = 1.2) +
  scale_x_continuous(breaks = THRESHOLDS) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  theme_bw(base_size = 16) +
  theme(
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "none",
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.75)
  )

ggsave("Total_transitions_vs_threshold.png", p_total_sens,
       width = 1500, height = 1500, dpi = 300, units = "px")

# ── 8d. Figure: proportion per mobility class vs threshold  →  Proportion_vs_threshold.png

p_mob_sens <- ggplot(df_mobility,
                     aes(x = threshold, y = prop,
                         colour = mobility_class, shape = mobility_class,
                         group = mobility_class)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3.5, stroke = 1.2) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  scale_x_continuous(breaks = THRESHOLDS) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.02, 0.08))) +
  scale_colour_manual(values = palette_mob,
                      name   = "Number of transitions",
                      labels = c("0", "1", "2", "3+")) +
  theme_bw(base_size = 16) +
  theme(
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "none",
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.75)
  )

ggsave("Proportion_vs_threshold.png", p_mob_sens,
       width = 1500, height = 1500, dpi = 300, units = "px")


# =============================================================================
# 9.  FILTERING FOR CLUSTERING
# =============================================================================
# Two sequential filters are applied before the MCA / k-means:
#   (a) Remove individuals with > 40 % of their transect labelled "unidentified"
#       (too little informative signal).
#   (b) Retain only the top 70 % by fork length (removes the smallest / youngest
#       fish whose otolith transects are truncated), then hard-cut at 12.9 cm
#       to ensure a comparable transect length across individuals.

# ── 9a. Filter by unidentified proportion ────────────────────────────────────
transect_filtered <- transect %>%
  group_by(Individual) %>%
  mutate(pct_unidentified = mean(Region_assigned == "unidentified")) %>%
  ungroup() %>%
  filter(pct_unidentified <= 0.40)

# ── 9b. Histogram of unidentified proportions  →  undetermined.png ───────────
undetermined_df <- transect_filtered %>%
  group_by(Individual) %>%
  summarise(pct_unid = mean(Region_assigned == "unidentified") * 100,
            .groups  = "drop")

undetermined_plot <- ggplot(undetermined_df, aes(x = pct_unid)) +
  geom_histogram(binwidth = 5, fill = "#5D97EB", color = "black") +
  scale_x_continuous(limits = c(0, 100)) +
  theme_bw(base_size = 16) +
  theme(
    axis.title       = element_blank(),
    axis.text        = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "none",
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.75)
  )

ggsave("undetermined.png", undetermined_plot,
       width = 1500, height = 1500, dpi = 300, units = "px")

# ── 9c. Add row order within each individual ─────────────────────────────────
transect_filtered <- transect_filtered %>%
  group_by(Individual) %>%
  mutate(Order = row_number() - 1) %>%
  ungroup()

# ── 9d. Keep top 70 % by fork length ─────────────────────────────────────────
length_info   <- transect_filtered %>%
  group_by(Individual) %>%
  summarise(Length_ind = first(Length), .groups = "drop")

length_thresh <- quantile(length_info$Length_ind, 0.3)

transect_top70 <- length_info %>%
  filter(Length_ind >= length_thresh) %>%
  dplyr::select(Individual) %>%
  inner_join(transect_filtered, by = "Individual")

# ── 9e. Hard length cut-off: retain only points ≤ LENGTH_CUTOFF cm ───────────
ind_above_cutoff <- transect_top70 %>%
  filter(Length >= LENGTH_CUTOFF) %>%
  pull(Individual) %>% unique()

transect_cut <- transect_top70 %>%
  filter(Individual %in% ind_above_cutoff) %>%
  group_by(Individual, Sex) %>%
  arrange(Length_at_dist) %>%
  filter(cumall(Length_at_dist <= LENGTH_CUTOFF)) %>%
  ungroup()

cat("Individuals retained for clustering:", n_distinct(transect_cut$Individual), "\n")


# =============================================================================
# 10.  MULTIPLE CORRESPONDENCE ANALYSIS (MCA)
# =============================================================================
# Each otolith transect is converted to a wide format where columns represent
# sequential positions (Order) and cells contain the assigned region. MCA
# reduces this high-dimensional categorical matrix to a 2D space used for
# clustering.

# ── 10a. Wide format: one column per sequential bin ──────────────────────────
transect_wide <- transect_cut %>%
  dplyr::select(Individual, Order, Region_assigned) %>%
  pivot_wider(names_from = Order, values_from = Region_assigned)

# Drop trailing columns with many NAs (last 4 positions tend to be incomplete)
transect_mca_input <- transect_wide %>%
  dplyr::select(Individual, where(~ any(!is.na(.)))) %>%
  dplyr::select(-tail(names(.), 4))

# ── 10b. Fit MCA ─────────────────────────────────────────────────────────────
df_mca <- transect_mca_input %>%
  dplyr::select(-Individual) %>%
  mutate(across(everything(), as.factor))

mca_res <- MCA(df_mca, graph = FALSE, ncp = 5)

cat("MCA eigenvalues (% variance):\n")
print(round(mca_res$eig, 2))

# ── 10c. Variance barplot  →  variances.png ───────────────────────────────────
eig_df <- data.frame(
  Dim      = paste0("Dim ", 1:10),
  Variance = mca_res$eig[1:10, 2]
)

variances_plot <- ggplot(eig_df, aes(x = reorder(Dim, -Variance), y = Variance)) +
  geom_bar(stat = "identity", fill = "grey80", colour = "black") +
  # Highlight the two retained dimensions in blue
  geom_bar(data = eig_df[1:2, ],
           aes(x = Dim, y = Variance),
           stat = "identity", fill = "#5D97EB", colour = "black") +
  theme_bw(base_size = 16) +
  theme(
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "none",
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.75)
  )

ggsave("variances.png", variances_plot,
       width = 1500, height = 1500, dpi = 300, units = "px")

# ── 10d. Extract 2D individual coordinates ───────────────────────────────────
mca_coords <- as.data.frame(mca_res$ind$coord) %>%
  dplyr::select(`Dim 1`, `Dim 2`) %>%
  mutate(Individual = transect_mca_input$Individual)

# Attach survey information for visual differentiation of data sources
survey_info <- transect_cut %>%
  group_by(Individual) %>%
  summarise(Survey = first(Survey), .groups = "drop")

mca_coords <- mca_coords %>%
  left_join(survey_info, by = "Individual")


# =============================================================================
# 11.  k-MEANS CLUSTERING ON MCA COORDINATES
# =============================================================================

# ── 11a. Diagnostic plots (choose k) ─────────────────────────────────────────
fviz_nbclust(mca_coords[, 1:2], kmeans, method = "wss")       # elbow
fviz_nbclust(mca_coords[, 1:2], kmeans, method = "silhouette") # silhouette
fviz_nbclust(mca_coords[, 1:2], kmeans, method = "gap_stat")   # gap statistic

# ── 11b. Final k-means (k = N_CLUSTERS) ──────────────────────────────────────
km_fit       <- kmeans(mca_coords[, 1:2], centers = N_CLUSTERS, nstart = 25)
mca_coords$cluster        <- factor(km_fit$cluster)

# Relabel clusters for consistency with the paper (mapping table)
cluster_relabel <- c("1" = 4, "2" = 3, "3" = 2, "4" = 1)
mca_coords$cluster_mapped <- factor(cluster_relabel[as.character(km_fit$cluster)])

var_explained <- round(mca_res$eig[1:2, 2], 1)
cat("Variance explained by Dim1:", var_explained[1], "%\n")
cat("Variance explained by Dim2:", var_explained[2], "%\n")

# ── 11c. Colour scheme: solid for Commercial, light for MPO ──────────────────
pal_cluster  <- c("#76B7B2", "#59A14F", "#E15759", "#B07AA1", "#76B7B2")
pal_light    <- lighten(pal_cluster, amount = 0.5)

mca_coords$cluster_survey <- paste0(mca_coords$cluster_mapped, "_",
                                    mca_coords$Survey)
cs_levels  <- unique(mca_coords$cluster_survey)
pal_cs     <- setNames(rep(NA, length(cs_levels)), cs_levels)
for (i in seq_along(pal_cluster)) {
  pal_cs[paste0(i, "_Commercial")] <- pal_cluster[i]
  pal_cs[paste0(i, "_MPO")]        <- pal_light[i]
}

# ── 11d. Cluster scatter plot  →  Clusters.png ───────────────────────────────
clusters_plot <- ggplot(mca_coords, aes(x = `Dim 1`, y = `Dim 2`)) +
  stat_ellipse(
    aes(fill = factor(cluster_mapped), group = factor(cluster_mapped)),
    type = "norm", geom = "polygon", alpha = 0.15, colour = NA
  ) +
  scale_fill_manual(
    values = setNames(pal_cluster,
                      as.character(sort(unique(mca_coords$cluster_mapped))))
  ) +
  stat_ellipse(
    aes(colour = factor(cluster_mapped), group = factor(cluster_mapped)),
    type = "norm", linewidth = 1.2, alpha = 0.85, show.legend = FALSE
  ) +
  new_scale_fill() +
  geom_point(
    aes(fill = cluster_survey, shape = factor(cluster_mapped)),
    colour = "black", size = 3, stroke = 0.8
  ) +
  scale_fill_manual(values  = pal_cs) +
  scale_colour_manual(
    values = setNames(pal_cluster,
                      as.character(sort(unique(mca_coords$cluster_mapped))))
  ) +
  scale_shape_manual(values = c(21, 24, 22, 23, 25)) +
  theme_bw(base_size = 16) +
  theme(
    axis.title       = element_blank(),
    axis.text        = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "none",
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.75)
  )

ggsave("Clusters.png", clusters_plot,
       width = 1500, height = 1500, dpi = 300, units = "px")

# ── 11e. Merge cluster assignments back into long-format data ─────────────────
transect_clustered <- transect_top70 %>%
  left_join(mca_coords %>% dplyr::select(Individual, cluster = cluster_mapped),
            by = "Individual") %>%
  filter(!is.na(cluster))

cat("Individuals per cluster:\n")
print(transect_clustered %>%
        distinct(Individual, cluster) %>%
        count(cluster, name = "n_individuals"))

# ── 11f. Chi-squared tests (cluster ~ Survey, cluster ~ Sex) ─────────────────
chisq.test(table(transect_clustered$cluster, transect_clustered$Survey))
chisq.test(table(transect_clustered$cluster, transect_clustered$Sex))

# Multinomial regression: cluster as a function of Survey
model_survey <- multinom(cluster ~ Survey, data = transect_clustered)
z_scores     <- summary(model_survey)$coefficients /
                summary(model_survey)$standard.errors
p_values     <- (1 - pnorm(abs(z_scores), 0, 1)) * 2
print(p_values)

# ── 11g. Region-per-cluster barplot  →  Barplot_clusters.png ─────────────────
# One row per individual; combine Region_caught × Survey into a single factor
# for a two-toned colour encoding (solid = Commercial, light = MPO).

cluster_region <- transect_clustered %>%
  distinct(Individual, cluster, Survey, Region_caught) %>%
  mutate(
    Region_Survey = case_when(
      Region_caught == "estuary" & Survey == "MPO"        ~ "Estuary MPO",
      Region_caught == "estuary" & Survey == "Commercial" ~ "Estuary Commercial",
      Region_caught == "north"   & Survey == "MPO"        ~ "North MPO",
      Region_caught == "north"   & Survey == "Commercial" ~ "North Commercial",
      Region_caught == "south"   & Survey == "MPO"        ~ "South MPO",
      Region_caught == "south"   & Survey == "Commercial" ~ "South Commercial",
      Region_caught == "strait"  & Survey == "MPO"        ~ "SBI MPO",
      Region_caught == "strait"  & Survey == "Commercial" ~ "SBI Commercial"
    ),
    Region_Survey = factor(Region_Survey,
                           levels = c("SBI MPO", "SBI Commercial",
                                      "South MPO", "South Commercial",
                                      "North MPO", "North Commercial",
                                      "Estuary MPO", "Estuary Commercial"))
  )

barplot_region_data <- cluster_region %>%
  count(cluster, Region_Survey, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

cols_region_survey <- c(
  "Estuary MPO"        = "#ffbb33", "Estuary Commercial" = "#e67e00",
  "North MPO"          = "#b8e07e", "North Commercial"   = "#66a21a",
  "South MPO"          = "#7bb9f0", "South Commercial"   = "#3b78b5",
  "SBI MPO"            = "#ff4d4d", "SBI Commercial"     = "#b20000"
)

barplot_clusters_plot <- ggplot(barplot_region_data,
                                aes(x = factor(cluster), y = prop,
                                    fill = Region_Survey)) +
  geom_bar(stat = "identity", colour = "black") +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5), size = 3, colour = "black") +
  scale_fill_manual(values = cols_region_survey) +
  theme_bw(base_size = 16) +
  theme(
    axis.title       = element_blank(),
    axis.text        = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "none",
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.75)
  )

ggsave("Barplot_clusters.png", barplot_clusters_plot,
       width = 1500, height = 1500, dpi = 300, units = "px")

# ── 11h. Sex-per-cluster barplot  →  Barplot_clusters_sex.png ────────────────
cluster_sex <- transect_clustered %>%
  distinct(Individual, cluster, Sex) %>%
  mutate(Sex = case_when(Sex %in% c("F", "Female") ~ "Female",
                         Sex %in% c("M", "Male")   ~ "Male",
                         TRUE                       ~ NA_character_),
         Sex = factor(Sex, levels = c("Female", "Male")))

barplot_sex_data <- cluster_sex %>%
  count(cluster, Sex, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

barplot_clusters_sex_plot <- ggplot(barplot_sex_data,
                                    aes(x = factor(cluster), y = prop, fill = Sex)) +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = c("Female" = "#ff6b6b", "Male" = "#4dabf7")) +
  theme_bw(base_size = 16) +
  theme(
    axis.title       = element_blank(),
    axis.text        = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "none",
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.75)
  )

ggsave("Barplot_clusters_sex.png", barplot_clusters_sex_plot,
       width = 1500, height = 1500, dpi = 300, units = "px")


# =============================================================================
# 12.  SYNTHETIC CLUSTER PROFILES — BUBBLE PLOT
# =============================================================================
# For each cluster and each valid fork-length bin, the number of individuals
# assigned to each region is shown as a bubble. The majority region is
# coloured; minority regions are shown in grey. Bubble size is scaled globally
# across clusters so sizes are directly comparable.

# ── 12a. Valid fork-length bins (≥ MIN_IND_PER_BIN individuals in any cluster)
valid_lengths <- transect_clustered %>%
  filter(!is.na(cluster)) %>%
  group_by(cluster, Length_at_dist) %>%
  summarise(n_ind_cl = n_distinct(Individual), .groups = "drop") %>%
  filter(n_ind_cl >= MIN_IND_PER_BIN) %>%
  pull(Length_at_dist) %>%
  unique()

# ── 12b. Count per cluster × length × region ─────────────────────────────────
bubble_data <- transect_clustered %>%
  filter(!is.na(cluster),
         Length_at_dist %in% valid_lengths) %>%
  mutate(Region_assigned = ifelse(Region_assigned == "north",
                                  "unidentified", Region_assigned)) %>%
  group_by(cluster, Length_at_dist, Region_assigned) %>%
  summarise(n_region = n(), .groups = "drop") %>%
  group_by(cluster, Length_at_dist) %>%
  mutate(
    is_majority  = n_region == max(n_region),
    color_to_use = ifelse(is_majority, as.character(Region_assigned), "other")
  ) %>%
  ungroup() %>%
  mutate(Region_assigned = factor(Region_assigned,
                                  levels = rev(c("estuary", "south",
                                                  "strait", "unidentified"))))

global_max_n <- max(bubble_data$n_region, na.rm = TRUE)

cols_bubble <- c(REGION_COLORS, "other" = "grey80")

# ── 12c. One PNG per cluster ──────────────────────────────────────────────────
bubble_dir <- "plots_cluster_majority_bis"
dir.create(bubble_dir, showWarnings = FALSE)

unique(bubble_data$cluster) %>%
  walk(function(clust) {
    p_bub <- ggplot(bubble_data %>% filter(cluster == clust),
                    aes(x = Length_at_dist, y = Region_assigned,
                        size = n_region, colour = color_to_use)) +
      geom_point(alpha = 0.85) +
      scale_colour_manual(values = cols_bubble) +
      scale_y_discrete(drop = FALSE) +
      scale_size_area(max_size = 10, limits = c(0, global_max_n)) +
      xlim(0, 17.1) +
      theme_bw(base_size = 16) +
      theme(
        axis.text        = element_blank(),
        axis.title       = element_blank(),
        panel.grid       = element_blank(),
        legend.position  = "none",
        panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.75)
      )

    ggsave(file.path(bubble_dir, paste0("cluster_", clust, "_bubbles.png")),
           plot = p_bub, width = 1500, height = 800, dpi = 300, units = "px")
  })


# =============================================================================
# 13.  MIGRATION SHIFT ANALYSIS (CLUSTER 3)
# =============================================================================
# Cluster 3 contains individuals that spent time in both the Strait of Belle
# Isle (SBI / "strait") and the southern Gulf ("south"). The proportion of
# individuals in each zone is tracked along the length axis to identify the
# fork length at which the stock shifts between the two areas.
# These proportions are then contextualised against observed fork-length
# distributions from northern and southern GSL trawl surveys.

# ── 13a. Wilson confidence interval helpers ───────────────────────────────────
wilson_lo <- function(x, n, z = 1.96) {
  p <- x / n
  (p + z^2/(2*n) - z * sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n) * 100
}
wilson_hi <- function(x, n, z = 1.96) {
  p <- x / n
  (p + z^2/(2*n) + z * sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n) * 100
}

# ── 13b. Subset cluster 3 & compute proportions ──────────────────────────────
transect_cl3 <- transect_clustered %>% filter(cluster == "3")

# Keep only individuals that visited at least one of the two focal zones
ind_in_cl3_zones <- transect_cl3 %>%
  filter(Region_assigned %in% c("strait", "south")) %>%
  pull(Individual) %>% unique()

prop_cl3 <- transect_cl3 %>%
  filter(Individual %in% ind_in_cl3_zones) %>%
  group_by(Length_at_dist) %>%
  summarise(
    n_total = n_distinct(Individual),
    n_sbi   = sum(Region_assigned == "strait"),
    n_south = sum(Region_assigned == "south"),
    .groups = "drop"
  ) %>%
  filter(n_total >= MIN_IND_PER_BIN) %>%
  mutate(
    prop_sbi    = n_sbi   / n_total * 100,
    prop_south  = n_south / n_total * 100,
    ci_sbi_lo   = wilson_lo(n_sbi,   n_total),
    ci_sbi_hi   = wilson_hi(n_sbi,   n_total),
    ci_south_lo = wilson_lo(n_south, n_total),
    ci_south_hi = wilson_hi(n_south, n_total)
  )

prop_cl3_long <- prop_cl3 %>%
  pivot_longer(cols = c(prop_sbi, prop_south),
               names_to = "Zone", values_to = "Proportion") %>%
  mutate(
    ci_lo = ifelse(Zone == "prop_sbi", ci_sbi_lo, ci_south_lo),
    ci_hi = ifelse(Zone == "prop_sbi", ci_sbi_hi, ci_south_hi),
    Zone  = case_when(Zone == "prop_sbi"   ~ "strait",
                      Zone == "prop_south" ~ "south")
  ) %>%
  dplyr::select(Length_at_dist, Zone, Proportion, ci_lo, ci_hi)

# ── 13c. Figure: SBI vs South proportions  →  shift_SBI_sud.png ──────────────
shift_plot <- ggplot(prop_cl3_long,
                     aes(x = Length_at_dist, y = Proportion,
                         colour = Zone, fill = Zone)) +
  geom_vline(data = SHIFT_VALUES, aes(xintercept = cross_x),
             linetype = "dashed", colour = "grey30", linewidth = 0.9) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = c("strait" = "#FE0000", "south" = "#5D97EB")) +
  scale_fill_manual(values   = c("strait" = "#FE0000", "south" = "#5D97EB")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     labels = paste0(seq(0, 100, 25), " %")) +
  scale_x_continuous(limits = c(0, 16.5)) +
  theme_bw(base_size = 16) +
  theme(
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "none",
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.75)
  )

ggsave("shift_SBI_sud.png", shift_plot,
       width = 1500, height = 1500, dpi = 300, units = "px")

# ── 13d. Prepare survey fork-length data ──────────────────────────────────────
# n-GSL: convert sexagesimal coordinates to decimal degrees; convert mm → cm;
#        restrict to the known capelin overwintering area.
n_gsl <- n_gsl_raw %>%
  filter(year > 2021, year < 2025) %>%
  mutate(
    lat_dd  =  floor(latitude  / 100) + (latitude  %% 100) / 60,
    lon_dd  = -(floor(longitude / 100) + (longitude %% 100) / 60),
    length_cm = length / 10
  ) %>%
  filter(lat_dd >= 49.67, lat_dd <= 51.92,
         lon_dd >= -58.85, lon_dd <= -55.65)

# s-GSL: already in cm; restrict to same time window
s_gsl <- s_gsl_raw %>% filter(year > 2021, year < 2025)

# ── 13e. Weighted summary statistics ─────────────────────────────────────────
weighted_summary <- function(df, length_col, weight_col) {
  x  <- df[[length_col]]; w <- df[[weight_col]]
  ok <- !is.na(x) & !is.na(w); x <- x[ok]; w <- w[ok]
  mu    <- weighted.mean(x, w)
  n_eff <- sum(w)
  sd_w  <- sqrt(sum(w * (x - mu)^2) / n_eff)
  se    <- sd_w / sqrt(n_eff)
  data.frame(mean_length = mu, sd_w = sd_w, n_eff = n_eff,
             se = se, ci_lower = mu - 1.96*se, ci_upper = mu + 1.96*se)
}

summary_n_gsl <- weighted_summary(n_gsl, "length_cm", "number.caught")
summary_s_gsl <- weighted_summary(s_gsl, "length",    "number.caught")
cat("n-GSL weighted mean fork length:", round(summary_n_gsl$mean_length, 2), "cm\n")
cat("s-GSL weighted mean fork length:", round(summary_s_gsl$mean_length, 2), "cm\n")

# ── 13f. Weighted kernel density estimates (same x-range as proportion plot) ──
x_range <- range(prop_cl3$Length_at_dist, na.rm = TRUE)

dens_n <- density(n_gsl$length_cm,
                  weights = n_gsl$number.caught / sum(n_gsl$number.caught, na.rm = TRUE),
                  from = x_range[1], to = x_range[2], n = 512)

dens_s <- density(s_gsl$length,
                  weights = s_gsl$number.caught / sum(s_gsl$number.caught, na.rm = TRUE),
                  from = x_range[1], to = x_range[2], n = 512)

dens_df <- bind_rows(
  data.frame(x = dens_n$x, y = dens_n$y * 100, Zone = "N-GSL"),
  data.frame(x = dens_s$x, y = dens_s$y * 100, Zone = "S-GSL")
)

# ── 13g. Figure: fork-length distributions  →  density_NGSL_SGSL.png ─────────
density_plot <- ggplot(dens_df, aes(x = x, y = y, colour = Zone, fill = Zone)) +
  geom_area(alpha = 0.15, position = "identity") +
  geom_line(linewidth = 1.2) +
  geom_vline(data = SHIFT_VALUES, aes(xintercept = cross_x),
             linetype = "dashed", colour = "grey30", linewidth = 0.9) +
  annotate("text", x = summary_n_gsl$mean_length, y = -Inf,
           label    = paste0(round(summary_n_gsl$mean_length, 1), " cm"),
           colour   = "#FE0000", vjust = 2.2, hjust = 0.5,
           size     = 4.5, fontface = "bold") +
  annotate("text", x = summary_s_gsl$mean_length, y = -Inf,
           label    = paste0(round(summary_s_gsl$mean_length, 1), " cm"),
           colour   = "#5D97EB", vjust = 2.2, hjust = 0.5,
           size     = 4.5, fontface = "bold") +
  scale_colour_manual(values = c("N-GSL" = "#FE0000", "S-GSL" = "#5D97EB")) +
  scale_fill_manual(values   = c("N-GSL" = "#FE0000", "S-GSL" = "#5D97EB")) +
  scale_x_continuous(limits  = c(0, 16.5)) +
  theme_bw(base_size = 16) +
  theme(
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "none",
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.75)
  )

ggsave("density_NGSL_SGSL.png", density_plot,
       width = 1500, height = 1000, dpi = 300, units = "px")

# =============================================================================
# END OF SCRIPT
# =============================================================================