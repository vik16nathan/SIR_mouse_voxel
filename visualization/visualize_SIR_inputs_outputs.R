library("pacman")
library("RMINC")
pacman::p_load(dplyr, readr, readxl, stringr, ggrepel, tidyr, ggplot2, ggalluvial, tidyr, patchwork)
setwd(".")
source("../ggslicer/plotting_functions/plotting_functions.R")
source("../ggslicer/plotting_functions/plotting_functions_labels.R")

allen_input_dir <- "../../preprocessed/templates/ccfv3_converted_RAS_MICe/"

allen_template_path_200um <- paste0(allen_input_dir, "average_template_200.mnc")
allen_mask_path_200um    <- paste0(allen_input_dir, "mask_200um.mnc")
allen_200um_template_df  <- prepare_masked_anatomy(allen_template_path_200um, allen_mask_path_200um, "y", seq(-7.5, 5, 1))[[2]]
allen_200um_template_df  <- allen_200um_template_df %>% filter(mask_value == 1)

# Helper: base anatomy layer shared by every row
base_anatomy <- function() {
  list(
    geom_raster(data = allen_200um_template_df,
                mapping = aes(x = x, y = z, fill = intensity),
                interpolate = TRUE),
    scale_fill_gradient(low = "black", high = "white",
                        oob = scales::squish, guide = "none"),
    ggnewscale::new_scale_fill()
  )
}

# Common theme / scale settings
common_scales <- function() {
  list(
    scale_x_continuous(expand = c(0, 0)),
    scale_y_continuous(expand = c(0, 0)),
    facet_wrap(~ slice_world, ncol = 13,
               labeller = labeller(slice_world = function(x) paste0("Slice: ", x))),
    coord_fixed(ratio = 1),
    theme_void(base_size = 32),
    theme(
      panel.spacing    = unit(0, "npc"),
      strip.text       = element_text(size = 20, face = "bold"),
      plot.title       = element_text(size = 28, face = "bold", hjust = 0.5),
      legend.title     = element_text(size = 22),
      legend.text      = element_text(size = 18),
      legend.position  = "right"
    )
  )
}

plots <- list()

# ── 1 & 2: Connectivity strength (cp, dg) ──────────────────────────────────
vol_conn <- "../../derivatives/SIR_inputs_200um/"

for (epi in c("cp", "dg")) {
  df <- prepare_masked_anatomy(
    paste0(vol_conn, epi, "_epi_knox_downsample_200um_RAS_MICe.mnc"),
    allen_mask_path_200um, "y", seq(-7.5, 5, 1))[[2]] %>%
    filter(intensity > 0)

  eps <- 1e-8
  p <- ggplot(allen_200um_template_df, aes(x = x, y = z)) +
    base_anatomy() +
    geom_raster(data = df, aes(fill = log(as.numeric(intensity) + eps))) +
    scale_fill_gradient(name = "Log conn. strength", low = "yellow", high = "red",
                        guide = "colourbar") +
    common_scales() +
    labs(title = paste0("Connectivity Strength → ", toupper(epi)))

  plots[[paste0("conn_", epi)]] <- p
}

# ── 3: Snca expression ─────────────────────────────────────────────────────
vol_gene <- "../../derivatives/gene_expression_kNN/RAS_MICe/"

df_snca <- prepare_masked_anatomy(
  paste0(vol_gene, "Snca_coronal.mnc"),
  allen_mask_path_200um, "y", seq(-7.5, 5, 1))[[2]] %>%
  filter(intensity > 0)

p_snca <- ggplot(allen_200um_template_df, aes(x = x, y = z)) +
  base_anatomy() +
  geom_raster(data = df_snca, aes(fill = scale(as.numeric(intensity))), alpha = 0.8) +
  scale_fill_gradient2(name = "Snca expression", low = "blue", mid = "white", high = "red", midpoint=0,
                      guide = "colourbar") +
  common_scales() +
  labs(title = "Snca Expression, 200 µm")

plots[["snca"]] <- p_snca

# ── 4 & 5: Atrophy t-maps (CP, DG) ────────────────────────────────────────
vol_atr <- "../../derivatives/voxel_atrophy_maps/"

pff_file <- list(
  CP = "M83MsPff_vs_M83PBS_t_stats_linear_allen200.mnc",
  DG = "M83HuPff_vs_M83PBS_t_stats_linear_allen200.mnc"
)

for (epi in c("CP", "DG")) {
  df <- prepare_masked_anatomy(
    paste0(vol_atr, epi, "/", pff_file[[epi]]),
    allen_mask_path_200um, "y", seq(-7.5, 5, 1))[[2]] %>%
    filter(mask_value == 1)

  p <- ggplot(allen_200um_template_df, aes(x = x, y = z)) +
    base_anatomy() +
    geom_tile(data = df, aes(fill = intensity), alpha = 0.8) +
    scale_fill_gradient2(name = "t value",
                         low = "cyan", mid = "white", high = "orange",
                         guide = "colourbar") +
    common_scales() +
    labs(title = paste0("Rel. Jac. t-values: PFF vs. PBS — ", epi))

  plots[[paste0("atr_", epi)]] <- p
}

# ── Stitch into one aligned figure ─────────────────────────────────────────
# patchwork stacks plots in the order they appear; plot_layout(ncol=1) makes
# each plot a full-width row, preserving the 13-column facet alignment.
combined <- (
  plots[["conn_cp"]] /
  plots[["conn_dg"]] /
  plots[["snca"]]    /
  plots[["atr_CP"]]  /
  plots[["atr_DG"]]
) +
  plot_layout(ncol = 1, guides = "keep") +
  plot_annotation(
    title   = "Brain-wide spatial maps — connectivity, gene expression & atrophy",
    theme   = theme(
      plot.title = element_text(size = 32, face = "bold", hjust = 0.5,
                                margin = margin(b = 12))
    )
  )

ggsave("../figures/combined_all_rows.png", combined,
       width = 30, height = 20,   # 4 units per row × 5 rows
       dpi   = 300)

message("Saved → figures/combined_all_rows.png")


####################SEPARATE SUBPLOTS##########################
snca_slices <- c(-1.5, 0)
allen_200um_template_df  <- prepare_masked_anatomy(allen_template_path_200um, allen_mask_path_200um, "y", snca_slices)[[2]]
allen_200um_template_df  <- allen_200um_template_df %>% filter(mask_value == 1)
df_snca <- prepare_masked_anatomy(
  paste0(vol_gene, "Snca_coronal.mnc"),
  allen_mask_path_200um, "y", snca_slices)[[2]] %>%
  filter(intensity > 0)

p_snca <- ggplot(allen_200um_template_df, aes(x = x, y = z)) +
  base_anatomy() +
  geom_raster(data = df_snca, aes(fill = scale(as.numeric(intensity))), alpha = 0.8) +
  scale_fill_gradient2(name = "Snca expression", low = "blue", mid = "white", high = "red", midpoint=0,
                      guide = "colourbar") +
  common_scales() +
  labs(title = "Snca Expression, 200 µm")

ggsave("../figures/Snca_plot.png", p_snca,
       width = 30, height = 20,   # 4 units per row × 5 rows
       dpi   = 300)



######Visualize striatum and hippocampus masks (200 um)########
slice_series_for_str <- seq(-2.5,   2.5,   0.5) 
str_mask_path <- "../../derivatives/downsampled_connectome/str_target_indices_filt_overlap_source_200um_RAS_MICe.mnc"
df_str <- prepare_masked_anatomy(
  str_mask_path,
  allen_mask_path_200um, "y", slice_series_for_str)[[2]] %>%
  filter(intensity > 0)

allen_200um_template_df  <- prepare_masked_anatomy(allen_template_path_200um, allen_mask_path_200um, "y", slice_series_for_str)[[2]]
allen_200um_template_df  <- allen_200um_template_df %>% filter(mask_value == 1)

p_str <- ggplot(allen_200um_template_df, aes(x = x, y = z)) +
  base_anatomy() +
  geom_raster(data = df_str, aes(fill = scale(as.numeric(intensity))), alpha = 0.8) +
  scale_fill_gradient(name = "Voxel label", low = "yellow", high = "red",
                      guide = "colourbar") +
  common_scales() +
  labs(title = "STR mask, 200 µm")

ggsave("../figures/str_mask_plot.png", p_str,
      width = 30, height = 20,   # 4 units per row × 5 rows
      dpi   = 300)

###HIP###
slice_series_for_hip <- seq(-4.5,   0.5,   0.5) 
hip_mask_path <- "../../derivatives/downsampled_connectome/hip_target_indices_filt_overlap_source_200um_RAS_MICe.mnc"
df_hip <- prepare_masked_anatomy(
  hip_mask_path,
  allen_mask_path_200um, "y", slice_series_for_hip)[[2]] %>%
  filter(intensity > 0)

allen_200um_template_df  <- prepare_masked_anatomy(allen_template_path_200um, allen_mask_path_200um, "y", slice_series_for_hip)[[2]]
allen_200um_template_df  <- allen_200um_template_df %>% filter(mask_value == 1)

p_hip <- ggplot(allen_200um_template_df, aes(x = x, y = z)) +
  base_anatomy() +
  geom_raster(data = df_hip, aes(fill = scale(as.numeric(intensity))), alpha = 0.8) +
  scale_fill_gradient(name = "Voxel label", low = "yellow", high = "red",
                      guide = "colourbar") +
  common_scales() +
  labs(title = "HIP mask, 200 µm")

ggsave("../figures/hip_mask_plot.png", p_hip,
      width = 30, height = 20,   # 4 units per row × 5 rows
      dpi   = 300)