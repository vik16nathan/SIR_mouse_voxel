library("pacman")
library("reticulate")
library("RMINC")
pacman::p_load(dplyr, readr, readxl, stringr, ggrepel, tidyr, ggplot2,
               ggnewscale, scales, purrr, patchwork)

setwd(".")
source("../ggslicer/plotting_functions/plotting_functions.R")
source("../ggslicer/plotting_functions/plotting_functions_labels.R")

pd            <- import("pandas")
cp_epi_index  <- 111045
dg_epi_index  <- 110656

eps <- 1e-10

# ── Region-specific slice sequences ───────────────────────────────────────────
# Define once here so brain plots, comparison plots, and movies are all consistent.
STR_SLICES <- seq(3,   8,   0.5)   # striatum / whole-brain
HIP_SLICES <- seq(6.5, 11.5, 0.5)  # hippocampus

### mask path###
mask_path <- "../../preprocessed/templates/mask_200um_PIR_CCFv3.mnc"

# ── Simulation file paths ──────────────────────────────────────────────────────

# ---- CP epicentre (full bilateral) -------------------------------------------
sim_file_str_retro <- paste0(
  "../simulations/retro_abm_spread_v.0.0017003872488696233.spread_rate.0.9971421696522822.dt.0.1.seed.111045.injection_amount.57.90868321082842.clearance_gene.None.k1.0.9345289210280755.k2.0.5583323404362099..csv")

sim_file_str_antero <- paste0(
  "../simulations/abm_spread_v.0.7130185134308772.spread_rate.0.8800327881355907.dt.0.1.seed.111045.injection_amount.56.16471050842456.clearance_gene.None.k1.0.8091228425272163.k2.0.9871793893663503..csv")

# ---- CP epicentre (ipsilateral only) -----------------------------------------
sim_file_str_retro_ipsi  <- "../simulations/retro_abm_spread_v.0.37707557027361915.spread_rate.0.7819432344877069.dt.0.1.seed.111045.injection_amount.55.92661754954737.clearance_gene.None.k1.0.09022400830018472.k2.0.6313247056289887._ipsi.csv"
sim_file_str_antero_ipsi <- NULL   # <-- INSERT CP ANTERO IPSI PATH HERE

# ---- CP epicentre (ipsilateral + clearance gene) -----------------------------
sim_file_cp_retro_cl_ipsi <- "../simulations/retro_abm_spread_v.0.9203497196853531.spread_rate.0.817754260835934.dt.0.1.seed.111045.injection_amount.93.28826921803088.clearance_gene.Sh3bgr.k1.0.10895190647204175.k2.0.799803880680589._ipsi.csv"

# ---- DG epicentre (full bilateral) -------------------------------------------
sim_file_hip_retro <- paste0(
  "../simulations/retro_abm_spread_v.0.011944797633727307.spread_rate.0.9832548416531878.dt.0.1.seed.110656.injection_amount.1.1848216069603854.clearance_gene.None.k1.0.08188919228370305.k2.0.8618764079000241..csv")

sim_file_hip_antero <- paste0(
  "../simulations/abm_spread_v.0.97071891251077.spread_rate.0.6670313303245704.dt.0.1.seed.110656.injection_amount.84.65707027363796.clearance_gene.None.k1.0.922920498771323.k2.0.42191190675374046..csv")

# ---- DG epicentre (ipsilateral only) -----------------------------------------
sim_file_hip_antero_ipsi <- paste0(
  "../simulations/abm_spread_v.0.008746234591760752.spread_rate.0.9666841451733642.dt.0.1.seed.110656.injection_amount.99.88133696594801.clearance_gene.None.k1.0.6870137681953692.k2.0.012801212760166067._ipsi.csv")

sim_file_hip_retro_ipsi <- NULL   # <-- INSERT DG RETRO IPSI PATH HERE

# ---- DG epicentre (ipsilateral + clearance gene) -----------------------------
sim_file_hip_antero_cl_ipsi <- "../simulations/abm_spread_v.0.8968132793840868.spread_rate.0.33268056212600106.dt.0.1.seed.110656.injection_amount.71.10174286797589.clearance_gene.Efna3.k1.0.9887085449902345.k2.0.0048143225873943105._ipsi.csv"

# ---- Params pkl paths --------------------------------------------------------
pkl_str_retro_params  <- "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/params/source_target_indices_filt_overlap_str_top30_retro.pkl"
pkl_str_antero_params <- "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/params/source_target_indices_filt_overlap_str_top30.pkl"

pkl_hip_retro_params  <- "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/params/source_target_indices_filt_overlap_hip_top30_retro.pkl"
pkl_hip_antero_params <- "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/params/source_target_indices_filt_overlap_hip_top30.pkl"

# ---- Atrophy pkl paths -------------------------------------------------------
pkl_str_atrophy <- "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/atrophy_pkl/str/M83PBS_vs_M83MsPff_t_stats_linear_allen200_PIR_CCFv3_target_filt.pkl"
pkl_hip_atrophy <- "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/atrophy_pkl/hip/M83PBS_vs_M83HuPff_hipp_t_stats_linear_allen200_PIR_CCFv3_target_filt.pkl"

# ── Shared ggplot themes ───────────────────────────────────────────────────────
theme_brain <- function(base_size = 28) {
  theme_void(base_size = base_size) +
    theme(
      panel.spacing   = unit(0, "npc"),
     
      strip.text = element_text(size = 20, face = "bold", colour = "grey20"),
      plot.title      = element_text(size = 34, face = "bold",   hjust = 0.5,
                                     colour = "black"),
      plot.subtitle   = element_text(size = 24, face = "italic", hjust = 0.5,
                                     colour = "grey20"),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = "white"),
      legend.title    = element_text(size = 17, colour = "black"),
      legend.text     = element_text(size = 15, colour = "black"),
      plot.margin = margin(0, 0, 0, 0, "pt"),
      legend.position = "right"
    )
}


theme_corr <- function() {
  theme_minimal(base_size = 20) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90"),
      plot.title       = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle    = element_text(hjust = 0.5, colour = "grey40"),
      axis.title       = element_text(face = "bold", size = 13),
      axis.text        = element_text(colour = "grey30", size = 12),
      legend.position  = "bottom",
      legend.title     = element_text(size = 14),
      legend.text      = element_text(size = 13)
    )
}

# ── Core transform: log then  ─────────────────────────────────────
# Applied to simulated intensities so colour always spans the full palette
# regardless of absolute magnitude (fixes the "all blue" problem when values
# are small and negative in log space).
log_transform <- function(x) {
  lx <- log(x + eps)
  lx
}

# ── Helper: load sim .mnc, apply log, optionally filter ipsi ─────
prepare_sim_df <- function(sim_mnc_path, slice_axis, slice_seq, ipsi_only,
                            template_z = NULL) {
  df <- prepare_masked_anatomy(sim_mnc_path, mask_path, slice_axis, slice_seq)[[2]] %>%
    filter(intensity > 0, mask_value == 1) %>%
    mutate(intensity = log_transform(intensity))

  if (ipsi_only && !is.null(template_z)) {
    z_mid <- (min(template_z) + max(template_z)) / 2
    df    <- df %>% filter(z > z_mid)
  }
  df
}

# ── Helper: overlay a sim volume on the template ──────────────────────────────
# Simulated values are log-transformed so the
# full colour palette is always used (avoids the all-blue problem).
# scale_limits: pass c(lo, hi) in [−1, 1] to keep movie frames consistent;
#               NULL = auto from data (peak plot).
brain_raster_plot <- function(sim_mnc_path,
                              title_str,
                              subtitle_str  = NULL,
                              legend_name   = "Simulated\natrophy\n(log)",
                              palette       = c("navy","dodgerblue","#FFFF00","orange","#CC0000"),
                              scale_limits  = c(-1, 1),
                              voxel_size    = 0.2,
                              slice_axis    = "x",
                              slice_seq     = STR_SLICES,
                              ipsi_only     = FALSE) {

tmpl_df <- prepare_masked_anatomy(
    "../../preprocessed/templates/average_template_200_pir.mnc", mask_path,
    slice_axis, slice_seq)[[2]] %>%
    filter(mask_value == 1)


sim_df <- prepare_sim_df(sim_mnc_path, slice_axis, slice_seq,
                          ipsi_only, template_z = tmpl_df$z)

  if (ipsi_only) {
    z_mid   <- (min(tmpl_df$z) + max(tmpl_df$z)) / 2
    tmpl_df <- tmpl_df %>% filter(z > z_mid)
  }

  ggplot(data = tmpl_df, aes(x = z, y = y)) +
    geom_raster(aes(fill = intensity), interpolate = TRUE) +
    scale_fill_gradient(low = "black", high = "white",
                        oob = squish, guide = "none") +
    ggnewscale::new_scale_fill() +
    geom_tile(data = sim_df, aes(fill = intensity),
              width = voxel_size, height = voxel_size) +
    scale_fill_gradientn(
      name    = legend_name,
      colours = palette,
      limits  = scale_limits,
      oob     = squish,
      guide   = guide_colourbar(barwidth = 0.45, barheight = 4.5,
                                title.position = "top", title.hjust=0.5)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), trans = "reverse") +
    facet_wrap(~ slice_world, ncol = 13,
               labeller = labeller(
                 slice_world = function(x) paste0(slice_axis, " = ", x))) +
    coord_fixed(ratio = 1) +
    labs(title = title_str, subtitle = subtitle_str) +
    theme_brain()
}

# ── Helper: write a numeric vector as a .mnc volume ───────────────────────────
write_vol <- function(vals, target_voxels, label_vol, mnc_like, out_path) {
  out_vol            <- numeric(length(label_vol))
  indices_to_replace <- lapply(target_voxels, function(v) which(label_vol == v))
  for (i in seq_along(target_voxels))
    out_vol[indices_to_replace[[i]]] <- vals[i]
  mincWriteVolume(out_vol, like = mnc_like, out_path)
  invisible(out_path)
}

# ── Helper: write empirical atrophy () as .mnc ───────────────────
make_empirical_vol <- function(atrophy_vec, target_voxels, label_vol, mnc_like, out_path) {
  atr_norm <- atrophy_vec
  write_vol(atr_norm, target_voxels, label_vol, mnc_like, out_path)
}

# ── Two-row comparison: log sim (top), empirical (bottom) ────────
# Each row uses its own colour limits; both share the same palette.
# Colorbars are small so both fit on the right margin.
peak_comparison_plot <- function(sim_mnc_path,
                                 emp_mnc_path,
                                 title_str,
                                 subtitle_str = NULL,
                                 palette      = c("navy","dodgerblue","#FFFF00","orange","#CC0000"),
                                 voxel_size   = 0.2,
                                 slice_axis   = "x",
                                 slice_seq    = STR_SLICES,
                                 ipsi_only    = FALSE) {

  tmpl_df <- prepare_masked_anatomy(
    "../../preprocessed/templates/average_template_200_pir.mnc", mask_path,
    slice_axis, slice_seq)[[2]] %>%
    filter(mask_value == 1)

  tmpl_df <- tmpl_df %>% filter(mask_value == 1)
  # Simulated: 
  sim_df <- prepare_sim_df(sim_mnc_path, slice_axis, slice_seq,
                           ipsi_only, template_z = tmpl_df$z)

  # Empirical: already  written to disk
  emp_df <- prepare_masked_anatomy(emp_mnc_path, mask_path, slice_axis, slice_seq)[[2]] %>%
    filter(intensity > 0, mask_value == 1)

  if (ipsi_only) {
    z_mid   <- (min(tmpl_df$z) + max(tmpl_df$z)) / 2
    tmpl_df <- tmpl_df %>% filter(z > z_mid)
    emp_df  <- emp_df  %>% filter(z > z_mid)
  }

  # Fixed limits: sim is normalised to [−1, 1]; empirical to [−1, 1]
  #sim_limits <- c(-1, 1)

  sim_limits <- NULL
  emp_limits <- range(emp_df$intensity, na.rm = TRUE)

  make_row <- function(overlay_df, tmpl, lims, row_title, fmt_fn) {
    ggplot(tmpl, aes(x = z, y = y)) +
      geom_raster(aes(fill = intensity), interpolate = TRUE) +
      scale_fill_gradient(low = "black", high = "white",
                          oob = squish, guide = "none") +
      ggnewscale::new_scale_fill() +
      geom_tile(data = overlay_df, aes(fill = intensity),
                width = voxel_size, height = voxel_size) +
      scale_fill_gradientn(
        name    = row_title,
        colours = palette,
        limits  = lims,
        oob     = squish,
        labels  = fmt_fn,
        guide   = guide_colourbar(
          barwidth       = 0.45, barheight = 4.5,
          title.position = "top", title.hjust = 0.5,
          title.theme    = element_text(size = 24, colour = "black",
                                        face = "bold", lineheight = 1.1),
          label.theme    = element_text(size = 22, colour = "black")
        )
      ) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0), trans = "reverse") +
      facet_wrap(~ slice_world, ncol = 13,
                 labeller = labeller(
                   slice_world = function(x) paste0(slice_axis, " = ", x))) +
      coord_fixed(ratio = 1) +
      theme_brain() +
      theme(
        strip.text      = element_text(size = 18, colour = "black"),
        legend.position = "right",
        plot.title      = element_text(size = 26, colour = "black",
                                       face = "bold", hjust = 0.5),
        plot.background = element_rect(fill = "white", colour = "white"),
        plot.margin = margin(0, 0, 0, 0, "pt")
      ) +
      labs(title = row_title)
  }

  p_sim <- make_row(sim_df,  tmpl_df, sim_limits,
                    "Simulated atrophy (log)",
                    function(x) sprintf("%.1f", x))
  p_emp <- make_row(emp_df, tmpl_df, emp_limits,
                    "Empirical atrophy (t)",
                    function(x) sprintf("%.2f", x))

  (p_sim / p_emp) +
    plot_annotation(
      title    = title_str,
      subtitle = subtitle_str,
      theme    = theme(
        plot.title      = element_text(size = 24, face = "bold",
                                       colour = "black", hjust = 0.5),
        plot.subtitle   = element_text(size = 18, colour = "black", hjust = 0.5),
        plot.background = element_rect(fill = "white", colour = "white"),
        plot.margin = margin(0, 0, 0, 0, "pt")
      )
    )
}

# ── Helper: Spearman correlation at every timestep ────────────────────────────
find_peak_timestep <- function(sim_matrix, empirical_vec) {
  cors <- vapply(seq_len(ncol(sim_matrix)), function(j) {
    cor(as.numeric(sim_matrix[, j]),
        as.numeric(empirical_vec),
        method = "spearman")
  }, numeric(1))
  list(peak_t   = which.max(cors),
       rho      = max(cors, na.rm = TRUE),
       all_cors = cors)
}

spearman <- function(x, y) cor(as.numeric(x), as.numeric(y), method = "spearman")

HEMI_COLS <- c("Right" = "#E8604C", "Left" = "#4CA8E8")

# ── Scatterplot helper ─────────────────────────────────────────────────────────
make_scatter <- function(x_vec, y_vec, colour, x_lab, title_str,
                         log = FALSE, gt0 = FALSE, hemi_mask = NULL) {
  if (gt0) {
    keep  <- x_vec > 0
    y_vec <- y_vec[keep]; x_vec <- x_vec[keep]
    if (!is.null(hemi_mask)) hemi_mask <- hemi_mask[keep]
  }
  if (log) x_vec <- log(x_vec + eps)

  rho_all <- spearman(x_vec, y_vec)

  if (!is.null(hemi_mask)) {
    df    <- data.frame(x = x_vec, y = y_vec,
                        hemi = factor(ifelse(hemi_mask, "Right", "Left"),
                                      levels = c("Right","Left")))
    rho_R <- spearman(x_vec[hemi_mask],  y_vec[hemi_mask])
    rho_L <- spearman(x_vec[!hemi_mask], y_vec[!hemi_mask])
    annot <- sprintf("ρ (R) = %.3f\nρ (L) = %.3f\nρ (all) = %.3f", rho_R, rho_L, rho_all)
    p <- ggplot(df, aes(x, y, colour = hemi)) +
      geom_point(alpha = 0.4, size = 1.2) +
      scale_colour_manual(name = "Hemisphere", values = HEMI_COLS,
                          guide = guide_legend(override.aes = list(alpha = 1, size = 3)))
  } else {
    df    <- data.frame(x = x_vec, y = y_vec)
    annot <- sprintf("ρ = %.3f", rho_all)
    p <- ggplot(df, aes(x, y)) + geom_point(alpha = 0.35, size = 1.2, colour = colour)
  }

  p +
    annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.4,
             label = annot, size = 4.2, fontface = "bold",
             colour = "grey20", lineheight = 1.3) +
    labs(title = title_str, x = x_lab, y = "Empirical atrophy (t-statistic)") +
    theme_corr() +
    theme(plot.title = element_text(size = 11))
}

# ══════════════════════════════════════════════════════════════════════════════
#  BARPLOT HELPERS
# ══════════════════════════════════════════════════════════════════════════════

theme_bar <- function(base_size = 16) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour = "grey88"),
      axis.text.y        = element_text(size = 14, colour = "grey20", face = "bold"),
      axis.text.x        = element_text(size = 13, colour = "grey30"),
      axis.title.x       = element_text(size = 14, face = "bold", margin = margin(t = 6)),
      axis.title.y       = element_blank(),
      plot.title         = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle      = element_text(size = 13, hjust = 0.5, colour = "grey20"),
      legend.position    = "bottom",
      legend.title       = element_text(size = 13, face = "bold"),
      legend.text        = element_text(size = 12),
      panel.spacing      = unit(0.8, "lines")
    )
}

SCOPE_COLS <- c("Overall"       = "#555555",
                "Ipsilateral"   = "#E8604C",
                "Contralateral" = "#4CA8E8")

compute_rhos <- function(x_vec, y_vec, hemi_mask, log = FALSE, gt0 = FALSE) {
  if (gt0) {
    keep      <- x_vec > 0
    x_vec     <- x_vec[keep]; y_vec <- y_vec[keep]; hemi_mask <- hemi_mask[keep]
  }
  if (log) x_vec <- log(x_vec + eps)
  list(
    overall = spearman(x_vec, y_vec),
    ipsi    = spearman(x_vec[hemi_mask],  y_vec[hemi_mask]),
    contra  = spearman(x_vec[!hemi_mask], y_vec[!hemi_mask])
  )
}

# ── Read sim CSV and extract the peak-timestep vector ─────────────────────────
read_sim_peak <- function(sim_csv, target_voxels_original, keep_mask, atrophy_vec) {
  sim_mat <- as.data.frame(read.csv(sim_csv, header = FALSE))
  rownames(sim_mat) <- target_voxels_original
  colnames(sim_mat) <- seq_len(ncol(sim_mat))
  sim_mat <- sim_mat[keep_mask, ]
  peak    <- find_peak_timestep(sim_mat, atrophy_vec)
  list(vec    = as.numeric(sim_mat[, peak$peak_t]),
       peak_t = peak$peak_t,
       rho    = peak$rho)
}

# ── Read ipsi sim CSV using params$sources indexing (mirrors run_condition) ───
# For ipsi simulations the CSV rows correspond to params$sources, not $targets.
# keep_mask here is already length(sources)-based (epicentre removed).
read_sim_peak_ipsi <- function(sim_csv, params_sources, keep_mask_ipsi, atrophy_vec_ipsi) {
  target_voxels_original <- as.integer(params_sources)
  sim_mat <- as.data.frame(read.csv(sim_csv, header = FALSE))
  rownames(sim_mat) <- target_voxels_original
  colnames(sim_mat) <- seq_len(ncol(sim_mat))
  sim_mat <- sim_mat[keep_mask_ipsi, ]
  peak    <- find_peak_timestep(sim_mat, atrophy_vec_ipsi)
  list(vec    = as.numeric(sim_mat[, peak$peak_t]),
       peak_t = peak$peak_t,
       rho    = peak$rho)
}

# ── Subset bilateral vector to ipsilateral voxels ─────────────────────────────
# Mirrors the run_condition ipsi logic: sources are the first
# length(params$sources) entries of the bilateral target list.
ipsi_subset_sources <- function(vec_bilateral, n_sources) {
  # The bilateral params vector is ordered [ipsi (sources) | contra (remaining)].
  # Take the first n_sources entries.
  vec_bilateral[seq_len(n_sources)]
}

make_barplot_bilateral <- function(predictors, atrophy_vec, hemi_mask,
                                   title_str, subtitle_str = NULL, out_path) {
  rows <- purrr::imap_dfr(predictors, function(p, nm) {
    rhos <- compute_rhos(p$vec, atrophy_vec, hemi_mask, log = p$log, gt0 = p$gt0)
    data.frame(
      predictor = p$label,
      scope     = c("Overall", "Ipsilateral", "Contralateral"),
      rho       = c(rhos$overall, rhos$ipsi, rhos$contra),
      stringsAsFactors = FALSE
    )
  })
  rows$predictor <- factor(rows$predictor, levels = rev(unique(rows$predictor)))
  rows$scope     <- factor(rows$scope, levels = c("Overall","Ipsilateral","Contralateral"))

  p <- ggplot(rows, aes(x = rho, y = predictor, fill = scope)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey40") +
    scale_fill_manual(name = "Scope", values = SCOPE_COLS) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2),
                       labels = function(x) sprintf("%.1f", x)) +
    labs(title = title_str, subtitle = subtitle_str,
         x = "Spearman ρ vs. empirical atrophy") +
    theme_bar()

  ggsave(out_path, p, width = 9, height = 3.5 + 0.45 * length(predictors), dpi = 300)
  invisible(p)
}

make_barplot_ipsi <- function(predictors, atrophy_vec_ipsi,
                               title_str, subtitle_str = NULL, out_path) {
  rows <- purrr::imap_dfr(predictors, function(p, nm) {
    x <- p$vec; y <- atrophy_vec_ipsi
    if (p$gt0) { keep <- x > 0; x <- x[keep]; y <- y[keep] }
    if (p$log)  x <- log(x + eps)
    data.frame(predictor = p$label, rho = spearman(x, y), stringsAsFactors = FALSE)
  })
  rows$predictor <- factor(rows$predictor, levels = rev(unique(rows$predictor)))

  p <- ggplot(rows, aes(x = rho, y = predictor)) +
    geom_col(width = 0.55, fill = SCOPE_COLS[["Ipsilateral"]]) +
    geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey40") +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2),
                       labels = function(x) sprintf("%.1f", x)) +
    labs(title = title_str, subtitle = subtitle_str,
         x = "Spearman ρ vs. empirical atrophy (ipsilateral)") +
    theme_bar()

  ggsave(out_path, p, width = 8, height = 3 + 0.45 * length(predictors), dpi = 300)
  invisible(p)
}

# ══════════════════════════════════════════════════════════════════════════════
#  Master helper: run all plots + volumes + movie for one condition
# ══════════════════════════════════════════════════════════════════════════════
run_condition <- function(label,
                          region_name,
                          sim_csv,
                          params,
                          epi_index,
                          keep_mask,
                          conn_vec,
                          dist_vec,
                          atrophy_vec,
                          mnc_like,
                          label_vol,
                          hemi_mask_original,
                          palette   = c("navy","dodgerblue","#FFFF00","orange","#CC0000"),
                          ipsi_only = FALSE,
                          slice_seq = STR_SLICES) {

  message(sprintf("--- %s ---", label))

  if (ipsi_only) {
    # Ipsi simulations index by params$sources (ipsilateral targets by construction)
    target_voxels_original <- as.integer(params$sources)
    n_src  <- length(target_voxels_original)
    # Trim all bilateral input vectors to the first n_src entries
    keep_mask  <- keep_mask[seq_len(n_src)]
    conn_vec   <- conn_vec[seq_len(n_src)]
    dist_vec   <- dist_vec[seq_len(n_src)]
    Snca       <- as.numeric(params$Snca)[seq_len(n_src)]
    # atrophy_vec is already keep_mask-filtered; trim to match (exclude epi row)
    atrophy_vec <- atrophy_vec[seq_len(n_src - 1)]
  } else {
    target_voxels_original <- as.integer(params$targets)
    Snca <- as.numeric(params$Snca)
  }

  sim_mat <- as.data.frame(read.csv(sim_csv, header = FALSE))
  rownames(sim_mat) <- target_voxels_original
  colnames(sim_mat) <- seq_len(ncol(sim_mat))

  sim_mat        <- sim_mat[keep_mask, ]
  target_voxels  <- target_voxels_original[keep_mask]
  conn           <- conn_vec[keep_mask]
  dist           <- dist_vec[keep_mask]
  Snca           <- Snca[keep_mask]
  hemi_mask_kept <- hemi_mask_original[keep_mask]

  cat(sprintf("%s  connectivity → atrophy  rho = %.3f\n", label, spearman(conn, atrophy_vec)))
  cat(sprintf("%s  Snca         → atrophy  rho = %.3f\n", label, spearman(Snca, atrophy_vec)))
  cat(sprintf("%s  distance     → atrophy  rho = %.3f\n", label, spearman(dist, atrophy_vec)))

  output_vol_dir <- sprintf("../../derivatives/sim_sir_vols/%s/", label)
  dir.create(output_vol_dir, showWarnings = FALSE, recursive = TRUE)

  indices_to_replace <- lapply(target_voxels, function(v) which(label_vol == v))

  # Write raw sim volumes — normalisation happens at plot time
  #for (j in seq_len(ncol(sim_mat))) {
  #  out_vol <- numeric(length(label_vol))
  #  for (i in seq_along(target_voxels))
  #    out_vol[indices_to_replace[[i]]] <- sim_mat[i, j]
  #  mincWriteVolume(out_vol, like = mnc_like,
  #                  paste0(output_vol_dir, "sir_t", j, ".mnc"))
  #}

  peak <- find_peak_timestep(sim_mat, atrophy_vec)
  cat(sprintf("%s  peak t = %d  rho = %.3f\n", label, peak$peak_t, peak$rho))

  # Empirical atrophy volume ()
  emp_vol_path <- paste0(output_vol_dir, "empirical_atrophy_norm.mnc")
  make_empirical_vol(atrophy_vec, target_voxels, label_vol, mnc_like, emp_vol_path)

  # Spearman-vs-time plot
  cor_df <- data.frame(timestep = seq_along(peak$all_cors), spearman_rho = peak$all_cors)
  p_cor <- ggplot(cor_df, aes(x = timestep, y = spearman_rho)) +
    geom_line(colour = "#E05C2D", linewidth = 0.9) +
    geom_point(data = cor_df[peak$peak_t, ], colour = "#E05C2D",
               size = 4, shape = 21, fill = "white", stroke = 1.5) +
    geom_vline(xintercept = peak$peak_t, linetype = "dashed",
               colour = "grey50", linewidth = 0.6) +
    annotate("text", x = peak$peak_t + nrow(cor_df) * 0.02, y = peak$rho - 0.03,
             label = sprintf("peak t = %d\nρ = %.3f", peak$peak_t, peak$rho),
             hjust = 0, size = 6, colour = "grey20") +
    labs(title    = sprintf("%s: SIR model fit over time", region_name),
         subtitle = "Spearman ρ between simulated and empirical atrophy",
         x = "Timestep", y = "Spearman ρ") +
    theme_corr()
  ggsave(sprintf("../figures/figure_%s_cor_vs_time.png", label),
         p_cor, width = 10, height = 5, dpi = 300)

  # ── Peak brain plot (log sim) ───────────────────────────────────
  p_brain <- brain_raster_plot(
    sim_mnc_path = paste0(output_vol_dir, "sir_t", peak$peak_t, ".mnc"),
    title_str    = sprintf("%s — SIR at peak (t = %d)", region_name, peak$peak_t),
    subtitle_str = sprintf("Spearman ρ = %.3f", peak$rho),
    palette      = palette,
    scale_limits = NULL,
    slice_seq    = slice_seq,
    ipsi_only    = ipsi_only)
  ggsave(sprintf("../figures/figure_%s_peak_sim.png", label),
         p_brain, width = 24, height = 3, dpi = 300)

  # ── Two-row comparison: sim (top) vs empirical (bottom), dual colorbars ──────
  p_compare <- peak_comparison_plot(
    sim_mnc_path = paste0(output_vol_dir, "sir_t", peak$peak_t, ".mnc"),
    emp_mnc_path = emp_vol_path,
    title_str    = sprintf("%s — simulated vs empirical atrophy at peak (t = %d)",
                           region_name, peak$peak_t),
    subtitle_str = sprintf("Spearman ρ = %.3f", peak$rho),
    palette      = palette,
    slice_seq    = slice_seq,
    ipsi_only    = ipsi_only)
  
  if(ipsi_only == FALSE) {
    ggsave(sprintf("../figures/figure_%s_peak_comparison.png", label),
         p_compare, width = 24, height = 6, dpi = 300)
  }

  else {
     ggsave(sprintf("../figures/figure_%s_peak_comparison.png", label),
         p_compare, width = 24, height = 9, dpi = 300)
  }

  # Scatterplots
  # Ipsi-only conditions are already restricted to one hemisphere, so the
  # R/L/all breakdown is meaningless there — show a single overall rho instead.
  scatter_hemi_mask <- if (ipsi_only) NULL else hemi_mask_kept

  p_sir  <- make_scatter(as.numeric(sim_mat[, peak$peak_t]), atrophy_vec,
                         colour = "#E05C2D",
                         x_lab = sprintf("Log simulated atrophy — SIR (t = %d)", peak$peak_t),
                         title_str = sprintf("%s: empirical atrophy vs. SIR", region_name),
                         log = TRUE, hemi_mask = scatter_hemi_mask)
  p_conn <- make_scatter(conn, atrophy_vec, colour = "#7B3F9E",
                         x_lab = "Log connectivity (> 0) to epicentre",
                         title_str = sprintf("%s: empirical atrophy vs. conn.", region_name),
                         log = TRUE, gt0 = TRUE, hemi_mask = scatter_hemi_mask)
  p_snca <- make_scatter(Snca, atrophy_vec, colour = "#2A7A4B",
                         x_lab = "Snca expression",
                         title_str = sprintf("%s: Snca vs. empirical atrophy", region_name),
                         log = FALSE, hemi_mask = scatter_hemi_mask)
  p_dist <- make_scatter(dist, atrophy_vec, colour = "#B05020",
                         x_lab = "Distance to epicentre",
                         title_str = sprintf("%s: empirical atrophy vs. distance", region_name),
                         log = FALSE, hemi_mask = scatter_hemi_mask)

  ggsave(sprintf("../figures/figure_%s_scatter_sir.png",  label), p_sir,  width = 7, height = 6, dpi = 300)
  ggsave(sprintf("../figures/figure_%s_scatter_conn.png", label), p_conn, width = 7, height = 6, dpi = 300)
  ggsave(sprintf("../figures/figure_%s_scatter_snca.png", label), p_snca, width = 7, height = 6, dpi = 300)
  ggsave(sprintf("../figures/figure_%s_scatter_dist.png", label), p_dist, width = 7, height = 6, dpi = 300)

  # ── Movie: log, fixed [−1, 1] limits, frames up to peak ────────
  message(sprintf("Generating %s movie frames (up to peak t = %d) …", label, peak$peak_t))
  frames_dir <- sprintf("../../derivatives/sim_sir_vols/%s/frames/", label)
  dir.create(frames_dir, showWarnings = FALSE, recursive = TRUE)

  #for (j in seq_len(peak$peak_t)) {
  for (j in c(1,2)) {
    p_frame <- brain_raster_plot(
      sim_mnc_path = paste0(output_vol_dir, "sir_t", j, ".mnc"),
      title_str    = sprintf("%s — SIR  t = %04d", region_name, j),
      subtitle_str = sprintf("Spearman ρ = %.3f", peak$all_cors[j]),
      palette      = palette,
      scale_limits = NULL,       # per-frame: [min, max] of that frame
      slice_seq    = slice_seq,
      ipsi_only    = ipsi_only)
    if(ipsi_only == FALSE) {
        ggsave(sprintf("%sframe_%05d.png", frames_dir, j),
          p_frame, width = 24, height = 3, dpi = 150)

    } else {
         ggsave(sprintf("%sframe_%05d.png", frames_dir, j),
          p_frame, width = 24, height = 4.5, dpi = 150)
    }
    if (j %% 50 == 0)
      message(sprintf("  %s frame %d / %d", label, j, peak$peak_t))
  }

  system(paste0(
    "ffmpeg -y -framerate 25 ",
   "-pattern_type glob -i '", frames_dir, "frame_*.png' ",
   "-vf 'scale=trunc(iw/2)*2:trunc(ih/2)*2' ",
    "-c:v libx264 -pix_fmt yuv420p -crf 18 ",
    sprintf("../figures/%s_sir_movie.mp4", label)))
    message(sprintf("%s movie → ../figures/%s_sir_movie.mp4", label, label))

  invisible(list(peak = peak, sim_mat = sim_mat))
}

# ── Scatter-only helper ────────────────────────────────────────────────────────
run_scatter_only <- function(label, region_name, params, epi_index,
                             keep_mask, conn_vec, dist_vec, atrophy_vec,
                             hemi_mask_original) {
  message(sprintf("--- %s (scatter only) ---", label))
  conn           <- conn_vec[keep_mask]
  dist           <- dist_vec[keep_mask]
  Snca           <- as.numeric(params$Snca)[keep_mask]
  hemi_mask_kept <- hemi_mask_original[keep_mask]

  cat(sprintf("%s  connectivity → atrophy  rho = %.3f\n", label, spearman(conn, atrophy_vec)))
  cat(sprintf("%s  Snca         → atrophy  rho = %.3f\n", label, spearman(Snca, atrophy_vec)))
  cat(sprintf("%s  distance     → atrophy  rho = %.3f\n", label, spearman(dist, atrophy_vec)))

  p_conn <- make_scatter(conn, atrophy_vec, colour = "#7B3F9E",
                         x_lab = "Log connectivity (> 0) to epicentre",
                         title_str = sprintf("%s: empirical atrophy vs. conn.", region_name),
                         log = TRUE, gt0 = TRUE, hemi_mask = hemi_mask_kept)
  p_snca <- make_scatter(Snca, atrophy_vec, colour = "#2A7A4B",
                         x_lab = "Snca expression",
                         title_str = sprintf("%s: Snca vs. empirical atrophy", region_name),
                         log = FALSE, hemi_mask = hemi_mask_kept)
  p_dist <- make_scatter(dist, atrophy_vec, colour = "#B05020",
                         x_lab = "Distance to epicentre",
                         title_str = sprintf("%s: empirical atrophy vs. distance", region_name),
                         log = FALSE, hemi_mask = hemi_mask_kept)

  ggsave(sprintf("../figures/figure_%s_scatter_conn.png", label), p_conn, width = 7, height = 6, dpi = 300)
  ggsave(sprintf("../figures/figure_%s_scatter_snca.png", label), p_snca, width = 7, height = 6, dpi = 300)
  ggsave(sprintf("../figures/figure_%s_scatter_dist.png", label), p_dist, width = 7, height = 6, dpi = 300)
  invisible(NULL)
}

# ══════════════════════════════════════════════════════════════════════════════
#  WHOLE-BRAIN  (scatter only)
# ══════════════════════════════════════════════════════════════════════════════
message("=== WHOLE-BRAIN CP retro (scatter only) ===")

wb_mnc_like  <- "../../derivatives/downsampled_connectome/target_indices_filt_overlap_source_200um.mnc"
wb_label_vol <- round(mincGetVolume(wb_mnc_like))

params_wb <- pd$read_pickle(
  "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/params/source_target_indices_filt_overlap_top30.pkl")

target_voxels_in_order_original_wb <- as.integer(params_wb$targets)
n_wb              <- length(target_voxels_in_order_original_wb)
hemi_mask_wb_orig <- seq_len(n_wb) <= floor(n_wb / 2)

keep_mask_wb_cp <- target_voxels_in_order_original_wb != cp_epi_index
conn_row_wb_cp  <- as.numeric(params_wb$weights[which(target_voxels_in_order_original_wb == cp_epi_index), ])
dist_row_wb_cp  <- as.numeric(params_wb$distance[which(target_voxels_in_order_original_wb == cp_epi_index), ])

atrophy_wb_cp <- as.numeric(as.data.frame(
  pd$read_pickle(
    "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/atrophy_pkl/M83PBS_vs_M83MsPff_t_stats_linear_allen200_PIR_CCFv3_target_filt.pkl"
  ))[keep_mask_wb_cp])

run_scatter_only(
  label = "wb_retro_cp", region_name = "Whole brain (retro)",
  params = params_wb, epi_index = cp_epi_index,
  keep_mask = keep_mask_wb_cp, conn_vec = conn_row_wb_cp, dist_vec = dist_row_wb_cp,
  atrophy_vec = atrophy_wb_cp, hemi_mask_original = hemi_mask_wb_orig)

message("=== WHOLE-BRAIN DG antero (scatter only) ===")

keep_mask_wb_dg <- target_voxels_in_order_original_wb != dg_epi_index
conn_row_wb_dg  <- as.numeric(params_wb$weights[which(target_voxels_in_order_original_wb == dg_epi_index), ])
dist_row_wb_dg  <- as.numeric(params_wb$distance[which(target_voxels_in_order_original_wb == dg_epi_index), ])

atrophy_wb_dg <- as.numeric(as.data.frame(
  pd$read_pickle(
    "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/atrophy_pkl/M83PBS_vs_M83HuPff_hipp_t_stats_linear_allen200_PIR_CCFv3_target_filt.pkl"
  ))[keep_mask_wb_dg])

run_scatter_only(
  label = "wb_antero_dg", region_name = "Whole brain (antero)",
  params = params_wb, epi_index = dg_epi_index,
  keep_mask = keep_mask_wb_dg, conn_vec = conn_row_wb_dg, dist_vec = dist_row_wb_dg,
  atrophy_vec = atrophy_wb_dg, hemi_mask_original = hemi_mask_wb_orig)

message("=== DONE ===")

# ══════════════════════════════════════════════════════════════════════════════
#  STRIATUM  (CP epicentre)
# ══════════════════════════════════════════════════════════════════════════════
message("=== STRIATUM ===")

str_mnc_like  <- "../../derivatives/downsampled_connectome/str_target_indices_filt_overlap_source_200um.mnc"
str_label_vol <- round(mincGetVolume(str_mnc_like))

params_str_retro  <- pd$read_pickle(pkl_str_retro_params)
params_str_antero <- pd$read_pickle(pkl_str_antero_params)

target_voxels_str_retro  <- as.integer(params_str_retro$targets)
target_voxels_str_antero <- as.integer(params_str_antero$targets)

keep_mask_str_retro  <- target_voxels_str_retro  != cp_epi_index
keep_mask_str_antero <- target_voxels_str_antero != cp_epi_index

n_str_retro               <- length(target_voxels_str_retro)
n_str_antero              <- length(target_voxels_str_antero)
hemi_mask_str_retro_orig  <- seq_len(n_str_retro)  <= floor(n_str_retro  / 2)
hemi_mask_str_antero_orig <- seq_len(n_str_antero) <= floor(n_str_antero / 2)

conn_row_str_retro  <- as.numeric(params_str_retro$weights[
  which(target_voxels_str_retro  == cp_epi_index), ])
conn_row_str_antero <- as.numeric(params_str_antero$weights[
  which(target_voxels_str_antero == cp_epi_index), ])

dist_row_str_retro  <- as.numeric(params_str_retro$distance[
  which(target_voxels_str_retro  == cp_epi_index), ])
dist_row_str_antero <- as.numeric(params_str_antero$distance[
  which(target_voxels_str_antero == cp_epi_index), ])

atrophy_str <- as.numeric(as.data.frame(pd$read_pickle(pkl_str_atrophy))[keep_mask_str_retro])

# ── STR retro (bilateral) ──────────────────────────────────────────────────────
 run_condition(
   label = "str_retro", region_name = "Striatum (retro)",
   sim_csv = sim_file_str_retro, params = params_str_retro,
   epi_index = cp_epi_index, keep_mask = keep_mask_str_retro,
   conn_vec = conn_row_str_retro, dist_vec = dist_row_str_retro,
   atrophy_vec = atrophy_str, mnc_like = str_mnc_like,
   label_vol = str_label_vol, hemi_mask_original = hemi_mask_str_retro_orig,
   ipsi_only = FALSE, slice_seq = STR_SLICES)

# ── STR antero (bilateral) ─────────────────────────────────────────────────────
 run_condition(
   label = "str_antero", region_name = "Striatum (antero)",
   sim_csv = sim_file_str_antero, params = params_str_antero,
   epi_index = cp_epi_index, keep_mask = keep_mask_str_antero,
   conn_vec = conn_row_str_antero, dist_vec = dist_row_str_antero,
   atrophy_vec = atrophy_str, mnc_like = str_mnc_like,
   label_vol = str_label_vol, hemi_mask_original = hemi_mask_str_antero_orig,
   ipsi_only = FALSE, slice_seq = STR_SLICES)

# ── STR retro ipsilateral ──────────────────────────────────────────────────────
 run_condition(
   label = "str_retro_ipsi", region_name = "Striatum (retro, ipsi)",
   sim_csv = sim_file_str_retro_ipsi, params = params_str_retro,
   epi_index = cp_epi_index, keep_mask = keep_mask_str_retro,
   conn_vec = conn_row_str_retro, dist_vec = dist_row_str_retro,
   atrophy_vec = atrophy_str, mnc_like = str_mnc_like,
   label_vol = str_label_vol, hemi_mask_original = hemi_mask_str_retro_orig,
   ipsi_only = TRUE, slice_seq = STR_SLICES)

# ── STR retro ipsilateral + clearance gene ─────────────────────────────────────
 run_condition(
   label = "str_retro_cl_ipsi", region_name = "Striatum (retro, ipsi, clearance)",
   sim_csv = sim_file_cp_retro_cl_ipsi, params = params_str_retro,
   epi_index = cp_epi_index, keep_mask = keep_mask_str_retro,
   conn_vec = conn_row_str_retro, dist_vec = dist_row_str_retro,
   atrophy_vec = atrophy_str, mnc_like = str_mnc_like,
   label_vol = str_label_vol, hemi_mask_original = hemi_mask_str_retro_orig,
   ipsi_only = TRUE, slice_seq = STR_SLICES)

# ══════════════════════════════════════════════════════════════════════════════
#  HIPPOCAMPUS  (DG epicentre)
# ══════════════════════════════════════════════════════════════════════════════
message("=== HIPPOCAMPUS ===")

hip_mnc_like  <- "../../derivatives/downsampled_connectome/hip_target_indices_filt_overlap_source_200um.mnc"
hip_label_vol <- round(mincGetVolume(hip_mnc_like))

params_hip_retro  <- pd$read_pickle(pkl_hip_retro_params)
params_hip_antero <- pd$read_pickle(pkl_hip_antero_params)

target_voxels_hip_retro  <- as.integer(params_hip_retro$targets)
target_voxels_hip_antero <- as.integer(params_hip_antero$targets)

keep_mask_hip_retro  <- target_voxels_hip_retro  != dg_epi_index
keep_mask_hip_antero <- target_voxels_hip_antero != dg_epi_index

n_hip_retro               <- length(target_voxels_hip_retro)
n_hip_antero              <- length(target_voxels_hip_antero)
hemi_mask_hip_retro_orig  <- seq_len(n_hip_retro)  <= floor(n_hip_retro  / 2)
hemi_mask_hip_antero_orig <- seq_len(n_hip_antero) <= floor(n_hip_antero / 2)

conn_row_hip_retro  <- as.numeric(params_hip_retro$weights[
  which(target_voxels_hip_retro  == dg_epi_index), ])
conn_row_hip_antero <- as.numeric(params_hip_antero$weights[
  which(target_voxels_hip_antero == dg_epi_index), ])

dist_row_hip_retro  <- as.numeric(params_hip_retro$distance[
  which(target_voxels_hip_retro  == dg_epi_index), ])
dist_row_hip_antero <- as.numeric(params_hip_antero$distance[
  which(target_voxels_hip_antero == dg_epi_index), ])

atrophy_hip <- as.numeric(as.data.frame(pd$read_pickle(pkl_hip_atrophy))[keep_mask_hip_antero])

hip_palette <- c("navy", "dodgerblue", "#FFFF00", "orange", "#CC0000")

# ── HIP retro (bilateral) ──────────────────────────────────────────────────────
run_condition(
  label = "hip_retro", region_name = "Hippocampus (retro)",
  sim_csv = sim_file_hip_retro, params = params_hip_retro,
  epi_index = dg_epi_index, keep_mask = keep_mask_hip_retro,
  conn_vec = conn_row_hip_retro, dist_vec = dist_row_hip_retro,
  atrophy_vec = atrophy_hip, mnc_like = hip_mnc_like,
  label_vol = hip_label_vol, hemi_mask_original = hemi_mask_hip_retro_orig,
  palette = hip_palette, ipsi_only = FALSE, slice_seq = HIP_SLICES)

# ── HIP antero (bilateral) ─────────────────────────────────────────────────────
run_condition(
  label = "hip_antero", region_name = "Hippocampus (antero)",
  sim_csv = sim_file_hip_antero, params = params_hip_antero,
  epi_index = dg_epi_index, keep_mask = keep_mask_hip_antero,
  conn_vec = conn_row_hip_antero, dist_vec = dist_row_hip_antero,
  atrophy_vec = atrophy_hip, mnc_like = hip_mnc_like,
  label_vol = hip_label_vol, hemi_mask_original = hemi_mask_hip_antero_orig,
  palette = hip_palette, ipsi_only = FALSE, slice_seq = HIP_SLICES)

# ── HIP antero ipsilateral ─────────────────────────────────────────────────────
run_condition(
  label = "hip_antero_ipsi", region_name = "Hippocampus (antero, ipsi)",
  sim_csv = sim_file_hip_antero_ipsi, params = params_hip_antero,
  epi_index = dg_epi_index, keep_mask = keep_mask_hip_antero,
  conn_vec = conn_row_hip_antero, dist_vec = dist_row_hip_antero,
  atrophy_vec = atrophy_hip, mnc_like = hip_mnc_like,
  label_vol = hip_label_vol, hemi_mask_original = hemi_mask_hip_antero_orig,
  palette = hip_palette, ipsi_only = TRUE, slice_seq = HIP_SLICES)

# ── HIP antero ipsilateral + clearance gene ────────────────────────────────────
run_condition(
  label = "hip_antero_cl_ipsi", region_name = "Hippocampus (antero, ipsi, clearance)",
  sim_csv = sim_file_hip_antero_cl_ipsi, params = params_hip_antero,
  epi_index = dg_epi_index, keep_mask = keep_mask_hip_antero,
  conn_vec = conn_row_hip_antero, dist_vec = dist_row_hip_antero,
  atrophy_vec = atrophy_hip, mnc_like = hip_mnc_like,
  label_vol = hip_label_vol, hemi_mask_original = hemi_mask_hip_antero_orig,
  palette = hip_palette, ipsi_only = TRUE, slice_seq = HIP_SLICES)

# ══════════════════════════════════════════════════════════════════════════════
#  BARPLOTS
#  BP1a  CP epicentre — bilateral  (3 bars: overall / ipsi / contra)
#  BP1b  DG epicentre — bilateral  (3 bars)
#  BP2a  CP epicentre — ipsilateral only (1 bar)
#  BP2b  DG epicentre — ipsilateral only (1 bar)
# ══════════════════════════════════════════════════════════════════════════════
message("=== BARPLOTS ===")
dir.create("../figures", showWarnings = FALSE)

# ── Ipsi source counts (for barplot ipsi subsetting) ──────────────────────────
n_str_retro_src  <- length(as.integer(params_str_retro$sources))
n_str_antero_src <- length(as.integer(params_str_antero$sources))
n_hip_retro_src  <- length(as.integer(params_hip_retro$sources))
n_hip_antero_src <- length(as.integer(params_hip_antero$sources))

# ── Bilateral atrophy for barplots ────────────────────────────────────────────
# (same as used in run_condition bilateral)
atrophy_str_bil <- atrophy_str   # already keep_mask_str_retro filtered
atrophy_hip_bil <- atrophy_hip   # already keep_mask_hip_antero filtered

# ── Ipsilateral atrophy for barplots ──────────────────────────────────────────
# Mirror run_condition ipsi: take first (n_src - 1) entries of the bilateral
# atrophy vector (which is itself keep_mask-filtered, so the epi is already gone;
# trim to match the n_src - 1 ipsi targets excluding the epi row).
atrophy_str_ipsi <- atrophy_str_bil[seq_len(n_str_retro_src - 1)]
atrophy_hip_ipsi <- atrophy_hip_bil[seq_len(n_hip_antero_src - 1)]

# ── Peak sim vectors (bilateral) ───────────────────────────────────────────────
cp_retro_bil  <- read_sim_peak(sim_file_str_retro,  target_voxels_str_retro,
                               keep_mask_str_retro,  atrophy_str_bil)
cp_antero_bil <- read_sim_peak(sim_file_str_antero, target_voxels_str_antero,
                               keep_mask_str_antero, atrophy_str_bil)
dg_retro_bil  <- read_sim_peak(sim_file_hip_retro,  target_voxels_hip_retro,
                               keep_mask_hip_retro,  atrophy_hip_bil)
dg_antero_bil <- read_sim_peak(sim_file_hip_antero, target_voxels_hip_antero,
                               keep_mask_hip_antero, atrophy_hip_bil)

# ── Peak sim vectors (ipsilateral) — use sources indexing ─────────────────────
cp_retro_ipsi_sim <- read_sim_peak_ipsi(
  sim_file_str_retro_ipsi,
  params_str_retro$sources,
  keep_mask_str_retro[seq_len(n_str_retro_src)],
  atrophy_str_ipsi)

# cp_antero_ipsi_sim <- read_sim_peak_ipsi(   # <-- uncomment when file available
#   sim_file_str_antero_ipsi,
#   params_str_antero$sources,
#   keep_mask_str_antero[seq_len(n_str_antero_src)],
#   atrophy_str_ipsi)

dg_antero_ipsi_sim <- read_sim_peak_ipsi(
  sim_file_hip_antero_ipsi,
  params_hip_antero$sources,
  keep_mask_hip_antero[seq_len(n_hip_antero_src)],
  atrophy_hip_ipsi)

# dg_retro_ipsi_sim <- read_sim_peak_ipsi(    # <-- uncomment when file available
#   sim_file_hip_retro_ipsi,
#   params_hip_retro$sources,
#   keep_mask_hip_retro[seq_len(n_hip_retro_src)],
#   atrophy_hip_ipsi)

# ── Peak sim vectors (ipsi + clearance gene) ──────────────────────────────────
cp_retro_cl_ipsi_sim <- read_sim_peak_ipsi(
  sim_file_cp_retro_cl_ipsi,
  params_str_retro$sources,
  keep_mask_str_retro[seq_len(n_str_retro_src)],
  atrophy_str_ipsi)

dg_antero_cl_ipsi_sim <- read_sim_peak_ipsi(
  sim_file_hip_antero_cl_ipsi,
  params_hip_antero$sources,
  keep_mask_hip_antero[seq_len(n_hip_antero_src)],
  atrophy_hip_ipsi)

# ── Ipsilateral parameter vectors (first n_src entries of bilateral vec) ───────
# This mirrors run_condition: bilateral params are ordered [ipsi | contra],
# so the first n_src entries are ipsilateral.
str_retro_snca_ipsi  <- ipsi_subset_sources(as.numeric(params_str_retro$Snca)[keep_mask_str_retro],   n_str_retro_src)
str_antero_snca_ipsi <- ipsi_subset_sources(as.numeric(params_str_antero$Snca)[keep_mask_str_antero], n_str_antero_src)

conn_str_retro_ipsi  <- ipsi_subset_sources(conn_row_str_retro[keep_mask_str_retro],   n_str_retro_src)
conn_str_antero_ipsi <- ipsi_subset_sources(conn_row_str_antero[keep_mask_str_antero], n_str_antero_src)

dist_str_retro_ipsi  <- ipsi_subset_sources(dist_row_str_retro[keep_mask_str_retro],   n_str_retro_src)

hip_antero_snca_ipsi <- ipsi_subset_sources(as.numeric(params_hip_antero$Snca)[keep_mask_hip_antero], n_hip_antero_src)
hip_retro_snca_ipsi  <- ipsi_subset_sources(as.numeric(params_hip_retro$Snca)[keep_mask_hip_retro],   n_hip_retro_src)

conn_hip_antero_ipsi <- ipsi_subset_sources(conn_row_hip_antero[keep_mask_hip_antero], n_hip_antero_src)
conn_hip_retro_ipsi  <- ipsi_subset_sources(conn_row_hip_retro[keep_mask_hip_retro],   n_hip_retro_src)

dist_hip_antero_ipsi <- ipsi_subset_sources(dist_row_hip_antero[keep_mask_hip_antero], n_hip_antero_src)

# ══════════════════════════════════════════════════════════════════════════════
#  BP1a — CP epicentre, bilateral
# ══════════════════════════════════════════════════════════════════════════════
bp1a_predictors <- list(
  snca = list(
    vec = as.numeric(params_str_retro$Snca)[keep_mask_str_retro],
    log = FALSE, gt0 = FALSE, label = "Snca expression"),
  conn_antero = list(
    vec = conn_row_str_antero[keep_mask_str_antero],
    log = TRUE, gt0 = TRUE, label = "Anterograde connectivity"),
  conn_retro = list(
    vec = conn_row_str_retro[keep_mask_str_retro],
    log = TRUE, gt0 = TRUE, label = "Retrograde connectivity"),
  distance = list(
    vec = dist_row_str_retro[keep_mask_str_retro],
    log = FALSE, gt0 = FALSE, label = "Distance to epicentre"),
  sim_retro = list(
    vec = cp_retro_bil$vec, log = TRUE, gt0 = FALSE,
    label = sprintf("Sim atrophy — retro (t=%d)", cp_retro_bil$peak_t)),
  sim_antero = list(
    vec = cp_antero_bil$vec, log = TRUE, gt0 = FALSE,
    label = sprintf("Sim atrophy — antero (t=%d)", cp_antero_bil$peak_t))
)

make_barplot_bilateral(
  predictors   = bp1a_predictors,
  atrophy_vec  = atrophy_str_bil,
  hemi_mask    = hemi_mask_str_retro_orig[keep_mask_str_retro],
  title_str    = "CP epicentre — bilateral predictors vs. empirical striatal atrophy",
  subtitle_str = "Three bars per predictor: overall / ipsilateral / contralateral",
  out_path     = "../figures/figure_barplot_bilateral_cp.png")

# ══════════════════════════════════════════════════════════════════════════════
#  BP1b — DG epicentre, bilateral
# ══════════════════════════════════════════════════════════════════════════════
bp1b_predictors <- list(
  snca = list(
    vec = as.numeric(params_hip_antero$Snca)[keep_mask_hip_antero],
    log = FALSE, gt0 = FALSE, label = "Snca expression"),
  conn_antero = list(
    vec = conn_row_hip_antero[keep_mask_hip_antero],
    log = TRUE, gt0 = TRUE, label = "Anterograde connectivity"),
  conn_retro = list(
    vec = conn_row_hip_retro[keep_mask_hip_retro],
    log = TRUE, gt0 = TRUE, label = "Retrograde connectivity"),
  distance = list(
    vec = dist_row_hip_antero[keep_mask_hip_antero],
    log = FALSE, gt0 = FALSE, label = "Distance to epicentre"),
  sim_retro = list(
    vec = dg_retro_bil$vec, log = TRUE, gt0 = FALSE,
    label = sprintf("Sim atrophy — retro (t=%d)", dg_retro_bil$peak_t)),
  sim_antero = list(
    vec = dg_antero_bil$vec, log = TRUE, gt0 = FALSE,
    label = sprintf("Sim atrophy — antero (t=%d)", dg_antero_bil$peak_t))
)

make_barplot_bilateral(
  predictors   = bp1b_predictors,
  atrophy_vec  = atrophy_hip_bil,
  hemi_mask    = hemi_mask_hip_antero_orig[keep_mask_hip_antero],
  title_str    = "DG epicentre — bilateral predictors vs. empirical hippocampal atrophy",
  subtitle_str = "Three bars per predictor: overall / ipsilateral / contralateral",
  out_path     = "../figures/figure_barplot_bilateral_dg.png")