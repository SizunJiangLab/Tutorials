

# Style -------------------------------------------------------------------

theme_col_count <- function() {
  list(
    # scale_y_continuous(expand = expansion(c(0, 0.05))), 
    theme_bw(), 
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
      legend.position = "top", 
      panel.grid.major.x = element_blank(), 
      panel.grid.minor.x = element_blank(), 
      axis.ticks.x = element_blank()
    )
  )
}

theme_col_percent <- function() {
  list(
    # scale_y_continuous(expand = expansion(c(0, 0))), 
    theme_bw(), 
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
      legend.position = "top", 
      panel.grid.major.x = element_blank(), 
      panel.grid.minor.x = element_blank(), 
      axis.ticks.x = element_blank()
    )
  )
}

theme_boxplot <- function() {
  list(
    # scale_y_continuous(expand = expansion(c(0.05, 0.05))), 
    theme_bw(), 
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
      panel.grid.major.x = element_blank(), 
      panel.grid.minor.x = element_blank(), 
      axis.ticks.x = element_blank()
    )
  )
}

COLOR <- c(
  "None" = "grey", 
  "Nuclear" = "#377EB8", 
  "Cytoplasm" = "#4DAF4A", 
  "Membrane" = "#E41A1C",
  
  "Unassigned" = "grey", 
  "Assigned" = "black", 
  
  "Gene" = "#E41A1C", 
  "Negative" = "#377EB8", 
  "SystemControl" = "#4DAF4A"
)


# Data Import -------------------------------------------------------------

#' Find paths of all CosMx output files for QC
#'
#' @param folder_input The path of the CosMx foldler following data standard (AtoMx folder).
#'
#' @return A list with the paths of all CosMx output files for QC. 
#' @export
#'
find_all_files <- function(folder_input){
  # pattern for flatfiles
  pattern_flatfile_list <- list(
    exprMat_file = ".*exprMat_file[.]csv[.]gz$", 
    metadata_file = ".*metadata_file[.]csv[.]gz$", 
    tx_file = ".*tx_file[.]csv[.]gz$", 
    polygons = ".*polygons[.]csv[.]gz$", 
    fov_positions_file = ".*fov_positions_file[.]csv[.]gz$"
  )
  # path for flatfiles
  flatfiles <- fs::dir_ls(fs::path(folder_input, "flatFiles/"), recurse = TRUE)
  flatfile_list <- pattern_flatfile_list %>% 
    purrr::map(~ stringr::str_subset(flatfiles, .x))
  
  # pattern for rawfiles
  pattern_rawfile_list <- list(
    fovs = "Logs/.*[.]fovs[.]csv$"
  )
  # path for rawfiles
  rawfiles <- fs::dir_ls(fs::path(folder_input, "RawFiles/"), recurse = TRUE)
  rawfile_list <- pattern_rawfile_list %>% 
    purrr::map(~ stringr::str_subset(rawfiles, .x))
  
  # path for all files
  path_list <- c(flatfile_list, rawfile_list)
  missing_files <- names(which(purrr::map_lgl(path_list, ~ length(.x) == 0)))
  if (length(missing_files) > 0) {
    cat("The following files are missing:\n")
    cat(paste(seq_along(missing_files), missing_files, sep = ". ", collapse = "\n"))
    cat("\nPlease make sure all required files exist and try again.\n")
    # Halt execution
    # q(save = "no",status = 1)
  }
  return(path_list)
}


#' Import all CosMx output files for QC
#' 
#' 1. Read all CosMx output files for QC. 
#' 2. Add information of FOV (locations and order).
#'
#' @param path_list A list with the paths of all CosMx output files for QC. 
#'
#' @return A list with the data of all CosMx output files for QC. 
#' @export
#'
load_all_files <- function(path_list){
  # Read data files
  exprMat_sub <- data.table::fread(path_list$exprMat_file) %>% distinct()
  metadata_sub <- data.table::fread(path_list$metadata_file) %>% distinct()
  tx_sub <- data.table::fread(path_list$tx_file) %>% distinct()
  polygons_sub <- data.table::fread(path_list$polygons) %>% distinct()
  positions <- data.table::fread(path_list$fov_positions_file) %>% distinct()
  positions <- positions %>% 
    rename_with(~ "X_mm", .cols = any_of(c("X_mm", "x_global_mm"))) %>% 
    rename_with(~ "Y_mm", .cols = any_of(c("Y_mm", "y_global_mm"))) 
  log.file <- data.table::fread(path_list$fovs) %>% distinct()
  
  SlideIndex <- unique(metadata_sub$slide_ID)
  # positions <- positions[positions$Slide==SlideIndex,]
  log.file <- log.file[log.file$Slide == SlideIndex,]
  
  # position of FOV
  metadata_sub <- metadata_sub %>% 
    left_join(
      positions[, c("FOV", "X_mm", "Y_mm")], 
      by = c("fov" = "FOV"), 
      relationship = "many-to-one"
    )
  metadata_sub$fov <- factor(metadata_sub$fov, levels = sort(unique(metadata_sub$fov)))
  tx_sub <- tx_sub %>% 
    left_join(
      positions[, c("FOV", "X_mm", "Y_mm")], 
      by = c("fov" = "FOV"), 
      relationship = "many-to-one"
    )
  tx_sub$fov <- factor(tx_sub$fov, levels = sort(unique(tx_sub$fov)))
  
  # factor: assignment, cell component, gene type
  target <- colnames(exprMat_sub)
  target <- target[!target %in% c("fov", "cell_ID", "cell")]
  target_system <- target[grepl("SystemControl", target)]
  target_negative <- target[grepl("Negative", target)]
  target_genetype <- setdiff(target, c(target_system, target_negative))
  df_target_genetype <- bind_rows(
    data.frame(target = target_genetype, target_type = "Gene"), 
    data.frame(target = target_system, target_type = "SystemControl"), 
    data.frame(target = target_negative, target_type = "Negative")
  )
  tx_sub <- tx_sub %>% 
    left_join(df_target_genetype, by = "target", relationship = "many-to-one") %>% 
    mutate(
      assignment = ifelse(cell_ID == 0, "Unassigned", "Assigned"), 
      across(assignment, ~ factor(.x, levels = c("Unassigned", "Assigned"))), 
      across(target_type, ~ factor(.x, levels = c("Negative", "SystemControl", "Gene"))), 
      across(CellComp, ~ factor(.x, levels = c("None", "Membrane", "Cytoplasm", "Nuclear")))
    )
  
  combined_data <- list(
    exprMat_sub = exprMat_sub,
    metadata_sub = metadata_sub,
    tx_sub = tx_sub,
    polygons_sub = polygons_sub,
    positions = positions,
    log = log.file
  )
  return(combined_data)
}


# General functions #####

## Range of pdata #####
get_pdata_range <- function(list_p, col = "n", quantile_cut = 1, vmin = 0, vmax = NULL) {
  v_max <- max(unlist(purrr::map(list_p, ~ quantile(.x$data[[col]], quantile_cut))))
  v_min <- min(unlist(purrr::map(list_p, ~ quantile(.x$data[[col]], 1 - quantile_cut))))
  if (!is.null(vmin)) {v_min <- min(v_min, vmin)}
  if (!is.null(vmax)) {v_min <- min(v_min, vmax)}
  v_range <- c(v_min, v_max)
  return(v_range)
}

## combined plot #####
plt_combined <- function(list_p, list_title = list_run_name) {
  p_combined <- purrr::map2(
    list_p, 
    list_title, 
    ~ .x + labs(subtitle = .y)
  ) %>% 
    wrap_plots(guides = "collect", widths = rep(1, length(list_title))) &
    theme(
      plot.subtitle = element_text(hjust = .5), 
      legend.position = "top"
    )
  return(p_combined)
}

plt_combined_grid <- function(list_p, list_title = list_run_name) {
  p_combined <- list_output[[1]]
  p_combined$data <- purrr::map2(
    list_output, 
    list_run_name, 
    ~ .x$data %>% mutate(tag_facet = .y)
  ) %>% 
    bind_rows() %>% 
    mutate(across(tag_facet, ~ factor(.x, levels = list_run_name)))
  p_combined <- p_combined +
    facet_grid(~ tag_facet, scales = "free")
  return(p_combined)
}


# Global-View Of Entire Run #####

#' Data for Global Summary
#'
#' @param combined_data A list with the data of all CosMx output files for QC. 
#'
#' @return A list with the data for global summary. 
#' @export 
#'
global_summarize_Fov <- function(combined_data){
  SlideIndex <- unique(combined_data$metadata_sub$slide_ID)
  RunName <- unique(combined_data$metadata_sub$Run_Tissue_name)
  # CellNumberTotal <- length(unique(combined_data$exprMat_sub$cell))
  FOVNumberTotal <- length(unique(combined_data$exprMat_sub$fov))
  TranscriptsNumberTotal <- nrow(combined_data$tx_sub)
  
  target_real <- sum(combined_data$tx_sub$target_type == "Gene")
  target_negative <- sum(combined_data $tx_sub$target_type == "Negative")
  target_false <- sum(combined_data $tx_sub$target_type == "SystemControl")
  # UnassignedTranscriptsNumberTotal <- sum(combined_data$tx_sub$cell_ID == 0)
  # CytoplasmTranscriptsNumberTotal <- sum(combined_data$tx_sub$CellComp == "Cytoplasm")
  # MembraneTranscriptsNumberTotal <- sum(combined_data$tx_sub$CellComp == "Membrane")
  # NuclearTranscriptsNumberTotal <- sum(combined_data$tx_sub$CellComp == "Nuclear")
  res_out <- list(
    "Run name" = RunName,
    "Slide Index" = SlideIndex,
    "Total FOV number" = FOVNumberTotal,
    # "Total Cell Number" = CellNumberTotal,
    "Total Transcripts" = TranscriptsNumberTotal,
    # "Unassigned Transcripts" = UnassignedTranscriptsNumberTotal,
    # "Cytoplasm Transcripts" = CytoplasmTranscriptsNumberTotal,
    # "Membrane Transcripts" = MembraneTranscriptsNumberTotal,
    # "Nuclear Transcripts" = NuclearTranscriptsNumberTotal
    "Real Gene Transcripts" = target_real, 
    "Negative Probe Transcripts" = target_negative, 
    "System Control Transcripts" = target_false
  )
  return(res_out)
}


#' Print Table for Global Summary
#'
#' @param Global_View A list with the data for global summary. 
#'
#' @return Table for global summary. 
#' @export
#'
GlobalViewTable_func <- function(Global_View){
  GlobalViewTable <- data.frame(
    Item = names(Global_View),
    Statistics = c(
      Global_View[[1]],
      Global_View[[2]],
      Global_View[[3]],
      Global_View[[4]],
      sprintf("%.0f (%.2f)%%", Global_View[[5]], round(Global_View[[5]] / Global_View[[4]] * 100, 2)), 
      sprintf("%.0f (%.2f)%%", Global_View[[6]], round(Global_View[[6]] / Global_View[[4]] * 100, 2)), 
      sprintf("%.0f (%.2f)%%", Global_View[[7]], round(Global_View[[7]] / Global_View[[4]] * 100, 2))
    )
  )
  knitr::kable(GlobalViewTable, align = "cc")
}


# Profiling position  (FOV index) #####

# function for Profiling position and order (FOV number/ Order)
plt_fov_position <- function(combined_data){
  p <- combined_data$positions %>%
    mutate(X_tag = X_mm + 0.25, Y_tag = Y_mm + 0.25) %>%
    arrange(FOV) %>% 
    ggplot() +
    geom_rect(aes(xmin = X_mm, xmax = X_mm + 0.5, 
                  ymin = Y_mm, ymax = Y_mm + 0.5),
              fill = "white", color = "grey50") +
    geom_text(aes(x = X_tag, y = Y_tag, label = FOV), color = "black", size = 3) +
    scale_x_reverse() +
    scale_y_reverse() +
    theme_void() + 
    theme(panel.border = element_rect(fill = NA)) +
    labs(x = NULL, y = NULL) +
    coord_equal()
  return(p)
}


# Unique targets detected per FOV #####

plt_n_unique_Fov <- function(combined_data, df_backgroupd = df_target_type) {
  target_n <- combined_data$tx_sub %>% 
    group_by(fov, target, target_type) %>%
    summarise(n = n(), .groups = "drop")
  p <- target_n %>% 
    ggplot(aes(x = fov, y = target, fill = log10(n))) +
    geom_tile(
      data = expand_grid(df_target_type, distinct(target_n, fov)), 
      fill = "white"
    ) +
    geom_tile() +
    facet_grid(~ target_type, scales = "free") +
    coord_flip() +
    labs(x = "FOV", y = "Target") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_bw() +
    theme(
      strip.clip = "off", 
      strip.background = element_blank(), 
      axis.text.x = element_blank(), 
      panel.grid = element_blank()
    )
  return(p)
}


# Number of unique targets detected per cell per FOV #####

plt_n_unique_CellFov <- function(combined_data) {
  p <- combined_data$metadata_sub %>% 
    ggplot(aes(x = fov, y = nFeature_RNA)) +
    geom_boxplot(outlier.alpha = 1/5, outlier.size = 0.5) +
    labs(x = "FOV", y = "Count") + 
    theme_boxplot()
  return(p)
}


# Count and percent of transcript

.get_tx_n_percent_by_group_Fov <- function(tx_sub, group.by = NULL) {
  if (is.null(group.by)) {
    n_percent_by_group <- tx_sub %>% 
      group_by(fov) %>% 
      summarise(n = n(), .groups = "drop")
  } else {
    group.by_tag <- setNames(
      c("Target Type", "Assignment", "Cell Compoent"), 
      c("target_type", "assignment", "CellComp")
    )[group.by]
    n_percent_by_group <- tx_sub %>% 
      select(setNames(c("fov", group.by), 
                      c("fov", "group.by"))) %>% 
      group_by(fov, group.by) %>% 
      summarise(n = n(), .groups = "drop") %>% 
      group_by(fov) %>% 
      mutate(
        N = sum(n), 
        percent = n / N * 100, 
        group.by_tag = group.by_tag
      )
  }
  return(n_percent_by_group)
}
get_tx_n_percent_by_group_Fov <- function(combined_data, group.by = NULL) {
  combined_data$tx_sub %>% 
    .get_tx_n_percent_by_group_Fov(group.by = group.by)
}
get_tx_real_n_percent_by_group_Fov <- function(combined_data, group.by = NULL) {
  combined_data$tx_sub %>% 
    filter(target_type == "Gene") %>% 
    .get_tx_n_percent_by_group_Fov(group.by = group.by)
}
get_tx_negative_n_percent_by_group_Fov <- function(combined_data, group.by = NULL) {
  combined_data$tx_sub %>% 
    filter(target_type == "Negative") %>% 
    .get_tx_n_percent_by_group_Fov(group.by = group.by)
}
get_tx_false_n_percent_by_group_Fov <- function(combined_data, group.by = NULL) {
  combined_data$tx_sub %>% 
    filter(target_type == "SystemControl") %>% 
    .get_tx_n_percent_by_group_Fov(group.by = group.by)
}

plt_n_by_group_Fov <- function(n_percent_by_group) {
  if (is.null(n_percent_by_group$group.by)) {
    n_percent_by_group$group.by <- "1"
  }
  group.by_tag <- unique(n_percent_by_group$group.by_tag)
  p <- n_percent_by_group %>% 
    ggplot(aes(x = fov, y = n, fill = group.by)) +
    geom_col(width = 0.6) +
    labs(x = "FOV", y = "Count", fill = group.by_tag) +
    scale_fill_manual(values = COLOR) +
    theme_col_count()
  return(p)
}

plt_percent_by_group_Fov <- function(n_percent_by_group) {
  if (is.null(n_percent_by_group$group.by)) {
    n_percent_by_group$group.by <- "1"
  }
  group.by_tag <- unique(n_percent_by_group$group.by_tag)
  p <- n_percent_by_group %>% 
    ggplot(aes(x = fov, y = percent, fill = group.by)) +
    geom_col(width = 0.6, position = "stack") +
    labs(x = "FOV", y = "Percent (%)", fill = group.by_tag) +
    scale_fill_manual(values = COLOR) +
    theme_col_percent()
  return(p)
}


.get_tx_mean_sd_by_group_Target <- function(tx_sub, group.by = NULL) {
  if (is.null(group.by)) {
    mean_sd_by_group <- tx_sub %>% 
      group_by(fov, target) %>% 
      summarise(n = n(), .groups = "drop") %>% 
      group_by(fov) %>% 
      summarise(
        mean = mean(n), 
        sd = sd(n), 
        group.by_tag = group.by_tag
      )
  } else {
    group.by_tag <- setNames(
      c("Target Type", "Assignment", "Cell Compoent"), 
      c("target_type", "assignment", "CellComp")
    )[group.by]
    mean_sd_by_group <- tx_sub %>% 
      select(setNames(c("fov", "target", group.by), 
                      c("fov", "target", "group.by"))) %>% 
      group_by(fov, target, group.by) %>% 
      summarise(n = n(), .groups = "drop") %>% 
      group_by(fov, group.by) %>% 
      summarise(
        mean = mean(n), 
        sd = sd(n), 
        cv = sd / mean, 
        mean_sd = mean + sd, 
        group.by_tag = group.by_tag
      )
  }
  return(mean_sd_by_group)
}
get_tx_mean_sd_by_group_Target <- function(combined_data, group.by = NULL) {
  combined_data$tx_sub %>% 
    .get_tx_mean_sd_by_group_Target(group.by = group.by)
}

plt_cv_by_group_Target <- function(mean_sd_by_group) {
  group.by_tag <- unique(mean_sd_by_group$group.by_tag)
  p <- mean_sd_by_group %>% 
    ggplot(aes(x = fov, y = cv, fill = group.by)) +
    geom_col(width = 0.6) +
    labs(x = "FOV", y = "Coefficient of Variation", fill = group.by_tag) +
    scale_fill_manual(values = COLOR) +
    theme_col_count()
  return(p)
}

plt_mean_sd_by_group_Target <- function(mean_sd_by_group) {
  group.by_tag <- unique(mean_sd_by_group$group.by_tag)
  p <- mean_sd_by_group %>% 
    ggplot() +
    geom_errorbar(aes(x = fov, ymin = mean, ymax = mean + sd), width = .5) +
    geom_col(aes(x = fov, y = mean, fill = group.by, color = group.by), width = 0.6) +
    labs(x = "FOV", y = "Mean + SD", fill = group.by_tag) +
    scale_fill_manual(values = COLOR) +
    scale_color_manual(values = COLOR) +
    guides(color = "none") +
    theme_col_count()
  return(p)
}


plt_n_target_CellFov <- function(tx_sub = combined_data$tx_sub, all_cell) {
  target_n <- tx_sub %>%
    filter(cell_ID != 0) %>% 
    group_by(fov, cell) %>% 
    summarise(n = n(), .groups = "drop")
  target_n <- all_cell %>% 
    left_join(target_n, by = c("cell", "fov")) %>%
    mutate(n = ifelse(is.na(n), 0, n))
  p <- target_n %>% 
    ggplot(aes(x = fov, y = n)) +
    geom_boxplot(outlier.alpha = 1/5, outlier.size = 0.5) +
    labs(x = "FOV", y = "Count") + 
    theme_boxplot()
  return(p)
}


# Total number of cells (per FOV) #####
plt_n_cell_Fov <- function(combined_data = combined_data) {
  cell_counts_per_fov <- combined_data$metadata_sub %>%
    group_by(fov) %>%
    summarise(n = n())
  p <- cell_counts_per_fov %>% 
    ggplot(aes(x = fov, y = n)) +
    geom_col() +
    labs(x = "FOV", y = "Count") +
    theme_col_count() 
  return(p)
}

# Size of cells per cell per FOV #####
plt_size_cell_CellFov <- function(combined_data = combined_data) {
  p <- combined_data$metadata_sub %>% 
    ggplot(aes(x = fov, y = Area.um2)) +
    geom_boxplot(outlier.alpha = 1/5, outlier.size = 0.5) +
    # geom_hline(yintercept = (1:2 * 10 / 2)^2 * pi, linetype = 2, color = "blue") +
    labs(x = "FOV") +
    theme_boxplot()
  return(p)
}

# Correlation #####
plt_cor_xy <- function(combined_data = combined_data, x, y) {
  df_p <- combined_data$metadata_sub %>% 
    select(setNames(c(x, y), c("x", "y")))
  res_cor <- cor.test(df_p$x, df_p$y)
  p_title <- paste(
    sprintf("cor = %.4f", round(res_cor$estimate, 4)), 
    ifelse(res_cor$p.value < 0.0001, "p < 0.0001", sprintf("p = %.4f", round(res_cor$p.value, 4))), 
    sep = ", "
  )
  p <- df_p %>% 
    ggplot(aes(x = x, y = y)) +
    geom_point() +
    geom_smooth(method = "lm") +
    labs(x = x, y = y, title = p_title) +
    theme_bw()
  return(p)
}




# Transcripts True Gene (Non-Control) vs Negative Control #####
get_negative_adjust_Fov <- function(combined_data) {
  count_target <- combined_data$tx_sub %>%
    group_by(fov, target_type) %>% 
    summarise(n = n(), .groups = "drop")
  unique_target <- combined_data$tx_sub %>%
    group_by(fov, target_type) %>% 
    summarise(n_unique = n_distinct(target), .groups = "drop")
  negative_adjust <- count_target %>% 
    left_join(unique_target, by = c("fov", "target_type")) %>% 
    mutate(mean = n / n_unique) %>% 
    pivot_wider(id_cols = fov, names_from = target_type, values_from = mean) %>% 
    mutate(
      division = (Gene + 1) / (Negative + 1), 
      subtract = Gene - Negative, 
    )
  return(negative_adjust)
}

plt_negative_adjust_Fov <- function(negative_adjust, plot.by) {
  y_lab <- c("division" = "(True Gene + 1) / (Negative Control + 1)",
             "subtract" = "True Gene - Negative Control")[plot.by]
  p <- negative_adjust %>% 
    select(setNames(c("fov", plot.by), 
                    c("fov", "plot.by"))) %>% 
    ggplot(aes(x = fov, y = plot.by)) +
    geom_col(width = 0.6) +
    labs(x = "FOV", y = y_lab) +
    theme_col_count() 
  return(p)
}


get_negative_adjust_CellFov <- function(combined_data) {
  count_target <- combined_data$tx_sub %>%
    filter(cell_ID != 0) %>% 
    group_by(fov, target_type, cell) %>% 
    summarise(n = n(), .groups = "drop")
  unique_target <- combined_data$tx_sub %>%
    filter(cell_ID != 0) %>% 
    group_by(fov, target_type, cell) %>% 
    summarise(n_unique = n_distinct(target), .groups = "drop")
  negative_adjust <- count_target %>% 
    left_join(unique_target, by = c("fov", "target_type", "cell")) %>% 
    mutate(mean = n / n_unique) %>% 
    pivot_wider(id_cols = c(fov, cell), names_from = target_type, values_from = mean) %>% 
    mutate(
      across(c(Negative, Gene), ~ ifelse(is.na(.x), 0, .x)), 
      division = (Gene + 1) / (Negative + 1), 
      subtract = Gene - Negative, 
    )
  return(negative_adjust)
}

plt_negative_adjust_CellFov <- function(negative_adjust, plot.by) {
  y_lab <- c("division" = "(True Gene + 1) / (Negative Control + 1)",
             "subtract" = "True Gene - Negative Control")[plot.by]
  p <- negative_adjust %>% 
    select(setNames(c("fov", plot.by), 
                    c("fov", "plot.by"))) %>% 
    ggplot(aes(x = fov, y = plot.by)) +
    geom_boxplot(outlier.alpha = 1/5, outlier.size = 0.5) +
    labs(x = "FOV", y = y_lab) +
    theme_boxplot() 
  return(p)
}

## Percent of Transcripts above LOD
plt_percent_above_lod <- function(combine_data, lod) {
  df_p <- combine_data$tx_sub %>% 
    group_by(fov, target, target_type) %>% 
    summarise(n = n()) %>% 
    left_join(lod, by = "fov") %>% 
    group_by(fov, target_type) %>% 
    summarise(percent = mean(n >= cut) * 100) %>% 
    mutate(percent = ifelse(target_type == "Gene", percent, -percent))
  p <- df_p %>% 
    ggplot(aes(x = fov, y = percent, fill = target_type)) +
    geom_col(data = ~ .x %>% filter(target_type == "Gene"), width = 0.6) +
    geom_col(data = ~ .x %>% filter(target_type != "Gene"), width = 0.6, position = "dodge") +
    scale_fill_manual(values = COLOR) +
    scale_y_continuous(labels = ~ abs(.x), limits = c(-100, 100)) +
    labs(x = "FOV", y = "Count of Calls Above Negative Probe LOD", fill = "Target Type") +
    theme_col_count()
  return(p)
}


## Segmentation and Transcript Density ##### 

plt_seg_density_Fov <- function(combined_data, FOV) {
  p1 <- polygon <- combined_data$polygons_sub %>% 
    mutate(y_local_px = -y_local_px) %>%
    filter(as.numeric(as.character(fov)) == FOV) %>% 
    st_as_sf(coords = c("x_local_px", "y_local_px")) %>% 
    group_by(cell) %>% 
    summarise(geometry = st_combine(geometry)) %>% 
    st_cast("POLYGON") %>% 
    ggplot() +
    geom_sf(fill = NA, color = "black") +
    labs(title = "Segmentation") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
  p2 <- combined_data$tx_sub %>% 
    filter(as.numeric(as.character(fov)) == FOV, cell_ID != 0) %>% 
    mutate(y_local_px = -y_local_px) %>% 
    ggplot(aes(x = x_local_px, y = y_local_px)) +
    geom_bin2d(bins = 700, geom = "tile", aes(fill = ..count..)) +
    labs(title = "Assigned") +
    coord_equal() + 
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
  p3 <- combined_data$tx_sub %>% 
    filter(as.numeric(as.character(fov)) == FOV, cell_ID == 0) %>% 
    mutate(y_local_px = -y_local_px) %>% 
    ggplot(aes(x = x_local_px, y = y_local_px)) +
    geom_bin2d(bins = 700, geom = "tile", aes(fill = ..count..)) +
    labs(title = "Unassigned") +
    coord_equal() + 
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
  return(list(p1, p2, p3))
}


# Other #####


# function for Total number of transcripts per FOV
plt_TranscriptsPerFov <- function(combined_data = combined_data, ylim_range = NULL) {
  Transcripts_per_fov <- combined_data$tx_sub %>%
    group_by(fov) %>%
    summarise(cell_count = n())
  
  p <- ggplot(Transcripts_per_fov, aes(x = fov, y = cell_count)) +
    theme_minimal() +
    geom_col(width = 0.5 * nrow(Transcripts_per_fov) / nrow(Transcripts_per_fov)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  if (!is.null(ylim_range)) {
    p <- p + ylim(ylim_range)
  } else {
    p <- p + ylim(0, max(Transcripts_per_fov$cell_count) + 10000)
  }
  print(p)
}

# function for Total number of unassigned transcripts per FOV
plt_UnassignedPerFOV <- function(combined_data = combined_data, ylim_range = NULL) {
  UnassignTranscripts_per_fov <- combined_data$tx_sub %>%
    filter(cell_ID == 0) %>%
    group_by(fov) %>%
    filter(CellComp == "None") %>%
    summarise(CellComp = n())
  
  p <- ggplot(UnassignTranscripts_per_fov, aes(x = fov, y = CellComp)) +
    theme_minimal() +
    geom_col(width = 0.5 * nrow(UnassignTranscripts_per_fov) / nrow(UnassignTranscripts_per_fov)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(y = "Unassigned Transcripts")
  if (!is.null(ylim_range)) {
    p <- p + ylim(ylim_range)
  } else {
    p <- p + ylim(0, max(UnassignTranscripts_per_fov$CellComp) + 10000)
  }
  print(p)
}

# Total number of nuclear transcripts per FOV
plt_NuclearPerFOV <- function(combined_data = combined_data, ylim_range = NULL) {
  NuclearTranscripts_per_fov <- combined_data$tx_sub %>%
    filter(cell_ID != 0) %>%
    group_by(fov) %>%
    filter(CellComp == "Nuclear") %>%
    summarise(CellComp = n())
  
  p <- ggplot(NuclearTranscripts_per_fov, aes(x = fov, y = CellComp)) +
    theme_minimal() +
    geom_col(width = 0.5 * nrow(NuclearTranscripts_per_fov) / nrow(NuclearTranscripts_per_fov)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(y = "Nuclear Transcripts")
  if (!is.null(ylim_range)) {
    p <- p + ylim(ylim_range)
  } else {
    p <- p + ylim(0, max(NuclearTranscripts_per_fov$CellComp) + 10000)
  }
  print(p)
}

# Total number of membrane transcripts per FOV
plt_MembraneTranscriptsPerFOV <- function(combined_data = combined_data, ylim_range = NULL) {
  MembraneTranscripts_per_fov <- combined_data$tx_sub %>% 
    filter(cell_ID != 0) %>%
    group_by(fov) %>%
    filter(CellComp == "Membrane") %>%
    summarise(CellComp = n())
  
  p <- ggplot(MembraneTranscripts_per_fov, aes(x = fov, y = CellComp)) +
    theme_minimal() +
    geom_col(width = 0.5 * nrow(MembraneTranscripts_per_fov) / nrow(MembraneTranscripts_per_fov)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(y = "Membrane Transcripts")
  if (!is.null(ylim_range)) {
    p <- p + ylim(ylim_range)
  } else {
    p <- p + ylim(0, max(MembraneTranscripts_per_fov$CellComp) + 10000)
  }
  print(p)
}

# Total number of cytoplasm transcripts per FOV
plt_CytoplasmTranscriptsPerFov <- function(combined_data = combined_data, ylim_range = NULL) {
  CytoplasmTranscripts_per_fov <- combined_data$tx_sub %>%
    group_by(fov) %>%
    filter(CellComp == "Cytoplasm") %>%
    summarise(CellComp = n())
  
  p <- ggplot(CytoplasmTranscripts_per_fov, aes(x = fov, y = CellComp)) +
    theme_minimal() +
    geom_col(width = 0.5 * nrow(CytoplasmTranscripts_per_fov) / nrow(CytoplasmTranscripts_per_fov)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(y = "Cytoplasm Transcripts")
  if (!is.null(ylim_range)) {
    p <- p + ylim(ylim_range)
  } else {
    p <- p + ylim(0, max(CytoplasmTranscripts_per_fov$CellComp) + 10000)
  }
  print(p)
}

## calculate cell summary
Summarize_CellComp <- function( combined_data = combined_data){
  # Calculate the count of transcripts in each localization for each cell
  localization_counts <- combined_data$tx_sub %>%
    group_by(cell,cell_ID, fov, CellComp,X_mm, Y_mm) %>%
    summarize(count = n()) %>%
    ungroup()
  
  # Calculate the total number of transcripts per cell
  total_transcripts_per_cell <- combined_data$tx_sub %>%
    group_by(cell,cell_ID, fov,X_mm, Y_mm) %>%
    summarize(total_transcripts = n()) %>%
    ungroup()
  
  # Merge the total counts with the localization counts
  merged_data <- merge(localization_counts, total_transcripts_per_cell, 
                       by = c("cell","cell_ID", "fov","X_mm", "Y_mm"))
  
  # Calculate the proportion for each localization
  merged_data <- merged_data %>%
    mutate(proportion = round(count / total_transcripts * 100,2))
  
  # Pivot the data to have a wide format, showing both count and proportion for each localization
  cell_comp_summary <- merged_data %>%
    pivot_wider(names_from = CellComp, 
                values_from = c(count, proportion), 
                names_sep = "_", 
                values_fill = list(count = 0, proportion = 0))
  return(cell_comp_summary)
}


## Transcripts (percentage) per cell per FOV

plt_TranscriptsPerCellPerFOVPct <- function(cell_comp_summary = cell_comp_summary, ylim_ranges = list(NULL, NULL, NULL, NULL)) {
  cell_comp_summary <- cell_comp_summary
  # plot 
  p0 <- ggplot(cell_comp_summary, aes(x = fov, y = proportion_None )) +
    theme_minimal() +
    geom_boxplot(width = 0.5 * length(unique(cell_comp_summary$fov)) / length(unique(cell_comp_summary$fov))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(y = "Unassigned RNA") +
    ylim(ylim_ranges[[1]]) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  
  p1 <- ggplot(cell_comp_summary, aes(x = fov, y = proportion_Membrane)) +
    theme_minimal() +
    geom_boxplot(width = 0.5 * length(unique(cell_comp_summary$fov)) / length(unique(cell_comp_summary$fov))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(y = "Membrane RNA") +
    ylim(ylim_ranges[[2]]) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  
  p2 <- ggplot(cell_comp_summary, aes(x = fov, y = proportion_Nuclear)) +
    theme_minimal() +
    geom_boxplot(width = 0.5 * length(unique(cell_comp_summary$fov)) / length(unique(cell_comp_summary$fov))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(y = "Nuclear RNA") +
    ylim(ylim_ranges[[3]]) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  
  p3 <- ggplot(cell_comp_summary, aes(x = fov, y = proportion_Cytoplasm)) +
    theme_minimal() +
    geom_boxplot(width = 0.5 * length(unique(cell_comp_summary$fov)) / length(unique(cell_comp_summary$fov))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(y = "Cytoplasm RNA") +
    ylim(ylim_ranges[[4]])
  
  psum <- ggpubr::ggarrange(p0, p1, p2, p3, ncol = 1, heights = c(3, 3, 3, 5))
  print(psum)
}

## Transcripts (count) per cell per FOV
plt_TranscriptsPerCellPerFOVCt <- function(cell_comp_summary = cell_comp_summary,ylim_ranges = list(NULL, NULL, NULL, NULL)){
  cell_comp_summary <- cell_comp_summary
  
  # plot 
  p0 <- ggplot(cell_comp_summary, aes(x = fov,y = count_None ))+theme_minimal()+
    geom_boxplot(width=0.5*length(unique(cell_comp_summary$fov))/length(unique(cell_comp_summary$fov)))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ylim(ylim_ranges[[1]])+labs(y = "Unassigned RNA")+theme(axis.title.x=element_blank(),
                                                            axis.text.x=element_blank())
  
  p1 <- ggplot(cell_comp_summary, aes(x = fov,y = count_Membrane))+theme_minimal()+
    geom_boxplot(width=0.5*length(unique(cell_comp_summary$fov))/length(unique(cell_comp_summary$fov)))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ylim(ylim_ranges[[2]])+labs(y = "Membrane RNA")+theme(axis.title.x=element_blank(),
                                                          axis.text.x=element_blank())
  
  p2 <- ggplot(cell_comp_summary, aes(x = fov,y = count_Nuclear))+theme_minimal()+
    geom_boxplot(width=0.5*length(unique(cell_comp_summary$fov))/length(unique(cell_comp_summary$fov)))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ylim(ylim_ranges[[3]])+labs(y = "Nuclear RNA")+theme(axis.title.x=element_blank(),
                                                         axis.text.x=element_blank())
  p3 <- ggplot(cell_comp_summary, aes(x = fov,y = count_Cytoplasm))+theme_minimal()+
    geom_boxplot(width=0.5*length(unique(cell_comp_summary$fov))/length(unique(cell_comp_summary$fov)))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ylim(ylim_ranges[[3]])+labs(y = "Cytoplasm RNA")
  psum <- ggpubr::ggarrange(p0,p1,p2,p3,ncol = 1,heights = c(3,3, 3,5))
  print(psum)
}

# find max value of probe
find_max_Nprobe <- function(combined_data = combined_data){
  NegativeProbeName <- unique(grep("Nega",(combined_data$tx_sub$target),value = T))
  NegativeProbePerFOV <- combined_data$tx_sub %>%
    filter(target %in% NegativeProbeName) %>%
    group_by(fov) %>%
    summarize(NegativeProbeCounts = n()) %>%
    ungroup()
  max_negative_probes <- max(NegativeProbePerFOV$NegativeProbeCounts, na.rm = TRUE)
  return(max_negative_probes)
}


## Negative probe per FOV

plt_NegativeProbePerFov <- function(combined_data = combined_data, ylim_range = NULL){
  NegativeProbeName <- unique(grep("Nega",(combined_data$tx_sub$target),value = T))
  NegativeProbePerFOV <- combined_data$tx_sub %>%
    filter(target %in% NegativeProbeName) %>%
    group_by(fov) %>%
    summarize(NegativeProbeCounts = n()) %>%
    ungroup()
  
  p <- ggplot(NegativeProbePerFOV, aes(x = fov,y = NegativeProbeCounts)) +
    geom_col(width = 0.5*nrow(NegativeProbePerFOV)/nrow(NegativeProbePerFOV)) + 
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(y = "Negative Probe counts")  
  if (!is.null(ylim_range)) {
    p <- p + ylim(ylim_range)
  }
  print(p)
}

## 
## The number of negative probe in unassigned transcripts per FOV
find_nax_NprobeC.0 <- function(combined_data = combined_data){
  NegativeProbeName <- unique(grep("Nega", combined_data$tx_sub$target, value = TRUE))
  # Summarize the number of negative probes per FOV for unassigned cells
  NegativeProbePerCellPerFOV <- combined_data$tx_sub %>%
    filter(cell_ID == 0, target %in% NegativeProbeName) %>%
    group_by(fov) %>%
    summarize(NegativeProbeCounts = n()) %>%
    ungroup()  # Always good practice to ungroup after summarization
  max_negative_probes_per_fov <- max(NegativeProbePerCellPerFOV$NegativeProbeCounts, na.rm = TRUE)
  return(max_negative_probes_per_fov)
}
find_nax_NprobeC.n0 <- function(combined_data = combined_data){
  NegativeProbeName <- unique(grep("Nega", combined_data$tx_sub$target, value = TRUE))
  # Summarize the number of negative probes per FOV for unassigned cells
  NegativeProbePerCellPerFOV <- combined_data$tx_sub %>%
    filter(cell_ID != 0, target %in% NegativeProbeName) %>%
    group_by(cell,fov) %>%
    summarize(NegativeProbeCounts = n()) %>%
    ungroup()  # Always good practice to ungroup after summarization
  max_negative_probes_per_fov <- max(NegativeProbePerCellPerFOV$NegativeProbeCounts, na.rm = TRUE)
  return(max_negative_probes_per_fov)
}

plt_UnassignedNegativeProbePerFov <- function(combined_data, ylim_range = NULL) {
  NegativeProbeName <- unique(grep("Nega", (combined_data$tx_sub$target), value = TRUE))
  NegativeProbePerCellPerFOV <- combined_data$tx_sub %>%
    filter(cell_ID == 0) %>%
    filter(target %in% NegativeProbeName) %>%
    group_by(fov) %>%
    summarize(NegativeProbeCounts = n())
  
  p <- ggplot(NegativeProbePerCellPerFOV, aes(x = fov, y = NegativeProbeCounts)) + theme_minimal() +
    geom_col() + labs(y = "Negative Probe counts") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  if (!is.null(ylim_range)) {
    p <- p + ylim(ylim_range)
  }
  print(p)
}



# Negative probe per cell per FOV 

plt_NegativeProbePerCellPerFov <- function(combined_data, ylim_range = NULL) {
  NegativeProbeName <- unique(grep("Nega", (combined_data$tx_sub$target), value = TRUE))
  
  df <- combined_data$tx_sub %>%
    mutate(NegativeProbe_expressed = ifelse(target %in% NegativeProbeName, 1, 0))
  
  expression_counts <- df %>%
    filter(cell_ID != 0) %>%
    group_by(cell,fov) %>%
    summarize(Negative_expressed = sum(NegativeProbe_expressed))
  
  p <- ggplot(expression_counts, aes(x = fov, y = Negative_expressed)) + theme_minimal() +
    geom_boxplot() + labs(y = "Negative Probe counts") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  if (!is.null(ylim_range)) {
    p <- p + ylim(ylim_range)
  }
  print(p)
}

# Function to calculate maximum density after using stat_bin2d
calculate_max_density <- function(data) {
  # Using stat_bin2d to calculate densities
  density_data <- ggplot(data, aes(x = x_global_px, y = y_global_px)) +
    stat_bin2d(bins = 700, geom = "tile", aes(fill = ..count..))
  # Use ggplot_build to extract the computed data
  built_plot <- ggplot_build(density_data)
  
  # Extract the calculated count values which represent the densities
  max_density <- max(built_plot$data[[1]]$count, na.rm = TRUE)
  return(max_density)
}

# Overall transcripts 2D density plot for one slide of whole flowcell
plt_DensityTranscripts <- function(combined_data, max_density) {
  p <- ggplot(combined_data$tx_sub, aes(x = x_global_px, y = y_global_px)) +
    geom_bin2d(bins = 700) +
    scale_fill_continuous(type = "viridis", limits = c(0, max_density)) +
    scale_y_reverse() +
    theme_bw() +
    # theme(legend.position = 'none') +
    labs(x = "X Global Position (px)", y = "Y Global Position (px)", 
         title = "Density of All Transcripts in 2D Space")
  print(p)
}

# Assigned transcripts 2D density plot for one slide of whole flowcell
plt_DensityAssignedTranscripts <- function(combined_data, max_density) {
  assigned_tx_sub <- combined_data$tx_sub %>% filter(cell_ID != 0)
  
  p <- ggplot(assigned_tx_sub, aes(x = x_global_px, y = y_global_px)) +
    geom_bin2d(bins = 700) +
    scale_fill_continuous(type = "viridis", limits = c(0, max_density)) +
    scale_y_reverse() +
    theme_bw() +
    # theme(legend.position = 'none') +
    labs(x = "X Global Position (px)", y = "Y Global Position (px)", 
         title = "Density of Assigned Transcripts in 2D Space")
  print(p)
}

# Unassigned transcripts 2D density plot for one slide of whole flowcell
plt_DensityUnassignedTranscripts <- function(combined_data, max_density) {
  unassigned_tx_sub <- combined_data$tx_sub %>% filter(cell_ID == 0)
  
  p <- ggplot(unassigned_tx_sub, aes(x = x_global_px, y = y_global_px)) +
    geom_bin2d(bins = 500) +
    scale_fill_continuous(type = "viridis", limits = c(0, max_density)) +
    theme_bw() +
    scale_y_reverse() +
    # theme(legend.position = 'none') +
    labs(x = "X Global Position (px)", y = "Y Global Position (px)", 
         title = "Density of Unassigned Transcripts in 2D Space")
  print(p)
}
