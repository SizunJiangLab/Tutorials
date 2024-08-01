#' A function for clustering and plot the marker heatmap 
#' 
#' This function clusters your single cell protein data using Phenograph
#' 
#' @param df A dataframe containing the markers and cell to be clustered
#' @param id_col A string specifying the column suitable to be the identifier of unique cells
#' @param  markers A vector of strings specifying the names of the columns that will be used for clustering
#' @param k_num A scalar specifying the k for Rphenoannoy's K-nearest neighborhood
#' @param z Boolean, if z == TRUE, the function will draw a Z score based heatmap. If z == FALSE, the function will draw a mean expression based heatmap
#' @param seed A scalar specifying the random seed


# create a z score matrix with rows being markers and columns being clusters

z_score_matrix <- function(x, cluster_col, markers){
  
  heatmap_df <- x %>% 
    dplyr::select(all_of(cluster_col), all_of(markers))
  
  # Calculate the sample mean and sample standard deviation
  
  pop_stat <- heatmap_df %>%
    # Pivot longer to get all the markers into the same column for later grouping
    pivot_longer(cols = -!!as.symbol(cluster_col), names_to = 'marker', values_to = 'value') %>% 
    # Group by marker
    group_by(marker) %>%
    # Calculate the mean and sd 
    summarise(pop_mean = mean(value),
              pop_sd = mean(value)) %>% 
    ungroup()
  
  mean_df <- heatmap_df %>% 
    pivot_longer(cols = -!!as.symbol(cluster_col), names_to = 'marker', values_to = 'value') %>% 
    group_by(marker, !!as.symbol(cluster_col)) %>% 
    summarise(mu = mean(value)) %>% 
    ungroup()
  
  match_idx <- match(mean_df$marker, pop_stat$marker)
  
  stratified_stat <- mean_df %>%
    # calculate z score
    mutate(z_score = (mu - pop_stat$pop_mean[match_idx])/pop_stat$pop_sd[match_idx])
  
  heatmap_df <- stratified_stat %>%
    dplyr::select(marker, z_score, all_of(cluster_col)) %>%
    pivot_wider(id_cols = 'marker', names_from = all_of(cluster_col), values_from = 'z_score')
  
  heatmap_mat <- heatmap_df %>%
    column_to_rownames(var = 'marker') %>%
    as.matrix()
  
  return(heatmap_mat)
  
}

# create a mean value matrix with rows being markers and columns being clusters

mean_matrix <- function(x, cluster_col, markers){
  
  heatmap_df <- x %>% 
    dplyr::select(all_of(cluster_col), all_of(markers))
  
  mean_df <- heatmap_df %>% 
    pivot_longer(cols = -!!as.symbol(cluster_col), names_to = 'marker', values_to = 'value') %>% 
    group_by(marker, !!as.symbol(cluster_col)) %>% 
    summarise(mu = mean(value))
  
  heatmap_df <- mean_df %>% 
    dplyr::select(marker, mu, all_of(cluster_col)) %>% 
    pivot_wider(id_cols = 'marker', names_from = all_of(cluster_col), values_from = 'mu')
  
  heatmap_mat <- heatmap_df %>%
    column_to_rownames(var = 'marker') %>%
    as.matrix()
  
  return(heatmap_mat)
  
}

# draw marker enrichment heatmap

marker_cluster_heatmap <- function(m, count_table, z){
  if (z) {
    library(circlize)
    col_fun = colorRamp2(c(-1, 0, 1), c('#4575b4', "white", '#d73027'))
    
    bar_vec <- count_table
    
    column_ha <- HeatmapAnnotation(count = anno_barplot(bar_vec$n))
    
    h <- Heatmap(m, cluster_rows = FALSE, cluster_columns  = FALSE, col = col_fun, name = 'Z Score', clustering_method_columns = 'complete',
                 row_names_side = 'left', top_annotation = column_ha, column_names_gp = grid::gpar(fontsize = 12),
                 row_names_gp = grid::gpar(fontsize = 12), rect_gp = gpar(col = "black", lwd = 2))
  } else {
    library(circlize)
    col_fun = colorRamp2(c(0, 0.5), c("white", '#d73027'))
    
    bar_vec <- count_table
    
    column_ha <- HeatmapAnnotation(count = anno_barplot(bar_vec$n))
    
    h <- Heatmap(m, cluster_rows = FALSE, cluster_columns  = FALSE, col = col_fun, name = 'Mean Expression', clustering_method_columns = 'complete',
                 row_names_side = 'left', top_annotation = column_ha, column_names_gp = grid::gpar(fontsize = 12),
                 row_names_gp = grid::gpar(fontsize = 12), rect_gp = gpar(col = "black", lwd = 2))
  }
  
  return(h)
}

to_cluster <- function(df, id_col, markers, k_num, seed, z = FALSE){
  
  # set seed
  set.seed(seed)
  
  # Subset the dataframe
  cluster_df <- df %>% 
    dplyr::select(id_col, markers)
  
  # Extract only the marker columns and convert dataframe into a matrix (Rphenoannoy takes matrix instead of dataframe)
  cluster_mat <- as.matrix(cluster_df[,2:ncol(cluster_df)])
  
  # Cluster
  pg <- Rphenoannoy::Rphenoannoy(cluster_mat, k_num)
  
  # Extract the clusters
  pg_cluster <- as.factor(membership(pg[[2]]))
  
  # Append the clusters to the dataframe with the id column
  df$pg_cluster <- pg_cluster
  
  # Count cluster
  cluster_count <- df %>% 
    group_by(pg_cluster) %>% 
    dplyr::count()
  
  if (z){
    mat <- z_score_matrix(df , 'pg_cluster', markers)
    h <- marker_cluster_heatmap(mat, cluster_count, z)
  } else {
    mat <- mean_matrix(df, 'pg_cluster', markers)
    h <- marker_cluster_heatmap(mat, cluster_count, z)
  }
  
  draw(h)
  return(df)
  
}





