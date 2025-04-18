#' Create a tree visualization with a variant heatmap
#'
#' This function generates a publication-quality tree visualization that includes
#' a heatmap showing variant states (REF/ALT) for each taxon. The visualization
#' can be customized with various parameters.
#'
#' @param tree A phylo object to visualize
#' @param data Dataframe with variant data to display as a heatmap
#' @param output_file Optional file path to save the plot
#' @return A ggtree plot object
#' @export
#' @importFrom ggtree ggtree geom_tiplab gheatmap theme_tree2
#' @importFrom ggplot2 labs scale_fill_manual 
create_variant_heatmap_tree <- function(tree, data, output_file,...) {
  # Create base tree
  p <- ggtree(tree, ladderize = FALSE) +
    geom_tiplab(size = 2.5, offset = 0.002, align = TRUE,linesize = 0.1) +
    ggtree::vexpand(.1, 1) +
    theme_tree2()
  
  # Define color palette for REF/ALT states
  ref_alt_colors <- c("REF" = "tomato", "ALT" = "royalblue")
  
  # Add heatmap alongside the tree
  # Each column represents a different variant or ancestry information
  p <- gheatmap(p, data, offset = 0.03, width = 0.12, 
                colnames = TRUE, colnames_angle = 45, 
                colnames_offset_y = 0.6, hjust = 0,
                colnames_position = "top") +
    scale_fill_manual(values = ref_alt_colors, name = "Allele") 
  
  # Return the plot object
  return(p)
}
