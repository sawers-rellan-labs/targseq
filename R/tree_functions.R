#' Reverse the tree plotting order by reversing edges
#'
#' @param tree A phylo object to be flipped
#' @return A phylo object with reversed edge order
#' @export
flip_tree <- function(tree) {
  # Create a copy of the input tree
  revtree <- tree
  
  # Get the current edges
  current_edges <- tree$edge
  
  # Reverse the edge order (this changes the plotting order)
  rev_edges <- current_edges[rev(1:nrow(current_edges)),]
  revtree$edge <- rev_edges
  
  # Also reverse the edge lengths to maintain correspondence
  revtree$edge.length <- rev(tree$edge.length)
  
  # Identify which edges connect to tips (nodes numbered <= n_tips)
  is_tip <- revtree$edge[,2] <= length(revtree$tip.label)
  
  # Get the ordered tips based on the reversed edges
  ordered_tips <- revtree$edge[is_tip, 2]
  
  # Check reverse order 
  print(ordered_tips)
  
  # Return the reversed tree
  return(revtree)
}

#' Rotate a tree so a specific tip appears at the top of the visualization
#'
#' @param tree A phylo object to be rotated
#' @param target_tip The name of the tip to be placed at the top
#' @return A rotated phylo object
#' @export
pivot_on <- function(tree, target_tip) {
  # Input validation: Check if the tree is a valid phylo object
  if (is.null(tree) || !inherits(tree, "phylo")) {
    stop("Input must be a valid phylogenetic tree (phylo object)")
  }
  
  # Input validation: Check if the target tip exists in the tree
  if (!target_tip %in% tree$tip.label) {
    stop("Target tip '", target_tip, "' not found in the tree.")
  }
  
  # Make a copy of the tree to modify
  rotated_tree <- tree
  
  # Get the total number of tips in the tree
  ntips <- Ntip(rotated_tree)
  
  # Find the node number corresponding to the target tip
  # In ape's phylo objects, tips are numbered 1:ntips
  target_tip_node <- which(rotated_tree$tip.label == target_tip)
  
  # Get the root node (typically ntips + 1 in ape trees)
  root_node <- ntips + 1
  
  # Try to find a path from the root to the target tip
  # This will be used to determine which nodes need rotation
  path_exists <- FALSE
  tryCatch({
    # Use ape's nodepath function to find the sequence of nodes from root to tip
    node_path <- nodepath(rotated_tree, from = root_node, to = target_tip_node)
    if (length(node_path) > 0) path_exists <- TRUE
  }, error = function(e) {
    # If path from root fails, we'll try a different approach
    path_exists <- FALSE
  })
  
  # If we couldn't find a path from the root, try using a reference tip
  if (!path_exists) {
    # Use the last tip as a reference point
    reference_node <- ntips
    
    # Ensure the reference node is different from the target
    if (reference_node == target_tip_node) {
      reference_node <- 1  # Use the first tip instead
    }
    
    # Try to find a path between reference and target tip
    tryCatch({
      node_path <- nodepath(rotated_tree, from = reference_node, to = target_tip_node)
      if (length(node_path) == 0) {
        stop("Cannot find a valid path between reference and target tip")
      }
    }, error = function(e) {
      stop("Failed to find a valid path in the tree: ", e$message)
    })
  }
  
  # We only need to consider rotating internal nodes on the path
  # Internal nodes have numbers greater than ntips
  internal_nodes_on_path <- node_path[node_path > ntips]
  
  # If no internal nodes to rotate, return the original tree
  if (length(internal_nodes_on_path) == 0) {
    return(rotated_tree)
  }
  
  # Iterate through each node and rotate if necessary
  # Process nodes in order from root towards tip for more predictable results
  nodes_to_check <- sort(internal_nodes_on_path)
  
  for (node_to_rotate in nodes_to_check) {
    # Find the direct children of this internal node
    children <- rotated_tree$edge[rotated_tree$edge[, 1] == node_to_rotate, 2]
    
    # Check if the node is bifurcating (has exactly two children)
    # Most phylogenetic trees are bifurcating
    if (length(children) != 2) {
      warning("Node ", node_to_rotate, " is not bifurcating. Rotation might not work as expected.")
      next # Skip rotation for non-bifurcating nodes
    }
    
    # Find which child is on the path towards the target tip
    current_node_index_in_path <- which(node_path == node_to_rotate)
    
    # Make sure we're not at the end of the path
    if (current_node_index_in_path < length(node_path)) {
      child_on_path <- node_path[current_node_index_in_path + 1]
      
      # If the child on the path is the first child, 
      # rotate the node to make it the right-hand child
      # This is because trees are typically plotted with right-side branches on top
      if (child_on_path == children[1]) {
        rotated_tree <- ape::rotate(rotated_tree, node = node_to_rotate)
      }
    }
  }
  
  # Ensure consistent ordering of tips using rotateConstr
  # This ensures tips are arranged in the desired order
  constrained_tree <- rotateConstr(rotated_tree, rotated_tree$tip.label)
  
  # Finally, flip the tree to get the desired orientation
  # This ensures the target tip appears at the top of the visualization
  return(flip_tree(constrained_tree))
}