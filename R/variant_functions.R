#' Extract variant information from an alignment
#'
#' Extract variant information from an alignment
#'
#' This function locates a specific variant in a sequence alignment and extracts
#' the nucleotide at that position for each sequence. It uses a pattern matching
#' approach to identify the exact position of the variant.
#'
#' @param alignment A DNAStringSet object containing the aligned sequences
#' @param variant_name Character string with the name of the variant (e.g., "I211V")
#' @param match_pattern Character string with the nucleotide pattern to locate the variant
#' @param reference_seq_index Integer index of the reference sequence (default: 1)
#' @return A factor with the nucleotide at the variant position for each sequence
#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom Biostrings vmatchPattern
extract_variant_line <- function(alignment, variant_name, match_pattern, reference_seq_index = 1) {
  # Find the position of the variant using the match pattern in the reference sequence
  variant_match <- vmatchPattern(match_pattern, alignment)[[reference_seq_index]]
  
  if (length(variant_match) == 0) {
    stop("Pattern '", match_pattern, "' not found in reference sequence")
  }
  
  # Get the starting position of the match
  variant_start <- start(variant_match)
  
  # Log information about the variant
  cat(variant_name, "variant found at position:", variant_start, "\n")
  
  # Number of sequences in the alignment
  n_taxa <- length(alignment)
  
  # Create genomic ranges to extract just the variant position from each sequence
  variant_range <- GRanges(
    seqnames = names(alignment),
    ranges = IRanges(
      start = rep(variant_start, n_taxa),
      end = rep(variant_start, n_taxa)
    )
  )
  
  # Extract the nucleotide at the variant position for each sequence
  variant_aln <- alignment[variant_range]
  variant_matrix <- as.matrix(variant_aln)
  
  # Replace gap characters with NA
  variant_matrix[variant_matrix == "-"] <- NA
  
  # Return the variant data as a factor
  return(factor(as.vector(variant_matrix)))
}

#' This function locates a specific variant in a sequence alignment and extracts
#' the nucleotide at that position for each sequence. It uses a pattern matching
#' approach to identify the exact position of the variant.
#' @param data A data.frame of variant names and  patterns to find
#' @param alignment A DNAStringSet object containing the aligned sequences
#' @param reference_seq_index Integer index of the reference sequence (default: 1)
#' @return A factor with the nucleotide at the variant position for each sequence
#' @export
get_variant_gt <- function(data, alignment, reference_seq_index=1) {
  variant_list <- split(data, factor(1:nrow(data)))
  
  result_list <- lapply(variant_list,
                        FUN = function(x) {
                          variant_result <- extract_variant_line(
                            alignment,
                            x$name[1],
                            x$pattern[1],
                            reference_seq_index)
                          return(variant_result)
                        })
  
  # Convert list to data frame
  result_df <- as.data.frame(result_list)
  
  # Use the variant names as column names
  colnames(result_df) <- data$name
  
  return(result_df)
}

