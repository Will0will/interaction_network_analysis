# Script to perform enrichment analysis

pathway_enrichment_fisher_test <- function(cluster_vec, pw_label_vec, tar_cluster, tar_pathway) {
  # Function to perform pathway enrichment analysis. 
  #
  # Args:
  #   cluster_vec: a vector containing which cluster a list of genes belong to
  #   pw_label_vec: a vector containing which pathway a list of genes are annotated with
  #   tar_cluster: the targeted cluster of genes  
  #   tar_pathway: indicate the pathway to see if the targeted cluster of genes are overrepreseneted in it
  #   using fisher's exact test.
  #
  # Return:
  #  the result from fisher's exact test. 
  
  is_in_cluster <- cluster_vec == tar_cluster
  is_in_pathway <- pw_label_vec == tar_pathway
  cross_table <- table(is_in_cluster, is_in_pathway)
  if(unique(dim(cross_table)[1] == dim(cross_table)[2])) {
    ft_res <- fisher.test(cross_table, alternative = "greater")$p.value
  } else {
    ft_res <- 1
  }
  return(ft_res)
}

map_a_to_b <- function(a, ind_a, ind_b){
  if (is.data.frame(a) | is.matrix(a)){
    result <- a[match(ind_b, ind_a), ]
  } else {
    result <- a[match(ind_b, ind_a)]
  }
  return(result)
}



# USAGE:
membership <- g_community$membership
universe_labels <- kegg_annotation_nodes$pw_des
mapped_membership <- map_a_to_b(a = membership, 
                                ind_a = vertex_names, 
                                ind_b = kegg_annotation_nodes$protein_ids)

mapped_membership <- mapped_membership[!is.na(universe_labels)]
universe_labels <- universe_labels[!is.na(universe_labels)]

pathways_for_test <- names(table(universe_labels[mapped_membership == 1]))

ft_res <- NULL
for (pw in pathways_for_test) {
  p_vals <- NULL
  for (cluster in 1:max(mapped_membership)) {
    p <- pathway_enrichment_fisher_test(mapped_membership, universe_labels, cluster, pw)
    p_vals <- c(p_vals, p)
  }
  print(length(p_vals))
  ft_res <- rbind(ft_res, p_vals)
}


row.names(ft_res) <- pathways_for_test
colnames(ft_res) <- paste0("cluster_", 1:max(mapped_membership))
colnames(ft_res)[is.null(ft_res)] <- 1
coodinates <- which(ft_res < 0.07, arr.ind = T)
coodinates[, 2]



