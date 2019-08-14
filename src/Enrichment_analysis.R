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

color_generator <- function(color_size) {
  R <- sample(1:225, size = color_size) / 225
  G <- sample(1:225, size = color_size) / 225
  B <- sample(1:225, size = color_size) / 225
  cols <- sapply(1:color_size, function(i) {
    rgb(red=R[i], green=G[i], blue=B[i]) 
  })
  return(cols)
}


# MAIN:

g_membership <- g_community$membership
universe_labels <- kegg_annotation_nodes$pw_des

# map the cluster to the equival length to the pathway list:

mapped_membership <- map_a_to_b(a = membership, 
                                ind_a = vertex_names, 
                                ind_b = kegg_annotation_nodes$protein_ids)

# remove the NAs in the list
mapped_membership <- mapped_membership[!is.na(universe_labels)]
universe_labels <- universe_labels[!is.na(universe_labels)]
table(universe_labels)[order(table(universe_labels))]
# list of patways for enrichment 
pathways_for_test <- names(table(universe_labels[mapped_membership == 1]))

# carry out enrichment analysis
ft_res <- NULL
for (pw in pathways_for_test) {
  p_vals <- NULL
  for (cluster in unique(mapped_membership)) {
    p <- pathway_enrichment_fisher_test(mapped_membership, universe_labels, cluster, pw)
    p_vals <- c(p_vals, p)
  }
  print(length(p_vals))
  ft_res <- rbind(ft_res, p_vals)
}

row.names(ft_res) <- pathways_for_test
colnames(ft_res) <- paste0("cluster_", 1:max(mapped_membership))

# the results:
cluster_results <- which(ft_res < 0.1, arr.ind = T)[, 2]
# remove meaningless annotations
cluster_results <- cluster_results[names(cluster_results) != "Metabolic pathways"] 

# plot:
#  generate the colors
color_size <- length(unique(membership))
col_for_cluster <- color_generator(color_size)


# shapes of the node:
shape_pos <- vertex_names %in% links$tar_proteins
v_shapes <- c("circle",  "square")[shape_pos+1] # target-square; potential-circle.

# visualization 
plot(g, vertex.shape = v_shapes, 
     vertex.color=col_for_cluster[membership],
     vertex.size = 7,
     vertex.label=NA)



legend("topleft", 
       legend= names(cluster_results), 
       col = col_for_cluster[cluster_results], 
       bty = "n", 
       pch=20 , 
       pt.cex = 1.5, 
       cex = 0.75, 
       text.col = col_for_cluster[cluster_results] , 
       horiz = FALSE, 
       inset = c(0.1, 0.1))

