library(parallel)

# Script to perform enrichment analysis

enrichment_fisher_test <- function(cluster_vec, annotation_vec, tar_cluster, tar_annotation) {
  # Function to perform pathway enrichment analysis. 
  #
  # Args:
  #   cluster_vec: a vector containing which cluster a list of genes belong to
  #   annotation_vec: a vector containing which functions a list of genes are annotated with
  #   tar_cluster: the targeted cluster of genes  
  #   tar_annotation: indicate the functions to see if the targeted cluster of genes are overrepreseneted in it
  #   using fisher's exact test.
  #
  # Return:
  #  the result from fisher's exact test. 
  
  is_in_cluster <- cluster_vec == tar_cluster
  is_in_pathway <- annotation_vec == tar_annotation
  cross_table <- table(is_in_cluster, is_in_pathway)
  if(dim(cross_table)[1] == dim(cross_table)[2]) {
    ft_res <- fisher.test(cross_table, alternative = "greater")$p.value
  } else {
    ft_res <- 1
  }
  return(ft_res)
}

enrichment_fisher_test_mass_apply <- function(cluster_vec, annotations_vec){
  enrichment_fisher_test <- function(cluster_vec, annotation_vec, tar_cluster, tar_annotation) {
    # Function to perform pathway enrichment analysis. 
    #
    # Args:
    #   cluster_vec: a vector containing which cluster a list of genes belong to
    #   annotation_vec: a vector containing which functions a list of genes are annotated with
    #   tar_cluster: the targeted cluster of genes  
    #   tar_annotation: indicate the functions to see if the targeted cluster of genes are overrepreseneted in it
    #   using fisher's exact test.
    #
    # Return:
    #  the result from fisher's exact test. 
    
    is_in_cluster <- cluster_vec == tar_cluster
    is_in_pathway <- annotation_vec == tar_annotation
    cross_table <- table(is_in_cluster, is_in_pathway)
    if(dim(cross_table)[1] == dim(cross_table)[2]) {
      ft_res <- fisher.test(cross_table, alternative = "greater")$p.value
    } else {
      ft_res <- 1
    }
    return(ft_res)
  }
  
  ft_res <- NULL
  for (cluster in unique(cluster_vec)) {
    cl <- makeCluster(4)
    row_segs <- parLapply(cl, unique(annotations_vec), function(annotation){
      enrichment_fisher_test(cluster_vec, annotations_vec, cluster, annotation)
    })
    row <- do.call(rbind, row_segs)
    ft_res <- cbind(ft_res, unlist(row))
    stopCluster(cl) 
  }
  row.names(ft_res) <- unique(annotations_vec)
  colnames(ft_res) <- paste0("cluster_", unique(cluster_vec))
  return(as.data.frame(ft_res))
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

color_range_match <- function(vals, range_vec, color_gradients){
  col_vec <- sapply(vals , function(x) {
    pos <- which(range_vec < x)[length(which(range_vec < x))]
    color_gradients[pos]
  })
  return(col_vec)
}

# MAIN:

g_membership <- g_community$membership
universe_labels <- kegg_annotation_nodes$pw_des

# map the cluster to the equival length to the pathway list:

mapped_membership <- map_a_to_b(a = g_membership, 
                                ind_a = vertex_names, 
                                ind_b = kegg_annotation_nodes$protein_ids)

# remove the NAs in the list
mapped_membership <- mapped_membership[!is.na(universe_labels)]
universe_labels <- universe_labels[!is.na(universe_labels)]
table(universe_labels)[order(table(universe_labels), decreasing = T)]

# carry out enrichment analysis
ft_res <- enrichment_fisher_test_mass_apply(mapped_membership, universe_labels)
p_adj_res <- ft_res
p_adj_res[,] <- p.adjust(unlist(p_adj_res), method = "BH")

# the results:
positions <- which(p_adj_res < 0.05, arr.ind = T)
cluster_results <- cbind(row.names(p_adj_res)[positions[, 1]], colnames(p_adj_res)[positions[, 2]])

# remove meaningless annotations

cluster_results <- cluster_results[names(cluster_results) != "Metabolic pathways"] 

# plot:
#  generate the colors
color_size <- length(unique(g_membership))
col_for_cluster <- color_generator(color_size)


# shapes of the node:
shape_pos <- vertex_names %in% links$tar_proteins
v_shapes <- c("circle",  "square")[shape_pos+1] # target-square; potential-circle.

# visualization 
plot(g, vertex.shape = v_shapes, 
     vertex.color=col_for_cluster[g_membership],
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

# ================== GO =========================

g_membership <- g_community$membership
universe_labels <- go_annotation_nodes$go_terms

# map the cluster to the equival length to the pathway list:

mapped_membership <- map_a_to_b(a = g_membership, 
                                ind_a = vertex_names, 
                                ind_b = go_annotation_nodes$protein_ids)

# remove the NAs in the list
mapped_membership <- mapped_membership[!is.na(universe_labels)]
universe_labels <- universe_labels[!is.na(universe_labels)]
table(universe_labels)[order(table(universe_labels), decreasing = T)]


# carry out enrichment analysis
start_t = Sys.time()
ft_res <- enrichment_fisher_test_mass_apply(mapped_membership, universe_labels)
end_t = Sys.time()
print(end_t - start_t)

# the results:


# adjust the p-values to the 
p_adj_res <- ft_res
p_adj_res[,] <- p.adjust(unlist(p_adj_res), method = "BH")

# the results:
cordi <- which(p_adj_res < 0.05, arr.ind = T)[, 2]
results <- colnames(p_adj_res)[cordi]
names(results) <- names(cordi)
descriptions <- go_annotation_nodes$go_descrips[match(names(results),  go_annotation_nodes$go_terms)]
View(cbind(descriptions, as.character(results), names(results)))
# go enrich visualiztion:
p_for_cluster_1 <- p_adj_res[,2]
top_20_p_vals <- p_for_cluster_1[order(p_for_cluster_1, decreasing = F)][1:50]
top_20_sig_gos <- row.names(p_adj_res)[order(p_for_cluster_1, decreasing = F)][1:50]

# query the hierachies: 
go_URIs <- paste0("oboGo:", gsub(":", "_" ,top_20_sig_gos))
go_hierachy_visualiztion <- data_require(GO_hierarchy_retrival_query_make, list(go_URIs), endpoint)

# make links to visualize the results in hierarchy

links_go_enrich <- cbind.data.frame(from = paste(go_hierachy_visualiztion[, 1], go_hierachy_visualiztion[, 2], sep = "\n" ), 
                 to = paste(go_hierachy_visualiztion[, 3], go_hierachy_visualiztion[, 4], sep = "\n" ))

g_go_hierachy <- graph_from_data_frame(d = links_go_enrich, directed = T)
V(g_go_hierachy)$label.cex <- 0.8
colors <- colorRampPalette(c("orange", "white"))(30)
ranges <- seq(0, 1, length.out = length(colors))
col_vec <- color_range_match(top_20_p_vals, ranges, colors)

plot(minimum.spanning.tree(g_go_hierachy),
     vertex.color = col_vec 
     # , vertex.size = 22
     )
