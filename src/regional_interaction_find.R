# or the interactions between 2 regions 
calls <- list(1, 23450000, 23650000, 5, 10450000, 10650000)
interactions <- data_require(interaction_two_regions_query_make, calls, endpoint)
patterns <- c("AT[0-9]G[0-9|.]+", "AT[0-9]G[0-9|.]+", "[0-9]+")
interactions <- regex_match_map(data_frame = interactions, regex_vec = patterns)

g <- graph_from_data_frame(d = interactions, directed = F) 
mst_net <- minimum.spanning.tree(g)
vertex_names <- V(mst_net)$name
shape_pos <- vertex_names %in% interactions$proteins_1
v_shapes <- c("circle",  "square")[shape_pos+1] 

plot(mst_net, 
     vertex.size = 5, 
     vertex.label= NA,  
     vertex.shape = v_shapes)

betw_centralities <- betweenness(mst_net, v = V(mst_net), directed = F, weights = NULL, normalized = T) betw_centralities <- betw_centralities[order(betw_centralities, decreasing = T)]
betw_cutoff <- quantile(betw_centralities, probs = 0.8)
nodes_with_higer_betw <- names(betw_centralities[betw_centralities > betw_cutoff])

protein_URIs <- URI_paste("http://rdf.ebi.ac.uk/resource/ensembl.protein/", nodes_with_higer_betw)
unip_anno <- data_require(uniprot_annotation_query_make, list(protein_URIs), endpoint)

# data clean
unip_anno <- regex_match_map(data_frame = unip_anno, regex_vec =  c(".*", "[A-Z|a-z]+_[A-Z|a-z|_]+", ".*")) #remember to use try here


# 'bind' key word:
bind_results <- unip_anno[grep(pattern = "bind", ignore.case = T, x = unip_anno$comments), ]
binding_pro_ids <- unique(bind_results$protein_ids)
bind_results$comments

# 'cataly.*' key words
cataly_results <- unip_anno[grep(pattern = "cataly", ignore.case = T, x = unip_anno$comments), ]
cataly_pro_ids <- unique(cataly_results$protein_ids)

# 'regulat.*' key words
regulat_results <- unip_anno[grep(pattern = "regulat", ignore.case = T, x = unip_anno$comments), ]
regulat_pro_ids <- unique(regulat_results$protein_ids)



# plot the results  
vernames_to_show <- vertex_names[match(vertex_names, regulat_pro_ids)]
plot(g, 
     vertex.size = 5, 
    # vertex.color=my_color, 
     vertex.label= vernames_to_show,  
     vertex.shape = v_shapes)

