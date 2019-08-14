library(networkD3)

unip_annotation_nodes
vertex_names
length(vertex_names)

table(unip_annotation_nodes$types)[order(table(unip_annotation_nodes$types))]
functional_anno <- unip_annotation_nodes[unip_annotation_nodes$types == "Function_Annotation", ]

annotation <- tapply(kegg_annotation_nodes$pw_des, kegg_annotation_nodes$protein_ids, function(x){
  pathes <- c()
  pathes <- c(pathes, x)
  paste(pathes, collapse = "; ")
})
btwness <- as.numeric(betw_centralities)*100
btwness[btwness < 1] <- 0.1
functional_anno <- map_a_to_b(a = functional_anno, ind_a = functional_anno$protein_ids, ind_b = vertex_names)
nodes <- cbind.data.frame(vertex_names,  group = g_membership, annotation = functional_anno$comments, btwness)


nodes$annotation <- sapply(strsplit(functional_anno$comments, "\\."), function(x) {x[1]})
# or 
nodes$annotation <- annotation 

# make the links
links <- interactions[, c(1, 2, 4)][interactions$weights > cutoff, ]
links$pot_proteins <- match(links$pot_proteins, vertex_names) - 1
links$tar_proteins <- match(links$tar_proteins, vertex_names) - 1
links$weights <- links$weights / 1000


calls <- list(
  # data frame inputs
  Links = links,  # linkage dataframe
  Nodes = nodes,  # the characters of the nodes
  # explanatory inputs:
  Source = "pot_proteins",  # source of the linkages
  Target = "tar_proteins",  #  targets of the linkages
  Value = "weights",  # thickness of the edges in the graph
  NodeID = "annotation",  # text shown in the nodes
  Group = "group",  # groups of the nodes
  Nodesize = "btwness" ,
  # aesthetics inputs:
  fontFamily="Arial", 
  fontSize = 10, 
  linkColour="black",  
  # colourScale, 
  # linkWidth, 
  charge = -13, # how strong the nodes should gather or repel to each other  
  opacity = 0.9,
  legend=T,
  arrows=T,
  bounded=F,  # whether turn on the boundaries for the plots
  opacityNoHover=200, # the degree of opacity when the mouse is not suspending on the nodes 
  zoom = T  # allow zoom(double click to )
  )

# change the default setting 
# calls["NodeID"] <- "vertex_names"

do.call(what = forceNetwork, args = calls)






