library(networkD3)
library(htmlwidgets)


tubular_paste <- function(lists_of_tubulars, separation = "; ") {
  paste_df <- as.data.frame(lists_of_tubulars)
  results <- do.call(paste, c(paste_df, sep= separation))
  return(results)
}


unip_annotation_nodes
vertex_names
length(vertex_names)

table(unip_annotation_nodes$types)[order(table(unip_annotation_nodes$types))]

# make the nodes

functional_anno <- unip_annotation_nodes[unip_annotation_nodes$types == "Function_Annotation", ]
betw_centralities <- map_a_to_b(betw_centralities, names(betw_centralities), vertex_names)
btwness <- as.numeric(betw_centralities)*100
btwness[btwness < 1] <- 0.1
functional_anno <- map_a_to_b(a = functional_anno, ind_a = functional_anno$protein_ids, ind_b = vertex_names)
annotation <- tubular_paste(list(vertex_names, functional_anno$comments))
nodes <- cbind.data.frame(vertex_names,  group = g_membership, annotation = annotation, btwness)

# make the links
links <- interactions[, c(1, 2, 4)][interactions$weights > cutoff, ]
links$pot_proteins <- match(links$pot_proteins, vertex_names) - 1
links$tar_proteins <- match(links$tar_proteins, vertex_names) - 1
links$weights <- links$weights / 200

# set the calls for forceNetwork
force_net <- forceNetwork(
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
  # opacityNoHover=10, # the degree of opacity when the mouse is not suspending on the nodes 
  zoom = T  # allow zoom(double click to )
  )
force_net
sankeyNetwork(Links = links,
              Nodes = nodes, 
              Source = "pot_proteins", 
              Target = "tar_proteins", 
              Value = "weights", 
              NodeID = "vertex_names",
              NodeGroup = "annotation", 
              fontSize = 10)


# onRender function enables d3.js commands to change the aesthetics: 
onRender(test, ' function(el,x) {
         d3.selectAll(".node text").style("font-size", "10px").attr("visibility", "hidden").on("mouseover", function(d){console.log(d.name); })
         } ' )

saveNetwork(network = test, file = 'test.html')

# d3.selectAll(".node text").remove()
# d3.selectAll(".node").append("foreignObject").attr("width", 200).attr("height", 50).attr("visibility", "hidden")
# .style("font-size", "5px").html(function(d) {
#   return d.name; }).on("mouseover", function(d) {})