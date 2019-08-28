library(networkD3)
library(visNetwork)
library(htmlwidgets)


tubular_paste <- function(lists_of_tubulars, separation = "; ") {
  paste_df <- as.data.frame(lists_of_tubulars)
  results <- do.call(paste, c(paste_df, sep= separation))
  return(results)
}



# if show the descriptions

annotation_list <- tapply(protein_descriptions_nodes$description, INDEX = protein_descriptions_nodes$protein_ids, unique)
rearranged_anno_list <- annotation_list[match(vertex_names, names(annotation_list))]
annotation_list_ready_to_paste <- lapply(rearranged_anno_list, function(vec) {
  vec_for_further_removal <- vec[!grepl("Generated via", x = vec)]
  res <- gsub(pattern = "\\[.*\\]", "",  vec_for_further_removal)
  res <- gsub(pattern = " $", "",  res)
  unique(res)
})
annotation <- paste(names(annotation_list_ready_to_paste), annotation_list_ready_to_paste, sep = ": ")

# if wanna see the annotation in details
# table(unip_annotation_nodes$types)
# uniprot_anno_category <- "Function_Annotation"
# functional_anno <- unip_annotation_nodes[unip_annotation_nodes$types == uniprot_anno_category, ]
# functional_anno <- map_a_to_b(a = functional_anno, ind_a = functional_anno$protein_ids, ind_b = vertex_names)
# annotation <- tubular_paste(list(vertex_names, functional_anno$comments))


betw_centralities <- map_a_to_b(betw_centralities, names(betw_centralities), vertex_names)
btwness <- as.numeric(betw_centralities)*100
btwness[btwness < 1] <- 0.1
nodes <- cbind.data.frame(vertex_names,  group = g_membership, annotation = annotation, btwness)


# make the links
links <- interactions[, c(1, 2, 4)][interactions$weights > cutoff, ]
links$pot_proteins <- match(links$pot_proteins, vertex_names) - 1
links$tar_proteins <- match(links$tar_proteins, vertex_names) - 1

# if highlight the betweenness
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

View(unip_annotation_nodes[unip_annotation_nodes$protein_ids == "AT5G63420.1", ])

sankeyNetwork(Links = links,
              Nodes = nodes, 
              Source = "pot_proteins", 
              Target = "tar_proteins", 
              Value = "weights", 
              NodeID = "vertex_names",
              NodeGroup = "annotation", 
              fontSize = 10)

saveNetwork(network = test, file = 'test.html')

# onRender function enables d3.js commands to change the aesthetics: 
# onRender(test, ' function(el,x) {
#          d3.selectAll(".node text").style("font-size", "10px").attr("visibility", "hidden").on("mouseover", function(d){console.log(d.name); })
#          } ' )



# d3.selectAll(".node text").remove()
# d3.selectAll(".node").append("foreignObject").attr("width", 200).attr("height", 50).attr("visibility", "hidden")
# .style("font-size", "5px").html(function(d) {
#   return d.name; }).on("mouseover", function(d) {})