# Scripts for finding the potential regulations via QTL analysis using semantic web

# ENV:

#  Default environment './example_data'
rm(list = ls())
default_env_dir <- file.path(dirname(getwd()), "example_data")
setwd(default_env_dir)
# setwd("~/thesis/interaction_network_analyssis/example_data")

# REQUIRED PACKAGES :

library(SPARQL)
library(igraph)
library(RColorBrewer)

# FUNCTIONS:

URI_paste <- function(URI_prefix, vec) {
  # Paste a list of strings into a list of URIs according to the input prefix
  # 
  # Args:
  #   URI_prefix: URI prefixes in str.
  #   vec: a vector containing targitted strings. 
  # 
  # Returns:
  #  A vector containing the results one on one to the the input vector
  
  URIs <- sapply(vec, function(strings) {
    paste0("<", URI_prefix, strings, ">")
  })
  names(URIs) <- NULL
  return(URIs)
}

# ================ Function to make the queries: ======================== 

potential_interaction_query_make <- function(target_genes_URIs, chr_num, start_pos, end_pos) {
  # Make the query of potential interactions between a list of genes and genes in a certain regiion
  #
  # Args: 
  #  target_genes: a vector containing genes in URI form. 
  #  chr: chromosom location of potential genes interacting with the input genes in numeric
  #  start_pos: the starting positions of a region where the potetial genes located on the chromosome
  #  end_pos: the starting positions of a region where the potetial genes located on the chromosome
  #  
  # Returns:
  #  A sparql query that obtain the potential interactions, together with the chomosome locations 
  #  and the scores of the interactions.  
  
  values <- paste(target_genes_URIs, collapse = "\n")  # make the vector as the input for the query
  query <- paste0("prefix obo: <http://purl.obolibrary.org/obo/> 
  prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> 
  prefix term: <http://rdf.ebi.ac.uk/terms/ensembl/> 
  prefix dc: <http://purl.org/dc/elements/1.1/> 
  prefix faldo: <http://biohackathon.org/resource/faldo#>
  prefix sio: <http://semanticscience.org/resource/>
  prefix core: <http://purl.uniprot.org/core/>
  prefix so: <http://purl.obolibrary.org/obo/so#>
  prefix gene: <http://rdf.ebi.ac.uk/resource/ensembl/>
  
  SELECT DISTINCT ?pot_proteins ?tar_proteins ?tar_loc_refs ?weights
  WHERE{
   VALUES ?target_genes {",values, "}
   VALUES ?loc_refs { <http://rdf.ebi.ac.uk/resource/ensembl/33/arabidopsis_thaliana/TAIR10/", as.character(chr_num), "> } 
   FILTER (?begin_positions >= " , as.character(start_pos), " &&
            ?end_positions <= " , as.character(end_pos), ") .
   
  #-- query the chromosomes where tar genes located
    ?target_genes a obo:SO_0001217 ;
                  faldo:location ?tar_locations .
    ?tar_locations rdfs:label ?tar_loc_refs .
  #-- query the protein encoded by tar genes
    ?tar_transcripts so:transcribed_from ?target_genes .
    ?tar_transcripts  so:translates_to ?tar_proteins .
  
  #-- query the potential genes and accoding to the locations
    ?pot_genes a obo:SO_0001217 ;
               faldo:location ?pot_genes_locations .
     ?pot_genes_locations faldo:begin ?begins ; 
                         faldo:end ?ends  ; 
                         faldo:reference ?loc_refs .
    ?begins faldo:position ?begin_positions .
    ?ends faldo:position ?end_positions .
  #-- query the protein encoded by pot genes
    ?pot_transcripts so:transcribed_from ?pot_genes .
    ?pot_transcripts  so:translates_to ?pot_proteins .
  
  #-- query the interactions between tar and pot genes
     ?pot_proteins (obo:RO_0002434|^obo:RO_0002434) ?tar_proteins . 
     ?pot_proteins sio:SIO_000062 ?interaction .
     ?tar_proteins sio:SIO_000062 ?interaction .
     ?interaction sio:SIO_000300 ?weights .
  }")
  return(query)
}

interaction_two_regions_query_make <- function(chr_num_1, start_pos_1, end_pos_1, chr_num_2, start_pos_2, end_pos_2) {
  query <- paste0("prefix obo: <http://purl.obolibrary.org/obo/> 
   prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> 
   prefix term: <http://rdf.ebi.ac.uk/terms/ensembl/> 
   prefix dc: <http://purl.org/dc/elements/1.1/> 
   prefix faldo: <http://biohackathon.org/resource/faldo#>
   prefix sio: <http://semanticscience.org/resource/>
   prefix core: <http://purl.uniprot.org/core/>
   prefix so: <http://purl.obolibrary.org/obo/so#>
   prefix gene: <http://rdf.ebi.ac.uk/resource/ensembl/>
   
   SELECT DISTINCT ?proteins_1 ?proteins_2 ?weights
   WHERE{
    VALUES ?genes_1_loc_refs { <http://rdf.ebi.ac.uk/resource/ensembl/33/arabidopsis_thaliana/TAIR10/", as.character(chr_num_1), "> } 
    FILTER (?genes_1_begin_positions >= ", as.character(start_pos_1), " &&
              ?genes_1_end_positions <= ", as.character(end_pos_1), ") .
    
    VALUES ?genes_2_loc_refs { <http://rdf.ebi.ac.uk/resource/ensembl/33/arabidopsis_thaliana/TAIR10/", as.character(chr_num_2), "> } 
    FILTER (?genes_2_begin_positions >= ", as.character(start_pos_2), " &&
              ?genes_2_end_positions <= ", as.character(end_pos_2), ") .
    
    # query the potential genes and accoding to the locations
    ?genes_1 a obo:SO_0001217 ;
    faldo:location ?genes_1_locations .
    ?genes_1_locations faldo:begin ?genes_1_begins ; 
    faldo:end ?genes_1_ends  ; 
    faldo:reference ?genes_1_loc_refs .
    ?genes_1_begins faldo:position ?genes_1_begin_positions .
    ?genes_1_ends faldo:position ?genes_1_end_positions .
    
    ?genes_2 a obo:SO_0001217 ;
    faldo:location ?genes_2_locations .
    ?genes_2_locations faldo:begin ?genes_2_begins ; 
    faldo:end ?genes_2_ends  ; 
    faldo:reference ?genes_2_loc_refs .
    ?genes_2_begins faldo:position ?genes_2_begin_positions .
    ?genes_2_ends faldo:position ?genes_2_end_positions .
    
    
    # query the protein encoded by genes located in region 1
    ?genes_1_transcripts so:transcribed_from ?genes_1 .
    ?genes_1_transcripts  so:translates_to ?proteins_1 .
    
    # query the protein encoded by pot genes located in region 2
    ?genes_2_transcripts so:transcribed_from ?genes_2 .
    ?genes_2_transcripts  so:translates_to ?proteins_2 .
    
    
    # query the interactions between proteins_1 and proteins_2
    ?proteins_1 (obo:RO_0002434|^obo:RO_0002434) ?proteins_2 . 
    ?proteins_1 sio:SIO_000062 ?interaction .
    ?proteins_2 sio:SIO_000062 ?interaction .
    ?interaction sio:SIO_000300 ?weights .
    }")
  return(query)
}

description_for_proteins_query_make <- function(protein_URIs){
  # Function to make the query for uniprot anotations of the input protein in URIs.
  # 
  # Args:
  #  protein_URIs: the URIs used for the proteining in the triplet store
  # 
  # Return:
  #  a query used for retrieving the uniprot anotations in strings
  
  values <- paste0(protein_URIs, collapse = " ")  # make the vector as the input for the query
  query <- paste0(
  "prefix obo: <http://purl.obolibrary.org/obo/> 
  prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> 
  prefix term: <http://rdf.ebi.ac.uk/terms/ensembl/> 
  prefix dc: <http://purl.org/dc/elements/1.1/> 
  prefix faldo: <http://biohackathon.org/resource/faldo#>
  prefix sio: <http://semanticscience.org/resource/>
  prefix core: <http://purl.uniprot.org/core/>
  prefix so: <http://purl.obolibrary.org/obo/so#>
  prefix protein: <http://rdf.ebi.ac.uk/resource/ensembl.protein/>
  
  
  SELECT DISTINCT ?protein_ids ?descriptions
  WHERE{
   VALUES ?proteins { ", values , " }
   ?proteins dc:identifier ?protein_ids . 
   OPTIONAL{
            {
            ?proteins (term:SEQUENCE_MATCH|term:CHECKSUM) ?uniprot_annos .
            ?uniprot_annos dc:description ?descriptions .
            ?uniprot_annos a core:Protein .
            } UNION {
             ?transcripts (core:translatedTo | so:translates_to) ?proteins .
             ?transcripts (core:transcribedFrom | so:transcribed_from) ?genes .
             ?genes dc:description ?descriptions .
            }
   }
 
  }")
  return(query)
}

KEGG_annotation_retrival_query_make <- function(protein_URIs){
  
  values <- paste0(protein_URIs, collapse = " ")
  query <- paste0("prefix obo: <http://purl.obolibrary.org/obo/> 
  prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> 
  prefix term: <http://rdf.ebi.ac.uk/terms/ensembl/> 
  prefix dc: <http://purl.org/dc/elements/1.1/> 
  prefix faldo: <http://biohackathon.org/resource/faldo#>
  prefix sio: <http://semanticscience.org/resource/>
  prefix core: <http://purl.uniprot.org/core/>
  prefix protein: <http://rdf.ebi.ac.uk/resource/ensembl.protein/>
  
  
  SELECT DISTINCT ?protein_ids ?pw_ids ?pw_des 
  WHERE{
    VALUES ?proteins { ", values , " }
    VALUES ?db { <http://purl.uniprot.org/database/KEGG> }
    ?proteins a term:protein .
    ?proteins dc:identifier ?protein_ids .
     OPTIONAL{ 
           ?proteins term:SEQUENCE_MATCH ?uniprot_annos .
           ?uniprot_annos rdfs:seeAlso ?xrefs .
           ?xrefs sio:SIO_000062 ?pathways .
           ?pathways dc:description ?pw_des ;
                     dc:identifier ?pw_ids .
          ?xrefs core:database ?db .     
     }
  }")
  return(query)
}

GO_annotation_retrival_query_make <- function(protein_URIs) {
  
  values <- paste0(protein_URIs, collapse = " ")
  query <- paste0("prefix obo: <http://purl.obolibrary.org/obo/> 
  prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> 
  prefix term: <http://rdf.ebi.ac.uk/terms/ensembl/> 
  prefix dc: <http://purl.org/dc/elements/1.1/> 
  prefix core: <http://purl.uniprot.org/core/>
  prefix protein: <http://rdf.ebi.ac.uk/resource/ensembl.protein/>
  prefix goFormat: <http://www.geneontology.org/formats/oboInOwl#>
  
  SELECT DISTINCT ?protein_ids ?go_terms ?go_descrips 
  WHERE{
    VALUES ?proteins { ", values , " }
    ?proteins a term:protein .
    ?proteins dc:identifier ?protein_ids .
     OPTIONAL{ 
           ?proteins term:SEQUENCE_MATCH ?uniprot_annos .
            ?uniprot_annos core:classifiedWith  ?classes .
            ?classes rdfs:subClassOf+ ?ancestors .
            ?ancestors goFormat:id ?go_terms ;
                              rdfs:label ?go_descrips .
     }
  }")
  return(query)
}

GO_hierarchy_retrival_query_make <- function(go_URIs) {
  
  values <- paste0(go_URIs, collapse = " ")
  query <- paste0("prefix obo: <http://purl.obolibrary.org/obo/> 
  prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> 
  prefix term: <http://rdf.ebi.ac.uk/terms/ensembl/> 
  prefix dc: <http://purl.org/dc/elements/1.1/> 
  prefix core: <http://purl.uniprot.org/core/>
  prefix protein: <http://rdf.ebi.ac.uk/resource/ensembl.protein/>
  prefix goFormat: <http://www.geneontology.org/formats/oboInOwl#>
  prefix oboGo: <http://purl.obolibrary.org/obo/>
  
  SELECT DISTINCT ?from_go_terms ?form_go_descrips ?to_go_terms ?to_go_descrips
  WHERE{
    VALUES ?go_from { ", values , " }
    VALUES ?go_to {",  values ,"}

    ?go_from rdfs:subClassOf ?go_offsprings .
    ?go_from goFormat:id ?from_go_terms ;
              rdfs:label ?form_go_descrips .
    ?go_to goFormat:id ?to_go_terms ; 
          rdfs:label ?to_go_descrips .
    filter(?go_to = ?go_offsprings)
  }")
  return(query)
}

# assemble the query bussiness: 
data_require <- function(query_function, calls, endpoint) {
  # Call the function query make function and get the results from the endpoints   
  # 
  # Args: 
  #  query_function: one of the query making functions defined above. 
  #  calls: a list containing the arguments for each parameters required by the 'query_function' .
  #  endpoint: the endpoint url in string where the triplets are stored.
  # 
  # Return:
  #  A matrix with queried table.
  
  require(SPARQL)
  
  query <- do.call(what = query_function, args = as.list(calls))
  results <- SPARQL(endpoint , query)
  res_df <- results$results
  return(res_df)
}

# ======================== Function for data cleaning ==============================

regex_match <- function(vec, pattern) {
  # Match the text with reg express. only text matched with the expression will be return.
  # 
  # Args: 
  #   vec: the vector containing strings for pattern matching
  #   pattern: the regular expression that will be appied to all the strings in the input vector.
  #
  # Returns:
  #  A vector containing the results one on one to the the input vector
  
  pos <- regexec(pattern = pattern, text = vec, ignore.case = T)
  matches <- regmatches(x = vec, m = pos)
  results <- sapply(1:length(matches), function(i) {
    matches[[i]][1]
  })
  return(results)
}

regex_match_map <- function(data_frame, regex_vec) {
  # Apply each column a regular expression formula to extract information out
  # 
  # Args:
  #  data_frame: a dataframe containing the strings
  #  regex_vec: the vector of the same length as the data_frame containig the regular
  #  expressions for each colum.
  # 
  # Returns:
  #  A dataframe containing the exatracted strings accoding to the reg ex applied to each column. 
  
  if (length(data_frame) != length(regex_vec)) {
    stop("Arguments data_frame and regex_vec have invalid lengths: ",
         length(data_frame), " and ", length(regex_vec), ".")
  }
  data_frame[,] <- sapply(1:length(data_frame), FUN = function(ind) {
    regex_match(vec = data_frame[, ind], pattern = regex_vec[ind])  # defined in the same environment
  })
  return(data_frame)
}

# MAIN:

# ====== 1) make the graph using according to the QTL data and STRING db ===============

# the list of gene is an example:
list_of_genes <- read.csv(file = "list_of_gene.txt", header = T, stringsAsFactors = F)
genes <- unlist(list_of_genes, use.names = F)
URIs <- paste0("gene:", genes)

# designate the endpoint and the calls query:
endpoint <- "http://localhost:8890/sparql"
calls <- list(URIs, 5, 15650000, 23850000)
interactions <- data_require(potential_interaction_query_make, calls, endpoint)

# data preprocessing:
patterns <- c("AT[0-9]G[0-9|.]+", "AT[0-9]G[0-9|.]+", "^chromosome [0-9]", "[0-9]+")
interactions <- regex_match_map(data_frame = interactions, regex_vec = patterns)
interactions[, 3] <- as.factor(interactions[, 3])
interactions[, 4] <- as.numeric(interactions[, 4])

# make the graph  

cutoff <- quantile(interactions$weights, probs = 0.25)  # use the third quantile as the cutoff
cutoff <- 0
links <- interactions[, c(1, 2, 4)][interactions$weights > cutoff, ]
g <- graph_from_data_frame(d = links, directed = F) 

# =================2) retrieve the annotations of Uniprot and KEGG for ================
# =================genes serving as vertexes in the graph =============================

# view the vertexes inside the network
vertex_names <- V(g)$name

# retrieve the annotations of proteomes and KEGG pw for the nodes

protein_URIs <- paste0("protein:", vertex_names)
protein_descriptions_nodes <- data_require(description_for_proteins_query_make, list(protein_URIs), endpoint)
kegg_annotation_nodes <- data_require(KEGG_annotation_retrival_query_make, list(protein_URIs), endpoint)
go_annotation_nodes <- data_require(GO_annotation_retrival_query_make, list(protein_URIs), endpoint)
View(protein_descriptions_nodes )
# data clean for annotation table

patterns <- c(".*", "map[0-9]+", "[A-Z][^\"]+")
kegg_annotation_nodes <- regex_match_map(data_frame = kegg_annotation_nodes, regex_vec = patterns)

# community detection using GN algorithm:

g_community <- edge.betweenness.community(g)

# Calculate the betweeness centrality
betw_centralities <- betweenness(g, v = V(g), directed = F, weights = NULL, normalized = T)
betw_centralities <- betw_centralities[order(betw_centralities, decreasing = T)]
betw_cutoff <- quantile(betw_centralities, probs = 0.8)
nodes_with_higer_betw <- names(betw_centralities[betw_centralities > betw_cutoff])

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
# specify the parameters for the network visulization:
#  shapes circle for target genes and sqare for potential genes:
shapes() #all shapes available
shape_pos <- vertex_names %in% links$tar_proteins
v_shapes <- c("circle",  "square")[shape_pos+1]  # diffrent shape for target and potential genes 
# colors stands for the pathway evolved:
color_vec <- kegg_annotation_nodes$pw_des
color_vec[is.na(color_vec)] <- "Unknown" 
color_list <- tapply(X = color_vec, INDEX = kegg_annotation_nodes$protein_ids, c)
#color_list <- sapply(color_list, as.factor)
color_list <- color_list[match(vertex_names, names(color_list))]
vex_values <- sapply(color_list, function(x) {
  rep(1,length(x)) 
})

chars <- sapply(color_list, function(x) {
  x[1]
})

chars <- as.factor(chars)

plot(g, vertex.shape = v_shapes, 
     vertex.color=chars,
     vertex.size = 7,
     vertex.label=NA)
col_legends <- as.numeric(as.factor(levels(as.factor(chars))))

legend("topleft", 
       legend=levels(as.factor(chars)), 
       col = col_legends, 
       bty = "n", 
       pch=20 , 
       pt.cex = 1.5, 
       cex = 0.75, 
       text.col=coul , 
       horiz = FALSE, 
       inset = c(0.1, 0.1))

# select the names to show
vernames_to_show <- links$tar_proteins[match(vertex_names, links$tar_proteins)]
chars <- interactions$tar_loc_refs[match(vertex_names, interactions$tar_proteins)]
chars[is.na(chars)] <- as.factor("chromosome 5")
col_pos <- as.numeric(chars)
coul <- brewer.pal(5, "Set1")
my_color <- coul[col_pos]
vertex_names[98]


# plot the graph
plot(g, 
     vertex.size = 5, 
     vertex.color=my_color, 
     vertex.label= vernames_to_show,  
     vertex.shape = v_shapes)


legend("topleft", 
       legend=levels(as.factor(chars)), 
       col = coul, 
       bty = "n", 
       pch=20 , 
       pt.cex = 1.5, 
       cex = 0.75, 
       text.col=coul , 
       horiz = FALSE, 
       inset = c(0.1, 0.1))