rm(list = ls())
REPOSITORY <- "https://cran.rstudio.com/"
REQUIRED_PACKS <- c("SPARQL", "igraph", "snowfall", "networkD3", "htmlwidgets")

installed_packs <- available.packages(repos = REPOSITORY)
package_names <- row.names(installed_packs)

for (required_pack in REQUIRED_PACKS){
  if(required_pack %in% installed_packs) {
    print(paste("The required packege:", required_pack, "is installed. "))
  } else {
    print(paste("This software requares:", required_pack, ". Installing"))
    install.packages(required_pack, repos = REPOSITORY)
  }
}
  