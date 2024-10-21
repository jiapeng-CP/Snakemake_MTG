# visualize taxon abundance with sankey plot by using Pavian https://github.com/fbreitwieser/pavian

if (!require(remotes)) { install.packages("remotes") }
remotes::install_github("fbreitwieser/pavian")

shiny::runGitHub("fbreitwieser/pavian", subdir = "inst/shinyapp")

# load files in the folder MetaPhlan, not merged.metaphlan.txt
