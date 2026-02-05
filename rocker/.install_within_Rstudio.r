# Install Packages for Class
# install.packages('igraph')

BiocManager::install("airway")
library(airway)
# BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)
BiocManager::install("DESeq2")
library(DESeq2)
BiocManager::install("DelayedArray", force = TRUE)
library(DelayedArray)
BiocManager::install("vsn")
library(vsn)
BiocManager::install("apeglm")
library(apeglm)
BiocManager::install("AnnotationDbi")
library(AnnotationDbi)
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)


# On terminal 
# apt-get install update 
# apt-get install r-cran-igraph
# apt-get install libglpk40
BiocManager::install("igraph")
library(igraph)



BiocManager::install("clusterProfiler")
library(clusterProfiler)
BiocManager::install("enrichplot")
library(enrichplot)



tempdir()
&& rm -rf /tmp/downloaded_packages \
&& strip /usr/local/lib/R/site-library/*/libs/*.so


