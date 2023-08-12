library(tasic2016data)
library(dendextend)
library(matrixStats)
library(Matrix)
library(scrattch.hicat)

source("D:/courses/cs6024_aacb/project/src/utils.R")


setwd("D:/courses/cs6024_aacb/project/datasets/zeisel_2015")

reads <- read.csv2("zeisel_2015_reads.csv", header=TRUE, row.names = 1, sep = ',')
anno <- read.csv2("zeisel_2015_clean_anno.csv", header=TRUE, row.names = 1, sep = ',')

### Select excitatory neurons

reads <- reads[anno$level1class == "pyramidal CA1" | anno$level1class == "pyramidal SS",]
anno <- anno[anno$level1class == "pyramidal CA1" | anno$level1class == "pyramidal SS",]

write.csv2(anno, file="zeisel_2015_clean_anno_eneurons.csv")

###

reads <- t(reads)
reads_matrix <- data.matrix(reads)
norm.dat <- log2(cpm(reads_matrix)+1)

write.csv2(norm.dat, file="zeisel_2015_logcpm_eneurons.csv")

# Convert to sparse matrix.
# norm.dat <- Matrix(cpm(tasic_2016_counts), sparse = TRUE)
# norm.dat@x <- log2(norm.dat@x+1)


# Parameter specification
de.param <- de_param(padj.th     = 0.05, 
                     lfc.th      = 1, 
                     low.th      = 1, 
                     q1.th       = 0.5, 
                     q.diff.th   = 0.7, 
                     de.score.th = 40)

filtered_dataset <- get_rd_dat(norm.dat, 
                                method = "kmeans",
                                dim.method = "pca", 
                                de.param = de_param(de.score.th=500),
                                verbose = TRUE)

filtered_dataset_df <- data.frame(filtered_dataset)

write.csv2(filtered_dataset_df, file="zeisel_2015_pc_eneurons.csv")

# Perform one-step clustering. Kmeans and PCA.
onestep.result <- onestep_clust(norm.dat, 
                                select.cells = select.cells, 
                                method = "kmeans",
                                dim.method = "pca", 
                                de.param = de_param(de.score.th=500))

# Perform iterative clustering.
kmeans.iter.clust <- iter_clust(norm.dat, 
                               select.cells = select.cells, 
                               method = "kmeans",
                               dim.method = "pca", 
                               de.param = de.param, 
                               result=onestep.result)

# Display
display.result = display_cl(onestep.result$cl, norm.dat, plot=TRUE, de.param=de.param)
