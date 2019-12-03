# http://rpubs.com/mathetal/clusterlabel
# \\Public dataset from the Allen Institute   https://celltypes.brain-map.org/rnaseq
#\\Labeling tutorial from Seurat https://satijalab.org/seurat/v3.0/vis...
#\\R code  http://rpubs.com/mathetal/clusterlabel
#install.packages('Seurat')  # https://satijalab.org/seurat/v3.0/visualization_vignette.html
library(Seurat)
library(dplyr)
#read in table
getwd()
setwd("~/Documents/STUDY/DATA_SCIENCE/DNA_human_brain/human_LGN_gene_expression_matrices_2018-06-14")

library(Seurat)
library(ggplot2)
#read in table
#setwd("~/Desktop/Classes/Bioinformatics/final_project/human")
r1 = read.table ("human_LGN_2018-06-14_exon-matrix.csv", header= TRUE, sep= ",", row.names = 1)
r1 %>% dim
#remove genes in less than 30 cells
rm_genes = which(rowSums(r1 > 0) <30)
r2 = r1 [-rm_genes,]

library(scran) # http://bioconductor.org/packages/release/bioc/html/scran.html
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# BiocManager::install("scran")
library(SingleCellExperiment)
library(SummarizedExperiment)
library(GenomicRanges)
library(stats4)
library(BiocGenerics)
library(parallel)
library(matrixStats)
library(S4Vectors)

## Loading required package: SingleCellExperiment
## Loading required package: SummarizedExperiment
## Loading required package: GenomicRanges
## Loading required package: stats4
## Loading required package: BiocGenerics
## Loading required package: parallel
##
## Attaching package: 'BiocGenerics'
## The following objects are masked from 'package:parallel':
##
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## The following objects are masked from 'package:stats':
##
##     IQR, mad, sd, var, xtabs
## The following objects are masked from 'package:base':
##
##     anyDuplicated, append, as.data.frame, basename, cbind,
##     colnames, dirname, do.call, duplicated, eval, evalq, Filter,
##     Find, get, grep, grepl, intersect, is.unsorted, lapply, Map,
##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
##     pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
##     setdiff, sort, table, tapply, union, unique, unsplit, which,
##     which.max, which.min
## Loading required package: S4Vectors
##
## Attaching package: 'S4Vectors'
## The following object is masked from 'package:base':
##
##     expand.grid
## Loading required package: IRanges
## Loading required package: GenomeInfoDb
## Loading required package: Biobase
## Welcome to Bioconductor
##
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## Loading required package: DelayedArray
## Loading required package: matrixStats
##
## Attaching package: 'matrixStats'
## The following objects are masked from 'package:Biobase':
##
##     anyMissing, rowMedians
## Loading required package: BiocParallel
##
##Attaching package: 'DelayedArray'
## The following objects are masked from 'package:matrixStats':
##
##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges
## The following objects are masked from 'package:base':
##
##     aperm, apply, rowsum
## Registered S3 methods overwritten by 'ggplot2':
##   method         from
##   [.quosures     rlang
##   c.quosures     rlang
##   print.quosures rlang
#create sce object
sce = SingleCellExperiment(list(counts = data.matrix(r2)))

#do clustering to reduce heterogeneity
clusters = quickCluster(sce, min.size=100)
clusters
## Warning: Setting 'use.ranks=TRUE' for the old defaults.
## Set 'use.ranks=FALSE' for the new defaults.
sce = computeSumFactors (sce, cluster = clusters)

#normalize, don't return log2
sce = normalize (sce, return_log = FALSE)

library(Seurat)
## Registered S3 method overwritten by 'R.oo':
##   method        from
##   throw.default R.methodsS3
#create Seurat object using scran object
s_obj= CreateSeuratObject(counts = log(counts(sce) +1))

#find variable features to reduce run time
s_obj=FindVariableFeatures(s_obj)

#regress out batch
s_obj_nobatch = ScaleData(s_obj)
## Centering and scaling data matrix
s_obj = RunPCA(s_obj_nobatch, features = VariableFeatures(s_obj))
## PC_ 1
## Positive:  4340, 4336, 5129, 4099, 51148, 57571, 115584, 9725, 55314, 60484
##     9705, 933, 745, 84735, 101929249, 2065, 5354, 94015, 1287, 57471
##     6035, 5376, 2628, 54443, 1307, 347, 79152, 10397, 5653, 1298
## Negative:  55384, 692218, 157627, 9699, 100130155, 1159, 440823, 285175, 6326, 6323
##     4133, 22854, 4897, 65009, 8745, 359822, 10451, 2561, 6857, 817
##     80144, 490, 54839, 6616, 151835, 9118, 2977, 3716, 2903, 115827
## PC_ 2
## Positive:  2572, 3815, 2571, 5649, 80309, 594855, 3479, 654790, 11122, 23236
##     6092, 1012, 2890, 10769, 132204, 53826, 9547, 57689, 148, 56937
##     116443, 30820, 84959, 11249, 10752, 400120, 2897, 129684, 3356, 6529
## Negative:  283131, 6474, 105369345, 57084, 6934, 128553, 1286, 253650, 6696, 153579
##     4741, 1285, 151835, 5774, 10777, 5816, 4744, 2620, 3736, 105377862
##     4747, 5168, 11113, 4124, 8404, 266722, 10451, 5121, 80144, 7018
## PC_ 3
## Positive:  477, 64344, 5625, 5803, 493, 2261, 361, 3371, 80036, 7026
##     254295, 2697, 5015, 165, 2152, 2202, 27151, 6563, 9628, 3280
##     50509, 2670, 54739, 6538, 10410, 59352, 57326, 3133, 718, 80162
## Negative:  285987, 117177, 64919, 5318, 54715, 3925, 24141, 85445, 6695, 11075
##     2554, 388121, 1404, 105375415, 11249, 56892, 27145, 115827, 885, 2561
##     230, 2786, 10368, 347689, 53942, 100287225, 100996645, 758, 1268, 22895
## PC_ 4
## Positive:  5015, 8403, 4212, 2625, 23414, 148281, 7026, 613212, 55515, 4880
##     100309464, 3356, 104355219, 644192, 56934, 105369212, 56243, 846, 132204, 5617
##     100131897, 401498, 6096, 777, 9628, 80059, 1010, 1008, 5334, 493
## Negative:  3398, 7164, 4781, 477, 285987, 64919, 24141, 5318, 64344, 10485
##     1462, 6925, 5625, 2261, 105375415, 1404, 2670, 885, 2066, 388121
##     361, 339789, 3937, 11281, 1268, 3730, 11170, 165, 27145, 30010
## PC_ 5
## Positive:  5625, 2261, 477, 64344, 361, 3371, 6563, 50509, 2697, 59352
##     5803, 2202, 27151, 3280, 165, 10840, 4502, 6385, 254295, 57326
##     2152, 10485, 80036, 54762, 4776, 123041, 3955, 6263, 7837, 3691
## Negative:  3127, 3687, 1524, 3128, 7305, 3123, 713, 1436, 7805, 718
##     3122, 2207, 81704, 963, 1535, 54518, 3903, 3394, 64581, 64092
##     1536, 4973, 6916, 3635, 6850, 3109, 6001, 10578, 10859, 8832
s_obj = FindNeighbors(s_obj)
## Computing nearest neighbor graph
## Computing SNN
s_obj = FindClusters(s_obj, resolution = 0.15)
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
##
## Number of nodes: 1576
## Number of edges: 48040
##
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9445
## Number of communities: 7
## Elapsed time: 0 seconds
table(Idents(s_obj))
##
##   0   1   2   3   4   5   6
## 918 185 139 120  94  92  28
s_obj = RunTSNE(s_obj)
DimPlot (s_obj, reduction = "tsne", label = TRUE)


#2571 = GAD1, GABAergic neuron; 4340 = MOG, oligodendrocyte; 361 = AQP4, astrocyte
#determine cell types
features <- c("2571")

RidgePlot(object = s_obj, features = features, ncol = 2)
## Picking joint bandwidth of 0.0866


VlnPlot(object = s_obj, features = features)


FeaturePlot(object = s_obj, features = features)


DotPlot(object = s_obj, features = features) + RotatedAxis()


#label clusters based on marker genes
current.cluster.ids = c(0, 1, 2, 3, 4, 5, 6)
new.cluster.ids = c("Glutamatergic neuron", "Oligodendrocyte", "GABAergic neuron", "Glutamatergic neuron",
                    "Oligodendrocyte precursor ", "GABAergic neuron", "\n Astrocyte")
names(x = new.cluster.ids) <- levels(x = s_obj)
pbmc <- RenameIdents(object = s_obj, new.cluster.ids)

#plot again
DimPlot(object = pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()

