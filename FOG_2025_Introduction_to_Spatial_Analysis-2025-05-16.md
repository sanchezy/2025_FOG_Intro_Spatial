# Introduction

In this tutorial you are going to learn some common analysis using
*Spatial Data*. Rather than exploring all the possible analysis
pipelines, this tutorial is intended to be an introduction to spatial
analysis. Additional analytical steps will depend on the technology
used, the tissue and the biological question you are addressing.

The spatial data in this tutorial is a toy dataset that was subsampled
from a real pancreatic cancer dataset provided by Adam Bryce. It is a
multiplex (mIF) experiment using the following markers to detect some
important cell types:

# CD4 and FOXP3: Tregs

# CD68: Macrophage

# PDL1: Functional / Immune-oncology

# aSMA: Fibroblast

# panCK: Epithelium

*IMPORTANT. Before the session*: Please download the folder (link
provided) and install the list of packages requiered in Rstudio (i.e
using the command: install.packages(‘rmarkdown’))

# Loading the requiered libraries and the data

    # load the libraries
    library(rmarkdown)
    library(reshape2)   # melt()
    library(ggplot2)    # ggplot()
    library(uwot)       # umap()
    library(RANN)       # nn2() - rapid nearest-neighbour computation
    library(patchwork) # combine plots
    library(Matrix) 

    # read the data in - You need to change the path 
    data = read.csv("~/Documents/Projects-analysis/FOG_2025_Introduction_Spatial_Analysis/data/workshop_data.csv",header=T)

    # Print the first 5 lines of the head of "data". What do the values represent? 
    head(data,5)

    ##       CellX     CellY        CD4      FOXP3       CD68      PDL1     panCK      aSMA
    ## 1 1139.2352 -10537.27 0.25609556 0.04672468 0.14541020 0.2209569 2.9953754 2.0016484
    ## 2  642.5411 -10842.73 0.42196316 0.12131447 0.08518437 0.2146205 4.1189364 1.5880690
    ## 3  588.9074 -11276.76 0.11724587 0.11188758 0.09883421 0.0977579 0.6153455 1.5510375
    ## 4  642.2683 -11345.83 0.09599299 0.08920708 0.79192211 0.3394328 1.8501472 2.0249361
    ## 5 1297.6909 -11001.18 0.21191688 0.15399958 0.47574753 0.1113520 0.9091848 0.4123984

# Exploratory visualisations

First, we are going to check the data (i.e. is the tissue we expect?).
<br>

    # view the tissue: plot the centroids of each cell
    plot_tissue <- ggplot(data, aes(CellX,CellY)) +
      geom_point()
    plot_tissue

![](/Users/ysanchez/Documents/Projects-analysis/FOG_2025_Introduction_Spatial_Analysis/notebooks-html/FOG_2025_Introduction_to_Spatial_Analysis-2025-05-16_files/figure-markdown_strict/viewthetissue-1.png)

    # colour-code by a marker of your choice
    # [insert you code here]

    # Check the intensity of different markers. What can you observe? 
      # boxplots  
      plot1 <- ggplot(melt(data[,3:8])) +
        geom_boxplot(aes(variable,value)) + ggtitle("Boxplot before scaling")

    ## No id variables; using all as measure variables

      # density plots
      plot2 <- ggplot(melt(data[,3:8])) +
        geom_density(aes(value,fill=variable),alpha=0.5) + ggtitle("Density plot before scaling")

    ## No id variables; using all as measure variables

    plot1 + plot2

![](/Users/ysanchez/Documents/Projects-analysis/FOG_2025_Introduction_Spatial_Analysis/notebooks-html/FOG_2025_Introduction_to_Spatial_Analysis-2025-05-16_files/figure-markdown_strict/before-after-scaling-1.png)

# Normalisation

It is often advisable to normalise and scale the data. The specifics
will depend on the type of data and method used. Here we do a Z-score
scaling.

    # recall the data structure. Print the head of the file:
    head(data,3)

    ##       CellX     CellY       CD4      FOXP3       CD68      PDL1     panCK     aSMA
    ## 1 1139.2352 -10537.27 0.2560956 0.04672468 0.14541020 0.2209569 2.9953754 2.001648
    ## 2  642.5411 -10842.73 0.4219632 0.12131447 0.08518437 0.2146205 4.1189364 1.588069
    ## 3  588.9074 -11276.76 0.1172459 0.11188758 0.09883421 0.0977579 0.6153455 1.551037

    # Use Z-score to normalised the intensity values
    data[,3:8] = scale(data[,3:8])

    # Print again the head of the file. Please notice the change on the values.
    head(data,3)

    ##       CellX     CellY        CD4      FOXP3       CD68       PDL1      panCK      aSMA
    ## 1 1139.2352 -10537.27  0.4795359 -0.3482475 -0.3959196  0.3070412  0.3477149 1.4612218
    ## 2  642.5411 -10842.73  1.4127056  0.3049579 -0.5893898  0.2489284  0.9204994 0.9406002
    ## 3  588.9074 -11276.76 -0.3016313  0.2224037 -0.5455408 -0.8228565 -0.8656099 0.8939842

    # Let's inspect again the boxplot and the density plots after normalisation. What do you notice?
    # boxplots 
    plot4 <-  ggplot(melt(data[,3:8])) +
        geom_boxplot(aes(variable,value)) + ggtitle("Boxplot after scaling")

    ## No id variables; using all as measure variables

    # densitity plot
    plot5 <-   ggplot(melt(data[,3:8])) +
        geom_density(aes(value,fill=variable),alpha=0.5) + ggtitle("Density plot after scaling")

    ## No id variables; using all as measure variables

    plot1 + plot4 + plot2 + plot5 + plot_layout(ncol = 2)

![](/Users/ysanchez/Documents/Projects-analysis/FOG_2025_Introduction_Spatial_Analysis/notebooks-html/FOG_2025_Introduction_to_Spatial_Analysis-2025-05-16_files/figure-markdown_strict/effects-on-scaling-1.png)

# Clustering to identify cell types

An common task in spatial analysis is to find the groups of cells
belonging to the same cell type or state. A typical workflow will
include a step of clustering, that groups cells depending on their
similarity of measured features (ie. transcriptome or proteome). Here,
we use a simple k-means here to illustrate this analysis.

    set.seed(123)    # set seed for reproducibility

    # calculate k means - here we set k = 4. Why do you think k = 4?
    kmeans = kmeans(scale(data[,3:8]),4)

    # Let's explore the results. View centres.
    kmeans$centers

    ##          CD4       FOXP3       CD68        PDL1      panCK       aSMA
    ## 1 -0.4166886 -0.36942857 -0.4037192 -0.03366421  0.8849618 -0.4566071
    ## 2 -0.3951473 -0.22532623  1.9765214  0.83358042 -0.4784028 -0.2439522
    ## 3  1.2467453  0.79124597 -0.3871663 -0.23211498 -0.4915175 -0.3670943
    ## 4 -0.3556531 -0.08491421 -0.1471558 -0.19801750 -0.5672805  1.1641307

    # Check cluster sizes
    table(kmeans$cluster)

    ## 
    ##   1   2   3   4 
    ## 739 282 479 500

<br> *How would you assign cell types identities to each group/clusters
of cells?* Hint: think about the definition of a marker cells. What is
the most likely identity of cluster 1?

    # looks like cluster 1 is Epi, 2 is Macrophage, 3 is Treg, 4 is CAF
    # Add the results from the k-mean to the 'data' 
    data$cluster = factor(kmeans$cluster,labels=c("Epi","Macrophage","Treg","CAF"))

    # check the data again
    head(data)

    ##       CellX     CellY        CD4       FOXP3       CD68       PDL1      panCK       aSMA    cluster
    ## 1 1139.2352 -10537.27  0.4795359 -0.34824751 -0.3959196  0.3070412  0.3477149  1.4612218        CAF
    ## 2  642.5411 -10842.73  1.4127056  0.30495786 -0.5893898  0.2489284  0.9204994  0.9406002       Treg
    ## 3  588.9074 -11276.76 -0.3016313  0.22240370 -0.5455408 -0.8228565 -0.8656099  0.8939842        CAF
    ## 4  642.2683 -11345.83 -0.4211997  0.02378371  1.6809429  1.3936226 -0.2361155  1.4905368 Macrophage
    ## 5 1297.6909 -11001.18  0.2309870  0.59119125  0.6652600 -0.6981804 -0.7158124 -0.5393562       Treg
    ## 6 1449.9161 -10774.09  3.5194795  2.44130290 -0.7146951 -0.1898074 -0.6714178 -0.8653561       Treg

    # Check the assigment of each cluster correspond to the expected cell type given the known markers
    plot6 <- ggplot(melt(data[,3:9])) +
      geom_boxplot(aes(cluster,value)) +
      facet_wrap(~variable,scales="free")

    ## Using cluster as id variables

    plot6

![](/Users/ysanchez/Documents/Projects-analysis/FOG_2025_Introduction_Spatial_Analysis/notebooks-html/FOG_2025_Introduction_to_Spatial_Analysis-2025-05-16_files/figure-markdown_strict/boxplots-markers-1.png)

*Are the markers corresponding to the expected cell types?*

    #  We can also visualise the assigments as a heatmap

    plot7 <- ggplot(melt(kmeans$centers)) +
      geom_tile(aes(Var1,Var2,fill=value)) +
      scale_fill_gradient2(low="#0072b2",mid="white",high="#d55e00",midpoint=0) +
      labs(x="cluster",y="marker",fill="intensity",title="heatmap of cluster centres")

    plot7

![](/Users/ysanchez/Documents/Projects-analysis/FOG_2025_Introduction_Spatial_Analysis/notebooks-html/FOG_2025_Introduction_to_Spatial_Analysis-2025-05-16_files/figure-markdown_strict/heatmap-of-clusters-1.png)

# Visualisation - dimensionality reduction

The number of features possible to profile in a spatial experiment is
increasing. We often apply a dimensionality reduction step to visalise
the data. Principal Component Analysis (PCA) and Uniform Manifold
Approximation and Projection (UMAP) are two common visualisation
techniques. Here we are using the later.

    # umap: get the data. 
    umap = cbind(data,as.data.frame(umap(data[,3:8],seed=123)))
    # # umap: plot (no colour)
    plot_umap <- ggplot(umap) +
      geom_point(aes(V1,V2)) +
      labs(x="umap_1",y="umap_2",title = "UMAP")

    # You can also coloured by something interesting. For instance, here we plot (coloured by panCK)
    plot_umap_feature <-   ggplot(umap) +
        geom_point(aes(V1,V2,colour=panCK)) +
        scale_colour_gradient(low="blue",high="red") +
        labs(x="umap_1",y="umap_2",title="UMAP (pancK)")

    plot_umap +  plot_umap_feature

![](/Users/ysanchez/Documents/Projects-analysis/FOG_2025_Introduction_Spatial_Analysis/notebooks-html/FOG_2025_Introduction_to_Spatial_Analysis-2025-05-16_files/figure-markdown_strict/umap-1.png)

    # overlay clusters on the umap and tissue
    plot_umap_celltypes <- ggplot(cbind(data,as.data.frame(umap(scale(data[,3:8]),seed=123)))) +
      geom_point(aes(V1,V2,colour=cluster)) + ggtitle("UMAP coloured by cell types")

    plot_tissue_celltypes <- ggplot(data,aes(CellX,CellY,colour=cluster)) +
      geom_point() + ggtitle("Tissue coloured by cell types")

    # visualise the plots
    plot_umap_celltypes + plot_tissue_celltypes

![](/Users/ysanchez/Documents/Projects-analysis/FOG_2025_Introduction_Spatial_Analysis/notebooks-html/FOG_2025_Introduction_to_Spatial_Analysis-2025-05-16_files/figure-markdown_strict/umap-tissue-by-celltype-1.png)

# Neighborhood analysis

A typical analysis used in spatial biology is to quantify whether a cell
type is often found next to another. Here, we use the nn2 package.

     # For every cell, return the row indices of its nearest two neighbours. The data[,1:2] corresponds to the CellX and CellY cell coordinates.
      nebs = nn2(data[,1:2],k=2)[[1]]
      # discard column 1 (these will almost always be the cell itself)
      nebs = nebs[,2]
      # supply these row indices to the dataset to return the phenotypes of the neighbours stored in the column called "cluster"
      nebs = data[nebs,"cluster"]
      # now we make a table where the rows are the cells and the columns are the phenotypes of their nearest neighbours
      nebs = as.matrix(table(data$cluster,nebs))
      # convert to proportion
      nebs = nebs/rowSums(nebs)
      # blank out the middle diagonal (since most cells will neighbour themselves - this information is not informative)
      diag(nebs) = NA

    # now we can melt and plot
    plot_nn <-  ggplot(melt(nebs)) +
        geom_tile(aes(Var1,nebs,fill=value)) +
        scale_fill_gradient(low="white",high="red") +
        labs(x="phenotype of cell",y="phenotype of nearest neighbour",fill="proportion",title="nearest neighbours")

    plot_nn

![](/Users/ysanchez/Documents/Projects-analysis/FOG_2025_Introduction_Spatial_Analysis/notebooks-html/FOG_2025_Introduction_to_Spatial_Analysis-2025-05-16_files/figure-markdown_strict/nearest-neigh-plot-1.png)
\# Change of cell type composition in space

Another common task in spatial biology is to know how a particular cell
type(s) abundance change as we move away from a boundary. Here we use
the epithelium “Epi” as such as boundary and check the changes of the
other cell types.

    # split data in to epithelial cells and non-epithelial cells
      epi = data[data$cluster == "Epi",]
      nonepi = data[data$cluster != "Epi",]
      
      # for every non-epithalial cell, get the distance to the nearest epithelial cell
      dist = nn2(epi,nonepi,k = 1)
      dist = dist[[2]]
      
      # add this information to the dataset to non-epithelial cells
      data$distance_to_nearest_epi = NA
      data$distance_to_nearest_epi[data$cluster != "Epi"] = dist

    # boxplot of distances. What do you observe?
     boxplot_distances <-  ggplot(data[data$cluster != "Epi",],aes(cluster,distance_to_nearest_epi)) +
        geom_boxplot()

     boxplot_distances 

![](/Users/ysanchez/Documents/Projects-analysis/FOG_2025_Introduction_Spatial_Analysis/notebooks-html/FOG_2025_Introduction_to_Spatial_Analysis-2025-05-16_files/figure-markdown_strict/distance-to-nn-1.png)
<br> The boxplot above is informative for median distance to Epi but how
do the numbers of each cell phenotype change as we move away from Epi?
<br>

      # bin distances in to 5 bins by 20% each. This allows to have discrete values, perhaps easier to see a trend.
      data$distance_to_nearest_epi = cut(data$distance_to_nearest_epi,quantile(data$distance_to_nearest_epi,seq(0,1,1/5),na.rm=T),include.lowest = T)
      # for every cell phenotype, how many of those cells are in particular bins?
      
      table = as.matrix(table(data$cluster,data$distance_to_nearest_epi))
      table  

    ##             
    ##              [7.03,16] (16,34.3] (34.3,75.3] (75.3,123] (123,329]
    ##   Epi                0         0           0          0         0
    ##   Macrophage        51       113          85         26         7
    ##   Treg               2        22         159        168       128
    ##   CAF              200       117           8         58       117

    # make the plot
    plot_distances <- ggplot(melt(table)) +
        geom_point(aes(Var2,value,colour=factor(Var1))) +
        geom_line(aes(as.numeric(Var2),value,colour=factor(Var1)))  +
        scale_colour_manual(values = c("transparent","#7CAE00","#00BFC4","#C77CFF")) +
        labs(x="distance from epithelium",y="number of cells",title="epithelium spatial relationships")


    plot_distances

![](/Users/ysanchez/Documents/Projects-analysis/FOG_2025_Introduction_Spatial_Analysis/notebooks-html/FOG_2025_Introduction_to_Spatial_Analysis-2025-05-16_files/figure-markdown_strict/epithelial-spatial-relationships-1.png)

*What can you observe?*

# Conclusion

Here we have confirmed what we thought after looking at the physical
tissue: different cell phenotypes have different distance relationships
with tumour epithelium. <br>

1.  Tregs (CD4 and FOXP3) are epithelium distal <br>
2.  Macrophages (CD68 and PDL1) are epithelium proximal <br>
3.  CAFs (aSMA) are epithelium adjacent <br>

<!-- -->

    sessionInfo()

    ## R version 4.3.2 (2023-10-31)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS Sonoma 14.3
    ## 
    ## Matrix products: default
    ## BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: Europe/London
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] Seurat_5.0.1       SeuratObject_5.0.1 sp_2.1-3           patchwork_1.2.0    RANN_2.6.1         uwot_0.1.16        Matrix_1.6-5      
    ##  [8] ggplot2_3.5.2      reshape2_1.4.4     rmarkdown_2.25    
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] deldir_2.0-2           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.6            magrittr_2.0.3        
    ##   [6] RcppAnnoy_0.0.22       spatstat.geom_3.2-8    matrixStats_1.2.0      ggridges_0.5.6         compiler_4.3.2        
    ##  [11] png_0.1-8              vctrs_0.6.5            stringr_1.5.1          pkgconfig_2.0.3        fastmap_1.1.1         
    ##  [16] ellipsis_0.3.2         labeling_0.4.3         utf8_1.2.4             promises_1.2.1         purrr_1.0.4           
    ##  [21] xfun_0.42              jsonlite_1.8.8         goftest_1.2-3          highr_0.10             later_1.3.2           
    ##  [26] spatstat.utils_3.0-4   irlba_2.3.5.1          parallel_4.3.2         cluster_2.1.6          R6_2.5.1              
    ##  [31] ica_1.0-3              stringi_1.8.4          RColorBrewer_1.1-3     spatstat.data_3.0-4    reticulate_1.35.0     
    ##  [36] parallelly_1.36.0      lmtest_0.9-40          scattermore_1.2        Rcpp_1.0.13            knitr_1.45            
    ##  [41] tensor_1.5             future.apply_1.11.1    zoo_1.8-12             sctransform_0.4.1      FNN_1.1.4             
    ##  [46] httpuv_1.6.14          splines_4.3.2          igraph_2.0.1.1         tidyselect_1.2.1       yaml_2.3.8            
    ##  [51] abind_1.4-5            rstudioapi_0.17.1      spatstat.random_3.2-2  codetools_0.2-19       miniUI_0.1.1.1        
    ##  [56] spatstat.explore_3.2-6 listenv_0.9.1          lattice_0.22-5         tibble_3.2.1           plyr_1.8.9            
    ##  [61] shiny_1.8.0            withr_3.0.1            ROCR_1.0-11            evaluate_0.23          Rtsne_0.17            
    ##  [66] future_1.33.1          fastDummies_1.7.3      survival_3.5-7         polyclip_1.10-6        fitdistrplus_1.1-11   
    ##  [71] pillar_1.9.0           KernSmooth_2.23-22     plotly_4.10.4          generics_0.1.3         RcppHNSW_0.6.0        
    ##  [76] munsell_0.5.0          scales_1.3.0           globals_0.16.2         xtable_1.8-4           glue_1.7.0            
    ##  [81] lazyeval_0.2.2         tools_4.3.2            data.table_1.16.0      RSpectra_0.16-1        leiden_0.4.3.1        
    ##  [86] dotCall64_1.1-1        cowplot_1.1.3          grid_4.3.2             tidyr_1.3.1            colorspace_2.1-0      
    ##  [91] nlme_3.1-164           cli_3.6.5              spatstat.sparse_3.0-3  spam_2.10-0            fansi_1.0.6           
    ##  [96] viridisLite_0.4.2      dplyr_1.1.4            gtable_0.3.4           digest_0.6.34          progressr_0.14.0      
    ## [101] ggrepel_0.9.5          htmlwidgets_1.6.4      farver_2.1.1           htmltools_0.5.7        lifecycle_1.0.4       
    ## [106] httr_1.4.7             mime_0.12              MASS_7.3-60.0.1
