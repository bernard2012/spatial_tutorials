## Dimensional reduction, clustering, and marker gene identification

After normalization of the dataset by library size and adjustment of gene expression for confounding factors, the next step is the selection of highly variable genes, dimensional reduction, and clustering.

<br>

### Usage

#### 1. Highly variable gene selection (HVG)

##### calculateHVG function
```R
calculateHVG <- function(gobject, expression_values = c('normalized', 'scaled', 'custom'), method = c('cov_groups','cov_loess'), reverse_log_scale = FALSE, logbase = 2, expression_threshold = 0, nr_expression_groups = 20, zscore_threshold = 1.5,HVGname = 'hvg',difference_in_cov = 0.1, show_plot = NA,return_plot = NA, save_plot = NA, save_param = list(),default_save_name = 'HVGplot', return_gobject = TRUE)
```

By default, we select the `expression_values="normalized"` to use the version of expression matrix normalized by library size, in log-normalized counts. 

There are two ways to perform HVG, as indicated by `method`:
- `method="cov_groups"` is covariance group. Genes are binned into expression groups (`nr_expression_groups`) based on the average expression of each gene. A zscore is computed for each expression group. Genes that exceed the zscore threshold (`zscore_threshold`) within each group are marked as HVG.
- `method="cov_loess"` is covariance from loess regression. A predicted COV is calculated for each gene using loess regression (COV ~ log(mean expression)). Genes that show a higher than predicted COV by a difference equal to `difference_in_cov` are considered highly variable (defined by `difference_in_cov`). 

This returns the HVG genes in the slot named `"hvg"` within the Giotto object.

* * *

<br>

#### 2. Dimensional reduction

The next main step is principle component analysis (PCA), achieved by `runPCA()` function.

##### runPCA function
```R
runPCA <- function(gobject, expression_values = c('normalized', 'scaled', 'custom'),reduction = c('cells', 'genes'), name = 'pca', genes_to_use = 'hvg', return_gobject = TRUE, center = F, scale_unit = F, ncp = 100, method = c('irlba','factominer'), rev = FALSE, verbose = TRUE, ...)
```
Depending on dataset size, one can use `irlba` package (`method="irlba"`) to perform PCA on very large dataset (such as proteomics data), or use the `factominer` package (`method="factominer"`) for small to medium size dataset. The `irlba` package uses approximate singular value decomposition (SVD) to speed up PCA.
The `reduction` is by default in cell-space. The genes that will be utilized for PCA analysis are highly variable genes (defined earlier by `calculateHVG()` function). The `ncp` is the number of principle components. 
Note the options `center` and `scale_unit`. It is recommended to perform centering and scaling the dataset before performing PCA.

Note that if `irlba` is chosen, it is absolutely necessary to center (set `center=TRUE`) before PCA. `irlba` package uses SVD to perform PCA, and if not centered, the result of noncentered PCA can be misleading. 

##### Reference:
- https://stats.stackexchange.com/questions/189822/how-does-centering-make-a-difference-in-pca-for-svd-and-eigen-decomposition

<br>

#### 2.1 tSNE or UMAP

One can use tSNE or UMAP to next visualize the dataset in 2D space.

##### runtSNE, runUMAP functions
```R
runtSNE <- function(gobject, expression_values = c('normalized', 'scaled', 'custom'),reduction = c('cells', 'genes'), dim_reduction_to_use = 'pca', dim_reduction_name = 'pca', dimensions_to_use = 1:10, name = 'tsne', genes_to_use = NULL, return_gobject = TRUE, dims = 2, perplexity = 30, theta = 0.5,do_PCA_first = F, set_seed = T, seed_number = 1234, verbose = TRUE)
                    
runUMAP <- function(gobject, expression_values = c('normalized', 'scaled', 'custom'), reduction = c('cells', 'genes'), dim_reduction_to_use = 'pca', dim_reduction_name = 'pca',dimensions_to_use = 1:10,name = 'umap',genes_to_use = NULL,return_gobject = TRUE,n_neighbors = 40, n_components = 2, n_epochs = 400, min_dist = 0.01, n_threads = 1, spread = 5, set_seed = T, seed_number = 1234, verbose = T,...) 
```

tSNE and UMAP can accept as input the original gene expression matrix (set `dim_reduction_to_use=NULL`) or the dimension reduced matrix from PCA (default) (`dim_reduction_to_use="pca"`).
If principle components are analyzed, then one specifies the `dimensions_to_use`. 

- For tSNE, one can further define `perplexity` and `theta` options (see tSNE guide here: https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html).
- For UMAP, one can further define the number of neighbors (`n_neighbors`), number of epochs (`n_epochs`), and `min_dist` (see UMAP parameter guide here: https://umap-learn.readthedocs.io/en/latest/parameters.html).

The options `set_seed` and `seed_number` are helpful to fix the random number generation seed so that the same result is returned each time the function is run.

A plot will be returned in the result.

* * *

<br>

#### 3. Clustering

##### Nearest neighbor network

Besides K-means clustering and hiercharchical clustering, Giotto contains two commonly used single-cell RNAseq based methods that are also suitable for spatial transcriptomic datasets: Leiden clustering (default and recommended) and Louvain clustering. There are two steps:

- Construct shared nearest neighbor network. 
- Perform community detection on nearest neighborhood network.

<br>

###### createNearestNetwork function
```R
createNearestNetwork <- function(gobject, type = c('sNN', 'kNN'), dim_reduction_to_use = 'pca', dim_reduction_name = 'pca', dimensions_to_use = 1:10,genes_to_use = NULL, expression_values = c('normalized', 'scaled', 'custom'), name = 'sNN.pca', return_gobject = TRUE, k = 30, minimum_shared = 5, top_shared = 3,verbose = T, ...)
```

Two modes are available:
- `type = sNN` (shared nearest neighbor network)
- `type = kNN` (K-nearest neighbor network)

One can perform clustering on gene expression matrix (by setting `dim_reduction_to_use = NULL`) or use dimension reduced matrix (`dim_reduction_to_use="pca"`). For latter, number of PCA components is specified (`dimensions_to_use`). 

For sNN, one defines the minimum number of shared neighbors needed (`minimum_shared`) and the number of top shared neighbors to keep per cell (`top_shared`). The `top_shared` acts as a filter on the sNN network. (see dbscan sNN).

For kNN, one defines the `k` the number of neighbors. (see dbscan kNN).

Returns the network in the `"sNN.pca"` slot of the Giotto object (slot name is defined by `name`).

<br>

##### Clustering (Leiden)

After creating nearest network, one can next perform either leiden or louvain clustering.

##### doLeidenCluster function
```R
doLeidenCluster = function(gobject, name = 'leiden_clus', nn_network_to_use = 'sNN', network_name = 'sNN.pca', python_path = NULL, resolution = 1, weight_col = 'weight', partition_type = c('RBConfigurationVertexPartition', 'ModularityVertexPartition'), init_membership = NULL, n_iterations = 1000, return_gobject = TRUE, set_seed = T, seed_number = 1234)
```

Slot name is `"leiden_clus"`. The `nn_network_to_use` can be either `"sNN"` or `"kNN"`. The `resolution` tunes the number of clusters to be detected. Usually this is set to a value from 0.5 to 1.0. The smaller the resolution, the fewer and larger the clusters are. When it is high (greater than 0.9), clusters become very refined. One can then display the result using `plotUMAP` to visualize the clusters. 

The edge weight (the `weight_col`) is set at 1/(1 + distance). The number of iterations (`n_iterations`) of Leiden is set at 1000. If the number of iterations is negative, the Leiden algorithm is run until an iteration in which there was no improvement (warning: this can take a long time in some cases). The `partition_type` defines the method of optimization:

- `"RBConfigurationVertexPartition"` - Implements Reichardt and Bornholdt’s Potts model. This quality function uses a linear resolution parameter
- `"ModularityVertexPartition"` - Implements modularity.

Since this function calls the `leidenalg` python package, the `python_path` defines where to find python binary (e.g. `/usr/bin/python`). However, if `python_path=NULL`, it will read the python path from the Giotto instructions. 


###### References:
- https://github.com/vtraag/leidenalg
- https://leidenalg.readthedocs.io/en/stable/index.html

<br>

##### Clustering (Louvain)

##### doLouvainCluster function
```R
doLouvainCluster_community <- function(gobject, name = 'louvain_clus',nn_network_to_use = 'sNN',network_name = 'sNN.pca', python_path = NULL, resolution = 1, weight_col = NULL, louv_random = F, return_gobject = TRUE, set_seed = F,seed_number = 1234)
```
Slot name is `"louvain_clus"`. The `resolution` tunes the number of clusters to be detected. Range of resolution is 0.5 - 1.0. A smaller resolution results in fewer and larger the clusters.
The `nn_network_to_use` is either `"sNN"` or `"kNN"`. The `weight_col` specifies the edge weight of sNN graph.

The `python_path` defines where to find python binary (e.g. `/usr/bin/python`). If `python_path=NULL`, it will read the python path from the Giotto instructions. 

This function is a wrapper for the Louvain algorithm implemented in Python.

###### Reference:
- https://python-louvain.readthedocs.io/en/latest/index.html

* * *

<br>

#### 4. Finding differentially expressed genes across clusters

Three choices are presented Gini, Scran or MAST. For more information and a comparison between the methods, see our paper. 

##### findMarkers_one_vs_all function
```R
findMarkers_one_vs_all <- function(gobject, expression_values = c('normalized', 'scaled', 'custom'), cluster_column, subset_clusters = NULL, method = c('scran','gini','mast'), pval = 0.01, logFC = 0.5, min_genes = 10, min_expr_gini_score = 0.5,  min_det_gini_score = 0.5, detection_threshold = 0, rank_score = 1, adjust_columns = NULL, verbose = TRUE, ...) 
```

This is one wrapper function for Gini, Scran, and MAST. Scran is fast and efficient, but Gini is able to better capture differenteial expressed genes for rare clusters. 

- For Gini (`method="gini"`), the relevant options are:
  - `min_expr_gini_score`: minimum gini coefficient on expression
  - `min_det_gini_score`: minimum gini coefficient on detection
  - `detection_threshold`: threshold for gene expression detection (any expression that is greater than threshold is considered detected"
  - `min_genes`: minimum number of genes to generate per cluster
  - `rank_score`: rank scores threshold for genes to include

  This will calculate the Gini-coefficients based on average expression, and detection fraction, and combine the two scores to identify marker genes for selected clusters. See `findGiniMarkers()` function for more details. Afterwards, it will next filter genes based on `min_expr_gini_score` and `min_det_gini_score`.

- For MAST (`method="mast"`), the options that are relevant are:
  - `adjust_columns`: column in `pDataDT` to adjust for, such as detection rate (default `NULL`)
  - `pval`: minimum p-value threshold
  - `logFC`: log-fold change threshold
  - `min_genes`: minimum number of genes to generate per cluster, overwrites pval and logFC

- For SCRAN (`method="scran"`), the following options apply:
  - `pval`: minimum p-value threshold
  - `logFC`: log-fold change threshold
  - `min_genes`: minimum number of genes to generate per cluster, overwrites pval and logFC

All methods require `cluster_column`, which is the cluster annotations to generate differentially expressed genes.

