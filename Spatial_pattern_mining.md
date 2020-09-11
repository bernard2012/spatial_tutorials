## Spatial pattern mining

Previously we have published a paper that describes the hidden markov random field (HMRF) method for identifying spatial patterns from spatial gene expression data.

Unlike typical clustering techniques, HMRF integrates both spatial locations and the gene expression to find groups of cells that are both spatially adjacent and which are similar to each other in expression space. 

The main steps are:

1. Create neighborhood network based on spatial positions
2. Define genes of interest
3. Call HMRF routine to find spatial clusters

<br>

### Usage

#### 1. Spatial network

To define spatial relationship between cells, we first create a spatial network using `createSpatialNetwork` function.

##### createSpatialNetwork function
```R
createSpatialNetwork <- function(gobject,name = NULL, dimensions = "all", method = c('Delaunay', 'kNN'), delaunay_method = c("deldir", "delaunayn_geometry", "RTriangle"), maximum_distance_delaunay = "auto", options = "Pp", Y = TRUE,j = TRUE,S = 0,minimum_k = 0,knn_method = "dbscan", k = 4,maximum_distance_knn = NULL, verbose = F, return_gobject = TRUE, ...)
```

Users can choose between two methods for the spatial network (Voronoi diagram by Delaunay triangulation or KNN network). 

- KNN network (`method="kNN"`) is simpler in concept and is faster to compute. The option `k` needs to be defined to specify the number of neighbors. When `maximum_distance_knn` is defined, no two cells are connected in the spatial network if they are further than `maximum_distance_knn` apart regardless of `k`. The `minimum_k` specifies how many minimum number of neighbors each cell will have in order to ensure that there are no singletons in the resulting network. For more details, see `createSpatialKNNnetwork`.
  - Note for 3D datasets, `dimensions="all"` will use all three dimensions to calculate the spatial distance. However it is sometimes desirable to use only x and y dimensions, but not z dimension. Pass `dimensions=c("sdimx", "sdimy")`. 
  - Note to create a network based on maximum distance only, one can set `k` to be a very high value (e.g. `k=100`)
  - Guidance for HMRF: we recommend approximately 5-10 neighbors per cell. 

- Voronoi network (`method="Delaunay"`).
  Voronoi diagrams are used to model a number of different biological structures, including cells and bone microarchitecture. It can be used to understand physical constraints in a tissue organization. 
  
  The default implementation of Delaunay triangulation is `"deldir"`. Besides the typical option for Delaunay, there are also a couple of other options to further filter the network. The `maximum_distance_delaunay` will cut-off the network based on a maximum distance. The `minimum_k` will ensure that there is a minimum number of neighbors to each cell. `Y`, `j`, and `S` parameters refer to the `RTriangle` settings (safely leave as default). The `dimension` option is also available for Voronoi. For details, see `createSpatialDelaunayNetwork()`.

* * *

<br>

#### 2. Define genes of interest

Genes of interest can be a list of spatial genes returned by a spatial gene detection algorithm (see tutorial).

* * *

<br>

#### 3. Call HMRF routine

##### doHMRF function
```R
doHMRF <- function(gobject, expression_values = c('normalized', 'scaled', 'custom'), spatial_network_name = 'Delaunay_network', spatial_genes = NULL, spatial_dimensions = c('sdimx', 'sdimy', 'sdimz'), dim_reduction_to_use = NULL, dim_reduction_name = 'pca', dimensions_to_use = 1:10, name = 'test', k = 10, betas = c(0, 2, 50), tolerance = 1e-10,zscore = c('none','rowcol', 'colrow'), numinit = 100, python_path = NULL,output_folder = NULL,overwrite_output = TRUE)
```

This is the function that performs the HMRF routine. 
- `expression_values`: the version of gene expression matrix to use. `"scaled"` is appropriate for HMRF as it is a row and column wise z-scored version of expression matrix.
- `spatial_network_name`: the version of spatial network to use. Either `Delaunay_network` or `kNN_network`.
- `dim_reduction_to_use`: whether to use dimension reduced expression matrix. Set it to `NULL` if we wish to rather use the original gene expression matrix. (Defaults to `NULL`)
- `dim_reduction_name`, `dimensions_to_use`: apply only if `dim_reduction_to_use` is set to `"pca"`
- `k`: number of spatial clusters
- `betas`: smoothness parameter with three values - (1) starting beta, (2) beta increment, (3) number of betas to run. E.g. `c(0, 2, 50)` means to run betas 50 times from beta=0, 2, 4, 6, 8, ..., to 100.
- `tolerance`: tolerance value
- `zscore`: set to `"none"` if the gene expression is z-scored already. Otherwise, perform row and col-wise z-scoring (`"rowcol"`). `"colrow"` means that the order is first col then row-wise zscoring. Note the longer dimension should be zscored first.
- `numinit`: number of initializations to try

A general guideline for choosing k is to use the gap-statistics by Tibshrani. The elbow point of the gap size vs. k plot usually indicates choices of k. 

To decide which beta one should use, the current state of algorithm repeats HMRF for a number of betas. Here are some guidelines:

- if the number of genes is from 10 to 50, the recommended range is 0 to 10 at increment of 0.5. E.g. `c(0,0.5,20)`
- if the number of genes is below 50, the recommended range is 0 to 15 at increment of 1. E.g. `c(0,1,15)`
- if the number of genes is between 50 to 100, the range is 0 to 50 at increment of 2. E.g. `c(0, 2, 25)`
- if the number of genes is between 100 and 500, the range is 0 to 100 at increment of 5. E.g. `c(0, 5, 20)`

Within the range of beta, we recommend selecting the best beta by the Bayes information criterion. This requires first performing randomization of spatial positions to generate the null distribution of log-likelihood scores for randomly distributed cells for the same range of betas. Then find the beta where the difference between the observed and the null log-likelihood is maximized. See the HMRF tutorial from the Chapter.

Returns a `list` type object with the following fields:
- `name`: name used
- `output_data`: output path
- `k`: k used,
- `betas`: betas used,
- `python_path`: python_path_used

After running HMRF, one can further view the HMRF spatial clusters (`viewHMRFresults2D`), and add HMRF annotation to `pDataDT` table (`addHMRF`).

##### Example:

```R
#from seqfish+ dataset
#100 spatial genes
HMRF_spatial_genes = doHMRF(gobject = VC_test, expression_values = 'scaled', spatial_genes = my_spatial_genes, spatial_network_name = 'Delaunay_network', k = 9, betas = c(0,1,50), output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_top100_k9_scaled'))
## view results of HMRF
for(i in seq(0, 50, by = 1)) {
  viewHMRFresults2D(gobject = VC_test,HMRFoutput = HMRF_spatial_genes,k = 9, betas_to_view = i,point_size = 2)
}
## add HMRF of interest to giotto object
VC_test = addHMRF(gobject = VC_test,HMRFoutput = HMRF_spatial_genes, k = 9, betas_to_add = c(28), hmrf_name = 'HMRF_2')

```



