## Guide to perform cell type enrichment analysis

Giotto implemented three algorithms for enrichment analysis of lower resolution of spatially expression datasets. It contains PAGE, RANK and hypergeometric. The aim of enrichment analysis is to use continuous values to represent the likelihood of the presence of a cell type of interest in specific spatial locations which may contain multiple cells.

* * *

<br>

### Introduction

#### PAGE
The method uses Parametric Analysis of Gene Set Enrichment (PAGE) method to evaluate cell type enrichment for each spatial location. Signature genes of interested cell types are used for enrichment analysis. Enrichment score was calculated based on signature gene expression for spatial locations. P-value could also be calculated via permutation test.

#### RANK
The method uses a rank method to calculate enrichment score for interested cell types. Single cell expression matrix as well as cell type labels are used for rank analysis. Rather than PAGE, RANK does not need signature gene selection. Based on the gene expression pattern of single cell RNA-seq, RANK could evaluate the cell type presence of spatial locations. P-value could also be calculated via permutation test.

#### Hypergeomatric
This method uses hypergeometric distribution test to evaluate cell type distribution of spatial locations based on signature genes of interested cell types. Enrichment score was calculated as -log10(p-value) of hypergeometric distribution test.

* * *

<br>

### Usage

#### 1. PAGE method

##### The makeSignMatrixPAGE function
```R
makeSignMatrixPAGE = function(sign_names, sign_list) 
```
This function converts a list of signature genes (e.g. for cell types or processes) into a binary matrix format that can be used with the PAGE enrichment option. Each cell type or process should have a vector of cell-type or process specific genes. These vectors need to be combined into a list (sign_list). The names of the cell types or processes that are provided in the list need to be given (sign_names).

###### Outputs:
Signature matrix is a binary (0/1) matrix. Rows are genes. Columns are cell types.

Once the signature matrix is created. The next step is to run runPAGEEnrich function.

##### The runPAGEEnrich function
```R
runPAGEEnrich <- function(gobject, sign_matrix, expression_values = c('normalized', 'scaled', 'custom'),reverse_log_scale = FALSE, logbase = 2, output_enrichment = c('original', 'zscore'), p_value = FALSE, n_times = 1000, name = NULL, return_gobject = TRUE)
```

| param | explanations |
| -------- | ---------------- |
| `sign_matrix` | The scRNAseq cell type gene signature matrix generated from `makeSignMatrixPAGE`. It should be  0/1 signature matrix for each cell type. |
| `expression_values=c('normalized', "raw", 'scaled', 'custom')` | The form of gene expression matrix to use for the spatial dataset. We recommend using the `normalized` form which is in log(normalized counts + 1) |
| `reverse_log_scale=FALSE` | In regards to calculating the mean expression per gene across all spots, whether to apply logbase^expr_values first. Applies only if `expression_values="normalized"`. Recommend to set to FALSE. |
| `logbase = 2` | Used in conjunction with reverse_log_scale |
| `p_value=TRUE` and `n_times` | P value could be calculated by setting `p_value=TRUE` using permutation test. `n_times` is the parameter for permutation test for p value calculation. |


##### Outputs:
A data frame with the enrichment score with each spatial location and cell type. In addition, if `p_value= TRUE` was specified, another data frame with p-value could also reported.

* * *

<br>

#### 2. RANK method

##### The makeSignMatrixRank function

```R
makeSignMatrixRank <- function(sc_matrix, sc_cluster_ids,gobject = NULL, ties.method=c("random", "max"))
```

This function will make a rank-based cell type gene signature matrix based on the scRNAseq dataset. In this signature matrix, there are N vectors where N is number of cell types. For each cell type, the vector is a rank-list of genes according to some criterion (in this case, according to `log2(mean_expr+1)-log2(av_expr+1)` where mean_expr is the cell type expression average, av_expr is the all cells' expression average). Where two values are sharing the same rank and thus creating a tie, the `ties.method` is used to break ties. 

The `sc_matrix` should be the gene expression matrix (in raw form). The `sc_cluster_ids` is the cluster annotation column. `ties.method` is tie breaking method for assigning ranks in case of ties. 

###### Example
```R
rank_matrix=makeSignMatrixRank(sc_matrix=cere_rnaseq2@raw_exprs, sc_cluster_ids=pDataDT(cere_rnaseq2)$leiden, ties.method="random")
```

The next step is to call runRankEnrich() function.

##### The runRankEnrich function

```R
runRankEnrich <- function(gobject,sign_matrix,expression_values = c('normalized', "raw", 'scaled', 'custom'), reverse_log_scale = FALSE, logbase = 2,output_enrichment = c('original', 'zscore'),ties.method = c("random", "max"),p_value = FALSE, n_times = 1000,name = NULL, return_gobject = TRUE, rbp_p = 0.99, num_agg=100 )
```

| param | explanations |
| -------- | ---------------- |
| `sign_matrix` | The scRNAseq cell type gene signature matrix generated from `makeSignMatrixRank` |
| `expression_values=c('normalized', "raw", 'scaled', 'custom')` | The form of gene expression matrix to use for the spatial dataset. We recommend using the `normalized` form which is in log(normalized counts + 1) |
| `reverse_log_scale=FALSE` | In regards to calculating the mean expression per gene across all spots, whether to apply logbase^expr_values first. Applies only if `expression_values="normalized"`. Recommend to set to FALSE. |
| `logbase = 2` | Used in conjunction with reverse_log_scale |
| `ties.method = c("random", "max")` | Breaking ties when ranking genes per spot. `"random"` means to assign ranks randomly among tied values. `"max"` means to assign the maximum rank. Recommend to set to "random" |
| `rbp_p = 0.99` | Local weighting on rank. Used when calculating an enrichment score per spot. Set to a value 0 - 1. Recommended range is 0.95 - 0.995. Lower means more emphasis to place on top rank. (Advanced setting) | 
| `num_agg=100` | Number of top rank values to aggregate in computing enrichment score. Recommend to leave as default. (Advanced setting) |
| `p_value=TRUE` and `n_times` | P value could be calculated by setting `p_value=TRUE` using permutation test. `n_times` is the parameter for permutation test for p value calculation. |

###### Example:

```R
Slide_test<-runRankEnrich(Slide_test, sign_matrix=rank_matrix, expression_values="norm", reverse_log_scale=F, logbase=2, output_enrichment="original", name="rank", rbp_p=0.99, num_agg=100, ties.method="random")
```

###### Outputs:
A data frame with the enrichment score with each spatial location and cell type. In addition, if `p_value= TRUE` was specified, another data frame with p-value could also reported.

* * *

<br>

#### 3. The Hypergeometric method
```R
createSpatialEnrich = function(gobject, enrich_method =  ' hypergeometric’, sign_matrix, p_value = FALSE, n_times = 1000 …)
```

The input is `sign_matrix` which is 0/1 signature matrix for each cell type. P value could be calculated by setting `p_value=TRUE`.

###### Outputs:
A data frame with the enrichment score with each spatial location and cell type. In addition, if `p_value= TRUE` was specified, another data frame with p-value could also reported.

