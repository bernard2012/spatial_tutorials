## Cell type interaction enrichment


##### cellProximityEnrichment function

```R
cellProximityEnrichment <- function(gobject, spatial_network_name = 'Delaunay_network', cluster_column,number_of_simulations = 1000,adjust_method = c("none", "fdr", "bonferroni","BH","holm", "hochberg", "hommel","BY")) {
```

This function computes the cell-type to cell-type interactions. It outputs the number of interactions (i.e. interaction frequency) per cell-type pair, and the depletion/enrichment score per cell-type pair, to be defined as the observed over expected interactions.

The input is a spatial network defined by `spatial_network_name`. The `cluster_column` refers to the column in the cell metadata table containing the cell types to use. The `number_of_simulations` is the number of times the spatial network is simulated. In each simulation, the labels of the cells are reshuffled in the spatial network and the interactions are recounted. `adjust_method` uses various ways to correct for multiple hyothesis testing. 

Returns a list of cell proximity scores in `data.table` format. The first
(`raw_sim_table`) shows the raw observations of both the original and simulated networks. The second (`enrichm_res`) shows the enrichment results.

It is recommended to use the functions `cellProximityBarplot`, `cellProximityHeatmap`, `cellProximityNetwork` to visualize the enrichment/depletion result.

###### Examples
```R
#seqFISH+ dataset
cell_proximities = cellProximityEnrichment(gobject = VC_test, cluster_column = 'cell_types',spatial_network_name = 'Delaunay_network',adjust_method = 'fdr',number_of_simulations = 2000)
#barplot
cellProximityBarplot(gobject = VC_test, CPscore = cell_proximities, min_orig_ints = 5, min_sim_ints = 5)
#network
cellProximityNetwork(gobject = VC_test, CPscore = cell_proximities, remove_self_edges = T, only_show_enrichment_edges = T)
```

* * *

<br>

##### findCPG function

```R
findCPG = function(gobject, expression_values = 'normalized',selected_genes = NULL, cluster_column, spatial_network_name = 'Delaunay_network',minimum_unique_cells = 1, minimum_unique_int_cells = 1, diff_test = c('permutation', 'limma', 't.test', 'wilcox'), mean_method = c('arithmic', 'geometric'), offset = 0.1,adjust_method = c("bonferroni","BH", "holm", "hochberg", "hommel","BY", "fdr", "none"),nr_permutations = 100,exclude_selected_cells_from_test = T,do_parallel = TRUE, cores = NA)
```
`findCPG` identifies genes that are differentially expressed due to proximity to other cell types. CPG stands for cell proximity genes. This analysis is particularly meaningful for whole transcriptome dataset such as seqFISH+.

The options `cluster_column`, `spatial_network_name`, `adjust_method` are explained previously (see `cellProximityEnrichment()`).

| param | explanations |
| -------- | ---------------- |
| `minimum_unique_cells` | minimum number of target cells required |
| `minimum_unique_int_cells` | minimum number of interacting cells required |
| `diff_test` | which differential statistical test to use, among "permutation", "limma", "t.test", or "wilcox" |
| `nr_permutations` | applies if `diff_test` is set to "permutation" |
| `mean_method` | whether to use arithmetic mean or geometric mean |
| `offset` | for calculating log2 ratio |
| `do_parallel` | whether or not to run in parallel using mcapply |
| `cores` | number of cores when `do_parallel=TRUE` |


Returns a `data.table` with the following information:

1. genes (tested genes)
2. target cell type
3. interacting cell type
4. number of cells for target cell type
5. number of cells for interacting cell type
6. average expression in interacting cells from the target cell type
7. average expression in the non-interacting cells from the target cell type
8. log2 fold change between (6) and (7)
9. difference of expression between (6) and (7)
10. pvalue 
11. p.adj

###### Examples
```R
gene_metadata = fDataDT(VC_test)
high_expressed_genes = gene_metadata[mean_expr_det > 1.31]$gene_ID
CPGscoresHighGenes =  findCPG(gobject = VC_test,selected_genes = high_expressed_genes,spatial_network_name = 'Delaunay_network', cluster_column = 'cell_types',diff_test = 'permutation',adjust_method = 'fdr',nr_permutations = 2000, do_parallel = T, cores = 2)
plotCellProximityGenes(VC_test, cpgObject = CPGscoresHighGenes, method = 'dotplot')
```

