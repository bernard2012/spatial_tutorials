## Virtual slicing in 3D Data

While majority of the current spatial omics technilogies are 2D in nature, real 3D technologies are emerging. In 2018, Wang et al. created a 3D spatial expression dataset consisting of 28 genes from 32,845 single cells in a visual cortex volume using the STARmap technology. 

To better understand 3D data and enable smooth swtich between 3D and 2D views, we created a set of functions in Giotto to generate virtual 2D slices (cross section) in 3D data. 

* * *

<br>

### Usage

#### Main steps

The main steps are:

- Create neighborhood network
- Choose a method of defining a 2D plane and identify the corresponding parameters
- Call CrossSection routine to generate the slices and visualize them

* * *

<br>

#### 1. Spatial network

To define spatial relationship between cells, we first create a spatial network using createSpatialNetwork function.

##### createSpatialNetwork function
```R
STAR_test <- createSpatialNetwork(gobject = STAR_test, dimensions = 'all', method='Delaunay',delaunay_method = 'delaunayn_geometry')
```

Here the spatial network is created using the Delaunay triangulation method. For 3D data, the option dimensions='all' will select all 3 dimensions to create the spatial network. The option `delaunay_method = 'delaunayn_geometry'`, calls the `"delaunay"` function from the `"geometry"` package to do the 3D triangulation. 

* * *

<br>

#### 2. Define a 2D plane with the right parameters

Giotto offers four different ways of defining the location and orientation of a 2D slice, i.e. `"equation"`, `"3 points"`, `"point and norm vector"` and `"point and two plane vectors"`. 

##### createCrossSection function

```R
createCrossSection <- function(gobject,name="cross_section",spatial_network_name = "Delaunay_network", thickness_unit = c("cell","natural"), slice_thickness = 2,cell_distance_estimate_method = "mean", extend_ratio = 0.2, method=c("equation","3 points","point and norm vector","point and two plane vectors"), equation=NULL, point1=NULL,point2=NULL,point3=NULL, normVector=NULL, planeVector1=NULL,planeVector2=NULL, mesh_grid_n = 20, return_gobject = TRUE)
```

##### The `method` parameter
- `method="equation"` choice:

Requires the parameter `equation` to be set. For example,  `equation=c(0,1,0,600)`. This command will create a 2D slice specified by the equation `0 * X + 1 * Y + 0 * Z = 600` (or simply `Y = 600`). In general, if the `"equation"` parameter is set to be `c("A","B","C","D")`, the plane defined by `Ax+By+Cz=D` will be created. This is the default option and is also the most flexible one to use.  

- `method="3 points"` choice:

Requires the paramters `point1`, `point2`, `point3` to be set. For example: `point1 = c(0,600,0), point2 = c(1,600,0), point3 = c(0,600,1)`. This command will create a 2D slice uniquely specified by three points `"point1"`, `"point2"` and `"point3"`. This option will be helpful if the equation of the slice plane is not easy to determine but the coordinates of a few anchor points for the slice in the tissue sample is available (or can be easily estimated). 

- `method="point and norm vector"` choice:

Requires the parameter `point1` and `normVector` to be set. For example: `point1 = c(0,600,0), normVector = c(0,1,0)`. This command will create a 2D slice uniquely specified by the point `"point1"` and a vector that is perpendicular (or "norm") to the plane. This option will be helpful if the orientation (direction) of the slice plane is predefined and one anchor point is available. 

- `methd="point and two plane vectors"` choice:

Requires the parameters `point1`, `planeVector1` and `planeVector2` to be set. For example: `point1 = c(0,600,0), planeVector1 = c(1,0,0), planeVector2 = c(0,0,1)`. This command will create a 2D slice uniquely specified by the point `"point1"` and two vectors that are parallel to the slice plane and not parallel with each other. This option will be helpful if the required parameters for the other options are not available. 

##### Other important parameters for defining a virtual slice

| param | explanations |
| -------- | ---------------- |
| `name = "cross_section"` | The default name for the virtual slice to be created. The user can pass customized names one at each time if multiple slices needs to be created. If the same name has been used already, the new result will overwrite the old one. |
| `spatial_network_name = "Delaunay_network"` | The default spatial network that is based on when estimating cell sizes. |
| `cell_distance_estimate_method = "mean"` | This parameter specifies how the average cell-cell distance in the 3D data is estimated. It can be "mean" (default) or "median" value across all the detected distances between neighboring cells. The resulting value will serve as an estimate of the average diameter (size) of the cells observed. |
| `thickness_unit = c("cell", "natural")` | The unit of the virtual section thickness. If "cell" (default), the average size of the observed cells is used as length unit. If "natural", the unit of cell location coordinates is used. |
| `slice_thickness = 2` | Thickness of slice. The default value is 2. If meanwhile the "thickness_unit" is specified as "cell" (default), a virtual slice with 2 layers of cells will be created. Alternatively, if the "slice_thickness" is specified as "10" and the "thickness_unit" is specified as "natural", a virtual slice with thickness of 10 units (e.g. um) will be created. |
| `extend_ratio = 0.2` and `mesh_grid_n = 20` | These two parameters specifies the range and density of the meshgrid when visualizing the virtual slice in 3D. "extend_ratio" specifies how much the slice mesh grid will extend outside the intersection between the slice plane and the tissue sample. "mesh_grid_n" specifies the total number of grid lines along each direction. |

* * *

<br>

#### 3. Call CrossSection routine

Here we walk through an example. First the preparation steps.

```R
# from STARmap dataset
## create consistent color code
mynames = unique(pDataDT(STAR_test)$cell_types)
mycolorcode = Giotto:::getDistinctColors(n = length(mynames))
names(mycolorcode) = mynames
## creat spatial network 
STAR_test <- createSpatialNetwork(gobject = STAR_test)
```

Next create the cross section.
```R
## creat virtual slice
STAR_test = createCrossSection(STAR_test,method="equation",equation=c(0,1,0,600),extend_ratio = 0.6)
```

Then create the cross section plots.
```R
## view the location and orientation of virtual slice in 3D (cell type)
insertCrossSectionSpatPlot3D(STAR_test, cell_color = 'cell_types', axis_scale = 'cube', point_size = 2, cell_color_code = mycolorcode)
## view the location and orientation of virtual slice in 3D (Gene expression)
insertCrossSectionGenePlot3D(STAR_test, expression_values = 'scaled',axis_scale = "cube", genes = "Slc17a7")
## view the cells captured by virtual slice (3D, cell type)
crossSectionPlot3D(STAR_test,point_size = 2, cell_color = "cell_types",  cell_color_code = mycolorcode,axis_scale = "cube")
## view the cells captured by virtual slice (2D, cell type)
crossSectionPlot(STAR_test, point_size = 2, cell_color = "cell_types",cell_color_code = mycolorcode)
## view the cells captured by virtual slice (3D, gene expression)
crossSectionGenePlot3D(STAR_test, point_size = 2,genes = c("Slc17a7"), expression_values = 'scaled')
## view the cells captured by virtual slice (2D, gene expression)
crossSectionGenePlot(STAR_test, genes = "Slc17a7",point_size = 2,cow_n_col = 1.5, expression_values = 'scaled')
```
