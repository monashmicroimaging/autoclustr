# Introduction
This repository contains a script (`autoclustr.R`) that attempts to automate the cluster analysis techniques for single molecule localisation microscopy first described by Owen et al. (2010)<sup>1</sup>. This includes, but is not limited to, assessing overall clustering using Ripley’s K-function analysis, generation of cluster maps and segmentation of clusters to extract various cluster measurements and statistics.

<sup>1</sup>Owen, D. M., Rentero, C., Rossy, J., Magenau, A., Williamson, D., Rodriguez, M., & Gaus, K. (2010). PALM imaging and cluster analysis of protein heterogeneity at the cell surface. _Journal of Biophotonics_, __3__(7), 446–454. [http://doi.org/10.1002/jbio.200900089](http://doi.org/10.1002/jbio.200900089)

## Installation instructions

### Requirements
* R
* Matlab

You will need to install *R* and some *R* package dependencies to run this script successfully. You will also need to have Matlab installed, although this won’t covered here.

1. Install *R*. Instructions for installing *R* for your platform can be found on the [R-project website](https://www.r-project.org/).
2. [Download](https://github.com/monashmicroimaging/autoclustr/archive/v0.1.0.zip) (or [git](https://git-scm.com/) clone) this repository to your local computer.

3.  `cd` into the downloaded directory and start R.
```console
$ cd /path/to/downloaded/directory
$ R
```

4. Start R and install the `packrat` package (note: `packrat` may be able to bootstrap itself when you start R — in that case you can skip the remainder of this step):
```r
install.packages(“packrat”)
q()
```


The next time you start R in this directory, `packrat` will setup a local _R_ library and install all the dependencies required to run this script.

Note: If you want to run this script frequently or develop it further, it is recommended that you install the [RStudio](https://www.rstudio.com/) IDE and open the `autoclustr.Rproj` project.

## Usage instructions
There are a couple of important conventions and settings required to run the script. Once you have setup the data as recommend in the [Conventions](#conventions) section and set all the [script parameters](#scriptsettings) appropriately for your data, you can simply `source` the `autoclustr.R` script from the R console and the analysis should get underway:
```r
source(“./autoclustr.R”)
```

Note: this assumes you started R in the downloaded autoclustr project and `packrat` has setup all the necessary dependencies (see [Installation instructions](#installationinstructions)).

Note: If you are using Mac OS X, you will probably need to specify the Matlab executable path manually using the `matlab_path` parameter in the [script settings](#scriptsettings).

### Conventions
First you need to supply the path to a project directory of data to be analysed using the `project_dir` variable (See [Script settings](#scriptsettings) below). The project directory should contain 2 subdirectories:
* Localizations - contains all lozalization files (e.g. 'localization1.txt’). This script currently only accepts localisation files exported from [rapidSTORM](http://www.super-resolution.biozentrum.uni-wuerzburg.de/research_topics/rapidstorm/) (although it should be trivial to support other formats).
* ROIs - Regions of Interest (ROIs) corresponding to each localization file. ROI files are simply a tab-delimited plain txt file that list the x and y coordinates of vertices of a polygon that specifies a region of interest to analyse in the corresponding localization file. These should be named the same as the corresponding localization file (not including the file extension) with an `_roiID` identifier appended (e.g. 'localization1_ROI1.txt'). We typically create ROIs using reconstructed super-resolution images and the ROI tools in [ImageJ](https://imagej.nih.gov/ij/) and export them via `File → Save As → XY Coordinates…`. Note that one needs to set the `roi_xy_pix_size` attribute to the pixel size of the image on which the ROI was drawn.

The following is an example of the overall input folder structure:
```
project_dir
  | ---- Localizations
  |          | ---- localizations1.txt
  |          | ---- localizations2.txt
  |
  | ---- ROIs
		  | ---- localizations1_ROI1.txt
		  | ---- localizations1_ROI2.txt
		  | ---- localizations2_ROI1.txt
		  | ---- localizations2_ROI2.txt
		  | ---- localizations2_ROI3.txt
```

## Script settings
Settings are located at the bottom of `autoclustr.R` script.

| Parameter | Description |
| :-------- | :---------- |
| **Input settings** ||
| `project_dir` | Path to the project directory |
| `roi_xy_pix_size` | Pixel size calibration for the ROIs i.e., the x/y pixel size of the image on which the ROIs were drawn. |
| **Output settings** ||
| `output_dir_name` | Name for the output directory. This will be created inside the `project_dir` if it does not exist. |
| **Local L and cluster heat-map settings** ||
| `localL_radius` | Spatial scale (r value) used for the local L analysis. |
| `use_global_k_peak` | Use the mean peak r value from the global Ripley’s K-function analysis as the spatial scale (i.e., r value) for the local L analysis. If this is set to `TRUE`, `localL_radius` is ignored. Default is `FALSE`. |
| `global_k_correction` | Correction factor to adjust mean peak r value from global Ripley’s K-function analysis that will be used in local L analysis. Adjusts the spatial scale (r value) either downward or upward where: a value of 1 means no adjustment, 0 - 1 reduce the spatial scale (detect smaller objects) and > 1 increases the spatial scale (detects larger objects). |
| `min_degree_clustering` | Threshold for cluster segmentation. Because scaling of L(r) is dependent on the value of r, it is difficult to set a manual threshold for intensity of the cluster heat-map. Rather, the script sets a threshold for the minimum degree of clustering (or local molecular density i.e., K value) above what would be expected in a random distribution with the same point density<sup>1</sup>. The actual threshold used for segmentation is calculated as `sqrt(min_degree_clustering) * localL_radius`. |
| `cluster_heatmap_xy_pix_size` | x & y pixel size of the output cluster heat-map.|
| **Watershed (splitting) parameters** ||
| `split_clumped_clusters` | Logical value determining whether to split cluster using watershed. Default is `TRUE` |
| `split_min_r_correction_factor` | Correction factor for `localL_radius` that defines the minimum distance between seeds for watershed |
| `split_min_threshold_correction_factor` | Correction factor for the cluster threshold that defines the minimum intensity that a watershed can have. |
| **Filtering parameters** ||
| `exclude_small_objects` | Logical value that specifies whether to exclude objects smaller than a given size. |
| `min_object_size` | If `exclude_small_objects = TRUE` then object smaller than this size will be excluded. |
| **Miscellaneous parameters** ||
| `matlab_path` | Path to Matlab executable. If this is not specified then the `R.matlab` package will attempt to determine this automatically, although this doesn’t seem to work on Mac OS X. |
| `matlab_port` | Port to run the Matlab server on. Default is 9999. |
| `parallel` | Logical value to specify whether the local L-function part of the analysis should be parallised<sup>2</sup>. |
| `ncores` | If `parallel = TRUE`, this specifies the number of cores to use. Defaults to the number of CPU cores in the computers - 2. If `parallel = FALSE` this setting is ignore.

<sup>1</sup> Owen, D. M., Magenau, A., Williamson, D. J., & Gaus, K. (2013). Super-resolution imaging by localization microscopy. *Methods in Molecular Biology* (Clifton, N.J.), **950**(Chapter 6), 81–93. http://doi.org/10.1007/978-1-62703-137-0_6

<sup>2</sup> This script has not been tested on Windows. The `doParallel` package may not run on Windows so you may need to use a sequential backend (or install and register a `doMC` backend).

## Notes
Many of the functions powering this script come from [supr](https://github.com/keithschulze/supr), which itself is built upon the extraordinary [spatstat](https://github.com/spatstat/spatstat) package. It is highly recommended that you get familiar with `spatstat` as it will help you to understand how this script and the `supr` package work.