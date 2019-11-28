# Australian rodent cranial ecomorphology (Chapter 4 of Thesis) - data and code
Code authors: Ariel E. Marcy, Dr Thomas Guillerme, Dr Vera Weisbecker

To cite the paper, data, and/or code:
> Coming soonish

As of December 2019, this is still a work in progress. Uses R (v. 3.6.1) and geomorph (v. 3.1.2).

After cloning this repo, remember to either set your working directory to the aus-rodent-skulls folder on your computer, or open an RStudio project from that folder.

## Data
**Landmarking data:**
* [MorphoSource Project 561](https://www.morphosource.org/MyProjects/Dashboard/dashboard/select_project_id/561) publically provides 3D meshes for all surface scanned crania in the study
* [Raw_Coordinates.csv](Data/Raw/Raw_Coord_Data.csv) provides the shape coordinates from landmarking 3D skulls in Viewbox 

**Ecological metadata:**
* [Trait data from Breed & Ford 2007](/Data/Processed/in_ex_traits.csv)

If you use these data, please cite the original authors:
> Breed B & Ford F. 2007. Native Mice and Rats. Australian Natural History Series, CSIRO Publishing: Colling-wood, Victoria, Australia, 185 pp. ISBN 978-0-6430-9166-5.

**Phylogenetic data:**
* [Fossil calibrated ultrametric tree from Smissen & Rowe 2018](/Data/Processed/Smissen-Rowe-2018-concat.tre)

If you use these data, please cite the original authors:
> Smissen PJ & Rowe KC. 2018. Repeated biome transitions in the evolution of Australian Rodents. Molecular Phylogenetics and Evolution. 128:182–191. doi: 10.1016/j.ympev.2018.07.015.
    
## Analyses
The analysis workflow is broken down into smaller scripts explained below. Each script loads data created by the script before, so this workflow requires you to run the scripts in order. The intermediate data -- stored as .rda files in the /Data/Processed folder -- are too large to upload to GitHub. 

All of the scripts below are in RMarkdown format (.Rmd), which can be opened in RStudio. There, you can edit and run code chunks as normal, or you can click the Knit button to create HTML versions with both code and output.

* **01-extract-data-for-analyses.Rmd** Extracts both 3D coordinate data as well as the metadata data from Viewbox and prepares it for analysis in `geomorph`. Separates coordinate data into big and small patch protocol datasets. Runs GPA with bilateral symmetry, merges the symmetric shape component with centroid size, and calculates asymmetric variation for both patch datasets.
* **02-calculate-user-error.Rmd** Allows users to view outliers and find major landmarking errors. Takes out replicated specimens from the shape data, finds their duplicates, and calculates user error based on repeatability for both patch datasets. 
* **03-compare-patch-protocols.Rmd** Plots PCAs colored by genus and runs shape ~ genus * centroid size ANOVAs for both patch protocols. Visualizes landmark variation differences with Dr Thomas Guillerme's new package, `landvR`, in color-coded heatmaps for both patch protocols and for both PC1 and PC2. 
* **04-plot-exploratory-PCAs.Rmd** Allows users to quickly plot PCAs with point colors and shapes according to taxa or trait information provided in the metadata table. Also loops through the 7 Australian regions to plot the PCA morphospace by region.
* **05-test-partitioning.Rmd** Tests for a phylogenetic signal difference between the 7 regional biomes and tests for partitioning by region when size-less shape (i.e. the residuals) are used instead; includes phylogenetic correction. 
* **06-visualize-residual-morphospace.Rmd** Plots morphospace differences in PCA of residuals by phylogenetic and guild (diet/locomotion) and generates landvR heatplots of shape changes over the PC axes.
* **07-test-LM-region-differences.Rmd** Sets modules by landmarks and creates 3D plot. Tests covariation via modularity and 2B-PLS using the modules defined for mammalian skulls in Goswami 2006 & 2007. Compares evolutionary rates of modules within crania.

### Custom functions
The analyses call custom functions, most of which are defined in the [..Data/Functions/utilities.R](/Data/Functions/utilities.R) file. A modified version of `geomorph`'s function `plotGMPhyloMorphoSpace` is defined in the [..Data/Functions/plotGMPhyloMorphoSpace_plotmod.R](/Data/Functions/plotGMPhyloMorphoSpace_plotmod.R) file.
