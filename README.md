# Australian rodent cranial ecomorphology (Chapter 4 of Thesis) - data and code
Code authors: Ariel E. Marcy, Dr Thomas Guillerme, Dr Vera Weisbecker

To cite the paper, data, and/or code:
> Coming soonish

As of April 2020, this is still a work in progress. Uses `R` (v. 3.6.1), `geomorph` (v. 3.1.3), and `landvR` (v. 0.4).

*All of the scripts are in RMarkdown format (.Rmd), which can be opened in RStudio. There, you can edit and run code chunks as normal, or you can click the Knit button to create HTML versions with both code and output. After cloning this repo, remember to either set your working directory to the allometry-rodents folder on your computer, or open an RStudio project from that folder.*

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
The analysis workflow is broken down into smaller scripts explained below. Each script loads data created by the script before, so this workflow requires you to run the scripts in order. The intermediate data -- stored as .rda files in the [..Data/Processed](/Data/Processed)  folder -- are too large to upload to GitHub. 

* **01-extract-data-for-analyses.Rmd** Extracts both 3D coordinate data as well as the metadata data from Viewbox and prepares it for analysis in `geomorph`. Runs GPA with bilateral symmetry, merges the symmetric shape component with centroid size, and calculates asymmetric variation.
* **02-calculate-user-error.Rmd** Allows users to view outliers and find major landmarking errors. Takes out replicated specimens from the shape data, finds their duplicates, and calculates user error based on repeatability for both patch datasets. 
* **03-prepare-data.Rmd** Preps datasets used throughout the later analyses for metadata, phylogenetic data, and graphics for plotting.
* **04-plot-allometry-PCAs.Rmd** Plots exploratory PCAs with point colors and shapes according to taxa or trait information provided in the metadata table.
* **05-plot-residual-PCAs.Rmd** Plots the residual PCAs by phylogenetic and by guild (diet/locomotion). **Generates Figure 1.**
* **06-heatmap-both-datasets.Rmd** Plots the `landvR` heatplots of shape changes over the PC axes for both allometric and residual datasets. **Generates Figures 2 and 3.**
* **07-test-LM-region-differences.Rmd** Sets modules by landmarks and creates 3D plot. Tests covariation via modularity using the modules defined for mammalian skulls in Goswami 2006 & 2007. Runs Mantel tests on the distance matrices of the PC scores for each module. Tests for global integration. **Generates Figure 4 and Table 1.**
* **08-visualize-phylo-distances** Plots phylo-morphological distance plots for both full shape and residual shape datasets. **Generates Figure 5.**

### Appendices
These contain related analyses, mostly ensuring the dataset and the methods are sound for the analyses above. The last script generates supplementary heatplots that would not fit in the manuscript.

* **App01-compare-sliding methods.Rmd** Double checks whether sliding in Viewbox before sliding in `geomorph` impacts the landmark variation patterns underlying our shape dataset. 
* **App02-compare-patch-protocols.Rmd** Plots PCAs colored by genus and runs shape ~ genus * centroid size ANOVAs for both big and small patch protocols. Visualizes landmark variation differences with Dr Thomas Guillerme's new package, `landvR`, in color-coded heatmaps for both patch protocols and for both PC1 and PC2. 
* **App03-test-fitted-allometry.Rmd** Plots `landvR` heatplots after performing 3 variations on methods to get the fitted allometric shapes from `procD.lm()`. **Generates Supplementary Figure TK.**
* **App04-residuals-without.Rmd** Plots `landvR` heatplots of shape changes for the residual dataset without Notomys PC2 and for the residual dataset without Notomys or the aquatic carnivores both PC1 and PC2. **Generates Supplementary Figure TK.**

### Custom functions
The analyses call custom functions, most of which are defined in the [..Data/Functions/utilities.R](/Data/Functions/utilities.R) file. A modified version of `geomorph`'s function `plotGMPhyloMorphoSpace` is defined in the [..Data/Functions/plotGMPhyloMorphoSpace_plotmod.R](/Data/Functions/plotGMPhyloMorphoSpace_plotmod.R) file.
