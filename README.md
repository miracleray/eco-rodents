# Australian rodent cranial ecomorphology (Chapter 4 of Thesis) - data and code
Code authors: Ariel E. Marcy, Dr Thomas Guillerme, Dr Vera Weisbecker

To cite the paper, data, and/or code:
> Coming soonish

As of June 2020, this is still a work in progress. Relies on `R` (v. 3.6.1), `geomorph` (v. 3.1.3), and `landvR` (v. 0.4).

*All scripts are in RMarkdown format (.Rmd), which can be opened in RStudio. There, you can edit and run code chunks as normal or use the Knit button to create HTML versions with both code and output. After cloning this repo, remember to either set your working directory to the eco-rodents folder on your computer or open an RStudio project from that folder.*

## Data
**Landmarking data:**
* [MorphoSource Project 561](https://www.morphosource.org/Detail/ProjectDetail/Show/project_id/561) publically provides 3D meshes for all surface scanned crania landmarked in the study.
* [Raw_Coordinates.csv](Data/Raw/Raw_Coord_Data.csv) provides the shape coordinates from landmarked 3D crania exported from Viewbox.

**Ecological metadata:**
* [Trait data from Breed & Ford 2007](/Data/Processed/in_ex_traits.csv)

If you use these data, please cite the original authors:
> Breed B & Ford F. 2007. Native Mice and Rats. Australian Natural History Series, CSIRO Publishing: Colling-wood, Victoria, Australia, 185 pp. ISBN 978-0-6430-9166-5.

**Phylogenetic data:**
* [Fossil calibrated ultrametric tree from Smissen & Rowe 2018 and Marcy et al. 2020](/Data/Processed/Marcy-BEAST01.con.tre)

If you use these data, please cite the original authors:
> Smissen PJ & Rowe KC. 2018. Repeated biome transitions in the evolution of Australian Rodents. Molecular Phylogenetics and Evolution. 128:182–191. doi: [10.1016/j.ympev.2018.07.015.](https://doi.org/10.1016/j.ympev.2018.07.015)

> Marcy AE, Guillerme T, Sherratt E, Rowe KC, Phillips MJ, and Weisbecker V. 2020. Australian rodents reveal conserved Cranial Evolutionary Allometry across 10 million years of murid evolution. bioRxiv. [https://www.biorxiv.org/content/early/2020/05/01/2020.04.30.071308](https://www.biorxiv.org/content/early/2020/05/01/2020.04.30.071308)
    
## Analyses
**The first three scripts prepare the data for analysis and plotting**, the intermediate data they generate are stored as .rda files in the [..Data/Processed](/Data/Processed) folder. These are too large to upload to GitHub so scripts must be run sequentially. 

* [**01-extract-data.Rmd**](/Analysis/01-extract-data.Rmd) Extracts both 3D coordinate data as well as the metadata data from Viewbox and prepares it for analysis with `geomorph`. Runs GPA with bilateral symmetry, merges the symmetric shape component with centroid size, and calculates asymmetric variation.
* [**02-calculate-user-error.Rmd**](/Analysis/02-calculate-user-error.Rmd) Allows users to view outliers, find major landmarking errors, and take out replicated specimens from the shape data. Also calculates user error to report in terms of repeatability. 
* [**03-prepare-data.Rmd**](/Analysis/03-prepare-data.Rmd) Prepares the shape datasets, the metadata, the phylogenetic data, and the vectors with graphics information for plotting. All of which are used throughout the later analyses and are saved in three different .rda files to improve efficiency.

**The next four scripts perform the analyses**, the tables and figures they generate are saved to the [..Data/Results](/Data/Results) folder.

* [**04-compare-PCAs.Rmd**](/Analysis/04-compare-PCAs.Rmd) Plots the full shape PCA versus shape residual PCAs colored by genus or by diet/locomotion. Tests for correlation between PC1 of the two datasets and runs Mantel tests comparing the entire morphospaces of both datasets. Plots screeplots of 3 different PCA datasets. **Generates Figure 1 and Supplementary Figure 1.**
* [**05-heatmap-both-datasets.Rmd**](/Analysis/05-heatmap-both-datasets.Rmd) Plots the `landvR` heatmaps of shape changes over the PC axes for both allometric and residual shape (allometry-free) datasets. **Generates Figures 2 and 3.**
* [**06-test-modularity.Rmd**](/Analysis/06-test-modularity.Rmd) Tests modularity and integration using the modules defined for mammalian skulls in Goswami 2006 & 2007. Modularity tests include `geomorph` function `phylo.modularity`and pairwise Mantel tests of each module. Tests for global integration with `geomorph` function `globalIntegration`. **Generates Figure 4, Table 1, and Supplementary Figure 3.**
* [**07-plot-phylomorph-dist.Rmd**](/Analysis/07-plot-phylomorph-dist.Rmd) Plots phylo-morphological distance plots for both the full shape and shape residual datasets. **Generates Figure 5.**

### Appendices
These contain related analyses, mostly ensuring the dataset and the methods are sound for the analyses above. The last two scripts generate supplementary heatmaps that would not fit in the main manuscript.

* [**App01-compare-sliding methods.Rmd**](/Analysis/App01-compare-sliding-methods.Rmd) Double checks whether sliding in Viewbox before sliding in `geomorph` impacts the landmark variation patterns underlying our shape dataset. 
* [**App02-compare-patch-protocols.Rmd**](/Analysis/App02-compare-patch-protocols.Rmd) Plots PCAs colored by genus and runs shape ~ genus * centroid size ANOVAs for both big and small patch protocols included in the data (the main analyses only use the "small patch" protocol. Visualizes landmark variation differences with `landvR`, in color-coded heatmaps for both patch protocols and for both PC1 and PC2. 
* [**App03-explore-PCAs.Rmd**](/Analysis/App03-explore-PCAs.Rmd) Plots exploratory PCAs with point colors and shapes according to taxa or trait information provided in the metadata table. Creates plots, including a phylomorphospace, better suited to powerpoint. 
* [**App04-test-fitted-allometry.Rmd**](/Analysis/App04-test-fitted-allometry.Rmd) Plots `landvR` heatmaps after performing 3 variations on methods to get the fitted allometric shapes from `procD.lm()`. Creates a methods figure not used on the manuscript.
* [**App05-supplementary-heatmaps.Rmd**](/Analysis/App05-supplementary-heatmaps.Rmd) Plots `landvR` heatmaps of shape changes for the residual dataset without Notomys PC2 and for the residual dataset without Notomys or the aquatic carnivores both PC1 and PC2. **Generates Supplementary Figure 2.**

### Custom functions
Some analyses call custom functions, most of which are defined in the [..Data/Functions/utilities.R](/Data/Functions/utilities.R) file. A modified version of `geomorph`'s function `plotGMPhyloMorphoSpace` is defined in the [..Data/Functions/plotGMPhyloMorphoSpace_plotmod.R](/Data/Functions/plotGMPhyloMorphoSpace_plotmod.R) file.
