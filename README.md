# Power Analyses

A collection of power analyses done for different sequencing experiments for researchers at Brown. Currently contains power analyses for RNAseq as well as 16s and wgs metagenomics experiments.

This repository is intended to serve as an internal resource for the [Computational Biology Core](cbc.brown.edu), although all data used is available from public sources. As such, no support is provided for external users who run into trouble with running the code, and notebooks may become out of date as various packages stop being maintained. Nevertheless we will do our best to keep things updated and self-contained.

# File Structure

Each folder in the repo is intended to have all the necessary pieces to run the associated power analysis notebook. In general there will be an RMarkdown notebook and pieces of data such as count matrices. To run the code all that needs to be done is to open the notebook and knit the code.

Top level folders are separated out by sequencing technology.

## Folders

 * **RNAseq**
 	* **mouse_muscle_2019:** A power analysis done using `RNASeqSampleSize` for an RNAseq experiment involving mouse muscle tissue.
	* **skin_2021:** A power analysis done using `RNASeqSampleSize` for an RNAseq experiment involving skin, lung tissue and looking into Reactome pathways. The data here is not in the repo due to files from GTeX that are ~2G in size, but can be provided if needed.
 * **microbiome**
	* **tutorial:** A tutorial on running `micropower` on data from HMP
 * **dscov**: Slides for a DSCoV workshop on power analyses given 3/3/2023.

# Other relevant information

The metagenomics power analysis was eventually published as [a blog post](https://medium.com/brown-compbiocore/power-analyses-for-microbiome-studies-with-micropower-8ff28b36dfe3).
