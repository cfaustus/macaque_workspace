This repository contains scripts for the analysis, processing and modelling of alpha globin in long-tailed macaques
Code is associated with the manuscript "Alpha globin variation in the long-tailed macaque indicates malaria selection"

There are three main parts of the manuscript and code is associated with each:
1. Analysis of geographic phenotypic variants of long-tailed macaque alpha globin
    hierarchical_models.py : script for hierarchical model of variants by malaria selection process (Python 2.7)
    mal_variant_hierarchical_model_results.R : visualising results of Bayesian models
    map_figures.R : visualizing alpha globin phenotypes in Southeast Asia
    primate_phylogeny.R : visualizing AQ variants on Macaca phylogeny
2. Genetic analysis of sequences
    processing_fastq_dada2.R: processing of Macaca fascicularis alpha globin Illumina reads
    summarizing_ngs_data.R : cleaning Macaca fascicularis alpha globin reads
    sequencing_viz.R : visualizing cleaned and processed Illumina reads
    MulattaAlpha.py : script to analyze alpha globin from published Macaca mulatta alpha globin genes
3. Population Genetic Model (folder /Population_genetic_model)
    AlphaGlobinPopulationGeneticModel.c : program to perform the population genetic simulations
    MatlabScriptToRunPopulationGeneticModel.m : facilitates running the MEX file in matlab

Cleaned data for sequencing is located in the dada2_data folder and raw sequences are located on the Sequence Read Archive (BioProject ID: PRJNA639946) under the accession numbers: SRR12404495- SRR12404572 (HBA1) and  SRR12404678- SRR12404755 (HBA2)).
Data for other scripts (i.e. shapefiles, phenotype by country, etc) are located in the 'data' folder
Output from scripts are located in the output folder
