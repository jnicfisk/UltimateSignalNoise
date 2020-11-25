# UltimateSignalNoise
This repositiory is a temporary repository until package is ready to be hosted on CRAN/Bioconductor.

Additionally, this is not yet a stable release. If you need a stable release, please use [PhyInformR](https://cran.r-project.org/src/contrib/Archive/PhyInformR/PhyInformR_1.0.tar.gz). Original [paper](https://bmcevolbiol.biomedcentral.com/articles/10.1186/s12862-016-0837-3) | Usage [manual](https://carolinafishes.github.io/software/phyinformR/) available.

## Related Prior Work
UlitmateSignalNoise is related to advances described in the following papers (among others):

Jeffrey P. Townsend et al | [Phylogenetic Signal and Noise: Predicting the Power of a Data Set to Resolve Phylogeny Systematic Biology OCTOBER 2012](https://www.jstor.org/stable/41677982)

Lopez-Giraldez et al | [Evaluating phylogenetic informativeness as a predictor of phylogenetic signal for metazoan, fungal, and mammalian phylogenomic data sets. BioMed Research International 2013](https://pubmed.ncbi.nlm.nih.gov/23878813/) 

Zuo et al | [The impact of incorporating molecular evolutionary model into predictions of phylogenetic signal and noise Front. Ecol. Evol. April 2014](https://doi.org/10.3389/fevo.2014.00011)

Alex Dornburg et al | [Optimal Rates for Phylogenetic Inference and Experimental Design in the Era of Genome-Scale Data Sets, Systematic Biology January 2019](https://doi.org/10.1093/sysbio/syy047)

## Overview
UltimateSignalNoise is a phylogenetic experimental design tool in R. An overview of the workflow is provided below.

![Workflow](https://github.com/jnickfisk/UltimateSignalNoise/blob/main/doc_images/Aim%204%402x.png)

## Purpose
Phylogenetic trees represent an amalgam of interconnected hypotheses, only some of which are of interest to particular research programs. Determining the optimal gene-taxa sampling schema prospectively allows for hypothesis-driven data collection and retrospective filtering to maximize the probability of achieving sufficient power to resolve specific hypotheses. UltimateSignalNoise aims to use phylogenetic signal/Q\*RP-based methods to provide users with an optimal sampling schema to best achieve power to address their specific hypotheses related to:
- differentiating the power of different genes and phylogenetic characters
- weighing the relative utility of increased taxonomic versus character sampling 
- designing taxonomically dense phylogenetic studies optimized by taxonomically sparse genome-scale data.
- ...and more!

*In short, UltimateSignalNoise is an implementation of the aforementioned theoretical tools to guide phylogenetic experimental design using advances to the phylogenetic signal framework to iteratively rank-order gene-taxa sampling schema and to ensure proposed sampling schema reach the desired power to answer specific phylogenetic hypotheses.*

This iterative process is represented in the following gifs where the grey cells represent currently examined taxa-gene pairs (convienently sampling the best next gene/taxa pair in this example), but each unresolved cell (white) will be calculated unless masked. 

Shannon Tree Collapse             |  Sampling Matrix
:-------------------------:|:-------------------------:
![collapse](https://github.com/jnickfisk/UltimateSignalNoise/blob/main/doc_images/collapse2.gif)  |  ![matrix](https://github.com/jnickfisk/UltimateSignalNoise/blob/main/doc_images/collapse1.gif)

## Guide 
A more robust user guide is forthcoming. For now, the program can be run by getting your data loaded into R in the following data types/structures:
- Base frequencies | freqs | numeric vector (length of 4) (T,C,A,G) 
- Substitution Parameters | subratevector | numberic vector (length of 6) (a,b,c,d,e,f)
- Site Rates | ratevector | numeric vector (length of # of input sites)
- Guide Tree | guideTree | ape tree 
- Emperical Tree | empiricalTree | ape tree
- Nucleotide Alignment | alnFile | Multifasta alignment
- internode specification | internode | numeric vector (internal node labels uniquely identifying interode of interest)

and invoking 

```R
TreeCollapse.signal.noise(freqs=c(0.34,0.16,0.32,0.18),subratevector=c(5.26,8.15,1,2.25,3.16,5.44),ratevector=c(rep(0.05,20),rep(0.1,10),rep(0.05,20)),guideTree=guideT, empiricalTree=empricialT, alnFile="aln.fasta", internode=c(9,10))
```

The stable package 1.0 version will be released upon completeion of the following:
- Integration of IQTREE for sitewise rate calculation
- Restructuring of data structure to rely on a custom S4 R object that can be converted between ape/phylo formated trees to prevent crashing due to quirks in those structures

Future versions will include:
- Parellelisation 
- Visualization
- More customization (e.g. ability to constrain not just on QIRP, but also on QIPP or QIHP). 
