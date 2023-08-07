
# Overview

This repo contains code and data used in the paper "Assessing the quality of 
comparative genomics data and results with the cogeqc R/Bioconductor package"


Reproducible reports for all the analyses we performed are available as 
a Quarto book at <https://almeidasilvaf.github.io/cogeqc_paper/>.


## Abstract

Comparative genomics has become an indispensable part of modern biology due to 
the advancements in high-throughput sequencing technologies and the accumulation 
of genomic data in public databases. However, the quality of genomic data and 
the choice of parameters used in software tools used for comparative genomics 
can greatly impact the accuracy of results. To address these issues, 
we present *cogeqc*, an R/Bioconductor package that provides researchers with 
a toolkit to assess genome assembly and annotation quality, 
orthogroup inference, and synteny detection. The package offers context-guided 
assessments of assembly and annotation statistics by comparing observed 
statistics to those of closely-related species on NCBI. To assess orthogroup 
inference, *cogeqc* calculates a protein domain-aware orthogroup score that aims 
at maximizing the number of shared protein domains within the same orthogroup. 
The assessment of synteny detection consists in representing anchor gene pairs 
as a synteny network and analyzing its graph properties, such as clustering 
coefficient, node count, and scale-free topology fit. The application of *cogeqc* 
to real data sets allowed for an evaluation of multiple parameter combinations 
for orthogroup inference and synteny detection, providing researchers with 
guidelines to aid in the selection of the most appropriate tools and parameters 
for their specific data.


