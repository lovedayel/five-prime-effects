# Five prime effects

This repository contains the source code for running the analysis described in the paper.


'FYP_key' contains the key to translate the figures from the paper to the figures in Loveday's final year project.

'genomes.csv' is a list of the set of 650 bacteria considered in the paper and their accession numbers. The corresponding EMBL files can be found in 'bacterial\_filtered_genomes.zip'

'Ecoli\_all_vectors.csv' contains all of the vectors of log odds ratios calculated for *E. coli* by codon.

'650_V5-core.csv' contains the values of V<sub>5-core</sub> for each of the 650 bacterial genomes.

'pax\_db\_genomes' and 'pax\_db_proteomes' contain the GenBank files and protein abundance data respectively for the five additional genomes whose proteomics data we consider in the paper


## Folders explained:

### The folders contain all of the scripts, input files and output files associated with each figure in the paper

- Fig1: calculate codon enrichment at the 5' ends of highly expressed *E. coli* transgenes and plot the log odds ratios for their enrichment (V<sub>edIO</sub>)
- Fig2: calculate codon enrichment at the 5' ends of *E. coli* genes compared to the cores and plot against V<sub>edIO</sub>
- Fig3: calculate codon enrichment at the 5' ends of genes compared to the cores across all 650 genomes (V<sub>5-core</sub>), correlate each vector with V<sub>edIO</sub> then plot a histogram of the resulting correlation coefficients
- Fig4a: plot the relationship between these correlation coefficients and GC3 content
- Fig4b: plot the relationship between these correlation coefficients and optimal growth temperature
- Fig5: calculate codon enrichment at the 5' ends of native *E. coli* genes with highly abundant proteins (V<sub>5-prot</sub>) and plot against V<sub>edIO</sub>
- Fig6: calculate codon enrichment at the cores of native *E. coli* genes with highly abundant proteins (V<sub>TO</sub>) and plot against V<sub>edIO</sub>
- Fig7: calculate codon enrichment at the 5' ends of *E. coli* genes with low noise measurements (V<sub>noise</sub>)and plot against V<sub>5-prot</sub>

- FigS1: plot V<sub>edIO</sub> against V<sub>TO</sub>
- FigS2: plot the relationship between GC3 content and the difference in predicted stability between 5' ends and cores across all 650 genomes
- FigS3: calculate codon enrichment at the 5' ends of genes with highly abundant proteins (V<sub>5-prot</sub>) for four additonal genomes and plot against V<sub>edIO</sub>
- FigS4a: calculate codon enrichment at the 5' ends of essential genes compared to non-essential (V<sub>ess</sub>) and plot against *E. coli* V<sub>5-prot</sub>
- FigS4b: plot V<sub>ess</sub> against V<sub>edIO</sub>

### A spreadsheet with detailed information on each script (including inputs and outputs) can be found in the document 'Script_details'
