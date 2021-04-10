# SVhound

SVhound is a framework to predict regions that harbour so far unidentified genotypes of Structural Variations. It uses a population size VCF file as input and reports the probabilities and regions across the population. 
SVhound was tested and applied to the 1000genomes VCF file and also other data sets. 


## Quickstart

First use the vcf_parser_for_svhound.py to extract the SV genotypes given a window size:

```
zcat my_Vcf.gz | python vcf_parser_for_svhound.py 10000 > my_Vcf_extracted10kbpwindows.txt 

```

You can apply filters or other commands to pipe the VCF file into the python script. 


Next, in R you are going to load the created table "my_Vcf_extracted10kbpwindows.txt": 

```
my_Vcf_10kbpwindows <- read.table("my_Vcf_extracted10kbpwindows.txt")

```

To then run SVhound:

```
library("SVhound")

svhound(SVallelesTable=my_Vcf_10kbpwindows, window_size=10000)

```

SVhound will write the results in the same folder a file named "results-svhound_analysis.RData". In order to change the name of the results file, use parameter "output_prefix"

```
library("SVhound")

svhound(SVallelesTable=my_Vcf_10kbpwindows, window_size=10000, output_prefix="my_file")

```


## Rhesus macaque SV Calling

We identified SV from 150 Illumina data sets based on bwa mem mappints to the v8 of the reference genome. The SV were identified by Manta and subsequently merged and filtered using SURVIVOR. 


## Citation


SVhound: Detection of future Structural Variation hotspots
Luis Felipe Paulin, Muthuswamy Raveendran, Ronald Alan Harris, Jeffrey Rogers, Arndt von Haeseler, Fritz J Sedlazeck
bioRxiv 2021.04.09.439237; doi: https://doi.org/10.1101/2021.04.09.439237
