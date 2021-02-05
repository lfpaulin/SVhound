# SVhound

SVhound is a framework to predict regions that harbour so far unidentified genotypes of Structural Variations. It uses a population size VCF file as input and reports the probabilities and regions across the population. 
SVhound was tested and applied to the 1000genomes VCF file and also other data sets. 


## Quickstart

First use the vcf_parser_for_svhound.py to extract the SV genotypes given a window size:

```
zcat my_Vcf.gz | python vcf_parser_for_svhound.py 10000 > my_Vcf_extracted10kbpwindows.txt 

```

You can apply filters or other commands to pipe the VCF file into the python script. 


Next, you going to run SVhound: 

```
Luis ?
```


## Rhesus macaque SV Calling

We identified SV from 150 Illumina data sets based on bwa mem mappints to the v8 of the reference genome. The SV were identified by Manta and subsequently merged and filtered using SURVIVOR. 


## Citation
