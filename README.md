# Ultrasens_DNAMethylation
This repository is a collection of code for the associated publication, "A Tunable, Ultrasensitive Threshold in Enzymatic Activity governs the DNA Methylation Landscape". It contains various terminal/console commands, Python scripts written in Jupyter Notebook, and MATLAB files which maniupulate the required genomic datasets, perform parameterization of ensuing models to experimentla data, and create the visualizations shown in the publication. 

## Goal
Taking .bed file of WGBS (or other methylation) dataset(s), perform CGI-Level analysis of methylation macrostates, Individual CpG analysis of methyl fraction dependence on local CpG landscape, perform Hill function fitting of Individual CpG mean methylation, and finally fit landscape using Master Equation.

## Files Needed
1. CpG Island Annotation (.bed File)
	1. Available directly from the UCSC Genome/Table Browser under "Regulation". Note that this is specifically created for the hg19 Assembly, but [This Page](https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=cpgIslandExt) describes the schema, along with how to manually perform the island annotation.
    
1. Whole-Genome Bisulfite Sequencing (WGBS) Datasets (.bed file(s) )
   1. Via their, GEO Accession numbers (shown in Table 1 within the Main Text), these files can be downloaded directly from the [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/)
   2. Note that some files may be in the ".wig" format, but can be converted within the terminal environment using the following command (available via the bedops package):
      
```
wig2bed < $WGBS_filename.wig > $WGBS_filename.bed
```
      
1. Entire Genome Sequence (hg19) (.fasta file)
   1. This may be downloaded also from the UCSC Genome Browser using the following command (Mac OS)

```
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ .
```
   2. Note that each chromosome is separated into its own file; the entire genome requires ~3.2gb of space.
   3. Using the "Seqkit" package, one can filter each chromosome for CpGs using the following command from the terminal environment; the example here performs the operation only on Chromosome 1:
```
cat chr1.fa | seqkit locate -P -p  cg > {OUTPUT_FILENAME}.csv
```

