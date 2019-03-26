# Assembly and Annotation Pipelines for Canine Influenza A - Vet Molecular Diagnostics  
This repo contains shell and python scripts used to assemble raw fastq sequencing files into full concatenated genomes.
All of these scripts are optimized to run in Cornell's HPC environment, but can easily be transfered to other HPCs or AWS.

The following scripts are including in this repo:
1. influenza_analysis.sh (De-novo assembly with Kraken binning and Trinity Assembly)
2. influenza_analysis_no-kraken.sh (De-novo assembly without Kraken binning and Trinity Assembly)
3. filter-genomic-segments.py (Filters out small contigs and keeps only the largest contig for each of the 8 genes)
4. influenza_analysis_reference_based.sh (Reference-based assembly, with Bowtie2 and Samtools)
5. influneza_analysis_reference_based_snippy.sh (Reference-based assembly, with the Snippy pipeline)


