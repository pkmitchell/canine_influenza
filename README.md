# Assembly and Annotation Pipelines for Canine Influenza A - Vet Molecular Diagnostics  
This repo contains shell and python scripts used to assemble raw fastq sequencing files into full concatenated genomes.
All of these scripts are optimized to run in Cornell's HPC environment, but can easily be transfered to other HPCs or AWS.

The following scripts are including in this repo:
1. influenza_analysis.sh (De-novo assembly with Kraken binning and Trinity Assembly)
2. influenza_analysis_no-kraken.sh (De-novo assembly without Kraken binning and Trinity Assembly)
3. filter-genomic-segments.py (Filters out small contigs and keeps only the largest contig for each of the 8 genes)
4. influenza_analysis_reference_based.sh (Reference-based assembly, with Bowtie2 and Samtools)
5. influneza_analysis_reference_based_snippy.sh (Reference-based assembly, with the Snippy pipeline)


# Reference-based Assembly Pipelines
I tend to use these pipelines more often than De-novo because they don't require as much sequencing depth to obtain complete and full genomes

1. Running Pipeline
```bash
./influenza_analysis_reference_based_snippy.sh
or
./influenza_analysis_reference_based.sh
```
- IMPORTANT - 
Analysis script needs to be in the same directory as the raw fastq files.
Raw fastq files need to be gzip.
No other files can be in the directory with the raw fastq files and analysis script.
There is no limit to the number of samples run at a time, the script scales to number of samples inputed 

2. User Provided Positional Arguments 
```bash
enter fasta genome PATH: <$PATH to reference fasta file>
```

3. Output files
```bash
fastqc/  - Fastqc output from raw reads and post Q/C reads

lighter-output/  - Error corrected fastq files 

trimmomatic-output/  - Error corrected and trimmed reads

trimmomatic-output/Sample_output/  - Snippy output files, contains vcf, bam, etc.

trimmomatic-output/consensus-fasta/final_consensus_fasta/  - Final fasta file containing called SNPs instantiated
```

