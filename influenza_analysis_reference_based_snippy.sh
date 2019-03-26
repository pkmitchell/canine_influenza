#!/bin/bash

## This script is for the assembly for influenza sequencing samples via reference-based assembly

## Author: Derek Rothenheber
## Date: 9-14-18
## Affliation: Cornell University - Animal Health and Diagnostic Center



# Exporting Program Paths to $PATH
 # Need to add all different programs to the command-line $PATH
 # so that they can be executed by calling the program name.



#Trinity Program (v2.6.6) - De-novo assembly program
export PATH=/programs/trinityrnaseq-Trinity-v2.6.6:$PATH
#Bowtie2 Program (v2.3.0) - Short read aligner
export PATH=/programs/bowtie2-2.3.0:$PATH
#Samtools Program (v1.8) - Suite of utilities for manipulating alingments in SAM format
export PATH=/programs/samtools-1.8/bin:$PATH
#Jellyfish Program (v2.2.3) - kmer counting program
export PATH=/programs/jellyfish-2.2.3/bin:$PATH
#Salmon Program (v0.10.0) - quanitifying gene expression (Used for quantifiying coverage of reads)
export PATH=/programs/salmon-0.10.0/bin:$PATH
# FastQC Program (v0.11.5) - Q/C  program checking a variety of metrics from sequencing reads
export PATH=/programs/FastQC-0.11.5:$PATH
# BBMap Program (v37.50) - Q/C program to obtain read length distribution for histogram
export PATH=/programs/bbmap-37.50:$PATH
# Lighter Program - Fast and memory efficient error corrector
# Lighter needs to be used with full PATH
# Khmer (v0.4) - Library and suite of tools to work with sequences
export PATH=/programs/khmer/scripts:$PATH
# Kraken2 (v2) - Read binning using hash tables and khmer frequency
export PATH=/programs/kraken2:$PATH
# HMMER (3.1b1) - Search sequence database for homologs of protein sequences uses hidden Markov models
export PATH=/programs/hmmer-3.1b2-linux-intel-x86_64/binaries:$PATH
# Prodigal (v2.6.3) - Program needed from Prokka
export PATH=/programs/prodigal-2.6.3:$PATH
# GNU Parallel
export PATH=/programs/parallel/bin:$PATH
# Prokka The genome annotation tool
source /workdir/prokka-1.12/source.sh
# Pigz - gzip in parallel
export PATH=/programs/pigz-2.4:$PATH
# Snippy - Reference-based assembly and SNP calling
export PATH=/programs/snippy/bin:/programs/snippy/binaries/linux/:$PATH
export PATH=/programs/seqtk:$PATH
export PATH=/programs/samclip:$PATH
export PATH=/programs/vt:$PATH




echo
echo =========================================================================
echo Script for Canine Influenza Reference-based genome assembly and analysis
echo Author: Derek Rothenheber
echo Affliation: Cornell University - Animal Health and Diagnostic Center
echo email: dr575@cornell.edu
echo github: @drothen15/canine_influenza
echo =========================================================================
echo
sleep 5s


echo -n 'enter fasta reference genome PATH: '
read reffasta


## Step 1 Initial FastQC on foward and reverse reads
echo
echo --------------------------------------------------------------------------------------------------------
echo Step 1: Running  FASTQC on concatanated fowards and reverse reads and obtaining read length distribution
echo --------------------------------------------------------------------------------------------------------
echo

rm -r -f fastqc/
mkdir fastqc

unpigz *R1* | unpigz *R2*
cat *R1* > fastqc/cat-r1.fastq | cat *R2* > fastqc/cat-r2.fastq
pigz *R1* | pigz *R2*

fastqc fastqc/cat-r1.fastq -o fastqc/ -t 40 --nogroup
fastqc fastqc/cat-r2.fastq -o fastqc/ -t 40 --nogroup


## awk command retriving read length distribution - in the future an R script will plot the results (two column data structure)
cd fastqc/
awk 'NR%4 == 2 {lengths[length($0)]++ ; counter++} END {for (l in lengths) {print l, lengths[l]}; print "total reads: " counter}' cat-r1.fastq > readlength-r1.csv
awk 'NR%4 == 2 {lengths[length($0)]++ ; counter++} END {for (l in lengths) {print l, lengths[l]}; print "total reads: " counter}' cat-r2.fastq > readlength-r2.csv
sort -t, -nk1 readlength-r1.csv > sorted.readlength-r1.csv
sort -t, -nk1 readlength-r2.csv > sorted.readlength-r2.csv
cd ..


## Step 2 Error Correction with Lighter
echo
echo ------------------------------------------------------------------------
echo Step 2: Error Correction wtih Lighter - program is for RNA based data
echo ------------------------------------------------------------------------
echo

rm -r -f lighter-output/
rm -f lighter.cmds
mkdir lighter-output


for R1 in *R1*
do
   R2=${R1//R1_001.fastq/R2_001.fastq}
#   echo $R1 $R2
   echo "/programs/Lighter/lighter -r $R1 -r $R2 -t 40 -k 25 13000 .1 -od lighter-output/" >> lighter.cmds
done

cat lighter.cmds
chmod u+x lighter.cmds

./lighter.cmds



## Step 3 Trimming Error Corrected Sequences with Trimmomatic - using sliding window 4:5 check for Nextera Adapaters
echo
echo ----------------------------------------------------------
echo Step 3 Trimming Error Corrected Sequences with Trimmomatic
echo ----------------------------------------------------------
echo

rm -f -r trimmomatic-output/
mkdir trimmomatic-output
rm -f trimmomatic.cmds

cd lighter-output/

for R1 in *R1*
do
   R2=${R1//R1_001.cor.fq/R2_001.cor.fq}
   R1paired=${R1//.fq/.paired.fq}
   R1unpaired=${R1//.fq/.unpaired.fq}
   R2paired=${R2//.fq/.paired.fq}
   R2unpaired=${R2//.fq/.unpaired.fq}
#   echo $R1 $R2 $R1paired $R1unpaired $R2paired $R2unpaired
   echo "java -jar /programs/trimmomatic/trimmomatic-0.36.jar PE lighter-output/$R1 lighter-output/$R2 trimmomatic-output/$R1paired trimmomatic-output/$R1unpaired trimmomatic-output/$R2paired trimmomatic-output/$R2unpaired ILLUMINACLIP:/programs/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:2 TRAILING:2 MINLEN:25 SLIDINGWINDOW:4:5" >> trimmomatic.cmds
done

cat trimmomatic.cmds
chmod u+x trimmomatic.cmds
mv trimmomatic.cmds ../

cd ../

./trimmomatic.cmds



## Step 4: Rerunning FastQC post error correction and trimming
echo
echo ------------------------------------------------------------------------------------
echo Step 4 Rerun FastQC - post error correction and trimmming, compare pre and post runs
echo ------------------------------------------------------------------------------------
echo

cd fastqc/
mkdir post-qc
cd ../trimmomatic-output

unpigz *R1* | unpigz *R2*
cat *R1_001.cor.paired* > ../fastqc/post-qc/postqc-cat-r1.fastq | cat *R2_001.cor.paired* > ../fastqc/post-qc/postqc-cat-r2.fastq
#pigz *R1* | pigz *R2*

cd ..

fastqc fastqc/post-qc/postqc-cat-r1.fastq -o fastqc/post-qc/ -t 40 --nogroup
fastqc fastqc/post-qc/postqc-cat-r2.fastq -o fastqc/post-qc/ -t 40 --nogroup



## Step 5: Reference-based Assembly with Bowtie2 Package (SAM/BAM file formating)

echo
echo ------------------------------------------------------------------------------------
echo Step 5 Reference-based Assembly with Snippy!!
echo ------------------------------------------------------------------------------------
echo



cd trimmomatic-output/

for file in *R1_001.cor.paired*
do
   R1=$file
   R2=${file//R1_001.cor.paired.fq/R2_001.cor.paired.fq}
   out=${file//L001_R1_001.cor.paired.fq/output}
   pre=${file//L001_R1_001.cor.paired.fq/snps}
   echo "snippy --outdir $out --R1 $R1 --R2 $R2 --reference $reffasta --mapqual 25 --prefix $pre --cpus 35 --ram 100" >> snippy.cmds
done


cat snippy.cmds
chmod u+x snippy.cmds
./snippy.cmds


rm -r -f consensus-fasta
mkdir consensus-fasta


for file in *output/*consensus.fa
do
   echo "cp $file consensus-fasta/" >> copy.cmds
done

cat copy.cmds
chmod u+x copy.cmds
./copy.cmds



#######################################################################
## Renaming fasta headers
echo
echo ------------------------------------------------------------------------------------
echo Step 6 Renaming Fasta Headers with Sample Name and Calculating Read Counts
echo ------------------------------------------------------------------------------------
echo
######################################################################

cd consensus-fasta

mkdir final_consensus_fasta

quote="'"

for file in *.fa
do
  output=${file//consensus.fa/consensus_FINAL.fa}
  echo awk $quote '/>/{sub(">","&"FILENAME" reference:");sub(/\.fa/,x)}1' $quote  $file ">" $output >> rename.cmds
done

cat rename.cmds
chmod u+x rename.cmds
./rename.cmds


mv *FINAL* final_consensus_fasta



echo Reference Assembly Complete
echo Final consensus fasta files are located in the final_consensus_fasta directory!!
