#!/bin/bash

## This script is for the assembly for influenza sequencing samples via de-novo and referenced based assembilies

## Author: Derek Rothenheber
## Date: 9-14-18
## Affliation: Cornell University - Animal Health and Diagnostic Center



# Exporting Program Paths to $PATH
 # Need to add all different programs to the command-line $PATH
 # so that they can be executed by calling the program name.



#Trinity Program (v2.6.6) - De-novo assembly program
export PATH=/programs/trinityrnaseq-Trinity-v2.6.6:$PATH
#Bowtie Program (v2.3.0) - Short read aligner
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



echo
echo =========================================================================
echo Analysis script for Canine Influenza De-novo genome assembly and analysis
echo Author: Derek Rothenheber
echo Affliation: Cornell University - Animal Health and Diagnostic Center
echo email: dr575@cornell.edu
echo =========================================================================
echo
sleep 5s



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
pigz *R1* | pigz *R2*

cd ..

fastqc fastqc/post-qc/postqc-cat-r1.fastq -o fastqc/post-qc/ -t 40 --nogroup
fastqc fastqc/post-qc/postqc-cat-r2.fastq -o fastqc/post-qc/ -t 40 --nogroup




















################################## Step 6: Trinity Assembly ######################################################
echo
echo ---------------------------------------------------------------------------------
echo Step 6: ~~~~~~~~~~~~~~~Trinity Assembly with digital normalization~~~~~~~~~~~~~~~
echo ---------------------------------------------------------------------------------
echo

rm -r -f trinity/
rm -f trinity.cmds
cd trimmomatic-output/


for R1 in *R1_001.cor.paired.fq.gz
do
   R2=${R1//R1_001.cor.paired.fq/R2_001.cor.paired.fq}
   wd=${R1//L001_R1_001.cor.paired.fq.gz/trinity-output}
   #echo $R1 $R2
   echo "Trinity --seqType fq --max_memory 200G --CPU 40 --left trimmomatic-output/$R1 --right trimmomatic-output/$R2 --output $wd/" >> trinity.cmds
done

cat trinity.cmds
chmod u+x trinity.cmds
mv trinity.cmds ../
cd ..
./trinity.cmds


mkdir trinity
mv *trinity-output/ trinity/



## Step 7 Estimate of adundance using RSEM within Trinity (TPM && FPKM)
echo
echo --------------------------------------------------------------------------------------------------
echo Step 7: ~~~~~~~~~~~~~~~Estimate of adundance using RSEM within Trinity TPM and FPKM~~~~~~~~~~~~~~~
echo --------------------------------------------------------------------------------------------------
echo

##################  Renaming all output files with sample specific prefixes and reformating fasta file to be blast-able


cd trinity/
rm -r -f trinity-output/

for file in */*
do
  fname=${file##*/} #This gives your base filename.
  fpath=${file%/*} # Your dir
  dname=${fpath##*/}
  mv $file ${fpath}/${dname}.${fname}
done

for fasta in */*Trinity.fasta
do
#   echo $fasta
   mv $fasta ./
done

for map in */*Trinity.fasta.gene_trans_map
do
#   echo $map
   mv $map ./
done

## Now starting RSEM analysis ##

mkdir trinity-output
mv *Trinity.fasta trinity-output/
mv *Trinity.fasta.gene_trans_map trinity-output/

rm -r -f resm-abundance/
rm -f resm.cmds

for file in trinity-output/*Trinity.fasta
do
   file2=${file//.fasta/.fasta.gene_trans_map}

   r1=${file//_trinity-output.Trinity.fasta/_L001_R1_001.cor.classified.fq}
   r1p=${r1//trinity-output/}

   r2=${file//_trinity-output.Trinity.fasta/_L001_R2_001.cor.classified.fq}
   r2p=${r2//trinity-output/}

   out=${file//_trinity-output.Trinity.fasta/_resm-abundance}
   outdir=${out//trinity-output/}
   echo $r1p $r2p
   echo $outdir

   echo "/programs/trinityrnaseq-Trinity-v2.6.6/util/align_and_estimate_abundance.pl --transcripts $file --seqType fq --left ../kraken-output$r1p --right ../kraken-output$r2p --est_method RSEM --aln_method bowtie2 --gene_trans_map $file2 --output_dir resm-adundance$outdir --thread_count 40 --prep_reference --include_rsem_bam" >> resm.cmds
done


cat resm.cmds
chmod u+x resm.cmds

./resm.cmds

###################### Reformating Trinity Output files for later BLASTING ###########################

cd trinity-output/

for fasta in *Trinity.fasta
do
   echo $fasta
   cat $fasta | sed "s/ path.*//gi" | sed "s/ /_/gi" | sed  "s/=//gi" | sed "s/TRINITY/$fasta/gi" > $fasta.temp
   awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $fasta.temp > $fasta.temp2
   tail -n +2 $fasta.temp2 > $fasta

   rm -f $fasta.temp
   rm -f $fasta.temp2
done

cd ../


echo
echo ----------------------------------------------------------------------------
echo Step 8 Finding the Longest ORFs searching for homologs with ProKKa with Pfam
echo ----------------------------------------------------------------------------
echo


rm -r -f prokka-output
mkdir prokka-output

for fasta in trinity-output/*Trinity.fasta
do
  input=$fasta
  temp=${fasta##*/}
  fname=${temp//_trinity-output.Trinity.fasta/}
  echo $input
  echo $fname
  echo "prokka --outdir prokka-output/$fname --force --compliant --kingdom Viruses --cpus 40 $input" >> prokka.cmds
done

cat prokka.cmds
chmod u+x prokka.cmds

./prokka.cmds

cd prokka-output/

for file in */*
do
  fname=${file##*/} #This gives your base filename.
  fpath=${file%/*} # Your dir
  dname=${fpath##*/}
  mv $file ${fpath}/${dname}.${fname}
done


#rm -r -f filtered.cmds

#for fasta in */*.ffn
#do
#  input=$fasta
#  temp=${fasta##*/}
#  output=${temp//.ffn/_final.fasta}
#  echo "python filter-genomic-segments.py /workdir/dr575/test-dir/trinity/prokka-output/$input $output" >> filtered.cmds
#done

#cat filtered.cmds
#chmod u+x filtered.cmds

#mv filtered.cmds /workdir/dr575/test-dir/

#cd /workdir/dr575/test-dir/
#rm -r -f annotated-fasta/
#mkdir annotated-fasta

#./filtered.cmds

#mv *final.fasta annotated-fasta/







