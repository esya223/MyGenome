#!/bin/bash

#SBATCH --time 48:00:00
#SBATCH --job-name=trim-spades
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal
#SBATCH --mem=500GB
#SBATCH --mail-type ALL
#SBATCH -A cea_farman_s25abt480
#SBATCH --mail-type ALL
#SBATCH --mail-user farman@uky.edu,szh231@uky.edu

echo "SLURM_NODELIST: "$SLURM_NODELIST

# Requires paired-end read datasets

#1. Place raw reads into RAW_READS directory

#2. Run assembly using command: $ sbatch trim-spades.sh RAW_READS SRRXXXXXXX

#3. Use for loop to submit assemblies in parallel: $ for f in `cat SRRlist.txt`; do sbatch trim-spades.sh RAW_READS $f; done


dir=$1

f=$2

mkdir $f 

mv $dir/$f*[_R]1\.f*q* $f/

mv $dir/$f*[_R]2\.f*q* $f/


cd $f

if [ $3 == 'yes' ]
then
  singularity run --app trimmomatic039 /share/singularity/images/ccs/conda/amd-conda2-centos8.sinf trimmomatic PE  -threads 16 -phred33 -trimlog ${f}_errorlog.txt \
  $f*[_R]1\.f*q* $f*[_R]2\.f*q* \
  ${f}_R1_paired.fq ${f}_R1_unpaired.fq \
  ${f}_R2_paired.fq ${f}_R2_unpaired.fq \
  ILLUMINACLIP:/project/farman_s25abt480/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:20:20 MINLEN:90;
fi

# now run spades

singularity run --app spades3155 /share/singularity/images/ccs/conda/amd-conda9-rocky8.sinf spades.py --pe1-1 ${f}_R1_paired.fq --pe1-2 ${f}_R2_paired.fq --pe1-s ${f}_R1_unpaired.fq --pe1-s ${f}_R2_unpaired.fq -o ${f}_assembly

mv ${f}_assembly/contigs.fasta ${f}_assembly/${f}.fasta

mv ${f}_assembly/scaffolds.fasta ${f}_assembly/${f}_scaffolds.fasta

perl /project/farman_s25abt480/SCRIPTs/SimpleFastaHeaders.pl ${f}_assembly/${f}.fasta $f

perl /project/farman_s25abt480/SCRIPTs/CullShortContigs.pl ${f}_nh.fasta

# values for --app and path to trimmomatic/spades images will need to be updated for your system

