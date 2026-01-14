#!/bin/bash

#SBATCH --time 48:00:00
#SBATCH --job-name=VelvetOptimiser
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=16
#SBATCH --partition=normal
#SBATCH --mem=180GB
#SBATCH --mail-type ALL
#SBATCH -A cea_farman_s25abt480
#SBATCH --mail-type ALL
#SBATCH --mail-user farman@uky.edu,szh231@uky.edu

echo "SLURM_NODELIST: "$SLURM_NODELIST

strainID=$1

lowK=$2

highK=$3

step=$4

mkdir $strainID 

  cp $strainID*1_paired*f*q* $strainID/

  cp $strainID*2_paired*f*q* $strainID/

  cd $strainID

# create hard-coded read names

  cp $strainID*1_paired*f*q* forward.fq

  cp $strainID*2_paired*f*q* reverse.fq

# load python2 module

#module load ccs/conda/python-2.7.15

# generate a dataset with interleaved R1 and R2 reads 

#  python2 /project/farman_s23cs485g/interleave-fastq.py *1_paired.fq *2_paired.fq > interleaved.fq

# unload python2 environment

#  module unload ccs/conda/python-2.7.15

# load singularity container

  module load ccs/singularity

# run velvetoptimiser in singularity

  singularity run --app perlvelvetoptimiser226 /share/singularity/images/ccs/conda/amd-conda2-centos8.sinf VelvetOptimiser.pl \
  -s $lowK -e $highK -x $step -d velvet_${strainID}_${lowK}_${highK}_${step}_noclean -o ' -clean no' -f ' -shortPaired -fastq -separate forward.fq reverse.fq'

# remove read files with hard-coded names

  rm forward.fq

  rm reverse.fq 

# rename optimal assembly from contigs.fa to <strainID>.fasta

  mv velvet_${strainID}_${lowK}_${highK}_${step}/contigs.fa velvet_${strainID}_${lowK}_${highK}_${step}/${strainID}_${lowK}_${highK}_${step}".fasta"

# copy assembly to a class ASSEMBLIES directory

  cp velvet_${lowK}_${highK}_${step}/${strainID}".fasta" /project/farman_s24cs485g/ASSEMBLIES/${strainID}_${lowK}_${highK}_${step}".fasta"

# copy logfile to class ASSEMBLIES directory and add strain identifier to the name

  logfile=`ls velvet_${lowK}_${highK}_${step}/*Logfile.txt`

  cp $logfile /project/farman_s23cs485g/ASSEMBLIES/${strainID}_${lowK}_${highK}_${step}_${logfile/*\//}

# backup assemblies and other output files to cloud based storage:

  #rclone copy ASSEMBLIES/${f}".fasta" GoogleDrive:LCC/ASSEMBLIES

  #rclone copy velvet_assembly GoogleDrive:LCC/${f}_velvet_assembly

# return to starting directory

  cd ..



