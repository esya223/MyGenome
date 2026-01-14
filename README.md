# MyGenome

_Pyricularia orzyzae_ or _Magnaporthe oryzae (2002)_, causes a devastating fungal disease that affects rice plants by creating lesions on leaves, stems and panicles. Through identifying and sequencing this fungus, information can be concluded about the pathogen's evolution, genetic variability, and host potential. _P. oryzae_ used for this data was collected from the host _Urochloa mutica_ in Nueve Ecija in the Philipines and was sequenced via ILLUMINA NovaSeq X. The assembled genome was submitted to NCBI (BioSample ID: SAMN46924042, SRA Accession: SRR32567630, BioProject: PRJNA1226220).

## Sequence Quality Assessment

<summary>Run forward and reverse sequences through FastQC for quality assessment.

```
fastqc Bm88312_1.fq.gz
fastqc Bm88312_2.fq.gz
```

## Trimming

<summary>Run relevant Trimmomatic command to trim forward and reverse sequences.

```
java -jar trimmomatic-0.38.jar PE -threads 2 \-phred33 -trimlog Br88312_errorlog.txt \ Br88312_1.fq.gz \ Br88312_2.fq.gz \ Br88312_1_paired.fq.gz \ Br88312_1_unpaired.fq.gz \Br88312_2_paired.fq.gz \ Br88312_2_unpaired.fq.gz \ ILLUMINACLIP:adaptors.fa:2:30:10 SLIDINGWINDOW:20:20 MINLEN:150
```

<summary>Run paired sequences through FASTQC to assess the quality of trimming.
#### Output

```
Bm88312_1_paired_fastqc.zip
Bm88312_2_paired_fastqc.zip
```
```
unzip Bm88312_1_paired_fastqc.zip and Bm88312_2_paired_fastqc.zip
```

#### Full FastQC Results 
<details>[Bm88312_1_paired_fastqc.html](https://github.com/esya223/MyGenome/blob/main/Bm88312_1_paired_fastqc.html) 
[Bm88312_2_paired_fastqc.html](https://github.com/esya223/MyGenome/blob/main/Bm88312_2_paired_fastqc.html)</details>

| Sequence | Reads |
| -------- | ----- |
|Raw Reads (single end) | 7,234,892|
|Cleaned Reads (single end)| 6,356,443|
|Total Bases in Cleaned Reads| 1,901,184,532|

## Assembly

<summary> Assembly was completed on MCC
<summary> scp was used to transfer trimmed sequences from the virtual machine to MCC.
<summary> Optimal k-mer values for the genome assembly was determined with Velvet Advisor (https://dna.med.monash.edu.au/~torsten/velvet_advisor/)
<summary> VelvetOptimiser ran at step 10 an step 2 to find the k-value. 


```
sbatch velvetoptimiser_noclean.sh Bm88312_step2 69 89 2
```
|Step|Genome Size|# of Contigs|N50|
|----|-----------|------------|---|
|2   |40,854,765 |2,210 |113,184| 

```
sbatch velvetoptimiser_noclean.sh Bm88312_step10 69 89 10
```

|Step|Genome Size|# of Contigs|N50|
|----|-----------|------------|---|
|10  |40,859,984 |2,197 |120,975| 

The highest fold coverage was 46.50 

OUTPUT:

```
slurm-30186556.out
slurm-30186446.out
contigs.fa
```

```
cp /project/farman_s25abt480/SCRIPTs/SimpleFastaHeaders.pl /project/farman_s25abt480/esya223
perl SimpleFastaHeaders.pl Bm88312.fasta Bm88312
```
OUTPUT:
Bm88312_nh.fasta

cd ./Bm88312_step2/velvet_Bm88312_step2_69_89_2_noclean/

```
cp /project/farman_s25abt480/SCRIPTs/CullShortContigs.pl ./ 
perl CullShortContigs.pl Bm88312_nh.fasta
```

OUTPUT:
Bm88312_final.fasta

```
cp /project/farman_s25abt480/SCRIPTs/SeqLen.pl ./
perl SeqLen.pl Bm88312_final.fasta
```

```
awk '/^>/ {if (seqlen >= 200) {print header; print seq}; header=$0; seq=""; seqlen=0; next} {seq=seq $0; seqlen+=length($0)} END {if (seqlen >= 200) {print header; print seq}}' Bm88312_nh.fasta > Bm88312_final.fasta
```

```
awk '{print $2}' seq_stats.txt | awk '
BEGIN {n=0; sum=0; max=0; min=1e9}
{n++; sum+=$1; if($1>max) max=$1; if($1<min) min=$1}
END {
  print "Total contigs:", n;
  print "Genome size (bp):", sum;
  print "Longest contig:", max;
  print "Shortest contig:", min;
  print "Average contig length:", int(sum/n)}'
```

```
cp /project/farman_s25abt480/SCRIPTs/BuscoSingularity.sh ./
sbatch /project/farman_s24cs485g/SLURM_SCRIPTS/BuscoSingularity.sh path/to/Bm88312_final.fasta
```
  
```
singularity exec busco_latest.sif \
busco -i Bm88312_final.fasta \
-l fungi_odb10 \
-o busco_Bm88312 \
-m genome \
--cpu 8
```

OUTPUT:
C:98.5%[S:98.2%,D:0.3%],F:0.8%,M:0.7%,n:758

myVM
scp ngs@10.163.183.71:/home/ngs/Desktop/MoMitochondrion.fasta ~/blast/

wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz

OUTPUT:
ncbi-blast-2.16.0+-x64-linux.tar.gz
tar -zxvpf ncbi-blast-2.16.0+-x64-linux.tar.gz
./blastn -version

cd ~/blast
blastn -query MoMitochondrion.fasta -subject Bm88312_nh.fasta -evalue 1e-50 -max_target_seqs 20000 -outfmt '6 qseqid sseqid slen length qstart qend sstart send btop' > B71v2sh.Bm88312.BLAST

OUTPUT:
B71v2sh.Bm88312.BLAST

awk '$4/$3 > 0.9 {print $2 ",mitochondrion"}' B71v2sh.Bm88312.BLAST > Bm88312_mitochondrion.csv

OUTPUT:
Bm88312_mitochondrion.csv

ssh esya223@mcc.uky.edu

cp /project/farman_s25abt480/FASTA/B71v2sh_masked.fasta /project/farman_s25abt480/esya223/

blastn -query B71v2sh_masked.fasta -subject Bm88312_final.fasta -evalue 1e-50 -max_target_seqs 20000 -outfmt '6 qseqid sseqid qstart qend sstart send btop' -out B71v2sh.Bm88312.BLAST

OUTPUT:
B71v2sh.Bm88312.BLAST

cp B71v2sh.Bm88312.BLAST /project/farman_s25abt480/CLASS_BLASTs/

scp ngs@10.163.183.71:/Users/ngs/Desktop/B71ref2.fasta .
scp ngs@10.163.183.71:/Users/ngs/Desktop/B71Ref2_a0.3.gff3 .

screen -S genes bash -l
echo '##FASTA' | cat B71Ref2_a0.3.gff3 - B71Ref2.fasta > B71Ref2.gff3
grep '##FASTA' -B 5 -A 5 B71Ref2.gff3
maker2zff B71Ref2.gff3
fathom genome.ann genome.dna -gene-stats
fathom genome.ann genome.dna -categorize 1000
fathom uni.ann uni.dna -export 1000 -plus
forge export.ann export.dna
hmm-assembler.pl Moryzae . > Moryzae.hmm
snap-hmm Moryzae.hmm Bm88312_final.fasta > Bm88312_final-snap.zff
fathom Bm88312_final-snap.zff Bm88312_final.fasta -gene-stats
snap-hmm Moryzae.hmm Bm88312_final.fasta -gff > Bm88312_final-snap.gff2

OUTPUT:
Bm88312_final-snap.zff
Bm88312_final-snap.gff2

snap-hmm Moryzae.hmm Bm88312_final.fasta > Bm88312_final-snap.zff
fathom Bm88312_final-snap.zff Bm88312_final.fasta -gene-stats
OUTPUT:
Bm88312_final-augustus.gff3
gff3_merge -d Bm88312_final.maker.output/Bm88312_final_master_datastore_index.log \
-o Bm88312_final-annotations.gff

OUTPUT:
Bm88312_final-annotations.gff 

maker -CTL

fasta_merge -d Bm88312_final.maker.output/Bm88312_final_master_datastore_index.log -o Bm88312_final.maker.proteins.fasta

OUTPUT:
Bm88312_final.maker.proteins.fasta
