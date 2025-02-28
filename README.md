# MyGenome
$env:DISPLAY = “localhost:0”

ssh -Y LinkBlue@LinkBlue.cs.uky.edu

cd MyGenome

fastqc &

java -jar trimmomatic-0.38.jar PE -threads 2 \-phred33 -trimlog Br88312_errorlog.txt \ Br88312_1.fq.gz \ Br88312_2.fq.gz \ Br88312_1_paired.fq.gz \ Br88312_1_unpaired.fq.gz \Br88312_2_paired.fq.gz \ Br88312_2_unpaired.fq.gz \ ILLUMINACLIP:adaptors.fa:2:30:10 SLIDINGWINDOW:20:20 MINLEN:150
