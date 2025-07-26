i=samplename

hisat2 --dta \
    -p 8 \
    -x /pathto/hisat2_index \
    -1 *_R1_001.fastq.gz \
    -2 *_R2_001.fastq.gz \
    -S $i.sam

samtools view -F 4 -Su $i.sam | samtools sort -T $i.accepted_hits -o $i.accepted_hits.bam

stringtie $i.accepted_hits.bam \
-G /pathto/gencode.v40.annotation1.gtf \
-e \
-o $i.gtf \
-p 8
