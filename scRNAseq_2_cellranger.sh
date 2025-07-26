#used for all three samples

id=sampleid
sample=myprefix

cellranger count \
  --id=$id \
  --fastqs=/pathto_fastqst \
  --transcriptome=/pathto/refdata-gex-GRCh38-2020-A \
  --sample=$sample
  
