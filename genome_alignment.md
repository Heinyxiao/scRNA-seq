### Download Cellranger
From 10X Genomics Website: https://www.10xgenomics.com/support/software/cell-ranger

Log in to remote computer

```
cd /N/slate/xuexiao/cellranger/
tar -xzvf cellranger-4.0.0.tar.gz
export PATH=/N/slate/xuexiao/cellranger/cellranger-4.0.0:$PATH
cellranger count --id=CAF --fastqs=/N/slate/xuexiao/scDataset/CAF --sample=CAF --transcriptome=/N/slate/xuexiao/Genome/GRCh38/refdata-gex-GRCh38-2020-A
```
