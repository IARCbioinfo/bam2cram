# bam2cram
A repository to convert BAM->CRAM files

## run examples

```
#using a tsv file, with custom working directory
nextflow run main.nf --bam_csv test/bams.txt --fasta test/ref.fa --output_folder csv2 -w tmp
#using a directory with bams
nextflow run main.nf --bams test --fasta test/ref.fa --output_folder dirs_bam
```
### tsv file
A tabular file with the following colums:

```
label   bam     index
S4      /path/S4.bam       /path/S4.bam.bai
S3      /path/S3.bam       /path/S3.bam.bai
```
