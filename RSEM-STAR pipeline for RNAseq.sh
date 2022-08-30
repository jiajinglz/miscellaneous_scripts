This allows you to align with STAR and quantify transcript.bam using RSEM
to do this instead of using star to build its own index cuz i have faced incompatibility problem


1. Assemble STAR reference index using RSEM, need to load both STAR and RSEM (same thing can be done with bowtie2):

#!/bin/bash
#SBATCH -p standard             ## partition/queue name
#SBATCH -A kwcho_lab            ## account to charge
#SBATCH --nodes=1               ## number of nodes the job will use
#SBATCH --ntasks=1              ## number of processes to launch
#SBATCH --cpus-per-task=8       ## number of OpenMP threads

mkdir ./star_build_by_RSEM
module load star/2.7.3a
module load rsem/1.3.3

rsem-prepare-reference -p 8 --gtf XEN10_gtf_modified_for_rsem_star_genome_assembly.gtf --star ./XENTR_10.0_genome.fasta ./star_build_by_RSEM/xen10_star_by_rsem


2. Align using RSEM but through STAR: Use STAR to align reads. Alignment parameters are from ENCODE3's STAR-RSEM pipeline. To save computational time and memory resources, STAR's Output BAM file is unsorted. It is stored in RSEM's temporary directory with name as 'sample_name.bam'. Each STAR job will have its own private copy of the genome in memory.

--------single-end-----------------
for i in ../rawdata/*.fastq.gz;do (echo "#! /bin/bash" &&
echo "#SBATCH -p standard" &&
echo "#SBATCH -A kwcho_lab" && 
echo "#SBATCH --nodes=1" &&          
echo "#SBATCH --ntasks=1" &&
echo "#SBATCH --cpus-per-task=8" &&
echo "module load rsem/1.3.3" &&
echo "module load star/2.7.3a" &&
echo "rsem-calculate-expression --star --star-gzipped-read-file -p 8 $i /dfs3b/mbt/jiajingz/genomeREF/xenbase10.0/star_build_by_RSEM/xen10_star_by_rsem --no-bam-output $(basename $i .fastq.gz)") > $(basename $i .fastq.gz)_align.sub;done

-------paired-end--------------------
for i in ../rawdata/*READ1.fastq.gz;do (echo "#! /bin/bash" &&
echo "#SBATCH -p standard" &&
echo "#SBATCH -A kwcho_lab" && 
echo "#SBATCH --nodes=1" &&
echo "#SBATCH --ntasks=1" &&
echo "#SBATCH --cpus-per-task=8" &&
echo "module load rsem/1.3.3" &&
echo "module load star/2.7.3a" &&
echo "rsem-calculate-expression --star --star-gzipped-read-file --paired-end -p 8 $i ../rawdata/$(basename $i READ1.fastq.gz)READ2.fastq.gz /dfs7/mbt/jiajingz/genomeREF/xenbase10.0/star_build_by_RSEM/xen10_star_by_rsem --no-bam-output $(basename $i -READ1.fastq.gz)") > $(basename $i -READ1.fastq.gz)_align.sub;done


3. counts and TPM are output in genes.result file

paste ... ... ... > all.txt

extract only TPM:
awk '{print $1 "\t" $6 "\t" $13 "\t" $20 "\t" $27 "\t" $34 "\t" $41 "\t" $48 "\t" $55 "\t" $62 "\t" $69 "\t" $76 "\t" $83 "\t" $90 "\t" $97 "\t" $104 "\t" $111 "\t" $118 "\t" $125 "\t" $132 "\t" $139 "\t" $146 "\t" $153 "\t" $160 "\t" $167}'

extract only counts:
awk '{print $1 "\t" $5 "\t" $12 "\t" $19 "\t" $26 "\t" $33 "\t" $40 "\t" $47 "\t" $54 "\t" $61 "\t" $68 "\t" $75 "\t" $82 "\t" $89 "\t" $96 "\t" $103 "\t" $110 "\t" $117 "\t" $124 "\t" $131 "\t" $138 "\t" $145 "\t" $152}'







