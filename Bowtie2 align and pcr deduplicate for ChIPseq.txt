#############alignment#################

genome index pathway: 

/dfs3b/mbt/jiajingz/genomeREF/xenbase10.0/bowtie2_build/xenbase10.0_bt2
/dfs3b/mbt/jiajingz/genomeREF/dMelanogaster_BDGP6/bowtie2_build/dm_ensembl_BDGP6_bt2
/dfs3b/mbt/jiajingz/genomeREF/hSapiens_hg38_ensembl/bowtie2-build/Hsapiens_hg38-bowtie2
/dfs3b/mbt/jiajingz/genomeREF/mMusculus_grcm38_ensembl/bowtie2_build/mMusculus38_bt2


###############single-end#####################
for i in ../rawdata/*.fastq.gz;do (echo "#! /bin/bash" &&
echo "#SBATCH -p standard" &&
echo "#SBATCH -A kwcho_lab" &&
echo "#SBATCH --nodes=1" &&
echo "#SBATCH --ntasks=1" &&
echo "#SBATCH --cpus-per-task=8" &&
echo "module load bowtie2/2.4.1" &&
echo "module load samtools/1.10" &&
echo "bowtie2 -p 8 -x /dfs7/mbt/jiajingz/genomeREF/xenbase10.0/bowtie2_build/xenbase10.0_bt2 $i -S $(basename $i .fastq.gz).sam" && 
echo "samtools view -bS $(basename $i .fastq.gz).sam > $(basename $i .fastq.gz).bam" &&
echo "samtools view -F 1804 -q 30 -b $(basename $i .fastq.gz).bam -o uniq_$(basename $i .fastq.gz).bam" &&
echo "samtools sort -n uniq_$(basename $i .fastq.gz).bam -o name.sorted.uniq_$(basename $i .fastq.gz).bam" &&
echo "samtools fixmate -m name.sorted.uniq_$(basename $i .fastq.gz).bam fixmate.name.sorted.uniq_$(basename $i .fastq.gz).bam" &&
echo "samtools sort fixmate.name.sorted.uniq_$(basename $i .fastq.gz).bam -o position.sorted.fixmate.name.sorted.uniq_$(basename $i .fastq.gz).bam" &&
echo "samtools markdup position.sorted.fixmate.name.sorted.uniq_$(basename $i .fastq.gz).bam -r pcrFree.$(basename $i .fastq.gz).bam")> $(basename $i .fastq.gz)_alignDedup.sub;done

#########################################################################################################
#########################################################################################################

#############alignment#################
###############paired-end#####################

for i in ../rawdata/*READ1.fastq.gz;do (echo "#! /bin/bash" &&
echo "#SBATCH -p standard" &&
echo "#SBATCH -A kwcho_lab" &&
echo "#SBATCH --nodes=1" &&
echo "#SBATCH --ntasks=1" &&
echo "#SBATCH --cpus-per-task=8" &&
echo "module load bowtie2/2.4.1" &&
echo "module load samtools/1.10" &&
echo "bowtie2 -p 8 -x /dfs7/mbt/jiajingz/genomeREF/xenbase10.0/bowtie2_build/xenbase10.0_bt2 -1 $i -2 ../rawdata/$(basename $i -READ1.fastq.gz)-READ2.fastq.gz -S $(basename $i -READ1.fastq.gz).sam" && 
echo "samtools view -bS $(basename $i -READ1.fastq.gz).sam > $(basename $i -READ1.fastq.gz).bam" &&
echo "samtools view -F 1804 -f 2 -q 30 -b $(basename $i -READ1.fastq.gz).bam -o uniq_$(basename $i -READ1.fastq.gz).bam" &&
echo "samtools sort -n uniq_$(basename $i -READ1.fastq.gz).bam -o name.sorted.uniq_$(basename $i -READ1.fastq.gz).bam" &&
echo "samtools fixmate -m name.sorted.uniq_$(basename $i -READ1.fastq.gz).bam fixmate.name.sorted.uniq_$(basename $i -READ1.fastq.gz).bam" &&
echo "samtools sort fixmate.name.sorted.uniq_$(basename $i -READ1.fastq.gz).bam -o position.sorted.fixmate.name.sorted.uniq_$(basename $i -READ1.fastq.gz).bam" &&
echo "samtools markdup position.sorted.fixmate.name.sorted.uniq_$(basename $i -READ1.fastq.gz).bam -r pcrFree.$(basename $i -READ1.fastq.gz).bam")> $(basename $i -READ1.fastq.gz)_alignDedup.sub;done


#############Extract unqiue reads single end reads##################
#############follow encode guide###################

use MAPQ score of 30

##################extract unique reads paired-end reads################
#############follow encode guide####################

############################################################################################################
#######################################clean up after alignment#############################################
############################################################################################################
nano clean.sh
####copy and paste following script####

#!/bin/bash
#check alignment % to make sure alignment is okay##
cat *.out > pcr_alignment_stat_summary.bamtxt
#remove files generated in intermediate steps##
#save only PCR deduplicated bam file ##
find . ! -name 'pcr*.bam*' -type f -exec rm -f {} +



