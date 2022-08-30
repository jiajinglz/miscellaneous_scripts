1. Prepare bam files for idr: rep1vsrep2, pool pseudo, rep1self, rep2self
################################################################################
##################USE BAM file that is PCR dup FREE#############################
################################################################################

#!/bin/bash

#SBATCH -p standard             ## partition/queue name
#SBATCH -A kwcho_lab            ## account to charge
#SBATCH --nodes=1               ## number of nodes the job will use
#SBATCH --ntasks=1              ## number of processes to launch
#SBATCH --cpus-per-task=16      ## number of OpenMP threads DO NOT CHANGE THIS FOR THIS PROCESS

module load samtools/1.10
SAMPLE=/dfs3b/mbt/jiajingz/chIPseq/foxh1/idr/bamprepare/pcrFree.st105_foxh1  ##prefix for sample format XX_rep1.bam/ XX_rep2.bam

#Merge BAMS for pool pseudo reps
samtools merge -u merged_$(basename $SAMPLE).bam $(basename $SAMPLE)_rep1.bam $(basename $SAMPLE)_rep2.bam
samtools view -H merged_$(basename $SAMPLE).bam > merged_$(basename $SAMPLE)_header.sam

#Split merged BAM
nlines=$(samtools view merged_$(basename $SAMPLE).bam | wc -l ) # Number of reads in the BAM file
nlines=$(( (nlines + 1) / 2 )) # half that number
samtools view merged_$(basename $SAMPLE).bam | shuf - | split -d -l ${nlines} - "merged_$(basename $SAMPLE)" # This will shuffle the lines in the file and split it into two SAM files named [whatever in ""]00 and 01 (follow the format) #
cat merged_$(basename $SAMPLE)_header.sam merged_$(basename $SAMPLE)00 | samtools view -bS - > merged_$(basename $SAMPLE)_00.bam
cat merged_$(basename $SAMPLE)_header.sam merged_$(basename $SAMPLE)01 | samtools view -bS - > merged_$(basename $SAMPLE)_01.bam

#Split for rep1 
samtools view -H $(basename $SAMPLE)_rep1.bam > $(basename $SAMPLE)_rep1_header.sam
nlines=$(samtools view $(basename $SAMPLE)_rep1.bam | wc -l ) # Number of reads in the BAM file
nlines=$(( (nlines + 1) / 2 )) # half that number
samtools view $(basename $SAMPLE)_rep1.bam | shuf - | split -d -l ${nlines} - "$(basename $SAMPLE)_rep1" # This will shuffle the lines in the file and split it into two SAM files
cat $(basename $SAMPLE)_rep1_header.sam $(basename $SAMPLE)_rep100 | samtools view -bS - > $(basename $SAMPLE)_rep1_00.bam
cat $(basename $SAMPLE)_rep1_header.sam $(basename $SAMPLE)_rep101 | samtools view -bS - > $(basename $SAMPLE)_rep1_01.bam

#Split for rep2
samtools view -H $(basename $SAMPLE)_rep2.bam > $(basename $SAMPLE)_rep2_header.sam
nlines=$(samtools view $(basename $SAMPLE)_rep2.bam | wc -l ) # Number of reads in the BAM file
nlines=$(( (nlines + 1) / 2 )) # half that number
samtools view $(basename $SAMPLE)_rep2.bam | shuf - | split -d -l ${nlines} - "$(basename $SAMPLE)_rep2" # This will shuffle the lines in the file and split it into two SAM files
cat $(basename $SAMPLE)_rep2_header.sam $(basename $SAMPLE)_rep200 | samtools view -bS - > $(basename $SAMPLE)_rep2_00.bam
cat $(basename $SAMPLE)_rep2_header.sam $(basename $SAMPLE)_rep201 | samtools view -bS - > $(basename $SAMPLE)_rep2_01.bam


2. Call peaks with each bam using MACS2: for good quality chip use default, so-so chip use p value of 0.001

for i in ../bamprepare/*.bam;do (echo "#! /bin/bash" &&
echo "#SBATCH -p standard" &&
echo "#SBATCH -A kwcho_lab " &&
echo "#SBATCH --nodes=1" &&
echo "#SBATCH --ntasks=1" &&
echo "#SBATCH --cpus-per-task=4" &&
echo "module load macs/2.2.7.1" &&
echo "macs2 callpeak -c /dfs3b/mbt/jiajingz/chIPseq/frog_chip_input_control/align/pcrFree_st8_input_foxh1.bam -t $i -f BAM -p 0.001 -g 1.1e+9 --outdir . -n $(basename $i .bam)_peak ") > $(basename $i .bam)_peakCall.sub;done


3. Sort narrowPeak file

for i in *st105*.narrowPeak;do (echo "#! /bin/bash" &&
echo "#SBATCH -p standard" &&
echo "#SBATCH -A kwcho_lab " &&
echo "#SBATCH --nodes=1" &&
echo "#SBATCH --ntasks=1" &&
echo "#SBATCH --cpus-per-task=1" &&
echo "module load bedtools2/2.29.2" &&
echo "bedtools sort -i $i > sorted.$(basename $i .narrowPeak).narrowPeak") > $(basename $i .narrowPeak)_peakSort.sub;done




4. IDR analysis


IDR local:

#!/bin/bash
## IDR thresholds default: 0.05 (similar to idea of FDR)

SAMPLE=pcrFree.st105_foxh1  ####change this everytime for new samples####

## Nt= rep1 vs rep2 true reps compare
idr --samples sorted.${SAMPLE}_rep1_peak_peaks.narrowPeak sorted.${SAMPLE}_rep2_peak_peaks.narrowPeak --input-file-type narrowPeak --rank signal.value --output-file ${SAMPLE}_rep1_vs_rep2_idr --plot --log-output-file ${SAMPLE}_rep1_vs_rep2_idr_log
##filter out peaks idr < 0.05
awk '{if($5 >= 540) print $0}' ${SAMPLE}_rep1_vs_rep2_idr > ${SAMPLE}_rep1_vs_rep2_idr_0.05_filtered.narrowPeak

## Np= pooled rep1 and rep2 then split for 2 pseudo reps to compare
idr --samples sorted.merged_${SAMPLE}_00_peak_peaks.narrowPeak sorted.merged_${SAMPLE}_01_peak_peaks.narrowPeak --input-file-type narrowPeak --rank signal.value --output-file ${SAMPLE}_pooled_pseudo00_vs_pseudo01_idr --plot --log-output-file ${SAMPLE}_pooled_pseudo00_vs_pseudo01_idr_log
##filter out peaks idr < 0.05
awk '{if($5 >= 540) print $0}' ${SAMPLE}_pooled_pseudo00_vs_pseudo01_idr > ${SAMPLE}_pooled_pseudo00_vs_pseudo01_idr_0.05_filtered.narrowPeak

## N1= rep1 splits into 2 pseudo samples to compare
idr --samples sorted.${SAMPLE}_rep1_00_peak_peaks.narrowPeak sorted.${SAMPLE}_rep1_01_peak_peaks.narrowPeak --input-file-type narrowPeak --rank signal.value --output-file ${SAMPLE}_rep1_pseudo00_vs_pseudo01_idr --plot --log-output-file ${SAMPLE}_rep1_pseudo00_vs_pseudo01_idr_log
##filter out peaks idr < 0.05
awk '{if($5 >= 540) print $0}' ${SAMPLE}_rep1_pseudo00_vs_pseudo01_idr > ${SAMPLE}_rep1_pseudo00_vs_pseudo01_idr_0.05_filtered.narrowPeak

## N2= rep2 splits into 2 pseudo samples to compare
idr --samples sorted.${SAMPLE}_rep2_00_peak_peaks.narrowPeak sorted.${SAMPLE}_rep2_01_peak_peaks.narrowPeak --input-file-type narrowPeak --rank signal.value --output-file ${SAMPLE}_rep2_pseudo00_vs_pseudo01_idr --plot --log-output-file ${SAMPLE}_rep2_pseudo00_vs_pseudo01_idr_log
##filter out peaks idr < 0.05
awk '{if($5 >= 540) print $0}' ${SAMPLE}_rep2_pseudo00_vs_pseudo01_idr > ${SAMPLE}_rep2_pseudo00_vs_pseudo01_idr_0.05_filtered.narrowPeak
#### quality control for reps
####(N1/N2 >= 2)& (Np/Nt >= 2) then it is low quality rep or low reproducibility




5. Find peak summit for IDR peaks:


1. sort 0.05 filtered IDR peaks by chromosome order

bedtools sort -i input.peak > output.peak


2. intersect with merged bam peak call result file: narrowPeak and summits.bed
!!!!if number of peaks do not match, play around with -f until number matches!!!!
eg. st9 hdac1 i used -f 0.511725
    st10.5 hdac1 i used -f 0.369

bedtools intersect -a merged.narrowPeak -b sorted.IDR.peaks -wa -f 0.2 -nonamecheck > result.narrowPeak

bedtools merge -i result.narrowPeak > final.narrowPeak

bedtools intersect -a merged.summits.bed -b final.narrowPeak -wa -nonamecheck > final.summits.bed












