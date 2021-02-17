#!/bin/bash


#reference:
#https://gist.github.com/darencard/72ddd9e6c08aaff5ff64ca512a04a6dd

cat ../bamfile_list | while read sample; do

echo "#$sample"
var1="$(echo $sample | cut -d '_' -f1)"
var2="$(echo $sample | cut -d '_' -f2-)"
echo "samtools sort -@ 15 -n ../${sample}.sorted.bam -o ${sample}.sorted.bam.qsort.bam"
echo "samtools sort -@ 15 -n ../../tophat_out/${var1}/${var2}/unmapped.bam -o ${sample}.unmapped.qsort.bam"
echo "samtools merge -@ 20 -n ${sample}.sorted.bam ${sample}.sorted.bam.qsort.bam ${sample}.unmapped.qsort.bam"

echo "samtools view -u -f 1 -F 12 ${sample}.sorted.bam > ${sample}.mapped.bam"
echo "samtools view -u -f 4 -F 264 ${sample}.sorted.bam >unmap_map.bam"
echo "samtools view -u -f 8 -F 260 ${sample}.sorted.bam >map_unmap.bam"
echo "samtools view -u -f 12 -F 256 ${sample}.sorted.bam >unmap_unmap.bam"

echo "samtools merge -u ${sample}.unmaped.bam unmap_map.bam map_unmap.bam unmap_unmap.bam"
echo "samtools sort -@ 20 -n ${sample}.mapped.bam -o ${sample}.mapped.sort.bam"
echo "samtools sort -@ 20 -n ${sample}.unmaped.bam -o ${sample}.unmapped.sort.bam"


# Get Mapped reads as FASTQ file
echo "bamToFastq -i ${sample}.mapped.sort.bam -fq ${sample}.mapped.1.fastq -fq2 ${sample}.mapped.2.fastq"

# Get unmapped reads as FASTQ file
echo "bamToFastq -i ${sample}.unmapped.sort.bam -fq ${sample}.unmapped.1.fastq -fq2 ${sample}.unmapped.2.fastq"

# Concatenate Reads pairs
echo "cat ${sample}.mapped.1.fastq ${sample}.unmapped.1.fastq >${sample}_1_dup.fastq"
echo "cat ${sample}.mapped.2.fastq ${sample}.unmapped.2.fastq >${sample}_2_dup.fastq"

echo "seqkit rmdup -j 20 -n ${sample}_1_dup.fastq >${sample}_1.fastq"
echo "seqkit rmdup -j 20 -n ${sample}_2_dup.fastq >${sample}_2.fastq"

# Compress FASTQ file
echo "pigz ${sample}_1.fastq"
echo "pigz ${sample}_2.fastq"


done
