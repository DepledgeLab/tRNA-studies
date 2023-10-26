#### This folder and its subfolder contain datasets and scripts associated with the [nano-tRNA Sequencing](https://github.com/novoalab/Nano-tRNAseq) used in this study.


#### R packages required
- ggplot2
- dplyr
- patchwork
- DESeq2
- EnhancedVolcano

#### Overview of data generation and processing

Nanopore datasets were generated using the [nano-tRNAseq](https://github.com/novoalab/Nano-tRNAseq) protocol and basecalling was performed without trimming using Guppy v6.5.7

```
guppy_basecaller -i nanopore_output_folder -s basecalled_data-hac.v6.5.7 -c rna_r9.4.1_70bps_hac.cfg -r --calib_detect --trim_strategy none --reverse_sequence true --u_substitution true -x auto
```

Sequence reads were subsequently aligned against a tRNA database in which adaptor sequences were added to the 5' and 3' of each tRNA sequence

```
# force/verify U to T substition in input fastq files
sed -i 's/U/T/g' basecalled-hac.v6.5.7.pass.notrim.fastq
sed -i 's/U/T/g' basecalled-hac.v6.5.7.fail.notrim.fastq

# align w/ bwa mem
bwa mem -W 13 -k 6 -T 20 -x ont2d Pa14_tRNAs_complete.fasta $BASE/"$i"-hac.v6.5.7.pass.notrim.fastq > aligned-pass.w13k6t20.all.sam
bwa mem -W 13 -k 6 -T 20 -x ont2d Pa14_tRNAs_complete.fasta $BASE/"$i"-hac.v6.5.7.fail.notrim.fastq > aligned-fail.w13k6t20.all.sam

# Parse passed reads and retain only primary alignments
samtools view -F 2324 -b -o aligned-pass.w13k6t20.all.bam aligned-pass.w13k6t20.all.sam
samtools sort -o aligned-pass.w13k6t20.all.sorted.bam aligned-pass.w13k6t20.all.bam
samtools index aligned-pass.w13k6t20.all.sorted.bam
bamToBed -i aligned-pass.w13k6t20.all.sorted.bam > aligned-pass.w13k6t20.all.sorted.bed

# Parse fail reads and retain only primary alignments
samtools view -F 2324 -b -o aligned-fail.w13k6t20.all.bam aligned-fail.w13k6t20.all.sam
samtools sort -o aligned-fail.w13k6t20.all.sorted.bam aligned-fail.w13k6t20.all.bam
samtools index aligned-fail.w13k6t20.all.sorted.bam
bamToBed -i aligned-fail.w13k6t20.all.sorted.bam > aligned-fail.w13k6t20.all.sorted.bed

# Merge pass and fail alignments
samtools merge -o aligned-merged.w13k6t20.all.bam aligned-pass.w13k6t20.all.sorted.bam aligned-fail.w13k6t20.all.sorted.bam
samtools sort -o aligned-merged.w13k6t20.all.sorted.bam aligned-merged.w13k6t20.all.bam
samtools index aligned-merged.w13k6t20.all.sorted.bam
bamToBed -i aligned-merged.w13k6t20.all.sorted.bam > aligned-merged.w13k6t20.all.sorted.bed

# Filter to remove all alignments with mapQ of 0 (multiple primary alignment)
samtools view -b -q 1 -o aligned-merged.w13k6t20.all.q1.sorted.bam aligned-merged.w13k6t20.all.sorted.bam
samtools index aligned-merged.w13k6t20.all.q1.sorted.bam
bamToBed -i aligned-merged.w13k6t20.all.q1.sorted.bam > aligned-merged.w13k6t20.all.q1.sorted.bed

# generate counts
samtools view aligned-merged.w13k6t20.all.q1.sorted.bam | cut -f3 | grep -v LN | sort | uniq -c | sed "s/^[ \t]*//" | sed "s/ / \t/g" > aligned-merged.counts.txt
```

#### Software used
- BWA v0.7.17
- SAMtools v1.15.1
- BEDTools v2.27.1
- Python v3.7.0
- GLib v2.72.1



