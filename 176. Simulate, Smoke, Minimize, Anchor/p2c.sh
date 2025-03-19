# index bam file
samtools index sim.bam
# index reference genome
samtools faidx sacCer3.fa
# inspect alignment pileup
samtools mpileup --reference sacCer3.fa sim.bam
