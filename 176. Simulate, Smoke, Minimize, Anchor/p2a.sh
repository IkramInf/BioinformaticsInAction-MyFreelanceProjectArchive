# sort sam file and store results into bam file
samtools view -bS sim.sam | samtools sort -o sim.bam
