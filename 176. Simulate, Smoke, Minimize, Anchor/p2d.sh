# all chromosomes with mean depth of coverage >= 19
samtools coverage sim.bam | awk '$7 >= 19 {print $0}'
