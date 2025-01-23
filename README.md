# Randomly_Useful_Commands
Mostly useful one-liners or commands in bash/unix

## Calculate the mean # of reads in all fastq files in a directory +/- SD
```bash
find fastp -type f -name '*R1.fastq' | parallel --bar "awk '(NR%4==1){count++} END {print \"{}\", count}' {}" | tee read_counts.txt | awk '{sum+=$2; sumsq+=$2*$2; n++} END {mean=sum/n; sd=sqrt(sumsq/n - mean^2); print "Mean:", mean, "SD:", sd}'
```
## Same as aove but also include median, and remove "Undetermined" files
```bash
find fastp -type f -name '*R1.fastq' | grep -Ev "(Undetermined-S20_R1.fastq|Undetermined-S20_R2.fastq)" | parallel -j 80 --bar "awk '(NR%4==1){count++} END {print \"{}\", count}' {}" | tee read_counts.txt | awk '
{counts[NR] = $2; sum+=$2; sumsq+=$2*$2; n++}
END {
    # Compute mean & standard deviation
    mean = sum / n;
    sd = sqrt(sumsq / n - mean^2);

    # Sort the counts array and compute the median
    asort(counts);
    if (n % 2 == 1) {
        median = counts[int(n/2) + 1];  # Odd number of values
    } else {
        median = (counts[n/2] + counts[n/2 + 1]) / 2;  # Even number of values
    }

    # Print results
    print "Mean:", mean, "SD:", sd, "Median:", median;
}'
```

## Transfer files between servers with rsync
```bash
 rsync -av -P /data/ccorbett/reads/ fhernan@narval.calculcanada.ca:/home/fhernan/scratch/stiltgrass/

```

## Calculate the mean +/- SD of coverage depth and number of variants for all .bam files in the current working directory

```bash
for bam in *.bam; do
    sample=$(basename "$bam" .bam)
    
    # Compute mean and SD for depth
    read mean_depth sd_depth <<< $(samtools depth "$bam" | 
        awk '{sum+=$3; sumsq+=$3*$3} 
             END { 
                 if (NR>0) { 
                     mean=sum/NR; 
                     sd=sqrt(sumsq/NR - mean^2); 
                     print mean, sd 
                 } else { 
                     print 0, 0 
                 } 
             }')

    # Count non-zero covered sites as a proxy for SNPs
    snp_count=$(samtools depth "$bam" | awk '$3>0 {count++} END {print count}')

    echo -e "${sample}\t${mean_depth}\t${sd_depth}\t${snp_count}"
done > ../mean_coverage_snps.txt
```

