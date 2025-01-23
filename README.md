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
### Calculates sites with depth >0, q20,30,40,50,60
### q42 = uniquely mapped to a single place in the genome

```bash

conda info --envs  ## list all conda environments

conda activate /usr/local/src/conda_envs/conda/envs/samtools  ## Activate samtools environment

for bam in *.bam; do
    sample=$(basename "$bam" .bam)

    # Compute mean and SD for depth using samtools depth
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

    # Count total variable sites (positions with >0 coverage)
    total_sites=$(samtools depth "$bam" | awk '$3>0 {count++} END {print count}')

    # Count sites with quality > 20 and > 30 using samtools mpileup
	# May need to adjust relative (or absolute) path to reference fasta file!
    sites_q20=$(samtools mpileup -Q 20 -f ../../reference/final_reference_assembly.fa "$bam" | awk '$6!~"^*$" {count++} END {print count}')
    sites_q30=$(samtools mpileup -Q 30 -f ../../reference/final_reference_assembly.fa "$bam" | awk '$6!~"^*$" {count++} END {print count}')
    sites_q40=$(samtools mpileup -Q 40 -f ../../reference/final_reference_assembly.fa "$bam" | awk '$6!~"^*$" {count++} END {print count}')

    echo -e "${sample}\t${mean_depth}\t${sd_depth}\t${total_sites}\t${sites_q20}\t${sites_q30}\t${sites_q40}"
done > ../mean_coverage_quality.txt
```

