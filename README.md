# Randomly_Useful_Commands
Mostly useful one-liners or commands in bash/unix


## Using 'screen' to set up long analyses that you can leave and come back to.

```bash

### The screen command will open a new terminal session. You can name the session anything. e.g.

# Open a new screen session called 'ml_analysis'
screen -S ml_analysis    # capital S

# Run your command 

# press 'ctrl+a+d' simultaneously (no plus signs) and then press enter, to detach from the session.

# you can remove all the previous text from the screen with:

clear      # then press enter

# reattach to the screen session, and you'll see everything since you left.

screen -r ml_analysis   # lower-case r

# press ctrl+a+d   # to detach again

# check all open screen sessions and their status

screen -ls

# output:
# (base) [cfb0001@LSB5211 phylonet]$ screen -ls
# There are screens on:
#         1306908.phylonetworks   (Detached)
#         892003.qtgild   (Detached)
# 2 Sockets in /run/screen/S-cfb0001.

# read the manual page for screen (or any other bash/unix command)

man screen   # press 'q' then press enter to exit manual.

```

## Download files directly from Illumina BaseSpace to a server
### Recommend using screen or nohup to run in background
```bash
mkdir -p /Data/cbarrett/stiltgrass_ref_panel/reads
cd /Data/cbarrett/stiltgrass_ref_panel/reads
bs download project -n "Barrett_Cam_Redemux_20260228" -o . --include="*.fastq.gz"

# Once the download finishes, pull all .fastq.gz file out of the subfolders and delete the empty subfolders
mv */*.fastq.gz .
rmdir */

```

## Calculate the mean # of reads in all fastq files in a directory +/- SD
```bash
find fastp -type f -name '*R1.fastq' | parallel --bar "awk '(NR%4==1){count++} END {print \"{}\", count}' {}" | tee read_counts.txt | awk '{sum+=$2; sumsq+=$2*$2; n++} END {mean=sum/n; sd=sqrt(sumsq/n - mean^2); print "Mean:", mean, "SD:", sd}'
```
## Same as above but also include median, and remove "Undetermined" files
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
# transfer from /data/craigdata/reads/ on current server to /home/craignewfolder/data/reads/ on remote server

 rsync -av -P /data/craigdata/reads/ craigbarrett@servername.wvu.edu:/home/craignewfolder/data/reads/

```

## Calculate the mean +/- SD of coverage depth and number of variants for all .bam files in the current working directory
### Calculates sites with depth >0, q20,30,40,50,60
### Capped at q42: likely uniquely mapped to a single place in the genome

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

## Count the number of reads for gzipped fastq files (R1 only) and create a tsv file with the output

```bash
#!/bin/bash

# Output file
output="read_counts.tsv"
echo -e "filename\tread_count" > "$output"

# Export a function to run in parallel
count_reads() {
    file="$1"
    count=$(zcat "$file" | wc -l)
    reads=$((count / 4))
    echo -e "$(basename "$file")\t$reads"
}

export -f count_reads

# Run in parallel over all _R1.fastq.gz files and append results to the output file
parallel count_reads ::: reads/*_R1.fastq.gz >> "$output"
```

# Upgrading R but keeping all packages from previous version
```{r}

# Save a list
installed_packages <- installed.packages()[, "Package"]
save(installed_packages, file = "installed_packages.RData")

# Install and update R. Do this from Rgui and not Rstudio!
install.packages("installr")
library(installr)
updateR()

# After updating R: Reload your packages. If R asks about missing libraries:

# Reinstall them easily 
load("installed_packages.RData")

# Reinstall only missing packages
new_packages <- setdiff(installed_packages, installed.packages()[,"Package"])
if(length(new_packages)) install.packages(new_packages) # Upgrade to newest package versions. This make take a while!
```

# Calculate the percent missing data for a fasta alignment

```bash
# Missing data using seqtk

total=$(seqkit seq -i my_alignment.fasta | grep -v '^>' | tr -d '\n' | wc -c)
missing=$(seqkit seq -i my_alignment.fasta | grep -v '^>' | tr -cd '\-\?N' | wc -c)

pct=$(echo "scale=2; 100 * $missing / $total" | bc)
echo "Missing data percentage: $pct%"
```

# Summarize .iqtree output with BIC and BIC weights analysis of a dataset under different models

```bash

# Code to extract information from several iqtree runs on the same dataset and compare them with BIC and BIC with

# For example, lets say we ran 5 analyses under different models in iqtree:
# iqtree2 -s accD_mafft_macse_NT.fasta -m GTR+R -bb 1000 -nt 50  --prefix gtr
# iqtree2 -s accD_mafft_macse_NT.fasta -m GTR+FO+H2 -bb 1000 -T 50  --prefix ghost_linked2
# iqtree2 -s accD_mafft_macse_NT.fasta -m GTR+FO+H3 -bb 1000 -T 50  --prefix ghost_linked3
# iqtree2 -s accD_mafft_macse_NT.fasta -m GTR+FO*H2 -bb 1000 -T 50  --prefix ghost_unlinked2
# iqtree2 -s accD_mafft_macse_NT.fasta -m GTR+FO*H3 -bb 1000 -T 50  --prefix ghost_unlinked3

# 1) Make the CSV
printf "file,n_taxa,n_sites,n_pars_inf,logL,n_params,AIC,AICc,BIC\n" > iqtree_model_summary.csv

# 2) Loop through all .iqtree output files in the current directory
for f in *.iqtree; do
  n_taxa=$(grep -m1 "Input data:" "$f" | sed -E 's/.*Input data: *([0-9]+) sequences.*/\1/')
  n_sites=$(grep -m1 "Input data:" "$f" | sed -E 's/.*with *([0-9]+) nucleotide sites.*/\1/')
  n_pars=$(grep -m1 "Number of parsimony informative sites:" "$f" | sed -E 's/.*: *([0-9]+).*/\1/')
  logL=$(grep -m1 "Log-likelihood of the tree:" "$f" | sed -E 's/.*: *([-0-9.eE+]+).*/\1/')
  n_par=$(grep -m1 "Number of free parameters" "$f" | sed -E 's/.*: *([0-9]+).*/\1/')
  AIC=$(grep -m1 "Akaike information criterion (AIC) score:" "$f" | sed -E 's/.*: *([-0-9.eE+]+).*/\1/')
  AICc=$(grep -m1 "Corrected Akaike information criterion (AICc) score:" "$f" | sed -E 's/.*: *([-0-9.eE+]+).*/\1/')
  BIC=$(grep -m1 "Bayesian information criterion (BIC) score:" "$f" | sed -E 's/.*: *([-0-9.eE+]+).*/\1/')

# 3) Print the captured text to the csv file
  printf "%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
    "$f" "$n_taxa" "$n_sites" "$n_pars" "$logL" "$n_par" "$AIC" "$AICc" "$BIC" \
    >> iqtree_model_summary.csv
done

# 4) Add wBIC column (BIC weights) -> new file, then replace
awk -F',' '
BEGIN { OFS="," }
NR==1 {
  print $0, "wBIC"
  next
}
{
  bic[NR] = $9 + 0
  line[NR] = $0
  if (NR==2 || bic[NR] < minbic) minbic = bic[NR]
  n = NR
}
END {
  # compute denominator
  denom = 0
  for (i=2; i<=n; i++) {
    w[i] = exp(-0.5 * (bic[i] - minbic))
    denom += w[i]
  }
  # print lines with normalized weights
  for (i=2; i<=n; i++) {
    printf "%s,%.10g\n", line[i], (w[i]/denom)
  }
}
' iqtree_model_summary.csv > iqtree_model_summary_wBIC.csv

mv iqtree_model_summary_wBIC.csv iqtree_model_summary.csv

```


