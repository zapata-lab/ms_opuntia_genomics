We ran Stacks with and without outgroup, specifying that a locus has to be in at least 80% of the dataset, and with either -M/-N of 1 or 2 (as determined in the optimization).

```
# With outgroup
denovo_map.pl --samples $reads_directory --popmap $popmap -o $out_dir -M 2 -N 2 --paired -T 8 -X "populations:--vcf -R 0.8" &> logfile.log

# Without outgroup
denovo_map.pl --samples $reads_dir --popmap $popmap -o $out_dir -M 1 -N 1 --paired -T 8 -X "populations:--vcf -R 0.8" &> logfile.log
```
