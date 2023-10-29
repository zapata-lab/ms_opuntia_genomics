[Cerca & Maurstad et al](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13562) have shown that removing 'bad apples' (i.e. individuals with high rates of missing data at the population-level) significantly improves the final dataset in terms of number of loci.

We followed the bad apples protocol by running stacks for each island separately:

```
# For example, for floreana this involved:
denovo_map.pl --samples $reads_dir --popmap floreana.tsv -o $out_dir -M 1 -N 1 --paired -T 8 -X "populations:--vcf -R 0.8" &> floreana_logfile.log
```

This led us to conclude:
floreana - no need to remove individuals
genovesa - excluding 531_O_helleri_UCBplate1 as it has 28% missing data when all others had <6% (with the exception of two specimens with 16% and 13%).
outgroup - no need to remove individuals
pinta - no need to remove individuals
santa_cruz - excluding 397_O_echios_UCBplate1 as it has 18% missing data when all others had <6% (with the exception of one specimens had 11%).
santiago - no need to remove individuals
wolf - excluding Op02_Wolf_O_sp_UCBplate1 as it had 24% missing data whereas the remaining had 4-6% missing data.
