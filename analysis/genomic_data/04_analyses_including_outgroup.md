The following analysis were done on the dataset including the outgroup:

- a fasta-based phylogeny
- PCA
- Dsuite

Fasta-based phylogeny:
```
## First, I re-ran populations to obtain fasta files (as opposed to a vcf file)
populations  -P ../../01_stacksRun/with_outgroup/ -O . -M ../../../00_popmap/popmap_cleaned_withOutgroup.tsv -R 0.8 --fasta-samples

## The code above spits out a single fasta file, which we parsed out in the following way.
in=populations.samples.fa
out=./out_folder/
script=processing.pl

# The script can be found [here](https://github.com/jcerca/Papers/blob/main/Stygocapitella_PeerJ/info/Split_STACKS_fasta_file.pl).

perl $script $in


02 - loci.sh
# First, we run loci.sh

# Then, we remove the files we do not need
rm *_[01].fas

# Then, we select loci to remove. We have 35 specimens, therefore 70 alleles. I'll remove loci with 5 or more missing. So, a locus has to be in at least 30 specimens.
grep -c ">" * | tr ":" "\t" | awk '$2<59' | awk '{print $1}' > ../loci_to_remove.tsv
while read locus; do echo $locus; rm $locus; done < ../loci_to_remove.tsv
# This removed 3672 loci


### 03 - I ran concatenate.sh
sbatch concatenate.sh

### 04 - Making a consensus matrix, see code inside

### 05 - iqtree
# I formatted a partitions files, and got the data
cat FcC_info.xls | grep "Locus" | awk '{print $1" = "$2"-"$3";"}' | sed 's/^/\tcharset /' | sed '1 i\begin sets;' | sed '1 i\#nexus' | sed '$ a\end;' > ../05_iqtree/opuntias.nex
ln  -s ../04_consensusMatrix/03_finalMatrix/opuntia.fasta  .
```

PCA
```
```

Dsuite
```
```
