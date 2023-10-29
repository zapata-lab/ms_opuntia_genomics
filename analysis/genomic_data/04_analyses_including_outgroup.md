The following analysis were done on the dataset including the outgroup:

- a fasta-based phylogeny
- PCA
- Dsuite

Fasta-based phylogeny:
The perl script (processing.pl) can be found [here](https://github.com/jcerca/Papers/blob/main/Stygocapitella_PeerJ/info/Split_STACKS_fasta_file.pl), where as FASconCAT-G_v1.04.pl can be found [here](https://github.com/PatrickKueck/FASconCAT-G).

```
## First, I re-ran populations to obtain fasta files (as opposed to a vcf file)
populations  -P ../../01_stacksRun/with_outgroup/ -O . -M ../../../00_popmap/popmap_cleaned_withOutgroup.tsv -R 0.8 --fasta-samples

## Second, The code above spits out a single fasta file, which we parsed out in the following way.
in=populations.samples.fa
out=./loci_folder/
script=processing.pl
perl $script $in

# This script will separate the data into loci, where every file will have the loci and individual information for each entry (>entry):
Loci_0000.fa
Loci_0001.fa
...
Loci_9999.fa

## Third, we removed the files we do not need:
rm *_[01].fas

# Forth, we select loci to remove. We have 35 specimens in this dataset, therefore a total of 70 alleles. I decided to remove loci with 5 or more missing. So, a locus has to be in at least 30 specimens.
grep -c ">" * | tr ":" "\t" | awk '$2<59' | awk '{print $1}' > ../loci_to_remove.tsv
while read locus; do echo $locus; rm $locus; done < ../loci_to_remove.tsv

# This removed 3672 loci, keeping a total of 10,484


## Forth, I concatenated the data using FASconCAT-G_v1.04.pl
script=FASconCAT-G_v1.04.pl
input_folder=./loci_folder/
output_f=./concatenated_matrix

perl $script -s -p -n

## Fifth, because FastConCat separates alleles, we make a consensus per individual.

# This invoves concatenating a specimen list ...
grep ">" ../../03_Concatenatedmatrix/FcC_supermatrix.fas | sed "s/_All.*//" | sort | sed "s/>//" > specimen_list

# ... and getting a file for each allele, each individual.
while read gene; do grep -A 1 --no-group-separator $gene ../../03_Concatenatedmatrix/FcC_supermatrix.fas > $gene.fas; done < specimen_list

# Now, we generate consensus sequences using EMBOSS v6.6.0
for i in *fas; do echo $i; consambig -sequence $i -outseq ${i%.fas}.consensus.fasta; sed -i "s/>.*/>${i%.fas}/" ${i%.fas}.consensus.fasta; mv ${i%.fas}.consensus.fasta ../02_consensus/; done

for i in *fasta; do echo $i; awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $i > ${i%.fasta}_deinterleaved.fasta; done

03_finalMatrix
cat *_deinterleaved.fasta > ../03_finalMatrix/opuntia.fasta

## There was some issue with the deinterleaved files - something weird with the file lines. I corrected it in this way:
sed "s/>/\n>/" opuntia.fasta  | sed '/^$/d' > tmp
mv tmp opuntia.tsv

## Seventh, and finally, we are ready to make a tree.
# This involves formating a partitions file, and gettingthe data:
cat FcC_info.xls | grep "Locus" | awk '{print $1" = "$2"-"$3";"}' | sed 's/^/\tcharset /' | sed '1 i\begin sets;' | sed '1 i\#nexus' | sed '$ a\end;' > ../05_iqtree/opuntias.nex
ln  -s ../04_consensusMatrix/03_finalMatrix/opuntia.fasta  .

iqtree2 -s opuntia.fasta -p opuntias.nex -T 10 -B 1000
```

PCA
```

```

Dsuite
```
```
