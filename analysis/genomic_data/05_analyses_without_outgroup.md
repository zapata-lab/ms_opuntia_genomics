These analyses include:

- PCA
- STRUCTURE
- Private alleles
- fineRAD
- FST
- splitsTree
- phylTree
- SNAPP

PCA & STRUCTURE analyses
```
# Both these analyses involved obtaining a vcf
populations  -P ../../01_stacksRun/without_outgroup/ -O . -M ../../../00_popmap/popmap_cleaned_withoutOutgroup.tsv -R 0.8 --vcf -t 8 --structure --write-random-snp

# ... and cleaning the vcf
ml VCFtools/0.1.16-intel-2018b-Perl-5.28.0
vcftools --vcf ../01_populationsRun/populations.snps.vcf --maxDP 200 --max-meanDP 200 --minDP 10 --min-meanDP 10 --max-missing 0.25 --recode --stdout > opuntia_without_outgroup_WriteRandomSNP_R80_MaxDP200_MinDP10_MaxMissing0_25.vcf

# For the PCA, I gave the R script in 04_analyses_including_outgroup.

# For the structure the code goes as follows

# We rename the vcf adding x's in front of the chr name
sed "/#/! s/^/x/" ../../02_vcfcleaning/opuntia_without_outgroup_WriteRandomSNP_R80_MaxDP200_MinDP10_MaxMissing0_25.vcf > renamed.vcf
plink --vcf  renamed.vcf --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--out structure --recode structure

# Then, I changed manually the second column to define populations.
# I also removed the line with -1's. I did it using control+K

# Notice, it didn't like population names, so I switched it for integers
sed -i "s/santa_cruz/1/; s/isabela/2/; s/genovesa/3/; s/pinta/4/; s/santiago/5/; s/wolf/6/; s/san_cristobal/7/; s/floreana/8/" structure.recode.strct_in

# For the structure analyses, we added number of specimens etc.
# For the number of loci, this works:
tail -1  ../01_converting/*in | tr " " "\n" | wc -l
23568 # remove 2 lines due to the first two being ind and pop., so 23566. Which, divided by two = 11783 loci

# For the structure analyses, I did:
for i in 1 2 3 4 5 6 7 8 9 10; do
  for j in 1 2 3 4 5; do
    echo sbatch --export=K=$i,run=$j structure.sh;
    sbatch --export=K=$i,run=$j structure.sh;
  done;
done

# structure.sh had:
structure -m ../../mainparams -e ../../extraparams -K $K -L 11783 -N 33 -i ../../../01_converting/structure.recode.strct_in -o ./run_${run} -D $RANDOM

```


Private alleles
```
# We obtained a vcf and cleaned it as follows:
populations -P ../../01_stacksRun/without_outgroup/ -O . -M ../../../00_popmap/popmap_cleaned_withoutOutgroup.tsv -R 0.8 --vcf -t 8
vcftools --vcf ../01_populations_run/populations.snps.vcf --maxDP 200 --max-meanDP 200 --minDP 10 --min-meanDP 10 --max-missing 0.25 --recode --stdout > opuntia_without_outgroup_WriteRandomSNP_R80_MaxDP200_MinDP10_MaxMissing0_25.vcf

library('vcfR')
library('tidyverse')
library(reshape2)

vcf <- read.vcfR("opuntia_without_outgroup_WriteRandomSNP_R80_MaxDP200_MinDP10_MaxMissing0_25.vcf")
vcf

## Gave it a list of populations, ordered as the individuals are in the vcf.
pop<-as.factor(c("santa_cruz", "santa_cruz", "santa_cruz",
      "isabela", "santa_cruz", "pinta",
      "genovesa","pinta","genovesa",
      "santa_cruz", "santa_cruz", "santa_cruz",
      "santiago", "santiago", "santa_cruz",
      "genovesa", "genovesa", "genovesa",
      "genovesa", "genovesa", "genovesa",
      "genovesa", "genovesa", "genovesa",
      "wolf", "wolf", "wolf",
      "wolf","wolf","wolf",
      "san_cristobal", "floreana", "floreana"))


myDiff <- genetic_diff(vcf, pops = pop, method = 'nei')
### This calculates heterozygotes and allelles per population

#manipulating.
df2<-myDiff [,1:10]

### Alleles private to floreana
# [1] "CHROM_POS"        "Hs_floreana"      "Hs_genovesa"      "Hs_isabela"       "Hs_pinta"
# [6] "Hs_san_cristobal" "Hs_santa_cruz"    "Hs_santiago"      "Hs_wolf"
floreana_private <- subset(df2, Hs_floreana > 0 & Hs_genovesa==0 & Hs_isabela==0 & Hs_pinta==0 & Hs_san_cristobal==0 & Hs_santa_cruz==0 & Hs_santiago==0 & Hs_wolf==0)
genovesa_private <- subset(df2, Hs_floreana==0 & Hs_genovesa > 0 & Hs_isabela==0 & Hs_pinta==0 & Hs_san_cristobal==0 & Hs_santa_cruz==0 & Hs_santiago==0 & Hs_wolf==0)
isabela_private <- subset(df2, Hs_floreana==0 & Hs_genovesa==0 & Hs_isabela > 0 & Hs_pinta==0 & Hs_san_cristobal==0 & Hs_santa_cruz==0 & Hs_santiago==0 & Hs_wolf==0)
pinta_private <- subset(df2, Hs_floreana==0 & Hs_genovesa==0 & Hs_isabela==0 & Hs_pinta > 0 & Hs_san_cristobal==0 & Hs_santa_cruz==0 & Hs_santiago==0 & Hs_wolf==0)
san_cristobal_private <- subset(df2, Hs_floreana==0 & Hs_genovesa==0 & Hs_isabela==0 & Hs_pinta==0 & Hs_san_cristobal > 0 & Hs_santa_cruz==0 & Hs_santiago==0 & Hs_wolf==0)
santa_cruz_private <- subset(df2, Hs_floreana==0 & Hs_genovesa==0 & Hs_isabela==0 & Hs_pinta==0 & Hs_san_cristobal==0 & Hs_santa_cruz > 0 & Hs_santiago==0 & Hs_wolf==0)
santiago_private <- subset(df2, Hs_floreana==0 & Hs_genovesa==0 & Hs_isabela==0 & Hs_pinta==0 & Hs_san_cristobal==0 & Hs_santa_cruz==0 & Hs_santiago > 0 & Hs_wolf==0)
wolf_private <- subset(df2, Hs_floreana==0 & Hs_genovesa==0 & Hs_isabela==0 & Hs_pinta==0 & Hs_san_cristobal==0 & Hs_santa_cruz==0 & Hs_santiago==0 & Hs_wolf > 0)


nrow(floreana_private)
nrow(genovesa_private)
nrow(isabela_private)
nrow(pinta_private)
nrow(san_cristobal_private)
nrow(santa_cruz_private)
nrow(santiago_private)
nrow(wolf_private)
```

fineRAD
```
# Getting the radpainter input
populations -P ../../01_stacksRun/without_outgroup/ -O . -M ../../../00_popmap/popmap_cleaned_withoutOutgroup.tsv -R 0.8 -t 8 --radpainter

# Now following (and taking plotting scripts from), https://github.com/millanek/fineRADstructure/
ml fineRADstructure/0.3.2r109-intel-2018b
ln -s ../01_populations/populations.haps.radpainter .

## Because it does not like specimens to begin with numbers, we will use sed.
sed -i "s/376_O_echios_UCBplate1/x376_O_echios_UCBplate1/; s/387_O_echios_UCBplate1/x387_O_echios_UCBplate1/; s/409_O_echios_UCBplate1/x409_O_echios_UCBplate1/; s/483_O_insularis_UCBplate1/x483_O_insularis_UCBplate1/; s/346_O_echios_UCBplate1/x346_O_echios_UCBplate1/; s/607_O_galapageia_UCBplate1/x607_O_galapageia_UCBplate1/; s/537_O_helleri_UCBplate1/x537_O_helleri_UCBplate1/; s/594_O_galapageia_UCBplate1/x594_O_galapageia_UCBplate1/; s/538_O_helleri_UCBplate1/x538_O_helleri_UCBplate1/; s/347_O_echios_UCBplate1/x347_O_echios_UCBplate1/; s/323_O_echios_UCBplate1/x323_O_echios_UCBplate1/; s/348_O_echios_UCBplate1/x348_O_echios_UCBplate1/; s/652_O_galapageia_UCBplate1/x652_O_galapageia_UCBplate1/; s/647_O_galapageia_UCBplate1/x647_O_galapageia_UCBplate1/; s/399_O_echios_UCBplate1/x399_O_echios_UCBplate1/;"  populations.haps.radpainter

RADpainter paint populations.haps.radpainter
# created populations.haps_chunks.out  populations.haps_missingnessMatrix.out  populations.haps_missingness.out
finestructure -x 100000 -y 100000 -z 1000 ./populations.haps_chunks.out ./populations.haps_chunks.mcmc.xml
finestructure -m T -x 10000
finestructure -m T -x 10000 populations.haps_chunks.out populations.haps_chunks.mcmc.xml populations.haps_chunks.mcmcTree.xml

# Rplots taken from the link above.
```

FST
```
## Getting a vcf, and cleaning it.
populations -P ../../01_stacksRun/without_outgroup/ -O . -M ../../../00_popmap/popmap_cleaned_withoutOutgroup.tsv -R 0.8 --vcf -t 8 --fasta-loci --batch-size 20000
vcftools --vcf ../01_populations/populations.snps.vcf --maxDP 200 --max-meanDP 200 --minDP 10 --min-meanDP 10 --max-missing 0.25 --recode --stdout > opuntia_without_outgroup_WriteRandomSNP_R80_MaxDP200_MinDP10_MaxMissing0_25.vcf

# FST on populations with >=5 individuals:
# santa_cruz.tsv
# wolf.tsv
# genovesa.tsv

# FST calculations:
vcftools --vcf ../02_vcfClean/opuntia_without_outgroup_WriteRandomSNP_R80_MaxDP200_MinDP10_MaxMissing0_25.vcf --fst-window-size 1000 --fst-window-step 1000 --weir-fst-pop santa_cruz.tsv --weir-fst-pop genovesa.tsv --out santa_cruz_genovesa
Weir and Cockerham mean Fst estimate: 0.10014
Weir and Cockerham weighted Fst estimate: 0.18869

vcftools --vcf ../02_vcfClean/opuntia_without_outgroup_WriteRandomSNP_R80_MaxDP200_MinDP10_MaxMissing0_25.vcf --fst-window-size 1000 --fst-window-step 1000 --weir-fst-pop santa_cruz.tsv --weir-fst-pop wolf.tsv --out santa_cruz_wolf
Weir and Cockerham mean Fst estimate: 0.13307
Weir and Cockerham weighted Fst estimate: 0.25564

vcftools --vcf ../02_vcfClean/opuntia_without_outgroup_WriteRandomSNP_R80_MaxDP200_MinDP10_MaxMissing0_25.vcf --fst-window-size 1000 --fst-window-step 1000 --weir-fst-pop genovesa.tsv --weir-fst-pop wolf.tsv --out genovesa_wolf
Weir and Cockerham mean Fst estimate: 0.15083
Weir and Cockerham weighted Fst estimate: 0.26648

# Selecting the top 5%. 
wc -l *fst
  8931 genovesa_wolf.windowed.weir.fst (5% is 446)
  9802 santa_cruz_genovesa.windowed.weir.fst (5% is 490)
  9321 santa_cruz_wolf.windowed.weir.fst (5% is 466)

# Adding an extra number to the tail so we remove the header
cat genovesa_wolf.windowed.weir.fst | grep -v "e" | grep -v "-" | sort -k5 | tail -n 447 | awk '{print $1}' | head -n 446 | sed "s/^/CLocus_/" > ../04_chrBLAST2go/genovesa_wolf_CHR_toExtract.tsv
cat santa_cruz_genovesa.windowed.weir.fst | grep -v "e" | grep -v "-" | sort -k5 | tail -n 491 | awk '{print $1}' | head -n 490 | sed "s/^/CLocus_/" > ../04_chrBLAST2go/santa_cruz_genovesa_CHR_toExtract.tsv
cat santa_cruz_wolf.windowed.weir.fst | grep -v "e" | grep -v "-" | sort -k5 | tail -n 467 | awk '{print $1}' | head -n 466 | sed "s/^/CLocus_/" > ../04_chrBLAST2go/santa_cruz_wolf_CHR_toExtract.tsv

# Code explained. I removed the e from entries such as 1e-10, and negative fst values. I sorted after weighted Fst, and get the top 447 fst values, printed only the chr name, removed the final line with the header, and added locus ID.


grep -A 1 --no-group-separator -F -w -f genovesa_wolf_CHR_toExtract.tsv ../01_populations/populations.loci.fa > genovesa_wolf.fa
grep -A 1 --no-group-separator -F -w -f santa_cruz_genovesa_CHR_toExtract.tsv ../01_populations/populations.loci.fa > santa_cruz_genovesa.fa
grep -A 1 --no-group-separator -F -w -f santa_cruz_wolf_CHR_toExtract.tsv ../01_populations/populations.loci.fa > santa_cruz_wolf.fa

# Downloaded the database: https://www.arabidopsis.org/download_files/Proteins/TAIR10_protein_lists/TAIR10_pep_20101214

makeblastdb -in TAIR10_pep_20101214.faa -parse_seqids -dbtype prot
blastx -out santa_cruz_genovesa.blast -db TAIR10_pep_20101214.faa -query ../04_chrBLAST2go/santa_cruz_genovesa.fa -outfmt "6 qseqid sseqid evalue sstart send sseq"
blastx -out santa_cruz_wolf.blast -db TAIR10_pep_20101214.faa -query ../04_chrBLAST2go/santa_cruz_wolf.fa -outfmt "6 qseqid sseqid evalue sstart send sseq"
blastx -out genovesa_wolf.blast -db TAIR10_pep_20101214.faa -query ../04_chrBLAST2go/genovesa_wolf.fa -outfmt "6 qseqid sseqid evalue sstart send sseq"


### conditions, e-value below 0.001 and length of aa sequence >20.
cat genovesa_wolf.blast | awk '$3<0.001 && length($6) >20 {print $0}'  > genovesa_wolf_tentativeArabidopsis.tsv
cat santa_cruz_genovesa.blast | awk '$3<0.001 && length($6) >20 {print $0}'  > santa_cruz_genovesa_tentativeArabidopsis.tsv
cat santa_cruz_wolf.blast | awk '$3<0.001 && length($6) >20 {print $0}'  > santa_cruz_wolf_tentativeArabidopsis.tsv
```

splitsTree
```
#Getting a vcf and cleaning it
populations -P ../../01_stacksRun/without_outgroup/ -O . -M ../../../00_popmap/popmap_cleaned_withoutOutgroup.tsv -R 0.8 -t 8 --vcf
vcftools --vcf ../01_populations/populations.snps.vcf --maxDP 200 --max-meanDP 200 --minDP 10 --min-meanDP 10 --max-missing 0.25 --recode --stdout > opuntia_without_outgroup_R80_MaxDP200_MinDP10_MaxMissing0_25.vcf

# Converting it using vcf2phyl.py
python2 vcf2phyl.py -i ../02_vcfClean/opuntia_without_outgroup_R80_MaxDP200_MinDP10_MaxMissing0_25.vcf --nexus

# Ran it on the GUI for splitstree.
```

phylogenetic tree
```
# The code is similar to the code used for the IQ-tree analyses with outgroup (see script 04_analyses_including_outgroup.md) :-)

```

SNAPP
```
# 00 - parsing out and reducing the vcf to 4,000 SNPs
ln -s ../../02_PCA_admixture/02_vcfcleaning/opuntia_without_outgroup_WriteRandomSNP_R80_MaxDP200_MinDP10_MaxMissing0_25.vcf  .

# First, we get 4000 snps.
for i in 1 2 3 4 5; do grep -v "^#" opuntia_without_outgroup_WriteRandomSNP_R80_MaxDP200_MinDP10_MaxMissing0_25.vcf | cut -f1-3 | shuf | head -n 4000 | cut -f 3 > random.$i.snps; done

for i in 1 2 3 4 5; do vcftools --vcf opuntia_without_outgroup_WriteRandomSNP_R80_MaxDP200_MinDP10_MaxMissing0_25.vcf --snps random.$i.snps --recode --out trimmed.$i.vcf; done

# Second, we use the vcf2nex.pl downloaded from the SNAPP website.
perl ../vcf2nex.pl trimmed.vcf > trimmed.nex

# Third, we need to ' missing=. ' manually to the nexus file so it knows how the missingdata should be processed.
nano trimmed.nex

# Forth, we format populations manually on beauti GUI (snapp window)
# Fifth, we uncheck "Include non-polymorphic sites"
# Sixth calculate mutation rates.

beast -seed $RANDOM -threads 10 -beagle trimmed.xls.xml
```
