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
```

splitsTree
```
```

phylTree
```
```

SNAPP
```
```
