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

PCA - preparation
```
# For the PCA, we re-ran populations, specifying we want a random snp per rad-tag (which is indirectly LD prunning for de novo datasets).
populations  -P ../../01_stacksRun/with_outgroup/ -O . -M ../../../00_popmap/popmap_cleaned_withOutgroup.tsv -R 0.8 --vcf -t 8 --structure --write-random-snp

# We then cleaned the vcf for depth and missing data 
ml VCFtools/0.1.16-intel-2018b-Perl-5.28.0
vcftools --vcf ../01_stacksPopulations/populations.snps.vcf --maxDP 200 --max-meanDP 200 --minDP 10 --min-meanDP 10 --max-missing 0.25 --recode --stdout > opuntia_with_outgroup_WriteRandomSNP_R80_MaxDP200_MinDP10_MaxMissing0_25.vcf
```

PCA - plotting
```
rm(list = ls())

library(vcfR)       #Data loading
library(adegenet)   #Data analysis
library(gridExtra)  #PCA
library(ape)        #Phylogenetic analyses
library(lemon)      #for the g_legend
library(tidyverse)
library(plotly)
setwd("~/Desktop/tmp/opuntia/")

# first read in and examine cross-validation error
info<- tbl_df(read.table("./pca_data.tsv", sep="\t", header=T))

# Loading the VCF
vcf_file<-read.vcfR("./opuntia_with_outgroup_WriteRandomSNP_R80_MaxDP200_MinDP10_MaxMissing0_25.vcf")

# Converting it to genelight
genlight_file<-vcfR2genlight(vcf_file)

# Plotting missing data
glPlot(genlight_file, posi="topleft")


# PCA now
pca <- glPca(genlight_file)
4 #axis

# Pretty plot
pca_scores<-rownames_to_column(as.data.frame(pca$scores),"samples_fromFelipe") #passing the row names to the first column
final_df<-merge(info, pca_scores)
# calculate PVE
pve <- (pca$eig/sum(pca$eig))*100

# make PCA scatterplot
a <- ggplot(final_df, aes(PC1, PC2, colour = popmap_island)) + geom_point(size=3)  + theme_bw() +
  xlab(paste0("PC1 (", signif(pve[1], 2), "%)")) + ylab(paste0("PC2 (", signif(pve[2], 2), "%)"))
b <- ggplot(final_df, aes(PC1, PC3, colour = popmap_island)) + geom_point(size=3)  + theme_bw() +
  xlab(paste0("PC1 (", signif(pve[2], 2), "%)")) + ylab(paste0("PC3 (", signif(pve[3], 2), "%)"))

# extract legend
c <- g_legend(a + theme(legend.position = "right"))
# make plot
grid.arrange(arrangeGrob(a + theme(legend.position = "none"),
                         b + theme(legend.position = "none"), ncol = 2),
             c, ncol = 2, widths = c(6, 1))


a + geom_text(aes(label = samples_fromFelipe))
plotly(a)
```

Dsuite
```
## For D-suite, we started with a vcf
ml Stacks/2.60-foss-2020a
populations  -P ../../01_stacksRun/with_outgroup/ -O . -M ../../../00_popmap/popmap_cleaned_withOutgroup.tsv -R 0.8 --vcf -t 8 --batch-size 20000

## Which was cleaned with vcf tools for depth and missing data
vcftools --vcf ../01_populationsRun/populations.snps.vcf --maxDP 200 --max-meanDP 200 --minDP 10 --min-meanDP 10 --max-missing 0.25 --recode --stdout > opuntia_with_outgroup_R80_MaxDP200_MinDP10_MaxMissing0_25.vcf

## For D-suite we followed https://github.com/millanek/tutorials/tree/master/analysis_of_introgression_with_snp_data
## Which involved specifying the tree, following the iq-tree run above.
# ((((((genovesa,pinta),santiago),(isabela,santa_cruz)),wolf),(floreana,san_cristobal)),Outgroup);
# Important, it needs to be Outgroup with a capital O
# Important, it does not handle  spaces.
# Important, it needs the ";" in the end, for the dtools.py code.

# and running Dsuite:
Dsuite Dtrios ../02_vcfCleaning/opuntia_with_outgroup_WriteRandomSNP_R80_MaxDP200_MinDP10_MaxMissing0_25.vcf ./popmap.tsv -o opuntia_trios -t tree.nwk

# We then did a rough plotting of the data to determine the max values for the final plots
### Exploring on R 
setwd("~/Desktop/tmp/opuntia/")

D_BBAA <- read.table("./opuntia_trios_BBAA.txt",as.is=T,header=T)
plot(D_BBAA$Dstatistic, ylab="D",xlab="trio number")

D_BBAA[which(D_BBAA$Dstatistic > 0.1),]

plot(D_BBAA$p.value, ylab="p value",xlab="trio number",ylim=c(0,0.05))

plot(p.adjust(D_BBAA$p.value,method="BH"), ylab="p value",xlab="trio number",ylim=c(0,0.05))

plot(D_BBAA$f4.ratio, ylab="f4-ratio",xlab="trio number", ylim=c(0,1))

# conclusion:
# 0.185116 is the max D value
# 0.189288 is the max f4 value

# Now, we get some nice svg plots, big shout-out to Micha Matschiner. I've paid him some beers, and you should too!

### Now getting some nice svg plots.
cut -f 2 ../03*/popmap.tsv | sort | uniq > plot_order.txt
# I changed the order.

# plot_d.rb came from
# https://github.com/millanek/tutorials/blob/master/analysis_of_introgression_with_snp_data/src/plot_d.rb
# plot_f4ratio.rb came from
# https://github.com/millanek/tutorials/blob/master/analysis_of_introgression_with_snp_data/src/plot_f4ratio.rb

#Plotting f4 and Dstat
# 0.185116 is the max D value - I saw this in the R code above.
module load Ruby/3.0.1-GCCcore-10.3.0
ruby plot_d.rb ../03_Dsuite_ABBABABA/opuntia_trios_tree.txt plot_order.txt 0.19 ABBA_BABA_BBAA.svg

# 0.189288 is the max f4 value - as determined above
ruby plot_f4ratio.rb ../03_Dsuite_ABBABABA/opuntia_trios_tree.txt plot_order.txt 0.19 f4.svg

# Now doing the F branch
ml Dsuite/20210309-GCC-9.3.0
# The install does not have the utils folder, so I downloaded dtools from :
# https://github.com/millanek/Dsuite/blob/master/utils/dtools.py
Dsuite Fbranch ../03_Dsuite_ABBABABA/tree.nwk ../03_Dsuite_ABBABABA/opuntia_trios_tree.txt > Fbranch.txt
./dtools.py Fbranch.txt ../03_Dsuite_ABBABABA/tree.nwk
```
