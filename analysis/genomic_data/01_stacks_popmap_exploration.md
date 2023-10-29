We began by exploring the data by defining four population maps:

- Islands where different specimens were collected (Floreana, Genovesa, Isabela, Pinta, San Cristóbal, Santa Cruz, Santiago, and Wolf)
- Taxonomic classification (echios, echios_mesaperma, galapageia, helleri, insularis, and megasperma)
- Islands where different specimens were collected, including outgroup taxa
- Taxonomic classification, including outgroup taxa

For all of these four, we ran the ['Paris optimization protocol'](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12775) including -M -N values of 1 to 10

```
for i in 1 2 3 4 5 6 7 8 9 10; do sbatch --export=M=$i stacks.sh ;done

# Where stacks.sh included:

# popmap= one of the four pop maps described above
# reads_dir= location of the reads
# out_dir= Output directory

denovo_map.pl --samples $reads_dir --popmap $popmap -o $out_dir -M $M -N $M --paired -T 8 -X "populations:--vcf -R 0.8" &> M_${M}_logfile.log
```

For each of the four analyses we did hockeystick plots (see Paris et al for details) and a PCA. The PCA analysis involved, cleaning the data in the following way:
```
# The denovo_map.pl command above output a vcf file, which was further cleaned by vcftools:
# --maf 0.05 - minimum allele frequency of 5 %
# --maxDP 150 and --max-meanDP 150 - max depth of 150 and max average depth of 150
# --minDP 150 and --min-meanDP 150 - min depth of 150 and min average depth of 150
# --max-missing 0.25 - Max missing data of 25 %
vcftools --vcf populations.snps.vcf --maf 0.05 --maxDP 150 --max-meanDP 150 --minDP 10 --min-meanDP 10 --max-missing 0.25 --recode --stdout
```
The PCA was done then by:

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

#we need the ind.species.pop.tsv file which is a "awk-modified" file from plink.fam
#the .Q files (they contain admixture proportions)
#the cv_error file :-)

# first read in and examine cross-validation error
info<- tbl_df(read.table("./pca_data.tsv", sep="\t", header=T))

# Loading the VCF
vcf_file<-read.vcfR("./opuntia_byIsland_M2.Maf0.05.maxDP150.minDP10.vcf")

# Converting it to genelight
genlight_file<-vcfR2genlight(vcf_file)

# Plotting missing data
glPlot(genlight_file, posi="topleft")


# PCA now
pca <- glPca(genlight_file)
4 #axis

scatter(pca1, posi="bottomright")
title("Opuntia\n axes 1-2")

# Pretty plot
pca_scores<-rownames_to_column(as.data.frame(pca$scores),"samples_fromFelipe") #passing the row names to the first column
final_df<-merge(info, pca_scores)
# calculate PVE
pve <- (pca$eig/sum(pca$eig))*100

# make PCA scatterplot
a <- ggplot(final_df, aes(PC1, PC2, pch = popmap_island, colour = popmap_species)) + geom_point(size=3)
b <- ggplot(final_df, aes(PC1, PC3, pch = popmap_island, colour = popmap_species)) + geom_point(size=3)
# add PVE labels
a <- a + xlab(paste0("PC1 (", signif(pve[1], 2), "%)")) + ylab(paste0("PC2 (", signif(pve[2], 2), "%)")) + theme_bw()
b <- b + xlab(paste0("PC1 (", signif(pve[2], 2), "%)")) + ylab(paste0("PC3 (", signif(pve[3], 2), "%)")) + theme_bw()
# extract legend
c <- g_legend(a + theme(legend.position = "right"))
# make plot
grid.arrange(arrangeGrob(a + theme(legend.position = "none"),
                         b + theme(legend.position = "none"), ncol = 2),
             c, ncol = 2, widths = c(6, 1))


a + geom_text(aes(label = samples_fromFelipe))
plotly(a)
```
