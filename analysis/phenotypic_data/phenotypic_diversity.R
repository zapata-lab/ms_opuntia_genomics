########################################################################
# Code used to estimate hypervolumes representing all taxa in Opuntia. #
# and overlap of taxa in phenotyoic space                              #
########################################################################


# This script is based on the code presented in: 
# Jacobs S. et al. 2021 Sci Rep. https://doi.org/10.1038/s41598-021-03419-0

#### Preliminaries

# Load libraries
library(tidyverse)
library(geometry)
library(reshape2)

# define function to calculate the intersection of hypercubes 
# for multiple species and store response variables: 
# volume of sp1, 
# volume of sp2, 
# volume of the intersection of sp1 and sp2

intersection_by_pair = function(i,j){
  inthull <- data.frame(NA)
  chullsp1 <- data.frame(NA)
  chullsp2 <- data.frame(NA)
  result = intersectn(i,j)
  chullsp1[1,1] = result$ch1$vol #store the volume of the first species n-cube
  chullsp2[1,1] = result$ch2$vol #store the volume of the second species n-cube
  inthull[1,1] = result$ch$vol #store the volume of the intersection hull
}


#### Load data

# read in data
indata_veg_raw <- read_csv("../data/VegetetativeTraits.csv")
indata_rep_raw <- read_csv("..data/ReproductiveTraits.csv")

# log transform data: traits at very different scales
indata_veg = indata_veg_raw  %>% mutate_at(c("min", "max"), log)
indata_rep = indata_rep_raw %>% mutate_at(c("min", "max"), log)

indata_all <- bind_rows(indata_veg, indata_rep)


#### ANALYSIS 


#########################
##### VEGETATIVE DATA
#########################
# Parse data, find VERTICES of hypercubes, and create 
# a list of vertices for each species based on those traits

all_sp_list_veg = lapply(split(indata_veg, indata_veg$Species), 
                         function(x) lapply(split(x[3:4], x$Morpho_trait), 
                                            unlist, use.names = FALSE))
list_of_all_vertices_veg = lapply(all_sp_list_veg, 
                                  function(x) as.matrix(expand.grid(x)))

# calculate the hypercubes for each species using the convhulln function in library geometry 
# - depending on the number of vertices, this can take a while to calculate (on the order 
# of tens of minutes)
chull_by_group_veg <- lapply(list_of_all_vertices_veg, 
                             function(x) convhulln(x))

# calculate the volume of the INTERSECTION of all pairwise combinations of 
# hypercubes - similarly, can take a while to calculate, depending on the number of
# vertices and the number of n-cubes considered

# uses the function defined at the beginning 
intersectionHulls_veg <- lapply(list_of_all_vertices_veg, 
                                function(x) lapply(list_of_all_vertices_veg, 
                                                   function(y) intersection_by_pair(x, y)))


#############################################################
# PROPORTION OF OVERLAP: Do hypercubes of different species overlap? 
# By how much?
############################################################# 
# convert the information about the intersection of all combinations of hypercubes 
# into a dataframe and store the diagonal element (the species specific
# n-cube volumes)
intersectnDF_veg <- data.frame(t(sapply(intersectionHulls_veg,c)))
DF_diag_veg <- unlist(diag((as.matrix(intersectnDF_veg)), names=FALSE))


# calculate the proportion of overlap per column, based on the diagonal element residing in that column 
# perform the column-wise calculation on the first column, using the first element in the diagonal list
overlapProp_veg <- data.frame()

for(i in 1:length(DF_diag_veg)){
  overlapProp_veg <- rbind(overlapProp_veg, lapply(intersectnDF_veg[,i], 
                                                   function(x) x/DF_diag_veg[[i]]))
}
rownames(overlapProp_veg) <- rownames(intersectnDF_veg)


# make the results tidy and long; as.matrix part is essential to get var1 and var2, as well as the value
melted_overlapProp_veg <- melt(as.matrix(overlapProp_veg))

melted_overlapProp_veg <- melted_overlapProp_veg %>% mutate(Var1 = str_replace_all(Var1, "_", ". "),
                                                            Var2 = str_replace_all(Var2, "_", ". "))

#### VISUALIZATION

# plot the overlap proportions of the n-cubes
propOverlap_veg <- ggplot(melted_overlapProp_veg, aes(x=Var1, y=Var2)) + 
  geom_tile(aes(fill=value)) +
  #geom_text(aes(label=round(value,digits=1)),color='darkgrey') +
  scale_fill_gradient("Proportion\noverlap", low='white', high="black") +
theme(axis.title = element_blank(),
      axis.text = element_text(size=8),
      axis.text.x = element_text(angle=90, hjust=0.1)) +
  theme(legend.text = element_text(size=6),
        legend.margin = margin(2,2,6,2),
        legend.title = element_text(size=8),
        legend.background = element_rect(color="#cccccc")) +
  guides(fill = guide_colorbar(barwidth = 0.5, barheight = 2))





#########################
##### REPRODUCTIVE DATA
#########################
# Parse data, find VERTICES of hypercubes, and create 
# a list of vertices for each species based on those traits

all_sp_list_rep = lapply(split(indata_rep, indata_rep$Species), 
                         function(x) lapply(split(x[3:4], x$Morpho_trait), 
                                            unlist, use.names = FALSE))
list_of_all_vertices_rep = lapply(all_sp_list_rep, 
                                  function(x) as.matrix(expand.grid(x)))

# calculate the hypercubes for each species using the convhulln function in library geometry 
# - depending on the number of vertices, this can take a while to calculate (on the order 
# of tens of minutes)
chull_by_group_rep <- lapply(list_of_all_vertices_rep, 
                             function(x) convhulln(x))

# calculate the volume of the INTERSECTION of all pairwise combinations of 
# hypercubes - similarly, can take a while to calculate, depending on the number of
# vertices and the number of n-cubes considered

# uses the function defined at the beginning 
intersectionHulls_rep <- lapply(list_of_all_vertices_rep, 
                                function(x) lapply(list_of_all_vertices_rep, 
                                                   function(y) intersection_by_pair(x, y)))

#### ANALYSIS 

#############################################################
# PROPORTION OF OVERLAP: Do hypercubes of different species overlap? 
# By how much?
############################################################# 
# convert the information about the intersection of all combinations of hypercubes 
# into a dataframe and store the diagonal element (the species specific
# n-cube volumes)
intersectnDF_rep <- data.frame(t(sapply(intersectionHulls_rep,c)))
DF_diag_rep <- unlist(diag((as.matrix(intersectnDF_rep)), names=FALSE))


# calculate the proportion of overlap per column, based on the diagonal element residing in that column 
# perform the column-wise calculation on the first column, using the first element in the diagonal list
overlapProp_rep <- data.frame()

for(i in 1:length(DF_diag_rep)){
  overlapProp_rep <- rbind(overlapProp_rep, lapply(intersectnDF_rep[,i], 
                                                   function(x) x/DF_diag_rep[[i]]))
}
rownames(overlapProp_rep) <- rownames(intersectnDF_rep)


# make the results tidy and long; as.matrix part is essential to get var1 and var2, as well as the value
melted_overlapProp_rep <- melt(as.matrix(overlapProp_rep))

melted_overlapProp_rep <- melted_overlapProp_rep %>% mutate(Var1 = str_replace_all(Var1, "_", ". "),
                                                            Var2 = str_replace_all(Var2, "_", ". "))

#### VISUALIZATION

# plot the overlap proportions of the n-cubes
propOverlap_rep <- ggplot(melted_overlapProp_rep, aes(x=Var1, y=Var2)) + 
  geom_tile(aes(fill=value)) +
  #geom_text(aes(label=round(value,digits=1)),color='darkgrey') +
  scale_fill_gradient("Proportion\noverlap", low='white', high="black") +
  theme(axis.title = element_blank(),
        axis.text = element_text(size=8),
        axis.text.x = element_text(angle=90, hjust=0.1)) +
  theme(legend.text = element_text(size=6),
        legend.margin = margin(2,2,6,2),
        legend.title = element_text(size=8),
        legend.background = element_rect(color="#cccccc")) +
  guides(fill = guide_colorbar(barwidth = 0.5, barheight = 2))




#########################
##### ALL DATA
#########################
# Parse data, find VERTICES of hypercubes, and create 
# a list of vertices for each species based on those traits

all_sp_list_all = lapply(split(indata_all, indata_all$Species), 
                         function(x) lapply(split(x[3:4], x$Morpho_trait), 
                                            unlist, use.names = FALSE))
list_of_all_vertices_all = lapply(all_sp_list_all, 
                                  function(x) as.matrix(expand.grid(x)))

# calculate the hypercubes for each species using the convhulln function in library geometry 
# - depending on the number of vertices, this can take a while to calculate (on the order 
# of tens of minutes)
chull_by_group_all <- lapply(list_of_all_vertices_all, 
                             function(x) convhulln(x))

# calculate the volume of the INTERSECTION of all pairwise combinations of 
# hypercubes - similarly, can take a while to calculate, depending on the number of
# vertices and the number of n-cubes considered

# uses the function defined at the beginning 
intersectionHulls_all <- lapply(list_of_all_vertices_all, 
                                function(x) lapply(list_of_all_vertices_all, 
                                                   function(y) intersection_by_pair(x, y)))

#### ANALYSIS 

#############################################################
# PROPORTION OF OVERLAP: Do hypercubes of different species overlap? 
# By how much?
############################################################# 
# convert the information about the intersection of all combinations of hypercubes 
# into a dataframe and store the diagonal element (the species specific
# n-cube volumes)
intersectnDF_all <- data.frame(t(sapply(intersectionHulls_all, c)))
DF_diag_all <- unlist(diag((as.matrix(intersectnDF_all)), names=FALSE))


# calculate the proportion of overlap per column, based on the diagonal element residing in that column 
# perform the column-wise calculation on the first column, using the first element in the diagonal list
overlapProp_all <- data.frame()

for(i in 1:length(DF_diag_all)){
  overlapProp_all <- rbind(overlapProp_all, lapply(intersectnDF_all[,i], 
                                                   function(x) x/DF_diag_all[[i]]))
}
rownames(overlapProp_all) <- rownames(intersectnDF_all)


# make the results tidy and long; as.matrix part is essential to get var1 and var2, as well as the value
melted_overlapProp_all <- melt(as.matrix(overlapProp_all))

melted_overlapProp_all <- melted_overlapProp_all %>% mutate(Var1 = str_replace_all(Var1, "_", ". "),
                                                            Var2 = str_replace_all(Var2, "_", ". "))

#### VISUALIZATION

# plot the overlap proportions of the n-cubes
propOverlap_all <- ggplot(melted_overlapProp_all, aes(x=Var1, y=Var2)) + 
  geom_tile(aes(fill=value)) +
  #geom_text(aes(label=round(value,digits=1)),color='darkgrey') +
  scale_fill_gradient("Proportion\nspecies\noverlap", low='white', high="black") +
  theme(axis.title = element_blank(),
        axis.text = element_text(size=8),
        axis.text.x = element_text(angle=90, hjust=0.1)) +
  theme(legend.text = element_text(size=6),
        legend.margin = margin(2,2,6,2),
        legend.title = element_text(size=8),
        legend.background = element_rect(color="#cccccc")) +
  guides(fill = guide_colorbar(barwidth = 0.5, barheight = 2))

