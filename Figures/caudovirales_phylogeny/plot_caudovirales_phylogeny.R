rm(list=ls())

library(ggplot2)
library(ggtree)
library(ape)
library(plyr)
library(phytools)
library(phangorn)
library(ggnewscale)
library(phylobase)
library(stringr)
library(tidyr)
library(ggpubr)
library(grid)
library(gridExtra)
library(dplyr)
library(aptheme)

# tree data --------------------------------------------------------------------

# load tree
caudovirales_tree <-  read.tree(file = paste0("Caudoviricetes.raxml.supportFBP"))
#caudovirales_tree$tip.label <- gsub("NODE_59_length_34079_cov_116.672731||full", "NODE_59", caudovirales_tree$tip.label)

# root it using an inovirus outgroup
rooted_caudovirales_tree <- root(caudovirales_tree, which(caudovirales_tree$tip.label=='V00604'))

# take a quick look in base R
p <- ggtree(rooted_caudovirales_tree) + 
  geom_nodelab(aes(label=label), size=1, nudge_x=-0.01, nudge_y=0.25) +
  geom_tiplab(align= FALSE, linetype="dotted", linesize = 0.1, size = 2.5) + geom_point(colour='red')

p

# drop tips that are distracting 
rooted_caudovirales_tree <- drop.tip(rooted_caudovirales_tree, c('NC_047818', 'NC_070915'))

# shorten outgroup branch length 
#outgroup_index <- which(final.rooted.kobu$tip.label == 'NC_026314')
#branch_index <- which(final.rooted.kobu$edge[,2] == outgroup_index)
#final.rooted.kobu$edge.length[branch_index] <- 8  

# tree metadata ----------------------------------------------------------------

# load manual tree csv data for caudovirales tree
# and fill in for NA values 
caudovirales_manual <- read.csv(file=paste0('sequences.csv'), 
                        header=T, stringsAsFactors = F)

caudovirales_manual <- read.csv(file=paste0('sequences.csv'), 
                        header=T, stringsAsFactors = F, na = "")

caudovirales_manual[is.na(caudovirales_manual)] = "NA"

# use dplyr to only include columns that will be in the final tip label

colnames(caudovirales_manual)

caudovirales_manual <- caudovirales_manual %>%
  dplyr::select(Accession,Organism_Name, Species, Geo_Location, Host, Collection_Date)

# tip labels on manual csv and raxml tree do not match, so fix that 
# first make a new df with the labels numbered by root tip on the raxml tree
caudovirales.dat<- data.frame(Accession=rooted_caudovirales_tree$tip.label, 
                         num =1:length(rooted_caudovirales_tree$tip.label))

# then create a 3rd df and right join the original manual and the 2nd df 
# by accession_number onto this 3rd df 
caudovirales_join <- join(caudovirales_manual, caudovirales.dat, by = "Accession", match = "all", 
                 type = "right")

# write this out to a csv 
#write.csv(caudovirales_join, "Caudovirales_NCBI_Tree_Data_Combined_FLG_09APR2025.csv", row.names = F)

# read in modified copy
caudovirales_manual_mod <- read.csv(file=paste0('Caudovirales_NCBI_Tree_Data_Combined_FLG_09APR2025 copy.csv'), 
                                header=T, stringsAsFactors = F)

# tree plot --------------------------------------------------------------------

# check unique hosts that will be used to color tip labels 
# and assign colors to them
unique(caudovirales_manual_mod$Organism_Name)

# alphabetize the order of the legend labels
C=c("Arthrobacter phage", "Bacillus phage", "Bifidobacterium phage", "Burkholderia phage", 
    "Clostridium phage", "Curtobacterium phage", "Enterococcus phage", "Escherichia phage",
    "Gordonia phage", "Haemophilus phage", "Helicobacter phage", "Lactobacillus phage",
    "Lactococcus phage", "Leuconostoc phage", "Lokiarchaeia virus", 
    "Microbacterium phage", "Propionibacterium phage", "Salmonella phage",
    "Staphylococcus phage", "Streptococcus phage", "Vibrio phage")

# collapsing clades ------------------------------------------------------------

# extract the colors used in the plot so i can reference them when collapsing clades
ggcolors <- scales::hue_pal()(length(C))
names(ggcolors) <- C

# print the colors and their assigned genus to match 
print(ggcolors)

# convert the phylo object to a ggtree object
ggtree_obj <- ggtree(rooted_caudovirales_tree)

# here I am defining the clades I want to collapse
# you'll have to re-type in clades of interest so not a great system, but it will work for now 
clades_to_collapse <- c("NC_049965",
                        "NC_049976")

# get the MRCA node
mrca_node <- getMRCA(rooted_caudovirales_tree, clades_to_collapse)

# visualize again 
p <- ggtree(rooted_caudovirales_tree) %<+% caudovirales_manual_mod + 
  geom_tippoint(aes(color=Host), size=2) + 
  geom_tiplab(size=3) 
  #geom_nodelab(size=1) +
  #scale_color_manual(values=colz, breaks=C) + 
  #theme(legend.position = c(.70, .70), legend.title = element_blank())
p

# normalize the node labels to a range between 0 and 1
rooted_caudovirales_tree$node.label <- as.numeric(rooted_caudovirales_tree$node.label)

min_val <- min(rooted_caudovirales_tree$node.label, na.rm = TRUE)
max_val <- max(rooted_caudovirales_tree$node.label, na.rm = TRUE)
rooted_caudovirales_tree$node.label <- (rooted_caudovirales_tree$node.label - min_val) / (max_val - min_val) * 100

# extract the node labels and ensure they are correctly aligned with the nodes in the tree
node_labels <- data.frame(node = (Ntip(rooted_caudovirales_tree) + 1):(Ntip(rooted_caudovirales_tree) + Nnode(rooted_caudovirales_tree)),
                          node_label = rooted_caudovirales_tree$node.label)

# convert the tree to a tibble and join with node labels
tree_data <- as_tibble(rooted_caudovirales_tree) %>%
  left_join(node_labels, by = "node")

# renaming tip labels ----------------------------------------------------------

caudovirales_manual_mod$new_label <- caudovirales_manual_mod$Accession

# create new tip labels


caudovirales_manual_mod$new_label <- paste(caudovirales_manual_mod$Accession,
                            ifelse(is.na(caudovirales_manual_mod$Organism_Name), "", paste(" | ", caudovirales_manual_mod$Organism_Name)),
                            ifelse(is.na(caudovirales_manual_mod$Species), "", paste(" | ", caudovirales_manual_mod$Species)),
                            ifelse(is.na(caudovirales_manual_mod$Host), "", paste(" | ", caudovirales_manual_mod$Host)),
                            ifelse(is.na(caudovirales_manual_mod$Geo_Location), "", paste(" | ", caudovirales_manual_mod$Geo_Location)),
                            sep = "")




caudovirales_manual_mod$Accession <- caudovirales_manual_mod$new_label

rooted_caudovirales_tree$tip.label <- caudovirales_manual_mod$Accession

rooted_caudovirales_tree$tip.label

# final tree 

# plot the tree
p <- ggtree(rooted_caudovirales_tree) %<+% tree_data %<+% caudovirales_manual_mod +
  geom_tippoint(aes(color=Organism_Name, fill=Organism_Name), size=2) +
  new_scale_fill() +
  geom_nodepoint(aes(fill = node_label), shape = 21, size = 1, color = "black") + 
  #scale_fill_manual(values=posfilz) + 
  scale_fill_gradient(low = "white", high = "black", breaks = seq(0, 100, by = 20), labels = seq(0, 100, by = 20),
                      guide = guide_colorbar(direction = "horizontal")) +   geom_treescale(y = 5, x = 10, fontsize = 4, offset = 1, color = "black", width = 0.5) + 
  geom_tiplab(geom="label", label.size = 0, alpha=.3, size=2.8, show.legend=F) +
  #xlim(c(0,30)) +
  theme(legend.position = c(0.7, 0.5), legend.title = element_blank())

p

# collapse clades 
p <- p %>%
  collapse(node = 52) + 
  geom_point2(aes(subset = (node == 52)), shape = 22, size = 1.5, fill = "#00B6EB") 
p <- p %>%
  collapse(node = 85) + 
  geom_point2(aes(subset = (node == 85)), shape = 22, size = 1.5, fill = "#F8766D") 
p <- p %>%
  collapse(node = 89) + 
  geom_point2(aes(subset = (node == 89)), shape = 22, size = 1.5, fill = "#00ABFD") 
p <- p %>%
  collapse(node = 94) + 
  geom_point2(aes(subset = (node == 94)), shape = 22, size = 1.5, fill = "#00B6EB") 
p <- p %>%
  collapse(node = 93) + 
  geom_point2(aes(subset = (node == 93)), shape = 22, size = 1.5, fill = "#FF63B9") 
p <- p %>%
  collapse(node = 95) + 
  geom_point2(aes(subset = (node == 95)), shape = 22, size = 1.5, fill = "#EB8335") 

# rename summarized clades and plot again
p <- p + 
  geom_text2(aes(subset = (node == 52), label = "Lactococcus Phage - Collapsed Clade"), 
             nudge_x = 0.73, nudge_y = 0.01, size = 2.7, color = "black") +
  geom_text2(aes(subset = (node == 85), label = "Arthrobacter Phage - Collapsed Clade"), 
             nudge_x = 0.72, nudge_y = 0.01, size = 2.7, color = "black") +
  geom_text2(aes(subset = (node == 89), label = "Leuconostoc Phage - Collapsed Clade"), 
             nudge_x = 0.73, nudge_y = 0.01, size = 2.7, color = "black") +
  geom_text2(aes(subset = (node == 94), label = "Lactococcus Phage - Collapsed Clade"), 
             nudge_x = 0.73, nudge_y = 0.01, size = 2.7, color = "black") + 
  geom_text2(aes(subset = (node == 93), label = "Streptococcus Phage - Collapsed Clade"), 
             nudge_x = 0.76, nudge_y = 0.01, size = 2.7, color = "black") +
  geom_text2(aes(subset = (node == 95), label = "Bacillus Phage - Collapsed Clade"), 
             nudge_x = 0.64, nudge_y = 0.01, size = 2.7, color = "black") 

p


# now export for visualization in adobe 
ggsave(file = paste0("Caudovirales_Phylo.png"),
       units="mm",  
       width=400, 
       height=250, 
       #limitsize = F,
       scale=1)#, 
