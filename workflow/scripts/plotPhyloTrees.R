### Define input ####

# input files
library(ape)
library(phytools)
library(stringr)
library(data.table)
library(dplyr)
library(glue)
f = glue

prefix    = "results/phylo/"
phy_list  = c("results/phylo/coeae1f_focal.fasta.treefile", "results/phylo/coeae1f_upstream.fasta.treefile", "results/phylo/coeae1f_downstream.fasta.treefile")
phy_name  = c("Focal haplotype","Upstream","Downstream")


name = 'coeae1f'


plot_phylo = function(phi, meta, var, region, name, min_edge_length = 5e-5){
  
  phi = f("results/phylo/{name}_{region}.fasta.treefile")
  meta_path = f("results/phylo/{name}_{region}.metadata.tsv")
  
  # tree to distance matrix
  phy    = read.tree(phi)
  dis    = cophenetic.phylo(phy)
  meta = fread(meta_path) %>% as.data.frame()
  phy_p = midpoint.root(phy)
  meta = meta %>% arrange(factor(hap, levels = phy_p$tip.label))
  
  phy_p$tip_label_sep  = gsub("_"," ",phy_p$tip.label)
  phy_p$tip_species    = meta$aim_species
  phy_p$tip_country = meta$country
  phy_p$tip_cluster   = ""
  phy_p$tip_karyotype  = ""
  #phy_p$tip_color     = brewe
  phy_p$tip.label      = rep("·", length(phy$tip.label))
  #phy_p$edge.length[phy_p$edge.length >  0.005] = 0.005
  phy_p$edge.length[phy_p$edge.length == 0]    = min_edge_length    #5e-5
  
  if (var == 'karyotype'){
    meta = meta %>% mutate("tip_color" = case_when(karyotype == '2l+a' ~ 'bisque2',
                                                   karyotype == '2la' ~ 'dodgerblue', TRUE ~ 'grey'))
  } else if (var == 'Sweep IDs'){
    meta$tip_color = rainbow(length(unique(meta[,var])))[as.factor(meta[, var])]
    meta = meta %>% mutate(tip_color = case_when(`Sweep IDs` == 'wt' ~ 'grey',
                                                 `Sweep IDs` %in% c("mela", "meru", "quad") ~ 'black',
                                                 TRUE ~ tip_color))
  } else if (var == 'aim_species'){
    meta = meta %>% mutate("tip_color" = case_when(aim_species == 'gambiae' ~ 'indianred',
                                                   aim_species == 'coluzzii' ~ 'dodgerblue',
                                                   aim_species == 'arabiensis' ~ 'aquamarine3',
                                                   aim_species == 'mela' ~ 'cornsilk',
                                                   aim_species == 'meru' ~ 'cornsilk4',
                                                   aim_species == 'quad' ~ 'darkolivegreen',
                                                   TRUE ~ 'grey'))
  }
  plot.phylo(phy_p, type="unr",
             use.edge.length=T, show.tip.label=T, show.node.label=F,
             tip.color = meta$tip_color, lab4ut = "axial",
             edge.color = "slategray3",
             font = 1, edge.width = 2, node.depth=1, cex=5,
             main=f("{name} {region}"), no.margin = TRUE)
}


plot_phylo(phi, meta, var='karyotype', region="focal", name=name)
plot_phylo(phi, meta, var='aim_species', region="focal", name=name)
plot_phylo(phi, meta, var='Sweep IDs', region="focal", name=name)

plot_phylo(phi, meta, var='karyotype', region="downstream", name=name)
plot_phylo(phi, meta, var='aim_species', region="downstream", name=name)
plot_phylo(phi, meta, var='Sweep IDs', region="downstream", name=name)

plot_phylo(phi, meta, var='karyotype', region="upstream", name=name)
plot_phylo(phi, meta, var='aim_species', region="upstream", name=name)
plot_phylo(phi, meta, var='Sweep IDs', region="upstream", name=name)

# 6 sweeps
# 1 arabiensis - CNV
# 3 gambiae one of which has CNV
# 1 coluzzii
# 1 mixed gambiae and coluzzii


