library(data.table)
library(tidyverse)
library(rtracklayer)










genomeTracksPlot = function(chr, mystart, myend, gffpath, gene2transcriptspath, highlight_genes=NULL){
  
  require(data.table)
  require(tidyverse)
  require(rtracklayer)
  
  #'
  #'
  #'
  #'
  
    #load gff with rtracklayer and filter to genes only 
  gff = rtracklayer::import(gffpath) %>% 
    as.data.frame() %>%
    dplyr::rename("chrom" = 'seqnames') %>% 
    filter(chrom == chr, type %in% c("exon", "gene")) %>% 
    arrange(chrom, start) %>% 
    select(chrom, start, end, strand, type, ID, Parent) %>% 
    as.data.table()
  
  gene_names = fread(gene2transcriptspath) %>% 
    dplyr::rename("ID" = "GeneID")
  
  gff = left_join(gff, gene_names)
  
  ### limit to loci of interest
  df = gff %>% 
    filter(start > mystart, start < myend) %>% 
    mutate("strand" = case_when(strand == "-" ~ 1,
                                strand == "+" ~ 2))     
  
  exons = df %>% filter(type == 'exon')
  genes = df %>% filter(type == 'gene') %>% select(start, end, GeneName, ID, strand)
  
  if (!is.null(highlight_genes)){
    genes = genes %>% 
      mutate("highlight" = case_when(!GeneName %in% highlight_genes ~ 0,
                                   GeneName %in% highlight_genes ~ 1))

    exons$ID = str_remove(exons$Parent, "-RA")
    exons$Parent = NULL  
    a = genes %>% select(ID, highlight)
    exons = left_join(exons, a)
  } else {  
    exons$ID = str_remove(exons$Parent, "-RA")
    exons$Parent = NULL  
    a = genes %>% select(ID)
    exons = left_join(exons, a)
  }

  #add items columns and gather to long format
  exons$item = 1:nrow(exons)
  exons = gather(exons, "state", "pos", c("start", "end"))
  
  #get midpoint for labelling genes
  labels = genes %>% mutate("midpoint" = (start+end)/2)
  labels$GeneName = c(str_to_title(labels$GeneName)) #c(str_to_title(labels$Gene_name[1:9]), as.character(labels$Gene_name[10]))
  
  return(list(exons, labels))
}


gffpath = "resources/reference/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gff3"
gene2transcripts = "resources/exampleGene2TranscriptMap.tsv"
mystart = 28524000
myend = 28565000

data = genomeTracksPlot('2L', mystart, myend, gffpath = gffpath, gene2transcriptspath = gene2transcripts, highlight_genes = c("COEAE2F"))
exons = data[[1]]
labels = data[[2]]


labels[labels$ID == 'AGAP006227', 'GeneName'] = 'Coeae1f'
labels[labels$ID == 'AGAP006227', 'highlight'] = 1
exons[exons$ID == 'AGAP006227', 'highlight'] = 1



ggplot(exons, aes(x=pos, y=as.factor(strand))) +
 # annotate("rect", xmin=28480190, xmax=28483470, ymin=-Inf, ymax=Inf, fill=rgb(71/256, 60/256, 200/256, 0.9)) +
  geom_segment(aes(x=(min(pos)-500), xend=(max(pos)+500), y=as.factor(strand), yend=as.factor(strand)), colour="black", size=1) + 
 # geom_vline(aes(xintercept=28480190), colour="slateblue4",linetype="dashed") + # cnv dup1 breakpoints 
#  geom_vline(aes(xintercept=28483470), colour="slateblue4", linetype="dashed") + 
 # geom_vline(aes(xintercept=28497967), colour="slategrey", linetype="dotted")+
  geom_line(aes(color=as.factor(highlight), group=item, size = 8)) + 
  theme_light() +
  labs(y=NULL, x=NULL) + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(colour="black", size=14),
        axis.ticks.x = element_line(colour="black"),
        axis.title.x = element_text(size=13),
        #    plot.title = element_text(hjust = 0.5),
        #   plot.subtitle = element_text(hjust = 0.5),
        legend.position = "None",  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.9),
        panel.background = element_blank()) + 
  scale_x_continuous(breaks=seq(mystart-10000, myend+10000, by=5000)) + 
  geom_text(data=labels, aes(x=midpoint, label=GeneName, color=as.factor(labels$highlight), y=strand, vjust=-2), #label genes using labels data frame
            size=4,
            fontface = "bold") + 
  scale_color_manual(values= c("black", "darkorange")) 
  
  #+ 
 # geom_text(aes(x=28480190, y=1, label="Dup1", vjust=1.4, hjust=1.5), # dup 1 text label
 #           colour="black", 
#            size=5, 
#            angle=90) + 
#  geom_text(aes(x=28497967, y=1, label="I236M", vjust=-0.75, hjust=1.3), # i236m text label
#            colour="slategrey", 
#            size=5, 
 #           angle=90)


