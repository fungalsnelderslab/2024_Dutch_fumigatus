setwd("/Users/ben/Desktop/")
library(ggplot2)
library(ggtree)
library(phangorn)
library(ape)
library(ggnewscale)
library(tidytree)
library(cowplot)
#weird glitches in ggtree, but this fixes it:
offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")
#this is the raw tree
whole_tree <- read.tree("/Users/ben/Downloads/global_combined_partial_filtered_nohets_smaller.min4.phy.treefile")
rooted_tree <- root(whole_tree,node=2112)

for (edge in 1:length(rooted_tree$edge.length)){
  if (rooted_tree$edge.length[edge] > 0.1){
    rooted_tree$edge.length[edge] <- rooted_tree$edge.length[edge]/25
    print(edge)
    print(rooted_tree$edge.length[edge])
    print(rooted_tree$edge[edge,])
    }
}
#these are the nodes affected

#this gets the node of the fumigatus strains
MRCA(rooted_tree,c("85C12","86C22","36B7","82C23","4A19","67A3"))

whole_plotting <- ggtree(rooted_tree,right=TRUE,lwd=0.1) + geom_rootedge()
open_tree(whole_plotting,400)
zoomed <- open_tree(scaleClade(rotate_tree(whole_plotting,30),1281,0.08),50) %>% collapse(1281,'max',fill="grey79",clade_name="fumigatus")


metadata <- read.delim("/Users/ben/Downloads/global_dataset_metadata.csv",sep="\t")

table(metadata$azole.resistance.sensitive)
metadata$azole.resistance.sensitive[metadata$azole.resistance.sensitive == ""] <- NA
metadata$azole.resistance.sensitive[metadata$azole.resistance.sensitive == "-"] <- NA
metadata$azole.resistance.sensitive[metadata$azole.resistance.sensitive == "Unknown"] <- NA
table(metadata$azole.resistance.sensitive)

metadata$azole_resistant <- FALSE
metadata$azole_susceptible <- FALSE

metadata$azole_resistant[metadata$azole.resistance.sensitive == "Resistant"] <- TRUE
metadata$azole_susceptible[metadata$azole.resistance.sensitive == "Susceptible"] <- TRUE

table(metadata$clinical.environmental)
metadata$clinical.environmental[metadata$clinical.environmental == "Other"] <- "Unknown"
metadata$clinical.environmental[metadata$clinical.environmental == "Unknown"] <- "Unknown"
metadata$clinical.environmental[grepl("Clinical",metadata$clinical.environmental)] <- "Clinical"
table(metadata$clinical.environmental)

mutations <- read.csv("/Users/ben/Downloads/antifungal_mutations.tsv",row.names = 1)
rownames(mutations) <- gsub("out.","",rownames(mutations),perl=T)

metadata_plotting <- metadata[,c("azole.resistance.sensitive","clinical.environmental")]
colnames(metadata_plotting) <- c("azole_resistant","clinical.environmental")
row.names(metadata_plotting) <- metadata$accession.id
metadata_plotting <- merge(mutations,metadata_plotting,by="row.names")
row.names(metadata_plotting) <- metadata_plotting$Row.names
metadata_plotting <- metadata_plotting[,-1]


resistance <- metadata_plotting[,"azole_resistant",drop=FALSE]
colnames(resistance) <- "Triazole\nresistant"
TR46 <- metadata_plotting[,"tr46", drop=FALSE]
colnames(TR46) <- "TR46"
TR34 <- metadata_plotting[,"tr34", drop=FALSE]
colnames(TR34) <- "TR34"
cytB <- metadata_plotting[,"cytB", drop=FALSE]
colnames(cytB) <- expression("cytB\nG143A")
benA <- metadata_plotting[,"benA", drop=FALSE]
colnames(benA) <- "benA\nF219Y"
sdhB <- metadata_plotting[,"sdhB", drop=FALSE]
colnames(sdhB) <- "sdhB\nH270Y"
msh6 <- metadata_plotting[,"msh6", drop=FALSE]

#to get clade of divergent group
getMRCA(rooted_tree,c("V254-50","SRR25183057"))

#to get clade of outgroup
getMRCA(rooted_tree,c("out.Afis_SRR11363404","SRR15010323"))

resistance_tree <- gheatmap(whole_plotting,resistance,offset=-0,width=0.15,font.size=3,colnames_position="top",colnames_offset_y=35)+
  theme(legend.position="none") +
  #geom_tiplab(size=0.5)+
  scale_fill_manual(breaks=c("Resistant","Susceptible","unknown"),values=c("black","grey90","white"),na.value = "white") +
  geom_cladelabel(node=2165,label="1") +
  geom_cladelabel(node=2113,label="2") +
  geom_cladelabel(node=1280,label="3") + 
  geom_cladelabel(node=1282,label="4") 

p1 <- gheatmap(resistance_tree + new_scale_fill(),TR34,offset=0.015,width=0.12,font.size=3,colnames_position="top",colnames_offset_y=35) +
  theme(legend.position="none") +
  scale_fill_manual(breaks=c("T","TTTCATTCGGCTCAGCACACATCCGGACCGCGTGA"),values=c("grey90","black"),na.value = "white")

p2 <- gheatmap(p1 + new_scale_fill(),TR46,offset=0.025,width=0.12,font.size=3,colnames_position="top",colnames_offset_y=35) + 
  theme(legend.position="none") +
  scale_fill_manual(breaks=c("G","G|G","GCAACTTTCATTCGGCTCAGCACACATCCGGACCGCGTGATTCTAGA","GCAACTTTCATTCGGCTCAGCACACATCCGGACCGCGTGATTCTAGACAACTTTCATTCGGCTCAGCACACATCCGGACCGCGTGATTCTAGA"),values=c("grey90","grey90","black","black"),na.value = "white")

p3 <- gheatmap(p2 + new_scale_fill(),cytB,offset=0.035,width=0.12,font.size=3,colnames_position="top",colnames_offset_y=35) +
  theme(legend.position="none") +
  scale_fill_manual(breaks=c("G","C"),values=c("grey90","black"),na.value="white")

p4 <- gheatmap(p3 + new_scale_fill(),benA,offset=0.045,width=0.12,font.size=3,colnames_position="top",colnames_offset_y=35) +
  theme(legend.position="none") +
  scale_fill_manual(breaks=c("A","A|T","T","T|T"),values=c("grey90","grey90","black","black"),na.value="white")

p5 <- gheatmap(p4 + new_scale_fill(),sdhB,offset=0.055,width=0.12,font.size=3,colnames_position="top",colnames_offset_y=35) +
  theme(legend.position="none") +
  scale_fill_manual(breaks=c("G","A"),values=c("grey90","black"),na.value="white")

p5 + scale_y_continuous(limits=c(0,1300))
#p6
#p5
ggsave("test4.svg",width=10,height=10)
ggsave("test4.png",width=8,height=8,dpi=300)
