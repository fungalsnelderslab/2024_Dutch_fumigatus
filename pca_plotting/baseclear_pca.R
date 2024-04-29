library(ggplot2)
library(dplyr)
library(cowplot)

# read in the eigenvectors, produced in PLINK
eigenvec <- read.table('/Users/ben/Desktop/baseclear/pca/dutch_data_clean_prune.eigenvec', header = FALSE, skip=0, sep = ' ')
rownames(eigenvec) <- eigenvec[,2]
eigenvec <- eigenvec[,3:ncol(eigenvec)]
colnames(eigenvec) <- paste('PC', c(1:20), sep = '')
eigenvec$Sample <- row.names(eigenvec)
class(eigenvec)

head(eigenvec)

meta <- read.delim("/Users/ben/Desktop/baseclear/pca/baseclear_metadata.csv",sep=";")

data <- left_join(eigenvec,meta,)

colors <- c("wt R"="purple1","46"="dodgerblue3","92"="dodgerblue3","34"="dodgerblue3","wt"="orange2")
shapes <- c("wt R"=18,"46"=15,"92"=17,"34"=16,"wt"=18)

clinvenv12 <- ggplot(data) + theme_classic(base_size=12) +
  geom_point(aes(x=PC1,y=PC2,color=Source),size=2) +
  labs(x="PC1 (22.1%)",y="PC2 (8.8%)")

data$status <- NA
data$status[data$Haplotype=="wt"] <- "Sensitive"
data$status[data$Haplotype %in% c("34","46","92","wt R")] <- "Resistant"
resvsens12 <- ggplot(data) + theme_classic(base_size=12) +
  geom_point(aes(x=PC1,y=PC2,color=status),size=2) +
  labs(x="PC1 (22.1%)",y="PC2 (8.8%)",color="Azole phentoype")

cyptype12 <- ggplot(data) + theme_classic(base_size=12) +
  geom_point(aes(x=PC1,y=PC2,color=Haplotype),size=2) +
  labs(x="PC1 (22.1%)",y="PC2 (8.8%)")

location12 <- ggplot(data) + theme_classic(base_size=12) +
  geom_point(aes(x=PC1,y=PC2,color=Location..province.),size=2) +
  labs(x="PC1 (22.1%)",y="PC2 (8.8%)",color="Province")


plot_grid(clinvenv12,
          resvsens12,
          cyptype12,
          location12,
          nrow=4,align="hv",
          labels="auto",label_fontface="plain")

ggsave("/Users/ben/Desktop/baseclear/pca/pca_plotting.svg",width=6,height=8)
