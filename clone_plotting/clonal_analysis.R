library(ggplot2)
library(cowplot)
library(UpSetR)
library(dplyr)
library(tidyr)
library(igraph)
#library(ggraph)
#library(tidygraph)

#now read in the global sample variant data----
global_diversity <- read.table("/Users/ben/Downloads/plink.genome",header=T)
global_metadata <- read.delim("/Users/ben/Downloads/global_dataset_metadata.csv",sep="\t")
global_metadata$new_name <- paste0(global_metadata$accession_id,"_",global_metadata$country,"_",global_metadata$clinical_environmental)

# first take alook at the raw data
ggplot(global_diversity) +
  geom_histogram(aes(x=DST),bins=400) +
  theme_classic() +
  labs(x="Genetic Distance",y="Pairs")

#now make plots
total <- ggplot(global_diversity) +
  geom_histogram(aes(x=DST),bins=400) + xlim(0.88,1) +
  theme_classic() +
  labs(x="Genetic Distance",y="Pairs")

close <- ggplot(global_diversity) +
  geom_histogram(aes(x=DST),bins=400) + xlim(0.995,1) +
  theme_classic() +
  geom_vline(aes(xintercept=0.998),lty=2)+
  labs(x="Genetic Distance",y="")

#read in the four locus analysis----
t <- read.delim("/Users/ben/Downloads/locus_typing_results.tsv")

#the file contains the read coverage for four control loci and 4 polymorphic genes
head(t)

# Compute averages for the specified control genes
target_regions <- c("NC_007196.1:1539720-1544676", "NC_007195.1:4660561-4662944","NC_007198.1:242866-244662","NC_007199.1:3491106-3496228")

avg_values <- t %>%
  filter(region %in% target_regions) %>%
  group_by(file) %>%
  summarise(Average = mean(depth, na.rm = TRUE))

#make the data long format, and then Recalculate third column based on the sample-specific average
t <- t %>%
  left_join(avg_values, by = "file") %>%
  mutate(Value = ifelse(!is.na(Average), Average, depth)) %>%
  select(-Average)

t$normalized_depth <- t$depth/t$Value

#after recalculating, make the data wide format again, and remove the now redundant control regions
wide <- pivot_wider(t[,c(1,2,5)],names_from=region,values_from=normalized_depth)
wide <- wide[,c(1,2,4,6,8)]
colnames(wide) <- c("Sample","mat","hetA","hetB","hetC")

#plot the histograms of read depth across the dataset
mat  <- ggplot(wide) + geom_histogram(aes(x=mat),bins=100) +  geom_vline(aes(xintercept=0.5),lty=2) + scale_x_continuous(limits = c(-0.3,1.5)) + theme_classic() + labs(x="Normalized Read Depth",y="Isolates")
hetA <- ggplot(wide) + geom_histogram(aes(x=hetA),bins=100) + geom_vline(aes(xintercept=0.5),lty=2) + scale_x_continuous(limits = c(-0.3,1.5)) + theme_classic() + labs(x="Normalized Read Depth",y="Isolates")
hetB <- ggplot(wide) + geom_histogram(aes(x=hetB),bins=100) + geom_vline(aes(xintercept=0.5),lty=2) + scale_x_continuous(limits = c(-0.3,1.5)) + theme_classic() + labs(x="Normalized Read Depth",y="Isolates")
hetC <- ggplot(wide) + geom_histogram(aes(x=hetC),bins=100) + geom_vline(aes(xintercept=0.5),lty=2) + scale_x_continuous(limits = c(-0.3,1.5)) + theme_classic() + labs(x="Normalized Read Depth",y="Isolates")

plot_grid(total,close,mat,hetA,hetB,hetC,ncol=2,labels=c("A)","B)", "C) mating-type",
                                             "D) hetA",
                                             "E) hetB",
                                             "F) hetC"),
                                              hjust=-0.5,
                                              rel_heights=c(1.7,1,1),
                                              label_fontface = "plain")
ggsave("/Users/ben/Downloads/Supp_S2.svg",width=8,height=12)

#now convert to binary barcode format
wide$mat[wide$mat <= 0.4] <- 0 #reference allele with Af293
wide$mat[wide$mat >  0.4] <- 1 #alternate allele compared to Af293

wide$hetA[wide$hetA <= 0.4] <- 0 #reference allele with Af293
wide$hetA[wide$hetA >  0.4] <- 1 #alternate allele compared to Af293

wide$hetB[wide$hetB <= 0.4] <- 0 #reference allele with Af293
wide$hetB[wide$hetB >  0.4] <- 1 #alternate allele compared to Af293

wide$hetC[wide$hetC <= 0.4] <- 0 #reference allele with Af293
wide$hetC[wide$hetC >  0.4] <- 1 #alternate allele compared to Af293

wide$combined <- paste0(wide$mat,wide$hetA,wide$hetB,wide$hetC)

#not used in the formal analysis, but count alleles
sum(wide$mat)
sum(wide$hetA)
sum(wide$hetB)
sum(wide$hetC)
nrow(wide)

# also unused, a summary of distribution of barcodes
table(wide$combined)

write.csv(wide,"/Users/ben/Downloads/locus_typing_inferred_results.csv",row.names = F)




#we want to merge 4 locus genotypes with clonal data
clonal_data <- subset(global_diversity,DST>0.998)
wide$Sample <- gsub("bams/.[0-9A-Za-z]+.","",wide$Sample,perl=T)
wide$Sample <- sub(".sorted.bam","",wide$Sample)

#this is the analysis for the pairs that are more than 99.8% similar
pairdata <- clonal_data %>% 
  left_join(wide, by=c("FID1" = "Sample")) %>%
  rename(barcode_FID1 = combined) %>%
  left_join(wide, by = c("FID2" = "Sample")) %>%
  rename(barcode_FID2 = combined)

sum(pairdata$barcode_FID1 == pairdata$barcode_FID2,na.rm=T)
sum(pairdata$barcode_FID1 != pairdata$barcode_FID2,na.rm=T)

#now we calculate for all samples
alldata <- global_diversity %>% 
  left_join(wide, by=c("FID1" = "Sample")) %>%
  rename(barcode_FID1 = combined) %>%
  left_join(wide, by = c("FID2" = "Sample")) %>%
  rename(barcode_FID2 = combined)

sum(alldata$barcode_FID1 == alldata$barcode_FID2,na.rm=T)
sum(alldata$barcode_FID1 != alldata$barcode_FID2,na.rm=T)



clonal_data <- subset(global_diversity,DST>0.998)
clonal_data <- clonal_data[order(clonal_data$FID1),]
global_metadata <- read.delim("/Users/ben/Downloads/global_dataset_metadata.tsv",sep="\t")
global_metadata$new_name <- paste0(global_metadata$accession_id,"_",global_metadata$country,"_",global_metadata$clinical_environmental)

#we want to add study names and countries to sample names
for (i in 1:nrow(clonal_data)){
  print(i)
  print(clonal_data$FID1[i])
  print(clonal_data$FID2[i])
  clonal_data$FID1[i] <- grep(clonal_data$FID1[i],global_metadata$new_name,value=T)
  clonal_data$FID2[i] <- grep(clonal_data$FID2[i],global_metadata$new_name,value=T)
}

#optional, if you want to only keep the Dutch clonal groups
clonal_data <- clonal_data[grepl("Netherlands",clonal_data$FID1),]
clonal_data <- clonal_data[grepl("Netherlands",clonal_data$FID2),]

#Now we want to get the list of individuals that are in this clonal dataset
#first we get all samples from the clonal datset
all_names <- unique(sort(c(global_diversity$FID1,global_diversity$FID2)))
clonal_names <- unique(sort(c(clonal_data$FID1,clonal_data$FID2)))
#but these are listed several times, so cannot be used
length(unique(all_names))
length(unique(clonal_names))
g <- graph_from_data_frame(clonal_data[,c(1,3)])

#this makes the connected clusters
clonal_groups <- clusters(g)
#looking at the membershp
clonal_groups$csize

#for the upset plot, we need to get which groups have UK, Netherlands, Germany, Japan, US members
GERMANY_clinical_clonal_group_members <- c()
GERMANY_environ_clonal_group_members <- c()
UK_clinical_clonal_group_members <- c()
UK_environ_clonal_group_members <- c()
USA_clinical_clonal_group_members <- c()
USA_environ_clonal_group_members <- c()
DUTCH_clinical_clonal_group_members <- c()
DUTCH_environ_clonal_group_members <- c()
JAPAN_clinical_clonal_group_members <- c()
#now we want to look at which countries included
for (i in 1:clonal_groups$no){
  curr_clone <- clonal_groups$membership[clonal_groups$membership == i]
  curr_names <- names(curr_clone)
  print(curr_names)
  if (sum(grepl("Netherlands_Clinical",curr_names))){
    DUTCH_clinical_clonal_group_members <- c(DUTCH_clinical_clonal_group_members,i)}
  if (sum(grepl("Netherlands_Environmental",curr_names))){
    DUTCH_environ_clonal_group_members <- c(DUTCH_environ_clonal_group_members,i)}
  
  if (sum(grepl("Germany_Clinical",curr_names))){
    GERMANY_clinical_clonal_group_members <- c(GERMANY_clinical_clonal_group_members,i)}
  if (sum(grepl("Germany_Environmental",curr_names))){
    GERMANY_environ_clonal_group_members <- c(GERMANY_environ_clonal_group_members,i)}
  
  if (sum(grepl("UK_Clinical",curr_names)) > 0){
    UK_clinical_clonal_group_members <- c(UK_clinical_clonal_group_members,i)}
  if (sum(grepl("UK_Environmental",curr_names))){
    UK_environ_clonal_group_members <- c(UK_environ_clonal_group_members,i)}
  
  if (sum(grepl("USA_Clinical",curr_names)) > 0){
    USA_clinical_clonal_group_members <- c(USA_clinical_clonal_group_members,i)}
  if (sum(grepl("USA_Environmental",curr_names)) > 0){
    USA_environ_clonal_group_members <- c(USA_environ_clonal_group_members,i)}
  if (sum(grepl("Japan",names(curr_clone))) > 0) {
    JAPAN_clinical_clonal_group_members <- c(JAPAN_clinical_clonal_group_members,i)}
}

listInput <- list(UK_clinical=UK_clinical_clonal_group_members,
                  UK_environmental=UK_environ_clonal_group_members,
                  USA_clinical=USA_clinical_clonal_group_members,
                  USA_environmental=USA_environ_clonal_group_members,
                  Germany_clinical=GERMANY_clinical_clonal_group_members,
                  Germany_environmental=GERMANY_environ_clonal_group_members,
                  Netherlands_clinical=DUTCH_clinical_clonal_group_members,
                  Netherlands_environmental=DUTCH_environ_clonal_group_members,
                  Japan_clinical=JAPAN_clinical_clonal_group_members)


# Combine both vectors to get all unique numbers
all_numbers <- sort(unique(c(listInput$UK_clinical, listInput$UK_environmental,
                             listInput$USA_clinical,listInput$USA_environmental,
                             listInput$Germany_clinical,listInput$Germany_environmental,
                             listInput$Netherlands_clinical,listInput$Netherlands_environmental,
                             listInput$Japan_clinical)))

# Create a matrix with zeros
table_matrix <- matrix(0, nrow = length(all_numbers), ncol = 9, dimnames = list(all_numbers, c("UK_clinical", "UK_environmental","USA_clinical","USA_environmental","Germany_clinical","Germany_environmental","Netherlands_clinical","Netherlands_environmental","Japan_clinical")))

# Fill in the matrix with values from the vectors
table_matrix[match(listInput$UK_clinical, all_numbers), "UK_clinical"] <- listInput$UK_clinical
table_matrix[match(listInput$UK_environmental, all_numbers), "UK_environmental"] <- listInput$UK_environmental
table_matrix[match(listInput$USA_clinical, all_numbers), "USA_clinical"] <- listInput$USA_clinical
table_matrix[match(listInput$USA_environmental, all_numbers), "USA_environmental"] <- listInput$USA_environmental
table_matrix[match(listInput$Germany_clinical, all_numbers), "Germany_clinical"] <- listInput$Germany_clinical
table_matrix[match(listInput$Germany_environmental, all_numbers), "Germany_environmental"] <- listInput$Germany_environmental
table_matrix[match(listInput$Netherlands_clinical, all_numbers), "Netherlands_clinical"] <- listInput$Netherlands_clinical
table_matrix[match(listInput$Netherlands_environmental, all_numbers), "Netherlands_environmental"] <- listInput$Netherlands_environmental
table_matrix[match(listInput$Japan_clinical, all_numbers), "Japan_clinical"] <- listInput$Japan_clinical



# Convert matrix to data frame
table_df <- t(as.data.frame(table_matrix))

# Print the data frame
print(table_df)
write.csv(table_df,"/Users/ben/Downloads/test.csv")

#if you want to calculate the number of 


#now to make the upset plot

pdf("/Users/ben/Desktop/baseclear/May7.upset_plot.pdf",width=6.6,height=4)
upset(fromList(listInput), order.by="freq",nsets=9,keep.order=T,
      sets=c("Netherlands_clinical","Netherlands_environmental","UK_clinical","UK_environmental","USA_clinical","USA_environmental","Germany_clinical","Germany_environmental","Japan_clinical"),
      mainbar.y.label = "# Clonal Groups",show.numbers="yes")
dev.off()
