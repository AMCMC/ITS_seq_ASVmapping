library(dada2);packageVersion("dada2")
library(phyloseq);packageVersion("phyloseq")
library(ggplot2);packageVersion("ggplot2")
source("../../MicrobiotaCentre/MiCA_ITS/Scripts/taxa_facet_barplot_asv.R")

#### regular bigdata approach ####

fqs <- list.files("./", pattern = "I.*fastq.gz")

dereps <- list()
for (file in fqs){
  dereps[[file]] <- derepFastq(file)
}
names(dereps) <- gsub("_.*","",names(dereps))

err8 <- readRDS("ITS_0008.R1.RDS")
err9 <- readRDS("ITS_0009.R1.RDS")

ddF <- list()
for (i in names(dereps)){
  if (substr(i,2,2)==8){
    ddF[[i]] <- dada(derep = dereps[[i]], err = err8)
  }else{
    ddF[[i]] <- dada(derep = dereps[[i]], err = err9)
  }
}

topx=10
topxseqs <- c()
for (i in dereps){
  topxseqs <- c(topxseqs,names(i$uniques[1:topx]))
}

length(unique(topxseqs)) # track a total of 13 sequences
table(table(topxseqs)) #6 sequences are in the top10 in all samples 


#build dataframe
df <- data.frame()
for(i in names(dereps)){
df <- rbind(df, data.frame(sequence=names(dereps[[i]]$uniques[1:topx]),
           Abundance=unname(dereps[[i]]$uniques[1:topx]),
           ASV=ddF[[i]]$sequence[ddF[[i]]$map[1:topx]],
           SampleID=i))
}

df$ASV2 <- df$ASV
levels(df$ASV2) <- paste0("ASV-",1:length(levels(df$ASV2)))
df$sequence2 <- df$sequence
levels(df$sequence2) <- paste0("sequence-",1:length(levels(df$sequence2)))
df$ASVrep <- as.character(df$sequence)==as.character(df$ASV)
df$ASVrep2[df$ASVrep] <- as.character(df$ASV2[df$ASVrep])

ggplot(df, aes(x = sequence2, y = Abundance, fill=ASV2)) + 
  geom_bar(stat="identity") + 
  facet_grid(SampleID ~ ., scale="free_y") + 
  theme(axis.text.x = element_text(angle=90, hjust=1))

ggplot(df, aes(x = sequence2, y = Abundance, fill=ASVrep2)) + 
  geom_bar(stat="identity") + 
  facet_grid(SampleID ~ ., scale="free_y") + 
  theme(axis.text.x = element_text(angle=90, hjust=1))

### identical sequences in the various samples get assigned to different ASVs
#### try to resolve this issue with priors ####

priors <- unique(topxseqs)
ddF2 <- list()
for (i in names(dereps)){
  if (substr(i,2,2)==8){
    ddF2[[i]] <- dada(derep = dereps[[i]], err = err8, priors = priors)
  }else{
    ddF2[[i]] <- dada(derep = dereps[[i]], err = err9, priors = priors)
  }
}

#build dataframe
df2 <- data.frame()
for(i in names(dereps)){
  df2 <- rbind(df2, data.frame(sequence=names(dereps[[i]]$uniques[1:topx]),
                             Abundance=unname(dereps[[i]]$uniques[1:topx]),
                             ASV=ddF2[[i]]$sequence[ddF2[[i]]$map[1:topx]],
                             SampleID=i))
}

df2$ASV2 <- df2$ASV
levels(df2$ASV2) <- paste0("ASV-",1:length(levels(df2$ASV2)))
df2$sequence2 <- df2$sequence
levels(df2$sequence2) <- paste0("sequence-",1:length(levels(df2$sequence2)))
df2$ASVrep <- as.character(df2$sequence)==as.character(df2$ASV)
df2$ASVrep2[df2$ASVrep] <- as.character(df2$ASV2[df2$ASVrep])

ggplot(df2, aes(x = sequence2, y = Abundance, fill=ASV2)) + 
  geom_bar(stat="identity") + 
  facet_grid(SampleID ~ ., scale="free_y") + 
  theme(axis.text.x = element_text(angle=90, hjust=1))

ggplot(df2, aes(x = sequence2, y = Abundance, fill=ASVrep2)) + 
  geom_bar(stat="identity") + 
  facet_grid(SampleID ~ ., scale="free_y") + 
  theme(axis.text.x = element_text(angle=90, hjust=1))

#### try to resolve this issue with 1 specific prior ####

priors <- as.character(df2$sequence[df2$sequence2=="sequence-4"][1])
ddF3 <- list()
for (i in names(dereps)){
  if (substr(i,2,2)==8){
    ddF3[[i]] <- dada(derep = dereps[[i]], err = err8, priors = priors)
  }else{
    ddF3[[i]] <- dada(derep = dereps[[i]], err = err9, priors = priors)
  }
}

#build dataframe
df3 <- data.frame()
for(i in names(dereps)){
  df3 <- rbind(df3, data.frame(sequence=names(dereps[[i]]$uniques[1:topx]),
                               Abundance=unname(dereps[[i]]$uniques[1:topx]),
                               ASV=ddF3[[i]]$sequence[ddF3[[i]]$map[1:topx]],
                               SampleID=i))
}

df3$ASV2 <- df3$ASV
levels(df3$ASV2) <- paste0("ASV-",1:length(levels(df3$ASV2)))
df3$sequence2 <- df3$sequence
levels(df3$sequence2) <- paste0("sequence-",1:length(levels(df3$sequence2)))
df3$ASVrep <- as.character(df3$sequence)==as.character(df3$ASV)
df3$ASVrep2[df3$ASVrep] <- as.character(df3$ASV2[df3$ASVrep])

ggplot(df3, aes(x = sequence2, y = Abundance, fill=ASV2)) + 
  geom_bar(stat="identity") + 
  facet_grid(SampleID ~ ., scale="free_y") + 
  theme(axis.text.x = element_text(angle=90, hjust=1))

ggplot(df3, aes(x = sequence2, y = Abundance, fill=ASVrep2)) + 
  geom_bar(stat="identity") + 
  facet_grid(SampleID ~ ., scale="free_y") + 
  theme(axis.text.x = element_text(angle=90, hjust=1))

#### how to the sequences relate ?
levels(df3$sequence)[c(1:7)]
# they clearly differ in the length of the homopolymer

#### is the issue due to the indel?


nwalign(names(dereps[[1]]$uniques)[2],names(dereps[[1]]$uniques)[3])
drmod <- dereps[[1]]

#make a substitution rather than an indel
names(drmod$uniques)[2] <- "AAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTAAAGAAATTTAATAATTTTGAAAATGGATTTTTTTTTTTAGTTTTGGCAAGAGCATGAGAGCTTTTACTGGGC"

ddFmod <- dada(derep = drmod, err = err8)
ddFref <- dada(derep = dereps[[1]], err = err8)

#### truncate homopolymers in forward read ####

fqs <- list.files("./", pattern = "I.*fastq.gz")

dereps <- list()
for (file in fqs){
  dereps[[file]] <- derepFastq(file)
}
names(dereps) <- gsub("_.*","",names(dereps))

err8 <- readRDS("ITS_0008.R1.RDS")
err9 <- readRDS("ITS_0009.R1.RDS")

ddF <- list()
for (i in names(dereps)){
  if (substr(i,2,2)==8){
    ddF[[i]] <- dada(derep = dereps[[i]], err = err8)
  }else{
    ddF[[i]] <- dada(derep = dereps[[i]], err = err9)
  }
}

for (i in names(ddF)){
  ddF[[i]]$sequence <- gsub("AAAAAA*","AAAAAA",ddF[[i]]$sequence)
  ddF[[i]]$sequence <- gsub("TTTTTT*","TTTTTT",ddF[[i]]$sequence)
  ddF[[i]]$sequence <- gsub("CCCCCC*","CCCCCC",ddF[[i]]$sequence)
  ddF[[i]]$sequence <- gsub("GGGGGG*","GGGGGG",ddF[[i]]$sequence)
}

topx=10
topxseqs <- c()
for (i in dereps){
  topxseqs <- c(topxseqs,names(i$uniques[1:topx]))
}

table(table(topxseqs)) #6 sequences are in the top10 in all samples 

#build dataframe
df <- data.frame()
for(i in names(dereps)){
  df <- rbind(df, data.frame(sequence=names(dereps[[i]]$uniques[1:topx]),
                             Abundance=unname(dereps[[i]]$uniques[1:topx]),
                             ASV=ddF[[i]]$sequence[ddF[[i]]$map[1:topx]],
                             SampleID=i))
}

df$ASV2 <- df$ASV
levels(df$ASV2) <- paste0("ASV-",1:length(levels(df$ASV2)))
df$sequence2 <- df$sequence
levels(df$sequence2) <- paste0("sequence-",1:length(levels(df$sequence2)))
df$ASVrep <- as.character(df$sequence)==as.character(df$ASV)
df$ASVrep2[df$ASVrep] <- as.character(df$ASV2[df$ASVrep])

ggplot(df, aes(x = sequence2, y = Abundance, fill=ASV2)) + 
  geom_bar(stat="identity") + 
  facet_grid(SampleID ~ ., scale="free_y") + 
  theme(axis.text.x = element_text(angle=90, hjust=1))

ggplot(df, aes(x = sequence2, y = Abundance, fill=ASVrep2)) + 
  geom_bar(stat="identity") + 
  facet_grid(SampleID ~ ., scale="free_y") + 
  theme(axis.text.x = element_text(angle=90, hjust=1))

st <- makeSequenceTable(ddF)
st2 <- st
colnames(st2) <- gsub("AAAAAA*","AAAAAA",colnames(st2))
colnames(st2) <- gsub("TTTTTT*","TTTTTT",colnames(st2))
colnames(st2) <- gsub("CCCCCC*","CCCCCC",colnames(st2))
colnames(st2) <- gsub("GGGGGG*","GGGGGG",colnames(st2))

st2 <- collapseNoMismatch(st2)

tt <- cbind(make.unique(substr(colnames(st),1,10)),colnames(st))
rownames(tt) <- tt[,2]

ps <- phyloseq(otu_table(st, taxa_are_rows = F), 
               sample_data(data.frame(row.names=names(ddF), 
                                      Sample_Namex=names(ddF),
                                      Lib=substr(names(ddF),1,2),
                                      Polymerase=c("T","P","T","P","T","P"))),
               tax_table(tt)
)

ps.rare <- rarefy_even_depth(ps)
ord <- ordinate(ps.rare, method = "PCoA", distance = "bray")
plot_ordination(ps, ord, label = "Sample_Namex", color="Polymerase")
ps.temp <- prune_taxa(colnames(ps.rare@otu_table) %in% colnames(ps.rare@otu_table)[1:10], ps.rare)
plot_bar(ps.temp, fill = "ta1")

ps.rare <- rarefy_even_depth(ps)
ord <- ordinate(ps.rare, method = "PCoA", distance = "bray")
plot_ordination(ps, ord, label = "Sample_Namex", color="Polymerase")

tt <- cbind(make.unique(substr(colnames(st2),1,10)),colnames(st2))
rownames(tt) <- tt[,2]

ps <- phyloseq(otu_table(st2, taxa_are_rows = F), 
               sample_data(data.frame(row.names=names(ddF), 
                                      Sample_Namex=names(ddF),
                                      Lib=substr(names(ddF),1,2),
                                      Polymerase=c("T","P","T","P","T","P"))),
               tax_table(tt)
               )
ps.rare <- rarefy_even_depth(ps)
ord <- ordinate(ps.rare, method = "PCoA", distance = "bray")
plot_ordination(ps, ord, label = "Sample_Namex", color="Polymerase")
ps.temp <- prune_taxa(colnames(ps.rare@otu_table) %in% colnames(ps.rare@otu_table)[1:10], ps.rare)
plot_bar(ps.temp, fill = "ta1")

#### where to collapse the homopolyer and nomismatch in the workflow? ####

fqs <- list.files("./", pattern = "I.*fastq.gz")

dereps <- list()
for (file in fqs){
  dereps[[file]] <- derepFastq(file)
}
names(dereps) <- gsub("_.*","",names(dereps))

err8 <- readRDS("ITS_0008.R1.RDS")
err9 <- readRDS("ITS_0009.R1.RDS")

ddF <- list()
for (i in names(dereps)){
  if (substr(i,2,2)==8){
    ddF[[i]] <- dada(derep = dereps[[i]], err = err8)
  }else{
    ddF[[i]] <- dada(derep = dereps[[i]], err = err9)
  }
}

for (i in names(ddF)){
  ddF[[i]]$sequence <- gsub("AAAAAA*","AAAAAA",ddF[[i]]$sequence)
  ddF[[i]]$sequence <- gsub("TTTTTT*","TTTTTT",ddF[[i]]$sequence)
  ddF[[i]]$sequence <- gsub("CCCCCC*","CCCCCC",ddF[[i]]$sequence)
  ddF[[i]]$sequence <- gsub("GGGGGG*","GGGGGG",ddF[[i]]$sequence)
}

table(ddF[[1]]$map)

aggregate(ddF[[i]]$denoised, by=list(ddF[[i]]$sequence), FUN=sum)
aggregate(ddF[[i]]$denoised, by=list(ddF[[i]]$sequence), FUN=sum)
aggregate(ddF[[i]]$denoised, by=list(ddF[[i]]$sequence), FUN=sum)
aggregate(ddF[[i]]$denoised, by=list(ddF[[i]]$sequence), FUN=sum)
aggregate(ddF[[i]]$denoised, by=list(ddF[[i]]$sequence), FUN=sum)

x <- c()
for (i in names(dereps)){
x[i] <- length(dereps[[i]]$uniques)/sum(dereps[[i]]$uniques)
}

ps@sam_data$seqcom <- 1-x
sequecne()