library(dada2);packageVersion("dada2")
library(phyloseq);packageVersion("phyloseq")
library(ggplot2);packageVersion("ggplot2")

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
### try to resolve this issue with priors ####

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

### try to resolve this issue with 1 specific prior ####

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
