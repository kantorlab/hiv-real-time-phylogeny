library(tidyverse)
library(ggplot2)
library(ggtree)
library(ape)
library(treeio)

args <- commandArgs(trailingOnly=TRUE)
study_id       <- args[1]
tree_file      <- args[2]
csv_file       <- args[3]
treatment_file <- args[4]
dataset        <- args[5]
pdf_file       <- args[6]
tips_file      <- args[7]

dataset <- as.integer(str_extract(dataset, "\\d+"))

treatment <- read_tsv(treatment_file) %>%
             mutate(Dataset=as.integer(str_extract(Dataset, "\\d+"))) %>%
             select(StudyID, Dataset)

metadata <- read_csv(csv_file) %>%
            left_join(treatment, by = c("StudyID")) %>%
            mutate(Index=as.factor(case_when(StudyID == study_id | Dataset == dataset | ThisMonth ~ "Current",
					     Dataset < dataset ~ "Previous",
					     TRUE ~ "")))
print(metadata)

tree <- read.raxml(tree_file)
tips <- tree@phylo$tip.label
droplist <- tips[!(tips %in% metadata$SequenceID)]
tree <- drop.tip(tree, droplist)
tips <- tree@phylo$tip.label
tree@phylo$tip.label <- sapply(tips, function(tip) { gsub("_.*", "", tip) })
print(tree@phylo$tip.label)
color <- metadata$Index
names(color) <- metadata$StudyID
print(color)

pdf(pdf_file, width=6, height=(0.2+nrow(metadata)+1)*0.2)
ggtree(tree, ladderize=FALSE) +
    geom_treescale(width=0.01, fontsize=2.4, x=0) +
    geom_tiplab(aes(color=color[label]), size=3.2, align=TRUE) +
    geom_text2(aes(subset=!isTip, label=bootstrap), hjust=-0.2, size=2.4) +
    scale_colour_manual(values=c("NA"="gray40", "Current"="red", "Previous"="forestgreen"), guide=FALSE) +
    xlim(0, 0.25)
dev.off()

write.csv(tips, file=tips_file, row.names=FALSE)

