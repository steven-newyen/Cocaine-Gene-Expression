# Names: Steven Nguyen
#
# Title: Identifying differentially expressed genes between cocaine addict deaths and, non cocaine addict deaths

# Importing libraries
library(affy)
library(dplyr)
library(ggplot2)
library(GEOquery)
library(limma)

# 1 & 2: Importing the data from directory
GSE54839 <- getGEO("GSE54839")

GSE54839.expr <- exprs(GSE54839[[1]])
GSE54839.p <- pData(GSE54839[[1]])


# 3: Generating a boxplot to see if processed data is already in log2 scale
GSE54839.exprTrunc <- head(GSE54839.expr, 10)
boxplot(GSE54839.exprTrunc, main = "processed data", xlab = "Sample Name", ylab = "Expression of Samples")

  # Since not in log2 scale, normalizing data with log2 scale
GSE54839.exprLog <- log2(GSE54839.expr)
GSE54839.exprLogTrunc <- head(GSE54839.exprLog, 10)
boxplot(GSE54839.exprLogTrunc, main = "log2 processed data", xlab = "Sample Name", ylab = "Normalized Expression of Samples")

# 4: Finding the number of probes and samples
ncol(GSE54839.expr)
nrow(GSE54839.expr)
  # There are 60 samples and 48761 probes in the dataset

# 5: Extracting the column that contains the categories that we would like to compare
addiction.status <- as.character(GSE54839.p$`disease state:ch1`)
table(addiction.status)
  # There are 30 individuals that were part of the group with cocaine addiction, and there are 30 individuals in the control group.

# 6: Finding the differentially expressed probes across the 2 groups

design <- model.matrix(~-1+addiction.status)
colnames(design) <- c("control", "cocaine.addiction")

fit <- lmFit(GSE54839.exprLog, design)
contrast.matrix <- makeContrasts(control - cocaine.addiction,levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 = eBayes(fit2)

tt20 <- topTable(fit2,sort.by = "p", p.value = 0.20, number = nrow(GSE54839.exprLog))

# 7: Constructing a boxplot for the probe with the lowest adjusted p-value

m <- match(rownames(tt20)[1], rownames(GSE54839.exprLog))
probe <- GSE54839.exprLog[m,]
df <- data.frame(expr = probe, addiction.status = addiction.status)

means <- df %>% group_by(addiction.status) %>% summarize(mean = mean(expr))
means

diff(means$mean)
tt20$logFC[1]

# convert to FC #
logFC <- tt20$logFC[1]
2**logFC

## visualize ##
FC <- paste0("FC = ", round(2**logFC, 2))
main <- paste0("Expression of ", rownames(tt20)[1], ", ", FC)

ggplot(df, aes(x = addiction.status, y = expr, fill = addiction.status)) + geom_boxplot() +
  ylab("log2 expression") + xlab("Addiction status") + ggtitle(main) +
  scale_fill_manual(values = c("red", "yellow")) +
  theme_classic() + theme(legend.position = "none")


# 8: Using R to output the annotation of the data
platform2 <- annotation(GSE54839[[1]])

# 9: Downloading the platform (GPL) for this data
pl2 <- getGEO(platform2)
pl2 <- Table(pl2)

probe <- rownames(tt20)[1]
m <- match(probe, pl2$ID)
pl2$`ILMN_Gene`[m]

# 10: Finding the gene names corresponding to all probes, and creating a table

match <- match(rownames(tt20), rownames(GSE54839.exprLog))
ILMNGene <- pl2$`ILMN_Gene`[match]

df20 <- data.frame(tt20)
df220 <- df20 %>% select(logFC,adj.P.Val)
df220["ILMN_Gene"] <- ILMNGene

head(df220, 5)

# 11: Using DAVID to identify Gene Ontology and KEGG pathways that are associated with the differentially expressed genes identified in (6) and (10)

probe <- rownames(tt20)
m <- match(probe, pl2$ID) 
genes<-pl2$ILMN_Gene[m] 
keep <- genes!=""
genes <- genes[keep]
genes <- strsplit(genes, " /// ")
genes <- unlist(genes)
genes <- unique(genes) 
write.table(genes, row.names = FALSE, quote = FALSE)

# 12:

# We analzyed the GSE54839 dataset, where this dataset compares 2 groups, one group being the control (non-addicted people) and the other group of individuals that
# have a cocaine addiciton. In this dataset, there are 60 samples and 48761 probes included. Each group was evenly split, with 30 samples in each group. With an FDR of 20%,
# we found there to be 1673 differentially expressed probes in GSE54839, and the names of the top 3 genes are RGL4, ITCB1BP1, and COX4I2 
