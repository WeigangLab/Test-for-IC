###### 04/03/2024 Test for IC post-doc
rm(list=ls())
# install.packages("mlbench")
library(mlbench)
library(tidyr)

# Set random seeds to ensure repeatable results
set.seed(123)

##### construct the data set
# Construct gene expression data set -> A-->WT;B--->KO
gene_dat <- data.frame(
                  A_rep1 = sample(0:1000,10000,replace = T),
                  A_rep2 = sample(0:1000,10000,replace = T),
                  A_rep3 = sample(0:1000,10000,replace = T),
                  B_rep1 = sample(0:1000,10000,replace = T),
                  B_rep2 = sample(0:1000,10000,replace = T),
                  B_rep3 = sample(0:1000,10000,replace = T))
row.names(gene_dat) <- paste("gene_",1:10000)

# Create a breast cancer dataset
breast_cancer_data <- na.omit(as.data.frame(lapply(BreastCancer,as.numeric))) 

# Create a individual dataset
individual_data <- data.frame(
  Gender = rep(c("A","B"), each = 25),
  Age = sample(20:60, 50, replace = T),
  Smoking = sample(c(0,1), 50, replace = T),
  Drinking = sample(c(0,1), 50, replace = T),
  Tumor = sample(c(0,1), 50, replace = T))

##### Differential expression analysis
library(sva)
library(pheatmap)
library(DESeq2)
library(ggplot2)
library(ggrepel)

# remove batch effect
batch <- c(rep(1:3,times = 2))
csif <- data.frame(group_name = colnames(gene_dat),batch = batch)
modcombat <- model.matrix(~1, data = csif)
exp <- ComBat(dat=gene_dat, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
head(exp)

matrix <- cor(exp,method="pearson") # pearson
pheatmap(matrix,legend = T,
         show_rownames = T, 
         show_colnames = T,
         cluster_rows =T,
         cluster_cols = T,
         clustering_method="complete",
         border_color = "gray",
         display_numbers=T,
         angle_col = 90,
         number_format = "%.3f",
         cellwidth=35,
         cellheight=35,
         number_color = "black",
         fontsize_number = 12,
         fontsize_row = 12,
         fontsize_col = 12)

# sample clustering
group_list <- factor(rep(c("WT","KO"),each =3),levels = c("WT","KO"))
colData <- data.frame(row.names = colnames(gene_dat),group_list = group_list)
head(colData)

dds <- DESeqDataSetFromMatrix(countData = gene_dat, colData=colData, design= ~group_list) # ~group_list+treatment

vsd <- varianceStabilizingTransformation(dds)

plotPCA(vsd,intgroup=c('group_list'))+theme_bw()
pcaData <- plotPCA(vsd, intgroup = c("group_list"), returnData = TRUE)
head(pcaData)

# Polish the graph
ggplot(pcaData,aes(x=PC1,y=PC2,color = group))+geom_point(size=4.5)+ # shape = group
  ylim(-35,35)+xlim(-35,35)+theme_bw()+
  labs(x="PC1:22% variance",y="PC2:21% variance")+
  geom_label_repel(data = pcaData, aes(label = rownames(pcaData)),
                   size = 4.5,box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.8, "lines"),
                   segment.color = "black",
                   show.legend = FALSE, max.overlaps = 10000)+
  theme(legend.position = c(.01,.01),
        legend.justification = c("left", "bottom"),
        legend.background = element_rect(fill=NA),
        legend.text = element_text(colour="black", size=10,face="bold"),
        legend.title = element_text(colour="black", size=10, face="bold"),
        plot.title = element_text(hjust = 0.5,size = 11,face = "bold"),
        axis.text.x=element_text(size=15,angle=0,face ="bold",vjust = 0.5,hjust = 0.5,colour="black"),
        axis.text.y=element_text(size=15,face ="bold",color="black"),
        axis.title=element_text(size=12,face ="bold"),
        axis.line.x=element_line(linetype=1,color="black",linewidth = 0),
        axis.line.y=element_line(linetype=1,color="black",linewidth = 0))

# differential analysis
exprSet <- gene_dat
exprSet <- dat[rowSums(dat) > 10,]
exprSet <- exprSet[which(rowSums(exprSet == 0) < 3),] # Genes with zero read count in more than three conditions were filtered out.

group_list <- factor(rep(c("WT", "KO"), each = 3),levels = c("WT","KO"))
colData <- data.frame(row.names = colnames(exprSet),group_list = group_list)
head(colData)

dds <- DESeqDataSetFromMatrix(countData = as.matrix(exprSet),colData = colData,design = ~group_list)

dds2 <- DESeq(dds)
res <- results(dds2)   # The differential analysis results are extracted
DEG <- res[order(res$padj),]
DEG <- as.data.frame(DEG)
DEG <- na.omit(DEG)
head(DEG)

# cutoff 
table(DEG$pvalue < 0.05 & (DEG$log2FoldChange > 0.58 | DEG$log2FoldChange < -0.58)) 

resdata <-  na.omit(merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE))
head(resdata)

DEGdata <- subset(resdata,pvalue < 0.05 & (log2FoldChange > 0.58 | log2FoldChange < -0.58))
head(DEGdata)
table(DEGdata$log2FoldChange>0)

# normalized count DEGs heatmap 
p <- t(scale(t(DEGdata[,8:13])))
pheatmap(p)
p[p > 1.5] <- 1.5
p[p < -1.5] <- -1.5

annotation_col <- data.frame(Replicate=c(rep(c("rep1","rep2","rep3"),times=2)),
                            CellType=c(rep("A",3),rep("B",3)))
rownames(annotation_col) <- colnames(p)

ann_colors <- list(CellType = c(A="blue",B="red")),
                  Regulation = c(Up="#F54EA2",Down="#00a03e")

# Polish the heatmap
pheatmap(p,legend = T,
         show_rownames = F, 
         show_colnames = T,
         angle_col=90,#cutree_row = 2,cutree_col = 2
         annotation_col=annotation_col,
         annotation_colors = ann_colors,
         annotation_legend=T,
         cluster_rows = T,
         cluster_cols = T,
         fontsize_col = 10,
         cellwidth=30,
         cellheight=0.5,
         color = colorRampPalette(c("#003399","#F7F7F7","#FF9933"))(100))

# DEGs vocano plot
table(DEGdata$log2FoldChange>0)
resdata$label <- ifelse(resdata$padj !=0 & resdata$pvalue < 0.0001 & abs(resdata$log2FoldChange) >= 1,resdata$Row.names,"")
resdata$regulation <- "Not"
resdata[which(resdata$log2FoldChange > 0.58),]$regulation <- "Up"
resdata[which(resdata$log2FoldChange < -0.58),]$regulation <- "Down"

ggplot(resdata,aes(x = log2FoldChange, y = -log10(pvalue), colour = regulation)) +  # padj;pvalue
  xlab("log2 FoldChange")+ylab("-log10 P-value")+ # adjusted
  geom_point(size = 1.5, alpha = 0.8) + 
  scale_color_manual(values = c("dodgerblue","gray50","orange")) + 
  theme_bw()+ xlim(-6,6)+ylim(0,12)+
  geom_hline(yintercept=-log10(0.05),linetype=2)+
  geom_vline(xintercept=c(-0.58,0.58),linetype=2)+
  theme(plot.title = element_text(hjust = 0.5,size = 12,face = "bold"),
        legend.text = element_text(colour="black", size=10,face="bold"),
        legend.title = element_text(colour="black", size=10, face="bold"),
        legend.position = "right",
        axis.text.x=element_text(size=11,angle=0,face ="bold",vjust = 0.5,hjust = 0.5,color="black"),
        axis.text.y=element_text(size=11,face ="bold",color="black"),
        axis.title.x=element_text(vjust=2, size=11,face = "bold",colour="black"),
        axis.title.y=element_text(vjust=2, size=11,face = "bold",colour="black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.x=element_line(linetype=1,color="black",linewidth=0),
        axis.line.y=element_line(linetype=1,color="black",linewidth=0))+
  annotate("text",x=-4.5,y=12,label="Down = 248",fontface="bold")+ 
  annotate("text",x=4.5,y=12,label="Up = 213",fontface="bold")+
  geom_label_repel(data = resdata, aes(label = label),
                   size = 3,
                   min.segment.length = 0, seed = 42,
                   box.padding = unit(0.5,"lines"),
                   point.padding = unit(0.5, "lines"),
                   segment.color = "black",
                   show.legend = FALSE, max.overlaps = 10000)



##### Model and predict whether a given sample (row of data set) is benign or malignan, based on 9 other cellular characteristics
dim(breast_cancer_data)
breast_cancer_data$type <- factor(breast_cancer_data$Class, levels=c(1,2), labels = c("A","B"))
train <- breast_cancer_data[1:400,c(-1,-11)]
test <- breast_cancer_data[401:nrow(breast_cancer_data),c(-1,-11)]
# logistic regression
cancer_glm <- glm(type~., data=train, family = binomial) # A logistic regression model is constructed using all the features
summary(cancer_glm)

cancer_glm_step <- step(cancer_glm)

# Then try to reduce the number of features from the model in last step and make another regression model.(use 0.1 as significant level)
summary(cancer_glm_step)

anova(cancer_glm, cancer_glm_step, test = "Chisq")
# predict
test$prob <- predict(cancer_glm_step, newdata = test, type = "response")
# he predict function return a value between 0 and 1, < 0.5 belong to the benign, >= 0.5 belong to the malignan
test$pred_type <- ifelse(test$prob < 0.5,"A","B") 
test$pred_type <- as.factor(test$pred_type)
accuracy <- mean(test$pred_type == test$type)
sprintf("accuracy = %.4f", accuracy)
# modest model accuracy



##### logistic regression
### The logistic regression models are used to calculate the model of age, smoking and drinking for cancer occurrence
dim(individual_data)
indi_glm <- glm(Tumor~., data=individual_data, family = binomial)
summary(indi_glm)

indi_glm_step <- step(indi_glm)

summary(indi_glm_step)

### So, age can affect the development of cancer (Limited to this test data)
