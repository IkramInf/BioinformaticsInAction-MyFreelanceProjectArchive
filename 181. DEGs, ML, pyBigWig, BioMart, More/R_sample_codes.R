### DESeq2
library(DESeq2)
library(airway)

data(airway) 
airway$dex <- relevel(airway$dex, ref="untrt")
dds <- DESeqDataSet(airway, design=~cell + dex)
dds <- DESeq(dds)
results <- results(dds)
summary(results)


### edgeR
library(edgeR)
library(airway)

data(airway)
counts <- assay(airway)
condition <- factor(colData(airway)$dex)
y <- DGEList(counts=counts, group=condition)
y <- calcNormFactors(y)
design <- model.matrix(~condition)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
topTags(lrt)


### GenomicRanges
library(GenomicRanges)
library(GenomicFeatures)

txdb <- makeTxDbFromUCSC(genome="hg19", tablename="knownGene")
genes <- genes(txdb)
query <- GRanges(seqnames="chr1", ranges=IRanges(start=1e6, end=2e6))
overlap <- findOverlaps(query, genes)
overlap


### AnnotationHub
library(AnnotationHub)

ah <- AnnotationHub()
gencode <- query(ah, "GENCODE")
human_annotations <- gencode[["AH75112"]]
human_annotations


### SummarizedExperiment
library(SummarizedExperiment)
library(airway)

data(airway)
assay_data <- assay(airway)
rowData <- rowData(airway)
colData <- colData(airway)
airway


### dplyr
library(dplyr)

filtered <- mtcars %>% filter(mpg > 20) %>% mutate(weight_kg=wt * 453.6)
filtered


### tidyr
library(tidyr)

data(billboard)
tidy_data <- billboard %>% pivot_longer(cols=starts_with("wk"), names_to="week", values_to="rank")
tidy_data


### ggplot2
library(ggplot2)

ggplot(iris, aes(x=Sepal.Length, y=Petal.Length, color=Species)) +
  geom_point() +
  theme_minimal()
  

### Lattice
library(lattice)

xyplot(mpg ~ wt | factor(cyl), data=mtcars, type="b", layout=c(1, 3))


### Plotly for R
library(plotly)

plot_ly(mtcars, x=~wt, y=~mpg, z=~disp, color=~factor(cyl), type="scatter3d", mode="markers")


### caret
library(caret)

model <- train(Species~., data=iris, method="rf", trControl=trainControl(method="cv", number=5))
model


### lme4
library(lme4)

model = lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
summary(model)
