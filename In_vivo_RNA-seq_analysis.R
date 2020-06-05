##R source code for in vivo RNA-seq analysis - In vivo transcriptome analysis provides insights into host-dependent expression of virulence factors by Yersinia entomophaga MH96, during infection of Galleria mellonella. Amber R. Paulson, Maureen O'Callagha, Xue-Xian Zhang, Paul Rainey and Mark R. H. Hurst. Pre-print 2020.
##All raw sequence data, raw count data and processed count data are made available on NCBI's Gene Expression Omnibus platform under "GSE142509_rockhopper_rawcounts.txt.gz" (https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE142509&format=file&file=GSE142509%5Frockhopper%5Frawcounts%2Etxt%2Egz)

##load libraries
library("GenomicFeatures")
library("stringi")
library("EDASeq")
library("RUVSeq")
library("lattice")
library("edgeR")
library("limma")
library("statmod")
library("EDASeq")
library("RColorBrewer")

##assign color from RColorBrewer
colors <- brewer.pal(12, "Paired")

##read in raw count data, row names are locus tags
seqdata <- read.delim("count_data.txt", header=TRUE, row.names=1)

##filter non-expressed genes (requires at least 5 reads in at least two samples to be included)
filter <- apply(seqdata, 1, function(x) length(x[x>5])>=2)
seqdata_filtered <- seqdata[filter, ]

##select all filtered rows
sel.rows <- row.names(seqdata_filtered) 

##assign treatment/collection batch as factors
myx <- as.factor(c("10^7 LB","10^7 LB","10^8 LB","10^8 LB","10^9 LB","10^9 LB","37 LB","37 LB","10^7 Galleria.A","10^7 Galleria.A","10^7 Galleria.A","10^7 Galleria.B","10^8 Galleria","10^8 Galleria","10^8 Galleria","10^8 Galleria","10^9 Galleria.A","10^9 Galleria.A","10^9 Galleria.A","10^9 Galleria.B","37 Galleria","37 Galleria","37 Galleria","37 Galleria"))
myset <- newSeqExpressionSet(as.matrix(seqdata_filtered),phenoData=data.frame(myx,row.names=colnames(seqdata_filtered)))

##apply upper quartile normalization
UQmyset <-betweenLaneNormalization(myset, which="upper")

##explore RLE plots of unnormalized and normalized data
plotRLE(myset, outline=FALSE,ylim=c(4,-4),col=colors[myx])
plotRLE(UQmyset, outline=FALSE,ylim=c(4,-4),col=colors[myx])

##explore PCA plots to look for outliers
plotPCA(myset, k=2, labels=FALSE, col=colors[myx])
plotPCA(UQmyset, k=2, labels=FALSE, col=colors[myx])

##subset to remove Late.3 outlier sample
myfilteredx <- seqdata_filtered[ , -which(names(seqdata_filtered) %in% c("Late_Galleria3"))]

##read in gene annotation information, locus tag is row names 
genes <- read.delim("annotations.txt",header=TRUE, row.names=1)
genes <- genes[which(row.names(genes) %in% sel.rows), ] ##Selects the filtered row names only

##assign treatment groups as factors
groups <- factor(c("early_LB","early_LB","mid_LB","mid_LB","late_LB","late_LB","X37_LB","X37_LB","early_host","early_host","early_host","early_host","mid_host","mid_host","mid_host","mid_host","late_host","late_host","late_host","X37_host","X37_host","X37_host","X37_host"),levels=c("early_host","mid_host","late_host","early_LB","mid_LB","late_LB","X37_host","X37_LB"))

##generate DGEList object and set design matrix based on treatment groups
dge <- DGEList(counts=myfilteredx, group=groups,genes=genes)
dge <- calcNormFactors(dge, method="upperquartile")
groups <- dge$samples$group
design <- model.matrix(~0+groups)

##apply voom mean-variance model
myvm <- voom(dge, design=design, plot=TRUE)

##print voom transformed LogCPM count data
write.table(myvm$E,file="LogCPM_count_data.txt")

#convert matrix to dataframe
myvm.df <- as.data.frame(myvm$E)

##fit linear models
fit <- lmFit(myvm,design)

##smooth standard errors with Empirical Bayes 
fit.eBayes <- eBayes(fit)

##specify host-specific contrast matrix
host_specific.contrast.matrix=makeContrasts(groupsearly_host-groupsearly_LB,groupsmid_host-groupsmid_LB, groupslate_host-groupslate_LB,levels=design)
colnames(host_specific.contrast.matrix)=c("Early_hostvsLB","Middle_hostvsLB","Late_hostvsLB")

##fit the host-specific contrast matrix
fit_host <- contrasts.fit(fit,host_specific.contrast.matrix)
fit_host <- eBayes(fit_host)

##identify DE genes
results_host <- decideTests(fit_host,adjust.method = "BY",p.value=0.05)

##print results
write.fit(fit_host, results_host,"results.txt")
vennDiagram(results_host)

##print all DE genes
AllDE.table <- topTable(fit_host, sort="non", n=Inf) 
write.table(AllDE.table, "allDE_table.txt")

##order significant results by F value
modFpvalue=fit_host$F.p.value
indx = p.adjust(modFpvalue, method="BY")
sig = modFpvalue[indx]
nsiggenes = length(sig)
modF = fit_host$F
modFordered = order(modF, decreasing = TRUE)

##order the genes in the dge object based on F value
ranked_genes = dge$genes$GeneID[modFordered[1:nsiggenes]]
clust.table <- topTable(fit_host, sort="non",n=Inf)
ranked <- clust.table[which(row.names(clust.table) %in% ranked_genes), ]

##select only Log2CPM value columns and re-name
keeps.clust <- c("Early_hostvsLB","Middle_hostvsLB","Late_hostvsLB")
clust.log2FC <- clust[ , keeps.clust, drop=FALSE]
colnames(clust.log2FC) = c("Early","Middle","Late")

##convert object to matrix
clust.log2FC <- data.matrix(clust.log2FC, rownames.force=NA)

##convert matrix to data frame
clust.log2FC.df <- as.data.frame(clust.log2FC, rownames.force=NA)
clust.log2FC.df.ordered <- clust.log2FC.df[order(row.names(clust.log2FC.df)),]

##combine annotation information
sel.rows.b <- row.names(clust.log2FC) 
genes.b <- genes[which(row.names(genes) %in% sel.rows.b), ] 
genes.b <- genes.b[order(row.names(genes.b)),]
genes_annotated_df <- new("AnnotatedDataFrame",data=genes.b)
clust.log2FC.matrix.ordered <- data.matrix(clust.log2FC.df.ordered,rownames.force = NA)

##generate an ExpressionSet using Biobase
library(Biobase)
DEGesetlog2FC <- ExpressionSet(assayData=clust.log2FC.matrix.ordered,featureData=genes_annotated_df)

##cluster the data using fuzzy clustering
library(Mfuzz)
fuzz <- mfuzz(DEGesetlog2FC, c = 20, m = 1.5)
mfuzz.plot2(DEGesetlog2FC,centre=TRUE,cl=fuzz)
overlap <- overlap(fuzz)
overlap.plot(fuzz, over=overlap,thres=0.3)
