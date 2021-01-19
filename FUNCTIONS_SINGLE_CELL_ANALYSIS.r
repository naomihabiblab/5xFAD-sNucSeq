library(Seurat)
library(dplyr)
library(sva)
library(statmod)
library(colorspace)

### Compute the group-wise mean of a dataset.
group.means <- function(counts, groups, fn=mean, use.data.table=F)
{
  counts <- aggregate(t(counts), by=list(groups), FUN=fn)
  rownames(counts) = counts$Group.1
  counts$Group.1 = NULL
  r = t(counts)
  return(r)
}

### Run ComBat batch correction from the SVA package
batch.normalise.comBat <- function(counts, batch.groups, max.val=6)
{
  batch.groups = factor(batch.groups) ## drop zero levels
  batch.id = 1:length(unique(batch.groups))
  names(batch.id) = unique(batch.groups)
  batch.ids = batch.id[batch.groups]
  
  correct.data = ComBat(counts,batch.ids, prior.plots=FALSE, par.prior=TRUE)
  correct.data[correct.data > max.val] = max.val
  as.data.frame(correct.data)
}


### Get variable genes. Code adapted from:
### | Brennecke et al, Accounting for technical noise in single-cell RNA-seq experiments
### | Nature Methods 10, 1093–1095 (2013), doi:10.1038/nmeth.2645
### 	See: https://images.nature.com/original/nature-assets/nmeth/journal/v10/n11/extref/nmeth.2645-S2.pdf
### 	and: http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
get.variable.genes <- function(ed, min.cv2=2, pdf=NULL, width=9, height=8, do.plot=T, p.thresh=0.05)
{
  means <- rowMeans(ed)
  vars <- apply(ed,1,var)
  cv2 <- vars/means^2
  minMeanForFit <- unname( quantile( means[ which( cv2 > min.cv2 ) ], .95 ) )
  useForFit <- means >= minMeanForFit # & spikeins
  info(sprintf("Fitting only the %s genes with mean expression > %s", sum(useForFit), minMeanForFit))
  fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ), cv2[useForFit] )
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"])
  if(do.plot){par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2))}
  xg <- exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
  vfit <- a1/xg + a0
  if(do.plot){lines( log(xg), log(vfit), col="black", lwd=3 )}
  
  df <- ncol(ed) - 1
  # add confidence interval
  if(do.plot){
    lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black")
    lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black")
  }
  afit <- a1/means+a0
  varFitRatio <- vars/(afit*means^2)
  varorder <- order(varFitRatio, decreasing=T)
  oed <- ed[varorder,]
  pval <- pchisq(varFitRatio*df,df=df,lower.tail=F)
  adj.pval <- p.adjust(pval,"fdr")
  r = data.frame(rownames(ed), varFitRatio, pval, adj.pval)
  colnames(r) = c("Gene", "VarianceFitRatio", "p", "p.adj")
  v = r[!is.na(r$p.adj),]
  n.sig = sum(v$p.adj<p.thresh)
  info(sprintf("Found %s variable genes (p<0.05)", n.sig))
  
  # add top 100 genes
  if(do.plot){
    points(log(means[varorder[1:n.sig]]),log(cv2[varorder[1:n.sig]]),col=2)
  }
  r = r[order(r$VarianceFitRatio, decreasing=T), ]
  r$Rank = 1:nrow(r)
  return(r)
}


num_cells_per_group = function(groups, total_cells=NULL, cells_per_group=NULL){
  num_cells = sort(table(groups))
  if(!is.null(cells_per_group)){
    num_cells[num_cells > cells_per_group] = cells_per_group
  } else {
    n = sort(table(groups))
    if(length(n) == 1){
      num_cells = total_cells
      names(num_cells) = names(n)
      u = c(0, cumsum(n)[1:(length(n)-1)])
      i = (total_cells - u)/seq(length(n), 1, -1) < n
      if(sum(i) > 0){
        num_cells[i] = as.integer(ceiling((total_cells - sum(n[!i]))/sum(i)))
      }
    }
  }
  num_cells
}
resample = function(x,...){if(length(x)==1) x else sample(x,...)}



########## 
# Init
readMat <- function(matname,dir,minG = 200, minC =3){
  matrix_dir <- paste0(dir,matname,"/")
  barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
  features.path <- paste0(matrix_dir, "genes.tsv.gz")
  matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
  mat <- readMM(file = matrix.path)
  #dim(mat)
  feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V1
  mat<- mat[,colSums(mat)>minG]
  mat <- mat[rowSums(mat>0)>minC,]
  #dim(mat)
  return(mat)
}

init.object <- function(expression_matrix,name.id,lowT=100,HighT=Inf){
  DatObj <- CreateSeuratObject(raw.data = Matrix(expression_matrix, sparse = T), min.cells = 5, project = name.id)
  mito.genes <- grep(pattern = "^MT", x = rownames(DatObj@raw.data), value = TRUE)
  if (length(mito.genes)==0) {
    mito.genes <- grep(pattern = "^Mt", x = rownames(x = DatObj@raw.data), value = TRUE)
  }
  percent.mito <- Matrix::colSums(DatObj@raw.data[mito.genes, ])/Matrix::colSums(DatObj@raw.data)
  DatObj <- AddMetaData(object = DatObj, metadata = percent.mito, col.name = "percent.mito")
  VlnPlot(object = DatObj, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3,point.size = 0,group.by = "orig.ident",do.sort=FALSE)
  DatObj <- FilterCells(DatObj, subset.names = "nGene", low.thresholds = lowT,  high.thresholds = HighT)
  DatObj <- FilterCells(DatObj, subset.names = "percent.mito", high.thresholds = 0.05)
  return(DatObj)
}

# Seurat V2 step 0
scale.object <- function(DatObj,reg.params=c("nUMI","percent.mito")){
  DatObj <- NormalizeData(DatObj)
  DatObj <- ScaleData(DatObj, vars.to.regress = reg.params, model.use = "linear",do.scale=TRUE,do.center=TRUE)
  #DatObj <- ScaleData(DatObj, vars.to.regress = c("nUMI"),)
  DatObj <- FindVariableGenes(DatObj, x.low.cutoff = 0.1,y.cutoff = 0.8)
  length(DatObj@var.genes)
  return(DatObj)
}

# Seurat V2 step 1
analyze.object <- function(DatObj,mn_i = 1,mx_i=25){
  P=100
  DatObj <- RunPCA(DatObj, pc.genes = DatObj@var.genes, do.print = TRUE, pcs.print = 1:20,pcs.compute = P, genes.print = 10)
  PCElbowPlot(DatObj, num.pc = P)
  DatObj <- ProjectPCA(object = DatObj, do.print = FALSE)
  DatObj <- RunTSNE(DatObj, dims.use = mn_i:mx_i, do.fast = T) #,genes.use = DatObj@var.genes, cells.use = NULL,reduction.use = "pca")
  DatObj <- FindClusters(object = DatObj, reduction.type = "pca", dims.use = mn_i:mx_i, resolution = 0.9, print.output = 0, save.SNN = TRUE)
  return(DatObj)
}

# Seurat V2 step 2
analyze.object1 <- function(DatObj,mn_i = 1,mx_i=25){
  P=100
  DatObj <- RunPCA(DatObj, pc.genes = DatObj@var.genes, do.print = TRUE, pcs.print = 1:20,pcs.compute = P, genes.print = 10)
  PCElbowPlot(DatObj, num.pc = P)
  VizPCA(DatObj, pcs.use = 1:6, num.genes = 20, use.full = FALSE, font.size = 0.5, nCol = NULL, do.balanced = TRUE)
  return(DatObj)
}

# Seurat V2 step 3
analyze.object2 <- function(DatObj,mn_i = 1,mx_i=25){
  DatObj <- RunTSNE(DatObj, dims.use = mn_i:mx_i, do.fast = T) #,genes.use = DatObj@var.genes, cells.use = NULL,reduction.use = "pca")
  TSNEPlot(DatObj, group.by = "orig.ident",pt.size = 0.5,do.label = F)
  DatObj <- FindClusters(object = DatObj, reduction.type = "pca", dims.use = mn_i:mx_i, resolution = 0.9, print.output = 0, save.SNN = TRUE)
  TSNEPlot(DatObj,pt.size = 0.5,do.label = T)
  return(DatObj)
}

# Seurat V3
NormAndClustV3 <- function(obj, maxPC = 20, nVar = 2000,ClusRes=0.8) {
  obj <- NormalizeData(object = obj,assay = "RNA",normalization.method = "LogNormalize")
  obj <- FindVariableFeatures(object = obj,assay = "RNA",nfeatures = nVar)
  obj <- ScaleData(object = obj, assay = "RNA",features = VariableFeatures(object = obj))
  obj <- RunPCA(object = obj,assay = "RNA",features = VariableFeatures(object = obj),npcs = 50,ndims.print = 6)
  ElbowPlot(object = obj,ndims = 50,reduction = "pca")
  obj <- RunTSNE(object = obj,reduction = "pca",dims = c(1:maxPC),features = VariableFeatures(obj),tsne.method = "Rtsne")
  DimPlot(object = obj,reduction = "tsne")
  obj <- FindNeighbors(object = obj,reduction = "pca",dims = c(1:maxPC),features = VariableFeatures(obj),k.param = 40)
  obj <- FindClusters(object = obj,resolution = ClusRes)
}



#############
# Pathways from MsigDB

Pathways.list <- list()
for (fname in c("h.all.v7.0.symbols.gmt.txt","c2.cp.kegg.v7.0.symbols.gmt.txt","c2.cp.pid.v7.0.symbols.gmt.txt","c2.cp.reactome.v7.0.symbols.gmt.txt","c2.cp.biocarta.v7.0.symbols.gmt.txt","c5.bp.v7.0.symbols.gmt.txt"))
{
  f <- paste0("~/Documents/Mouse_AD/Differential_genes_mASC/Pathways/",fname)
  inlist <- strsplit(readLines(f), "[[:space:]]+")
  pathways <- lapply(inlist, tail, n = -2)
  names(pathways) <- lapply(inlist, head, n = 1)
  Pathways.list <- append(Pathways.list ,pathways)
}
library(bc3net)
tab.hypg=enrichment(candidate, reference, pathways, verbose=TRUE)



#############
# Plots


# Volcano Plots
volc.plot <- function(m,tname,ngene1=10,ngene2=10) {
  eps <- 1e-300
  volc = ggplot(m,aes(x=avg_logFC,y=-log(p_val_adj+eps)))+geom_point(size = 0.5)+ggtitle(label = "",subtitle = tname)
  m<-cbind(m,gene<-rownames(m))
  colnames(m) <- c(colnames(m)[1:length(colnames(m))-1],"gene")
  P <- volc+geom_text_repel(data=m[m$avg_logFC>0,][1:ngene1,], aes(label=gene),size=3,color="steelblue4")+geom_text_repel(data=m[m$avg_logFC<0,][1:ngene2,], aes(label=gene),size=3,color="plum4")
  return(P)
}

# Volcano Plots
for (id in c(1:11)) {
  i1 <- which(Markers$cluster==id & Markers$avg_logFC>0 & Markers$pct.2<0.3 )
  i2 <- which(Markers$cluster==id & Markers$avg_logFC<0 & Markers$pct.1<0.3 )
  i <- union(i1,i2)
  volc = ggplot(Markers[i,],aes(x=avg_logFC,y=-log(p_val_adj)))+geom_point(size = 0.1)+ggtitle(paste("cluster", id))
  plot6 <- volc+geom_text_repel(data=Markers[union(i1[1:8],i2[1:4]),], aes(label=gene))
  dev.copy(pdf,paste(Dpath, nameA,'.Volcano.cluster.All',id,'.pdf' , sep = "") ); dev.off()
  ii <- which(Markers$cluster==id)
  volc = ggplot(Markers[ii,],aes(x=avg_logFC,y=-log(p_val_adj)))+geom_point(size = 0.1)+ggtitle(paste("(negbin) cluster", id))
  volc+geom_text_repel(data=Markers[ii[1:20],], aes(label=gene))
  dev.copy(pdf,paste(Dpath, nameA,'.Volcano.negbin.cluster',id,'.pdf' , sep = "") ); dev.off()
}


# Dot Plots
plot.box <- function(dat,name){
  qplot( x=Condition , y=Data , data=dat, geom=c("boxplot","jitter") , fill=Condition,xlab = paste("cluster",name),ylab = "Frequency")
}

# Dot Plots
dot.plot <- function(data,genes.of.interest,group.vector,order.genes=TRUE,group.levels=NULL,do.return.gene.order=F){
  
  alpha <- function(x) {mean(x > 0)}
  exp.median <- function(x) {median(x[x > 0])}
  
  genes.of.interest <- genes.of.interest[genes.of.interest %in% rownames(data)]
  data.summarized <- merge(melt(data.frame(data.frame(t(data[genes.of.interest,]),groups = as.factor(group.vector)) %>% 
                                             group_by(groups) %>% 
                                             summarise_all(funs(alpha)))),
                           melt(data.frame(data.frame(t(data[genes.of.interest,]),groups = as.factor(group.vector)) %>% 
                                             group_by(groups) %>% 
                                             summarise_all(funs(exp.median)))),by=c('groups','variable'),suffixes=c('.alpha','.exp.median'))
  data.summarized$variable <- make.names(data.summarized$variable)
  # genes.of.interest <- make.names(genes.of.interest)

  if (order.genes){
    hclustering <- hclust(as.dist(1-cor(t(data[genes.of.interest,]))),method='ward.D2')
    data.summarized$variable <- factor(data.summarized$variable,levels=make.names(genes.of.interest[hclustering$order]))
  } else {
    data.summarized$variable <- factor(data.summarized$variable,levels=make.names(genes.of.interest))
  }
  data.summarized <- data.summarized[order(data.summarized$variable),]
  
  if (!is.null(group.levels)){
    data.summarized$groups <- factor(data.summarized$groups,levels=group.levels)
  }
  
  P <- ggplot(data.summarized,aes(groups,variable)) + 
    geom_point(aes(size=value.alpha,color=value.exp.median),stroke=0) + 
    scale_color_viridis() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(size='% cells expressed',color='median x > 0') + coord_flip()
  
  if (do.return.gene.order){return(genes.of.interest[hclustering$order])}
  return(P)
}


# BoxPlots
colors4 <- c("steelblue4","steelblue3","steelblue1","skyblue2","turquoise2","palevioletred1","pink2","plum1","plum3","mediumpurple1","mediumpurple3")
category <- 'Condition'
T <- summary(as.factor(Obj@meta.data[,category]))
N <- names(T)
C <- length(summary(as.factor(Obj@ident)))
s <- matrix(,C[1],length(N))
for (i in 1:length(N)) { s[,i] <- t(summary(Obj@ident[which(Obj@meta.data[,category]==N[i])])) }
d <- mapply(`/`, data.frame(s), colSums(s))
colnames(d) <- N
rownames(d) <- names(summary(as.factor(Obj@ident)))
barplot(t(d), beside=TRUE,  xlab='Clusters', ylab='Fraction Nuclei',main="Nuclei per ident",col=colors4)
legend("topright", colnames(d), cex=0.8, pch=1, pt.cex = 1,col=colors4)
dev.copy(pdf,paste(Dpath, nameA,'.BarPlot.By',category,'.pdf' , sep = "") ); dev.off()
#colnames(d) <- stringr::str_replace(colnames(d),".*WT.*","Wt")
#colnames(d) <- stringr::str_replace(colnames(d),".*Wt.*","Wt")
#colnames(d) <- stringr::str_replace(colnames(d),".*AD.*","AD")
#colnames(d) <- stringr::str_replace(colnames(d),".*Untreated.*","AD")
for (id in c(1:11)) {
  dd <- data.frame(Data<-d[id,],Condition<-colnames(d))
  qplot( x=Condition , y=Data , data=dd, geom=c("boxplot","jitter") , fill=Condition,xlab = paste("cluster",rownames(d)[id]),ylab = "Frequency")
  dev.copy(pdf,paste(Dpath, nameA,'.Boxplot.cluster-',id,'.pdf' , sep = "") ); dev.off()
}




#############
# Add meta data

ADmetadatAdd <- function(DatObj){
  
  oi <- DatObj@meta.data$orig.ident
  d <- DatObj@meta.data$orig.ident
  Condition = c("Wt","AD","Untreated")
  for (pt in Condition) {d[grep(pt, oi)] = pt}
  d[grep("Untreated", oi)] = "AD"
  d[grep("WT", d)] = "Wt"
  DatObj <- AddMetaData(DatObj,metadata = d,col.name = "Condition")
  DatObj@meta.data[,'Condition'] <- d
  
  oi <- DatObj@meta.data$orig.ident
  d <- DatObj@meta.data$orig.ident
  d[] <- "7m"
  d[grep("13m", oi)] = "13m"
  d[grep("14m", oi)] = "14m"
  d[grep("8w", oi)] = "8w"
  d[grep("4m", oi)] = "4m"
  d[grep("10m", oi)] = "10m"
  d[grep("6w", oi)] = "6w"
  DatObj <- AddMetaData(DatObj,metadata = d,col.name = "Age")
  DatObj@meta.data[,'Age'] <- d
  
  d <- paste0(DatObj@meta.data$Condition,'.',DatObj@meta.data$Age)
  DatObj <- AddMetaData(DatObj,metadata = d,col.name = "Condition.Age")
  DatObj@meta.data[,'Condition.Age'] <- d
  
  d <- DatObj@meta.data$orig.ident
  region = c( "Hip", "Crtx")
  for (pt in region) {d[grep(pt, oi)] = pt}
  d[grep("CorA", oi)] = "Crtx"
  DatObj <- AddMetaData(DatObj,metadata = d,col.name = "Region")
  DatObj@meta.data[,'Region'] <- d
  
  d <- DatObj@meta.data$orig.ident
  batch = c("Hip-S1-L", "Hip-S2-L")
  for (pt in batch) {d[grep(pt, oi)] = "NP40"}
  batch = c("Crtx-S1","Crtx-S2","Hip-S1-R", "Hip-S2-R")
  for (pt in batch) {d[grep(pt, oi)] ="NP40"}
  batch = c("G1-4w","G3-2w","G2-4w","G2-2w")
  for (pt in batch) {d[grep(pt, oi)] = "EZ"}
  d[grep("CorA", oi)] = "EZ"
  batch = c("Crtx-7m","Crtx-10m","Crtx-4m","Crtx-6w")
  for (pt in batch) {d[grep(pt, oi)] = "Crtx-TC"}
  batch = c("Hip-7m","Hip-10m","Hip-4m","Hip-7m","Hip-6w")
  for (pt in batch) {d[grep(pt, oi)] = "Hip-TC"}
  batch = c("G1-13m","G1-8w","G2-8w","G2-14m")
  for (pt in batch) {d[grep(pt, oi)] = "Hip-TC.2"}
  DatObj <- AddMetaData(DatObj,metadata = d,col.name = "Batch")
  DatObj@meta.data[,'Batch'] <- d
  
  return(DatObj)
}
