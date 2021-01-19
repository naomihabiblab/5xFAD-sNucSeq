library(destiny)
library(GGally)
library(ggrepel)
library(rgexf)
library(cccd)


DC.gene.correlation <- function(DC,data,dc.genes){
  curr.DC.gene.cor <- apply(data[dc.genes,],1,function(x) {cor(x,DC)})
  data.to.plot <- data.frame(melt(curr.DC.gene.cor))
  data.to.plot <- data.to.plot[order(data.to.plot$value),,drop=F]
  data.to.plot$index <- 1:nrow(data.to.plot)
  data.to.name <- data.to.plot[c(1:10,(nrow(data.to.plot)-10):nrow(data.to.plot)),]
  P <- ggplot(data.to.plot,aes(index,value)) + geom_point() + geom_text_repel(data=data.to.name,label=rownames(data.to.name))
  return(P)
}

# dataframe genes/cells (@data or scale.data), gene-array (var.genes? PC-loading?), ceel.info=metadata,n_eigen=difusion-componenet(30-50) num, K=KNN, sigma=kernal-width(number)or 'local'
compute.diffusion.map <- function(data,dc.genes,cell.info,n_eigen=30,k=40,sigma,out.dir)
{
  dmap <- DiffusionMap(data = t(data[rownames(data) %in% dc.genes,]),n_eigs = n_eigen,k = k,sigma = sigma)
  rownames(dmap@eigenvectors) <- colnames(data)
  saveRDS(object = dmap,file = paste0(out.dir,'dmap.rds'))
  
  pdf(paste0(out.dir,'Dmap_Info.pdf'),width = 10,height = 4)
  par(mfrow=c(1,2))
  plot(dmap@eigenvalues,pch=20) 
  plot(dmap@sigmas)
  dev.off()
  
  # visualize first few DC
  DCs.to.visualize <- c(1:8)
  data.to.plot <- data.frame(cell.info,dmap@eigenvectors[,DCs.to.visualize])
  P <- ggpairs(data = data.to.plot,columns = which(colnames(data.to.plot) %in% paste0('DC',DCs.to.visualize))) #mapping = ggplot2::aes(color=cell.id)
  ggsave(plot = P,filename = paste0(out.dir,'DC_top_by_cell_id.png'),width = 18,height = 18)
  # P <- ggpairs(data = data.to.plot,columns = which(colnames(data.to.plot) %in% paste0('DC',DCs.to.visualize)),mapping = ggplot2::aes(color=tmp))
  # ggsave(plot = P,filename = paste0(out.dir,'DC_top_by_tmp.png'),width = 18,height = 18)
  # 
  DC.w.nGene <- melt(data.frame(nGene=cell.info$nGene,dmap@eigenvectors[,DCs.to.visualize]),id.vars='nGene')
  P <- ggplot(DC.w.nGene,aes(nGene,value)) + geom_point(size=1) + facet_wrap(~variable)
  ggsave(plot = P,filename = paste0(out.dir,'DC_versus_nGene.png'))
  
  plot.list <- apply(dmap@eigenvectors[,1:8],2,function(x) {DC.gene.correlation(x,data,dc.genes)})
  ggsave(plot = plot_grid(plotlist = plot.list,ncol = 4),
         filename = paste0(out.dir,'DC_top_cor_with_genes.png'),
         width=15,height = 12)
  return(dmap)
}