library(rmarkdown)
library(tinytex)
library(knitr)
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(edgeR) # well known package for differential expression analysis, but we only use for the DGEList object and for normalization methods
library(matrixStats) # let's us easily calculate stats on rows or columns of a data matrix
library(cowplot) # allows you to combine multiple plots in one figure
library(DT)
library(gt)
library(plotly)
library(ggrepel)
library(pheatmap)
library(limma) #linear modeling for microarray data
library(Biobase) #for expressionset objects
library(RColorBrewer)
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(gplots)
library(ggpubr)
library(uwot)

#setwd----
setwd("/Users/elizabethli/Box Sync/Castleman")

#read in metadata and qc information, then make complete study design file----
study.design = readxl::read_xlsx("AI panel annotation with other info.xlsx")
qc = readxl::read_xlsx("Copy of VLP00554_IR_PLATE1_19FEB2021_Parsed-forReveal_QualityControlled.xlsx", sheet = "QC Summary", range = "A8:i79")

colnames(qc)[1] = "sample.name"
study.design = study.design %>% left_join(qc)
study.design$group[study.design$group == "React- C"] = "ReactC"
study.design$disease[study.design$disease == "React- C"] = "ReactC"
study.design$group[study.design$group == "React -k"] = "ReactK"
study.design$disease[study.design$disease == "React -k"] = "ReactK"

#only keep kikuchi and ReactK samples that passed QC1 in study design (Kikuchi = KFD, reactK = control)
study.design = study.design %>% dplyr::filter(disease %in% c("Kikuchi", "ReactK"), `QC Status` != "QC1 Fail")

#read in counts, and clean up matrix----
counts = readxl::read_xlsx("Copy of VLP00554_IR_PLATE1_19FEB2021_Parsed-forReveal_QualityControlled.xlsx", sheet = "Raw", range = "A10:BT2031")

#drop total counts (first row)
if(counts[1,1] == "Total Counts"){
  total.counts.from.htg = counts[1,]
  counts = counts[-1,]
}

#rename column
colnames(counts)[1] = "gene.name"

#subset counts matrix only for KFD/control samples that passed QC1 (from study design)
raw.counts = counts[,colnames(counts) %in% c("gene.name", study.design$sample.name)]

#make matrix + clean
gene.rows = raw.counts$gene.name
raw.counts = raw.counts %>% dplyr::select(-gene.name) %>% as.matrix()
rownames(raw.counts) = gene.rows

#order study design samples the same as raw.counts to get sample labels----
study.design = study.design[match( raw.counts %>% colnames(), study.design$sample.name),]
samplelabels <- study.design$sample.name
#make DGElist----
myDGEList <- DGEList(raw.counts, samples = study.design)
log2.cpm <- cpm(myDGEList, log = T)
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = samplelabels, # column names to be stored as a SINGLE variable
                                  names_to = "sample_id", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values

#plot log2 counts
p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=sample_id, y=expression, fill=sample_id) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 20, 
               size = 2, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()

#filter DGEList
cpm <- cpm(myDGEList, log = F)
keepers <- rowSums(cpm>15)>=3
table(keepers)
myDGEList.filtered <- myDGEList[keepers,]

#plot filtered log2 counts
log2.cpm.filtered.df <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered.df, rownames = "geneID")
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = samplelabels, # column names to be stored as a SINGLE variable
                                           names_to = "sample_id", # name of that new variable (column)
                                           values_to = "expression") # name of new variable

p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=sample_id, y=expression, fill=sample_id) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 20, 
               size = 2, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()

#normalize log2 counts
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)


log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = samplelabels, # column names to be stored as a SINGLE variable
                                                names_to = "sample_id", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=sample_id, y=expression, fill=sample_id) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 20, 
               size = 2, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()

#plot all figures
p4 = plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12, nrow = 1)

figure.directory = paste0(getwd(),'/Figures_KikuchiOnly')
dir.create(figure.directory)
ggsave(p4, filename = "ExpressionQC.pdf", height =10, width = 30, path = figure.directory)

#doPCA----
pca.res <- prcomp(t(log2.cpm.filtered.norm), center=TRUE, scale.=F, retx=T)

pc.var<-pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x)

summary(pca.res) # Prints variance summary for all principal components.
screeplot(pca.res, npcs = 30)

#make PCA plots
study_design = study.design

pdf(paste0(figure.directory, "/PCA-Analysis.pdf"), height = 10, width = 10)
ggplot(pca.res.df) + # pcaplot (1,2) colored by treatment, shape by cell line
  aes(x=PC1, y=PC2) +
  geom_point(aes(shape = study_design$`QC Status`, color=study_design$disease), size=3, show.legend = T) +
  geom_text_repel(aes(label = paste0(study_design$sample.name, ":", study_design$comments)), size = 2.5) +
  xlab(paste0("PC1 (", pc.per[1],"%)")) + 
  ylab(paste0("PC2 (", pc.per[2],"%)")) +
  labs(title="PCA plot, Kikuchi vs ReactK",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggplot(pca.res.df) + # pcaplot (3,4) colored by treatment, shape by cell line
  aes(x=PC3, y=PC4) +
  geom_point(aes(shape = study_design$`QC Status`, color=study_design$disease), size=3, show.legend = T) +
  geom_text_repel(aes(label = paste0(study_design$sample.name, ":", study_design$comments)), size = 2) +
  xlab(paste0("PC3 (", pc.per[3],"%)")) + 
  ylab(paste0("PC4 (", pc.per[4],"%)")) +
  labs(title="PCA plot, Kikuchi vs ReactK",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggplot(pca.res.df) + # pcaplot (5,6) colored by treatment, shape by cell line
  aes(x=PC5, y=PC6) +
  geom_point(aes(shape = study_design$`QC Status`, color=study_design$disease), size=3, show.legend = T) +
  geom_text_repel(aes(label = paste0(study_design$sample.name, ":", study_design$comments)), size = 2) +
  xlab(paste0("PC5 (", pc.per[5],"%)")) + 
  ylab(paste0("PC6 (", pc.per[6],"%)")) +
  labs(title="PCA plot, Kikuchi vs ReactK",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggplot(pca.res.df) + # pcaplot (7,8) colored by treatment, shape by cell line
  aes(x=PC7, y=PC8) +
  geom_point(aes(shape = study_design$`QC Status`, color=study_design$disease), size=3, show.legend = T) +
  geom_text_repel(aes(label = paste0(study_design$sample.name, ":", study_design$comments)), size = 2) +
  xlab(paste0("PC7 (", pc.per[7],"%)")) + 
  ylab(paste0("PC8 (", pc.per[8],"%)")) +
  labs(title="PCA plot, Kikuchi vs ReactK",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()
dev.off()

#make UMAP----
npc = 10 #number of PCs from PCA

umap.reduction = uwot::umap(pca.res.df[,1:npc]) %>% as_tibble()


colnames(umap.reduction) = c("UMAP_1", "UMAP_2")

umap.reduction = cbind(umap.reduction, study_design)

ggplot(umap.reduction, aes(x= UMAP_1, y = UMAP_2)) +
  geom_point(aes(shape = `QC Status`, color= disease), size=3) +theme_classic() + ggtitle(paste("UMAP:", npc, "PCs"))

plot.list = list()
for(npc in c(2:20)){
  
  umap.reduction = uwot::umap(pca.res.df[,1:npc])
  
  umap.reduction = umap.reduction %>% as_tibble()
  colnames(umap.reduction) = c("UMAP_1", "UMAP_2")
  
  umap.reduction = cbind(umap.reduction, study_design)
  
  plot.list[[npc -1 ]]= ggplot(umap.reduction, aes(x= UMAP_1, y = UMAP_2)) +
    geom_point(aes(shape = `QC Status`, color= disease), size=3) +theme_classic() + ggtitle(paste("UMAP:", npc, "PCs"))
  
  
  
}

pdf(paste0(figure.directory,"/screeplot.pca.pdf"))
screeplot(pca.res, npcs = npc, type = "lines")
dev.off()

umap.plot = cowplot::plot_grid(plotlist = plot.list)

plot.list[[8]]+ geom_text(aes(label = sample.name == weird))

ggsave(umap.plot, filename = "umap.plots.pdf", height =20, width = 25, path = figure.directory)


#perform DE testing----


log2.filtered.norm.counts <- myDGEList.filtered.norm$counts #assayData for eset



metadata <- as.data.frame(study_design) #phenoData for eset
rownames(metadata) <- samplelabels #same rownames

genes <- rownames(log2.filtered.norm.counts) %>% as.data.frame() #featureData for eset
rownames(genes) <- genes[1:nrow(genes), 1]
colnames(genes) <- "geneID"

eset <- ExpressionSet(assayData = log2.filtered.norm.counts, #eset object
                      phenoData = AnnotatedDataFrame(metadata),
                      featureData = AnnotatedDataFrame(genes))

group <- with(pData(eset), paste0(disease)) #cat-ing groups for design matrix
group <- factor(group) %>% print()

design <- model.matrix(~0 + group) %>% print() #design matrix
colnames(design) <- levels(group)
colSums(design) #checking design matrix

v.eset <- voom(eset, design, plot = TRUE) #running voom on eset --> log2cpm counts and weights.

cm <- makeContrasts(Kikuchi.vs.Control = Kikuchi -ReactK,
                    levels = design)


fit <- lmFit(v.eset, design)
fit_all <- contrasts.fit(fit, contrasts = cm)
fit_all <- eBayes(fit_all)

summary(decideTests(fit_all, method="separate", adjust.method="BH", p.value=0.25))

results_all_p025 <- decideTests(fit_all, adjust.method="BH", p.value=0.25, method = "separate")

#2 or more, same direction, down
dges <- results_all_p025 %>% as.data.frame.matrix() %>% 
  as.data.frame(rownames = "geneID")

#make volcano plots----

library(EnhancedVolcano)

volcano.fig.dir = paste0(figure.directory, "/VolcanoPlots")
dir.create(volcano.fig.dir)

limma.list = list()

for(i in 1:length(colnames(fit_all))){
  
  comparison = colnames(fit_all)[i] %>% print()
  
  
  limma.list[[i]] = topTable(fit_all, coef = i, number=nrow(genes), adjust ="BH", p.value = 0.5) %>%
    as_tibble() %>% print()
  
  if(nrow(limma.list[[i]]) ==0){
    print("no sig genes")
    
  }else{
    
    #volcanoplot(fit_all, coef = i, style = "p-value", highlight = 20, names = fit$genes$geneID, hl.col="blue",
    #xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)
    #ggsave(path = volcano.fig.dir, paste0("limma_default_", comparison, ".png"), height = 10, width = 10)
    
    
    
    
    top.genes = limma.list[[i]]$geneID[limma.list[[i]]$logFC>0] %>% head(25)
    top.genes = c(top.genes, limma.list[[i]]$geneID[limma.list[[i]]$logFC<0] %>% head(25))
    
    ylim = limma.list[[i]]$adj.P.Val %>% min()
    
    
    EnhancedVolcano(limma.list[[i]],
                    lab = limma.list[[i]]$geneID,
                    x = 'logFC',
                    y = 'adj.P.Val',
                    title = comparison,
                    #selectLab = top.genes,
                    pCutoff = 0.05,
                    FCcutoff = 1,
                    pointSize = 3.0,
                    #ylim = c(0, -log10(10e-5)),
                    ylim = c(0, -log10(ylim / 1.1)),
                    labSize = 3, drawConnectors = T) 
    
    ggsave(path = volcano.fig.dir, paste0(comparison, "_ALL.pdf"), height = 10, width = 10)
    
    
    EnhancedVolcano(limma.list[[i]],
                    lab = limma.list[[i]]$geneID,
                    x = 'logFC',
                    y = 'adj.P.Val',
                    title = comparison,
                    #selectLab = top.genes,
                    pCutoff = 0.05,
                    FCcutoff = 0.5,
                    pointSize = 3.0,
                    labSize = 3,
                    drawConnectors = F, 
                    ylim = c(0, -log10(ylim / 1.1))) 
    
    ggsave(path = volcano.fig.dir, paste0(comparison, "_Default.pdf"), height = 10, width = 10)
    
    
    EnhancedVolcano(limma.list[[i]],
                    lab = limma.list[[i]]$geneID,
                    x = 'logFC',
                    y = 'adj.P.Val',
                    title = paste(comparison, "relaxed cutoffs"),
                    #selectLab = top.genes,
                    pCutoff = 0.1,
                    FCcutoff = 0.25,
                    pointSize = 3.0,
                    labSize = 3,
                    drawConnectors = F, 
                    ylim = c(0, -log10(ylim / 1.1))) 
    
    ggsave(path = volcano.fig.dir, paste0(comparison, "_0.25FC-0.1Pval.pdf"), height = 10, width = 10)
    
    
  }
  names(limma.list)[[i]] = comparison
  
}
#do heatmap with genes > 0.25 and p<0.1
heatmap.genes = limma.list[[i]] %>% filter(abs(logFC) > 0.25, adj.P.Val < 0.1) %>% arrange(logFC)
metadata = metadata %>% arrange(disease)
heatmap.plot = log2.cpm.filtered.norm.df %>% filter(geneID %in% heatmap.genes$geneID) 
rownames(heatmap.plot) = heatmap.plot$geneID
heatmap.plot = heatmap.plot %>% as.data.frame()
heatmap.plot = heatmap.plot[, -1]
heatmap.plot = heatmap.plot[heatmap.genes$geneID,metadata$sample.name]

heatmap.plot.test = heatmap.plot[heatmap.genes$geneID,metadata$sample.name]

heat_colors <- brewer.pal(6, "YlOrRd")

row.anno = data.frame(Enriched = ifelse(heatmap.genes$logFC > 0,yes = "Kikuchi.Enriched", no = "Reactive-LN.Enriched"))
rownames(row.anno) = heatmap.genes$geneID            

pdf(paste0(figure.directory, "/heatmap_prelim_options.pdf"))
#scale by row
pheatmap(heatmap.plot,
         #color = heat_colors,
         cluster_rows = F,
         cluster_cols  = F,
         show_rownames = T,
         annotation = dplyr::select(metadata, disease),
         annotation_row = row.anno,
         scale = "row",
         #cale = "column",
         fontsize_row = 5
)

#scale by column
pheatmap(heatmap.plot,
         #color = heat_colors,
         cluster_rows = F,
         cluster_cols  = F,
         show_rownames = T,
         annotation = dplyr::select(metadata, disease),
         annotation_row = row.anno,
         #scale = "row",
         scale = "column",
         fontsize_row = 5
)

#clustered genes
pheatmap(heatmap.plot,
         #color = heat_colors,
         cluster_rows = T,
         cluster_cols  = F,
         show_rownames = T,
         annotation = dplyr::select(metadata, disease),
         annotation_row = row.anno,
         scale = "row",fontsize_row = 5
)

#clustered genes and samples
pheatmap(heatmap.plot,
         #color = heat_colors,
         cluster_rows = T,
         cluster_cols  = T,
         show_rownames = T,
         annotation = dplyr::select(metadata, disease),
         annotation_row = row.anno,
         scale = "row",fontsize_row = 5
)
dev.off()

#now do a strict cutoff
heatmap.genes = limma.list[[i]] %>% filter(abs(logFC) > 0.5, adj.P.Val < 0.05) %>% arrange(logFC)
metadata = metadata %>% arrange(disease)
heatmap.plot = log2.cpm.filtered.norm.df %>% filter(geneID %in% heatmap.genes$geneID) 
rownames(heatmap.plot) = heatmap.plot$geneID

heatmap.plot = heatmap.plot %>% as.data.frame()
heatmap.plot = heatmap.plot[, -1]
heatmap.plot = heatmap.plot[heatmap.genes$geneID,metadata$sample.name]
#heat_colors <- brewer.pal(6, "YlOrRd")

row.anno = data.frame(Enriched = ifelse(heatmap.genes$logFC > 0,yes = "Kikuchi.Enriched", no = "Reactive-LN.Enriched"))
rownames(row.anno) = heatmap.genes$geneID            

pdf(paste0(figure.directory, "/heatmap_prelim_options_strict_cutoff.pdf"))

#scale by row
pheatmap(heatmap.plot,
         #color = heat_colors,
         cluster_rows = F,
         cluster_cols  = F,
         show_rownames = T,
         annotation = dplyr::select(metadata, disease),
         annotation_row = row.anno,
         scale = "row",
         #cale = "column",
         fontsize_row = 5
)

#scale by column
pheatmap(heatmap.plot,
         #color = heat_colors,
         cluster_rows = F,
         cluster_cols  = F,
         show_rownames = T,
         annotation = dplyr::select(metadata, disease),
         annotation_row = row.anno,
         #scale = "row",
         scale = "column",
         fontsize_row = 5
)

#clustered genes
pheatmap(heatmap.plot,
         #color = heat_colors,
         cluster_rows = T,
         cluster_cols  = F,
         show_rownames = T,
         annotation = dplyr::select(metadata, disease),
         annotation_row = row.anno,
         scale = "row",fontsize_row = 5
)

#clustered genes and samples
pheatmap(heatmap.plot,
         #color = heat_colors,
         cluster_rows = T,
         cluster_cols  = T,
         show_rownames = T,
         annotation = dplyr::select(metadata, disease),
         annotation_row = row.anno,
         scale = "row",fontsize_row = 5
)
dev.off()

#now only plot IFI genes
heatmap.genes = limma.list[[i]] %>% filter(abs(logFC) > 0.25, adj.P.Val < 0.1) %>% arrange(logFC)

ifi.genes = grep(heatmap.genes$geneID, pattern = "IFI|ISG|IRF", value = T) %>% print()

heatmap.genes = heatmap.genes %>% filter(geneID %in% ifi.genes)
metadata = metadata %>% arrange(disease)
heatmap.plot = log2.cpm.filtered.norm.df %>% filter(geneID %in% heatmap.genes$geneID) 
rownames(heatmap.plot) = heatmap.plot$geneID
heatmap.plot = heatmap.plot %>% as.data.frame()
heatmap.plot = heatmap.plot[, -1]
heatmap.plot = heatmap.plot[heatmap.genes$geneID,metadata$sample.name]
heat_colors <- brewer.pal(6, "YlOrRd")

row.anno = data.frame(Enriched = ifelse(heatmap.genes$logFC > 0,yes = "Kikuchi.Enriched", no = "Reactive-LN.Enriched"))
rownames(row.anno) = heatmap.genes$geneID            

pdf(paste0(figure.directory, "/heatmap_IFI_genes_only.pdf"))
#scale by row
pheatmap(heatmap.plot,
         #color = heat_colors,
         cluster_rows = F,
         cluster_cols  = F,
         show_rownames = T,
         annotation = dplyr::select(metadata, disease),
         annotation_row = row.anno,
         scale = "row",
         #cale = "column",
         fontsize_row = 5
)

#scale by column
pheatmap(heatmap.plot,
         #color = heat_colors,
         cluster_rows = F,
         cluster_cols  = F,
         show_rownames = T,
         annotation = dplyr::select(metadata, disease),
         annotation_row = row.anno,
         #scale = "row",
         scale = "column",
         fontsize_row = 5
)

#clustered genes
pheatmap(heatmap.plot,
         #color = heat_colors,
         cluster_rows = T,
         cluster_cols  = F,
         show_rownames = T,
         annotation = dplyr::select(metadata, disease),
         annotation_row = row.anno,
         scale = "row",fontsize_row = 5
)

#clustered genes and samples
pheatmap(heatmap.plot,
         #color = heat_colors,
         cluster_rows = T,
         cluster_cols  = T,
         show_rownames = T,
         annotation = dplyr::select(metadata, disease),
         annotation_row = row.anno,
         scale = "row",fontsize_row = 5
)
dev.off()

#try scaled heatmap
scaled.heatmap.plot = scale(heatmap.plot, center = T)

scaled.heatmap.plot[scaled.heatmap.plot > 2] = 2
scaled.heatmap.plot[scaled.heatmap.plot < -2] = -2

pheatmap(scaled.heatmap.plot,
         #color = heat_colors,
         cluster_rows = F,
         cluster_cols  = F,
         show_rownames = T,
         annotation = dplyr::select(metadata, disease),
         annotation_row = row.anno,
         scale = "none",fontsize_row = 5
)




up.genes = limma.list[[i]] %>% filter(logFC > 0.25, adj.P.Val < 0.1)
down.genes = limma.list[[i]] %>% filter(logFC < -0.25, adj.P.Val < 0.1)

write_csv(up.genes, file = paste0(volcano.fig.dir, "/", names(limma.list)[i], "_UPregulated.csv"))
write_csv(down.genes, file = paste0(volcano.fig.dir,"/", names(limma.list)[i], "_DOWNregulated.csv"))



comparisons = list(c("Kikuchi", "ReactK"))

titles = list(c("Kikuchi", "Reactive LN"))

colors = list(c("#e87d71", "#52b64a"))

for (i in 1:length(limma.list)){
  
  up.genes = limma.list[[i]][limma.list[[i]]$logFC>0,] %>% arrange(adj.P.Val) %>% head(15) %>% print()
  
  up.genes$geneID = factor(up.genes$geneID, levels = up.genes$geneID)
  
  down.genes = limma.list[[i]][limma.list[[i]]$logFC<0,] %>% arrange(adj.P.Val) %>% head(15) %>% print()
  
  down.genes$geneID = factor(down.genes$geneID, levels = down.genes$geneID)
  
  log2.cpm.filtered.norm.df.pivot.study.design = log2.cpm.filtered.norm.df.pivot %>% left_join(study.design, by = c("sample_id" = "sample.name"))
  
  to.graph = log2.cpm.filtered.norm.df.pivot.study.design %>% filter(geneID %in% up.genes$geneID, disease %in% comparisons[[i]])
  
  to.graph$geneID = factor(to.graph$geneID, levels = up.genes$geneID)
  
  y.pos = to.graph %>% group_by(geneID) %>% summarize(max.gene = max(expression))
  
  label.string = round(up.genes$adj.P.Val,3)
  label.string[up.genes$adj.P.Val < 0.1 & up.genes$adj.P.Val > 0.01] = "*"
  label.string[up.genes$adj.P.Val < 0.05 & up.genes$adj.P.Val > 0.01] = "**"
  label.string[up.genes$adj.P.Val <= 0.01] = "***"
  
  #with rounding
  p1 = ggplot(to.graph, aes(x = geneID, y = expression)) + 
    geom_boxplot(aes(fill = disease, color = disease), alpha = 0.5) + theme_classic() + 
    geom_point(aes(fill = disease, color =  disease), position = position_jitterdodge(jitter.width = 0.01))+ 
    geom_bracket(
      xmin=seq(from = 0.75, to = 14.75, by = 1), 
      xmax=seq(from = 1.25, to = 15.25, by = 1), 
      y.position = y.pos$max.gene * 1.05,
      label=label.string, 
      tip.length=0.01, label.size=3.1, vjust=-0.1) + labs(x = "Probe", y = expression("-Log"[2]*"(CPM)"))+
    ggtitle(paste("Top 15 Probes Upregulated in", titles[[i]][1], "vs", titles[[i]][2])) +
    theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values=c(colors[[i]][1], colors[[i]][2])) + 
    scale_fill_manual(values=c(colors[[i]][1], colors[[i]][2]))
  
  boxplot.dir = paste0(figure.directory, "/BoxPlots_final")
  dir.create(boxplot.dir)
  
  comparison = paste0(comparisons[[i]][1], "_vs_", comparisons[[i]][2]) %>% print()
  
  ggsave(p1, path = boxplot.dir, filename = paste0(comparison, "_UP.pdf"), height = 10, width = 10)
  
  to.graph = log2.cpm.filtered.norm.df.pivot.study.design %>% filter(geneID %in% down.genes$geneID, disease %in% comparisons[[i]])
  
  to.graph$geneID = factor(to.graph$geneID, levels = down.genes$geneID)
  
  y.pos = to.graph %>% group_by(geneID) %>% summarize(max.gene = max(expression))
  
  
  
  label.string = round(down.genes$adj.P.Val,3)
  label.string[down.genes$adj.P.Val < 0.1 & down.genes$adj.P.Val > 0.01] = "*"
  label.string[down.genes$adj.P.Val < 0.05 & down.genes$adj.P.Val > 0.01] = "**"
  label.string[down.genes$adj.P.Val <= 0.01] = "***"
  
  p2 = ggplot(to.graph, aes(x = geneID, y = expression)) + 
    geom_boxplot(aes(fill = disease, color = disease), alpha = 0.5) + theme_classic() + 
    geom_point(aes(fill = disease, color = disease), position = position_jitterdodge(jitter.width = 0.01))+ 
    geom_bracket(
      xmin=seq(from = 0.75, to = 14.75, by = 1), 
      xmax=seq(from = 1.25, to = 15.25, by = 1), 
      y.position = y.pos$max.gene * 1.05,
      label= label.string, tip.length=0.01, label.size=3.1, vjust=-0.1) + labs(x = "Probe", y = expression("-Log"[2]*"(CPM)"))+
    ggtitle(paste("Top 10 Probes Downregulated in", titles[[i]][1], "vs", titles[[i]][2])) +
    theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values=c(colors[[i]][1], colors[[i]][2])) + 
    scale_fill_manual(values=c(colors[[i]][1], colors[[i]][2]))
  
  ggsave(p2, path = boxplot.dir, filename = paste0(comparison, "_DOWN.pdf"), height = 10, width = 10)
  
}

#save PCA and UMAP

pca.plot = ggplot(pca.res.df) + # pcaplot colored by treatment, shape by cell line
  aes(x=PC1, y=PC2) +
  geom_point(aes(color=study_design$disease, shape = study_design$disease), size=3, show.legend = T) +
  #geom_text_repel(aes(label = paste0(study_design$sample.name, ":", study_design$comments)), size = 2.5) +
  #stat_ellipse(aes(x=PC1, y=PC2, shape=study_design$timepoint)) +
  xlab(paste0("PC1 (", pc.per[1],"%)")) + 
  ylab(paste0("PC2 (", pc.per[2],"%)")) +
  labs(title="PCA plot, Kikuchi vs ReactK") +
  coord_fixed() +
  theme_classic() 
ggsave(pca.plot, path = figure.directory, filename = paste0("finalPCA.pdf"), height = 10, width = 10)


npc = 9

umap.reduction = uwot::umap(pca.res.df[,1:npc])

umap.reduction = umap.reduction %>% as_tibble()
colnames(umap.reduction) = c("UMAP_1", "UMAP_2")

umap.reduction = cbind(umap.reduction, study_design)

umap.plot = ggplot(umap.reduction, aes(x= UMAP_1, y = UMAP_2)) +
  geom_point(aes(shape = disease, color= disease), size=3) +theme_classic() + ggtitle(paste0("UMAP: ", "Kikuchi vs ReactK (", npc, " PCs)"))

ggsave(umap.plot, path = figure.directory, filename = paste0("finalUMAP.pdf"), height = 10, width = 10)

