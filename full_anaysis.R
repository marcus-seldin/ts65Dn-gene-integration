library(ggVennDiagram)
library(qgraph)
library(factoextra)
library(FactoMineR)
library(qvalue)
library(reshape2)
library(WGCNA)
library(ggplot2)
library(MetBrewer)

load('GTEx NA included env.RData')
working_dataset=GTEx_subfiltered
row.names(working_dataset) = working_dataset$gene_tissue
working_dataset$gene_tissue=NULL
working_dataset = as.data.frame(t(working_dataset))
working_dataset[1:5,1:5]
test1 = working_dataset[,grepl('ITIH5', colnames(working_dataset))]
colnames(test1)

liver = read.csv('mouse DEGS/liver degs.csv')
liver$tissue = paste0('liver')
gwat = read.csv('mouse DEGS/GWAT degs.csv')
gwat$tissue = paste0('gwat')
iwat = read.csv('mouse DEGS/iwat degs.csv')
iwat$tissue = paste0('iwat')
muscle = read.csv('mouse DEGS/sk musc degs.csv')
muscle$tissue = paste0('muscle')


gene_list <- list(Liver= liver$Symbol[liver$Liver.P.Value<0.01],
                  Iwat = iwat$Symbol[iwat$iWAT.P.Value<0.01],
                  Gwat = gwat$Symbol[gwat$gWAT.P.Value<0.01],
                  Muscle = muscle$Symbol[muscle$SkMusc.P.Value<0.01])
pdf(file = 'Venn Diagram Significant Genes across tissues P less 0.01.pdf')
ggVennDiagram(gene_list, label_alpha = 0) + scale_fill_distiller(palette = "RdYlBu") + ggtitle('Tissue-sharing of DEGs P <0.01')
dev.off()



gene_list <- list(Liver= liver$Symbol[liver$Liver.P.Value<0.05],
                  Iwat = iwat$Symbol[iwat$iWAT.P.Value<0.05],
                  Gwat = gwat$Symbol[gwat$gWAT.P.Value<0.05],
                  Muscle = muscle$Symbol[muscle$SkMusc.P.Value<0.05])
pdf(file = 'Venn Diagram Significant Genes across tissues P less 0.05.pdf')
ggVennDiagram(gene_list, label_alpha = 0) + scale_fill_distiller(palette = "RdYlBu") + ggtitle('Tissue-sharing of DEGs P <0.05')
dev.off()


gene_list <- list(Liver= liver$Symbol[liver$Liver.P.Value<0.001],
                  Iwat = iwat$Symbol[iwat$iWAT.P.Value<0.001],
                  Gwat = gwat$Symbol[gwat$gWAT.P.Value<0.001],
                  Muscle = muscle$Symbol[muscle$SkMusc.P.Value<0.001])
pdf(file = 'Venn Diagram Significant Genes across tissues P less 0.001.pdf')
ggVennDiagram(gene_list, label_alpha = 0) + scale_fill_distiller(palette = "RdYlBu") + ggtitle('Tissue-sharing of DEGs P <0.001')
dev.off()



df1 = as.data.frame(liver$Symbol[liver$Liver.P.Value<0.001])
colnames(df1) = 'gene_symbol'
df1$tissue = paste0('Liver')

df2 = as.data.frame(iwat$Symbol[iwat$iWAT.P.Value<0.001])
colnames(df2) = 'gene_symbol'
df2$tissue = paste0('iwat')

df3 = as.data.frame(gwat$Symbol[gwat$gWAT.P.Value<0.001])
colnames(df3) = 'gene_symbol'
df3$tissue = paste0('gwat')

df4 = as.data.frame(muscle$Symbol[muscle$SkMusc.P.Value<0.001])
colnames(df4) = 'gene_symbol'
df4$tissue = paste0('muscle')

full_sig_degs = as.data.frame(rbind(df1, df2, df3, df4))
full_sig_degs = na.omit(full_sig_degs)
table(full_sig_degs$tissue)

orth_table = read.delim('Mouse Gene info with Human Orthologues.txt')
full_sig_degs$human_orth = orth_table$human_orth[match(full_sig_degs$gene_symbol, orth_table$Symbol)]


####################################################################


full_sig_degs$human_tissue = ifelse(full_sig_degs$tissue=='Liver', 'Liver', '')
full_sig_degs$human_tissue = ifelse(full_sig_degs$tissue=='iwat', 'Adipose - Subcutaneous', paste0(full_sig_degs$human_tissue))
full_sig_degs$human_tissue = ifelse(full_sig_degs$tissue=='gwat', 'Adipose - Visceral (Omentum)', paste0(full_sig_degs$human_tissue))
full_sig_degs$human_tissue = ifelse(full_sig_degs$tissue=='muscle', 'Muscle - Skeletal', paste0(full_sig_degs$human_tissue))
full_sig_degs = na.omit(full_sig_degs)
full_sig_degs$gene_tissueH = paste0(full_sig_degs$human_orth, '_', full_sig_degs$human_tissue)
gene_set = full_sig_degs$gene_tissueH
isogenes = working_dataset[,colnames(working_dataset) %in% gene_set]
ii = na.omit(isogenes)

cc1 = bicorAndPvalue(ii, ii, use = 'p')
cc3 = cc1$bicor
cc3[is.na(cc3)] = 0
cc4 = ifelse(cc1$p<0.01, '*', '')
 
colnames(ii)[1:10]

anno = data.frame(row.names(cc3), Group=gsub(".*_","", row.names(cc3)))
row.names(anno) = row.names(cc3)
library(RColorBrewer)
anno$row.names.cc3.=NULL
pdf(file = 'global cor structure of DEGS.pdf')
breaksList = seq(-1, 1, by = .1)
pheatmap::pheatmap(cc3, annotation_row = anno, display_numbers = F, labels_row = F, fontsize_number = 0, labels_col = F, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdGy")))(length(breaksList)) )
dev.off()


cors_df = reshape2::melt(as.matrix(cc1$bicor))  
colnames(cors_df) = c('gene_tissue_1', 'gene_tissue_2', 'bicor')
nn1 = reshape2::melt(as.matrix(cc1$p))  
cors_df$pvalue = nn1$value
qq1 = qvalue::qvalue(cors_df$pvalue)
cors_df$qvalue = qq1$qvalues

cors_df = cors_df[order(cors_df$qvalue, decreasing = F),]
cors_df$gene_symbol_1 = gsub("\\_.*","",cors_df$gene_tissue_1)
cors_df$gene_symbol_2 = gsub("\\_.*","",cors_df$gene_tissue_2)
cors_df$tissue_1 = gsub(".*_","",cors_df$gene_tissue_1)
cors_df$tissue_2 = gsub(".*_","",cors_df$gene_tissue_2)
head(cors_df)
tt1= cors_df[!cors_df$gene_tissue_1==cors_df$gene_tissue_2,]

sec_prots = read.delim('human secreted proteins.tab')
#untag to use secrete proteins only
#tt1 = tt1[tt1$gene_symbol_1 %in% sec_prots$Gene.names...primary..,]
tt2 = tt1 %>% dplyr::select(gene_tissue_1, bicor) %>% dplyr::group_by(gene_tissue_1) %>% dplyr::summarise(mean_abs=mean(abs(bicor), na.rm=T))

top_genes = tt2[order(tt2$mean_abs, decreasing = T),]
topgset = as.vector(top_genes$gene_tissue_1[1:25])
nn1 = tt1[tt1$gene_tissue_1 %in% topgset,]


write.csv(nn1, file = 'gene enrichments with cumulative DEGs.csv', row.names = F)


tissue_list = c('Adipose - Subcutaneous', 'Liver')
new_working = working_dataset[, grepl('Adipose - Subcutaneous', colnames(working_dataset)) | grepl('Liver', colnames(working_dataset)) | grepl('Adipose - Visceral (Omentum)', colnames(working_dataset), fixed=T) | grepl('Muscle - Skeletal', colnames(working_dataset))]
################################################################################################################
#Now run on top gene

pathways_annots = read.delim('uniprot-human-genes and goterms mapping.tab')
pathw_set_adip = pathways_annots[grepl('endoplasmic reticulum', pathways_annots$Gene.ontology..biological.process.),]
pathw_set_adip$gene_tissue = paste0(pathw_set_adip$Gene.names...primary.., '_', 'Adipose - Subcutaneous')
pathw_set_liv = pathways_annots[grepl('interferon-alpha', pathways_annots$Gene.ontology..biological.process.),]
pathw_set_liv$gene_tissue = paste0(pathw_set_liv$Gene.names...primary.., '_', 'Liver')

full_paths = c(pathw_set_adip$gene_tissue[1:200], pathw_set_liv$gene_tissue[1:200])
tt3 = tt1[tt1$gene_tissue_2 %in% full_paths,]
tt3 = na.omit(tt3)
tt3 = tt3[tt3$gene_symbol_1 %in% sec_prots$Gene.names...primary..,]
tt2 = tt3 %>% dplyr::select(gene_tissue_1, bicor) %>% dplyr::group_by(gene_tissue_1) %>% dplyr::summarise(mean_abs=mean(abs(bicor), na.rm=T))

top_genes = tt2[order(tt2$mean_abs, decreasing = T),]
topgset = as.vector(top_genes$gene_tissue_1[1:20])
nn1 = tt1[tt1$gene_tissue_1 %in% topgset,]
nn1$logq = -log(nn1$qvalue)
library(forcats)
nn1$meanbics = tt2$mean_abs[match(nn1$gene_tissue_1, tt2$gene_tissue_1)]
nn1 = na.omit(nn1)
length(unique(nn1$gene_tissue_1))
pdf(file = 'top-ranked candidates mean connectivity.pdf')
ggplot(nn1, aes(x=fct_reorder(gene_tissue_1, meanbics, .desc = TRUE), y=logq, fill=tissue_1)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + geom_violin(width=0.6) + geom_boxplot(width=0.2, position = position_dodge(width=0.6), alpha=0.3, color='grey') + xlab('') + ylab('gene-gene connectivity -log(qvalue)')
dev.off()

set1 = c(paste(nn1$gene_tissue_1[1]), full_paths)


all_tog = working_dataset[,colnames(working_dataset) %in% set1]

cc1 = bicor(all_tog, all_tog, use = 'p')
map1 = as.data.frame(cc1)
map1 = reshape2::melt(as.matrix(map1))
map1$tissue1 = gsub(".*_","",map1$Var1)
map1$tissue2 = gsub(".*_","",map1$Var2)
head(map1)
map1$tissue_col = ifelse(grepl('ATRN', map1$Var1), 'seagreen1', 'darkorange2')
map1$tissue_col = ifelse(grepl('Liver', map1$tissue1), 'darkorchid2', paste0(map1$tissue_col))
#map1$value = ifelse(map1$tissue1=='Adipose - Subcutaneous' & map1$tissue2=='Adipose - Subcutaneous', paste0(map1$value*0.3), paste0(map1$value))
map1$value = as.numeric(map1$value)
map1$value[map1$value > 0.999999] <- 0
map2 = reshape2::dcast(map1, Var1 ~ Var2, value.var = 'value')
row.names(map2) = map2$Var1
map2$Var1 = NULL
table(map1$tissue1)
colkey1 = colnames(map2)
names(colkey1) = map1$tissue_col[match(colkey1, map1$Var1)]
pdf(file = paste('Undirected network qvalLess 1e-3 - ', "ATRN", '.pdf'))
qgraph(map2, minimum = 0.25, cut = 0.6, vsize = 2, color=names(colkey1), legend = F, borders = TRUE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3", label.cex=F, directed=F, labels = F) + ggtitle('')
dev.off()


names(colkey1) = map1$tissue_col[match(colkey1, map1$Var1)]
pdf(file = paste('Undirected network qvalLess 1e-3 - ', "ATRN", ' with labels.pdf'))
qgraph(map2, minimum = 0.2, cut = 0.6, vsize = 2, color=names(colkey1), legend = F, borders = TRUE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3", label.cex=3, directed=F, labels = gsub("\\_.*","",colnames(map2))) + ggtitle('')
dev.off()

