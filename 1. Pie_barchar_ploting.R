library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggthemes)
library(cowplot)
library(data.table)
library(parallel)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(reshape2)
library(ggsci)
library(threed)
library(ggthreed)
library(tibble)
library(dendextend)
library(scales)
library(DescTools)
library(scales)
library(rlang)
library(GSEABase)
library(GSVA)
library(future)
library(parallel)

### the annotation metadata of 84363 cells of AHCA1.0 could be used as demo file. For further use please contact HeShuai! 
meta.dat <- read.table(file = "clipboard",
                       sep = "\t",
                       row.names = NULL,
                       header = T,
                       stringsAsFactors = F,
                       comment.char = "", fileEncoding = "UTF-8")
meta.dat_for_plot <-  meta.dat

DEGs <- read.table("Myeloid_top50_20_1.csv",
                   header = T,
                   row.names = 1,
                   stringsAsFactors = F,
                   sep = ",",
                   comment.char = "") 

levels <- DEGs$cell.type %>% unique()

meta.dat_for_plot$B_subtype <- factor(meta.dat_for_plot$annotation, levels = levels)
meta.dat_for_plot$B_P <- 0
meta.dat_for_plot$B_P[grep("Mon", meta.dat_for_plot$B_subtype, value = F)] <- "Mon"
meta.dat_for_plot$B_P[grep("Mac", meta.dat_for_plot$B_subtype, value = F)] <- "Mac"
meta.dat_for_plot$B_P[grep("DC|Langerhans", meta.dat_for_plot$B_subtype, value = F)] <- "DC"

tmp1 <- meta.dat_for_plot %>% group_by(orig.ident, B_P) %>% summarise(number = n()) 
tmp2 <- tmp1 %>% as.data.frame %>% group_by(orig.ident) %>% summarise_at(.vars = "number", .funs = sum)

tmp1$pre <- tmp1$number/rep(tmp2$number, time = tmp1$orig.ident %>% table)

pdf("pie_plot_according_tissue_legend.pdf", width = 15, height = 7.5)
ggplot(meta.dat_for_plot, aes(x = Tissue, y = Number, fill = Cell.type)) +
  geom_bar(width = 0.8, stat = "identity", color = "white", position="fill") +
  scale_fill_manual(values = color_used) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  coord_polar("y", start = 0) +
  facet_wrap(~ orig.ident, labeller = label_both, nrow = 1) +
  theme_void() +
  theme(legend.position = 'bottom',
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x.top = element_blank(),
        strip.text.x = element_blank())
dev.off()       

pdf("pie_plot_according_tissue.pdf", width = 15, height = 7.5)
ggplot(meta.dat_for_plot) + 
  geom_threedpie(aes(as.factor(B_P))) + 
  facet_wrap(~ orig.ident, labeller = label_both, nrow = 1) +
  scale_fill_brewer(name = "Cut", palette = 'Set2') +
  theme_void() + 
  theme(legend.position = 'bottom',
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x.top = element_blank(),
        strip.text.x = element_blank())
dev.off()

pdf("bar_plot_according_tissue.pdf", width = 15, height = 7.5)
ggplot(data = meta.dat_for_plot, mapping = aes(x = orig.ident, fill = B_subtype)) + 
  geom_bar(position = "fill") +
  scale_fill_manual(values = color_used) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = guide_legend(ncol= 1,
                             label.theme = element_text(size = 5),
                             title.theme = element_text(size = 5))) +
  theme(panel.background = element_rect(fill = NULL, colour = "black", 
                                        size = 1, linetype = "solid"),
        axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.2))
dev.off()

dat <- read.table("merged_annotations.csv", sep = ",", header = T, row.names = NULL, stringsAsFactors = F)
adult_dat <- dat %>% subset(Development_stage == "Adult")

adult_dat_immune <- adult_dat[adult_dat$Celltype %>% grepl(pattern = "B cell|T cell|Dendritic|Eosinophil|Kuppfer|Macrophage|Mast|Myeloid|NK|Plasm"), ]

Celltype <- adult_dat_immune$Celltype 
Celltype[Celltype %>% grepl(pattern = "B cell")] <- "B Cell"
Celltype[Celltype %>% grepl(pattern = "T cell")] <- "T Cell"
Celltype[Celltype %>% grepl(pattern = "Mast")] <- "Myeloid"
Celltype[Celltype %>% grepl(pattern = "Mast")] <- "Myeloid"
Celltype[Celltype %>% grepl(pattern = "Macrophage")] <- "Myeloid"
Celltype[Celltype %>% grepl(pattern = "Myeloid")] <- "Myeloid"
Celltype[Celltype %>% grepl(pattern = "Dendritic|dendritic")] <- "Myeloid"
Celltype[Celltype %>% grepl(pattern = "Kuppfer")] <- "Myeloid"
Celltype[Celltype %>% grepl(pattern = "Eosinophil")] <- "Myeloid"
Celltype[Celltype %>% grepl(pattern = "Plasmocyte")] <- "Plasmocyte"
Celltype[Celltype %>% grepl(pattern = 'Intercalated cell_SPINK1 high')] <- "disc"

adult_dat_immune$Celltype <- Celltype
adult_dat_immune <- adult_dat_immune %>% subset(Celltype != "disc")

pdf("pie_plot_according_tissue_legend.pdf", width = 15, height = 7.5)
ggplot(adult_dat_immune, aes(Sample)) +
  geom_bar(width = 0.8, color = "white", position="fill", aes(fill = Celltype)) +
  scale_fill_manual(values = color_used) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
coord_polar("y", start = 0) +
  facet_wrap(~ orig.ident, labeller = label_both, nrow = 1) +
  theme_void() +
  theme(legend.position = 'bottom',
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x.top = element_blank(),
        strip.text.x = element_blank())
dev.off()       

