#!usr/bin/env Rscript

# R version 4.2
set.seed(34)

## Import 
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggh4x)
library(conflicted)
conflict_prefer("select", "dplyr") ## required in %>% dplyr
conflict_prefer("filter", "dplyr") ## required in %>% dplyr

set.seed(100)
#setwd("D:/backup files/getegbm/GBM-TE-scRNAseq/results/GBM-GSC neuroblastoma TE/2024-08-06_metabolicgenes-genesetenrichment")
print(getwd())

annotations <- data.frame(
    cluster = factor(rep(0:11), levels=rep(0:11)),
    cell.type = factor(c(
        "AC", "OL", "INTER", "LYM", "P.NPC", "MG", "MACRO",
        "P.NPC", "OPC", "D.OPC", "D.OPC", "ENDO"
        ),levels = c(
            "AC", "OL", "OPC","D.OPC","P.NPC","MG", 
            "MACRO","LYM","ENDO","INTER"
        ),#labels = c(
            # "Astrocytes", "Oligodendrocytes", "OPC", "Diff. OPC", 
            # "Prolif. NPC","Microglia","Macrophage",
            # "Lymphoid", "Endoderm","Intermediate"))
    ))
# cycling.anno <- data.frame(
#     cluster = factor(rep(0:11), levels=rep(0:11)),
#     cycling = factor(c(
#         "G1/M","G2/M" # todo
#     ), levels = c("G1/M", "G2/M")),
# )
annotations$cell.type2 <- annotations$cell.type

cell.type.col = c("#FF0066FF", "#328C97FF", "#D1AAC2FF", "#A5506DFF", 
                  "#B3E0BFFF", "#2A9D3DFF", "#EDF181ff", "#db7003ff",
                  "#fba600ff", "#f8c1a6ff", "#a30000ff", "#ff3200ff",
                  "#011a51ff", "#97d1d9ff", "#916c37ff")

gsea.celltypes <- read.csv("gec8-celltypes-summary.csv", header = TRUE)
gsea.celltypes$cluster <- factor(gsea.celltypes$cluster,
                                 levels = levels(annotations$cluster))
gsea.celltypes$NAME <- factor(gsea.celltypes$NAME, 
                              levels = rev(unique(gsea.celltypes$NAME)))
gsea.celltypes$nlog.FDRq <- -log10(gsea.celltypes$FDR.q.val+.00001)
# cell.types <- pivot_wider(names_from = teamA, values_from = win) %>%
#     column_to_rownames('teamB') %>%
#     as.matrix


#simple dot plot with cluster by GSEA terms
# ggplot(gsea.celltypes, aes(x=cluster, y=NAME, color=nlog.FDRq, size=NES)) + 
#     geom_point() + 
#     coord_fixed() + # squared tiles
#     # cowplot::theme_cowplot() +
#     theme(axis.line = element_blank()) +
#     theme(axis.ticks=element_blank()) +
#     ylab(NULL) +
#     viridis::scale_color_viridis() +
#     scale_y_discrete(position="right")

strip_cell.type <- strip_themed(
  # Horizontal strips
  background_x = elem_list_rect(fill = cell.type.col, colour=cell.type.col),
  text_x = elem_list_text(face = "bold", size=6),
  by_layer_x = FALSE,
  # Vertical strips
  background_y = element_blank(),
  text_y = elem_list_text(face = "bold", size=6),
  by_layer_y = FALSE,
)

#simple dot plot with cell.type annotations by GSEA terms
p <- gsea.celltypes %>% left_join(annotations) %>% 
ggplot(aes(x=cluster, y=NAME, color=nlog.FDRq, size=NES)) + 
    geom_point() + 
    # coord_fixed() + # squared tiles
    # cowplot::theme_cowplot() +
    theme_bw()+
    theme(axis.line = element_blank()) +
    theme(axis.ticks=element_blank()) +
    theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0, size=6.6))+
    labs(
      x=NULL,
      y=NULL,
      col="-log(FDR q-val)"
    )+
    scale_y_discrete(position="right")+
    scale_x_discrete(position="top")+
    facet_grid2(cell.type2~cell.type, 
                scales="free", space="free", 
                strip=strip_cell.type)
size <- 5
ggsave("figs_ge_Dotplot_genesetc8-all.tiff", plot = p, units="in", width=size*2, height=size*2.5, dpi=300, compression = 'lzw')


# subset into 3 different plots (normal brain cells)
p <- gsea.celltypes %>% left_join(annotations) %>% 
  filter(cell.type %in% c("AC", "OL", "MG", "OPC", "D.OPC")) %>%
  ggplot(aes(x=cluster, y=NAME, color=nlog.FDRq, size=NES)) + 
      geom_point() + 
      # coord_fixed() + # squared tiles
      # cowplot::theme_cowplot() +
      theme_bw()+
      theme(axis.line = element_blank()) +
      theme(axis.ticks=element_blank()) +
      theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0, size=6.6))+
      labs(
        x=NULL,
        y=NULL,
        col="-log(FDR q-val)"
      )+
      scale_y_discrete(position="right")+
      scale_x_discrete(position="top")+
      facet_grid2(~cell.type, 
                  scales="free", space="free", 
                  strip=strip_cell.type)

ggsave("figs_ge_Dotplot_genesetc8-braincells.tiff", plot = p, units="in", width=size*2, height=size*1.5, dpi=300, compression = 'lzw')

# subset into 3 different plots (Immune)
strip_cell.type <- strip_themed(
  # Horizontal strips
  background_x = elem_list_rect(fill = cell.type.col[4:6], colour=cell.type.col[4:6]),
  text_x = elem_list_text(face = "bold", size=6),
  by_layer_x = FALSE,
  # Vertical strips
  background_y = element_blank(),
  text_y = elem_list_text(face = "bold", size=6),
  by_layer_y = FALSE,
)

p <- gsea.celltypes %>% left_join(annotations) %>% 
  filter(cell.type %in% c("P.NPC")) %>%
  ggplot(aes(x=cluster, y=NAME, color=nlog.FDRq, size=NES)) + 
      geom_point() + 
      # coord_fixed() + # squared tiles
      # cowplot::theme_cowplot() +
      theme_bw()+
      theme(axis.line = element_blank()) +
      theme(axis.ticks=element_blank()) +
      theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0, size=6.6))+
      labs(
        x=NULL,
        y=NULL,
        col="-log(FDR q-val)"
      )+
      scale_y_discrete(position="right")+
      scale_x_discrete(position="top")+
      facet_grid2(~cell.type, 
                  scales="free", space="free", 
                  strip=strip_cell.type)

ggsave("figs_ge_Dotplot_genesetc8-stem.tiff", plot = p, units="in", width=size*2, height=size*1.5, dpi=300, compression = 'lzw')

strip_cell.type <- strip_themed(
  # Horizontal strips
  background_x = elem_list_rect(fill = cell.type.col[6:8], colour=cell.type.col[6:8]),
  text_x = elem_list_text(face = "bold", size=6),
  by_layer_x = FALSE,
  # Vertical strips
  background_y = element_blank(),
  text_y = elem_list_text(face = "bold", size=6),
  by_layer_y = FALSE,
)

p <- gsea.celltypes %>% left_join(annotations) %>% 
  filter(cell.type %in% c("MG","MACRO","LYM")) %>%
  ggplot(aes(x=cluster, y=NAME, color=nlog.FDRq, size=NES)) + 
      geom_point() + 
      # coord_fixed() + # squared tiles
      # cowplot::theme_cowplot() +
      theme_bw()+
      theme(axis.line = element_blank()) +
      theme(axis.ticks=element_blank()) +
      theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0, size=6.6))+
      labs(
        x=NULL,
        y=NULL,
        col="-log(FDR q-val)"
      )+
      scale_y_discrete(position="right")+
      scale_x_discrete(position="top")+
      facet_grid2(~cell.type, 
                  scales="free", space="free", 
                  strip=strip_cell.type)

ggsave("figs_ge_Dotplot_genesetc8-immune.tiff", plot = p, units="in", width=size*1.5, height=size*0.7, dpi=300, compression = 'lzw')

# Single Barplots
p <- gsea.celltypes %>% left_join(annotations) %>% 
  filter(cell.type %in% c("P.NPC")) %>%
  ggplot(aes(y=NAME, x=nlog.FDRq)) + 
    geom_col(aes(fill = NES))+
    theme_minimal()+
    theme(legend.position="top")+
    labs(
      x="-log10(FDR q.val)",
      y="",
      fill = "NES"
    )+
    scale_x_continuous(breaks=c(0,2.5, 5))+
    scale_y_discrete(position="right")+
    scale_fill_continuous(trans="reverse")+
    facet_grid2(~cell.type+cluster, 
                strip = strip_nested(bleed=TRUE, 
                                     background_x = elem_list_rect(fill = c(cell.type.col[5], 
                                                                            "white",
                                                                            "white"), 
                                                                   colour=c(cell.type.col[5], 
                                                                            "white",
                                                                            "white")),
                                     text_x = elem_list_text(face = "bold", size=6),
                                     by_layer_x = TRUE))+
    theme(panel.spacing.x = unit(4, "mm"))

ggsave("figs_ge_barplot_genesetc8-npc.tiff", plot = p, units="in", width=size*1, height=size*0.9, dpi=300, compression = 'lzw')

# For loop on all barplot cell types
ind.plots <- vector(mode="list", length = length(levels(annotations$cell.type)))
p.titles <- vector(mode="string", length = length(levels(annotations$cell.type)))

for (i in 1:length(levels(annotations$cell.type))) {
  ind.plots[[i]] <- gsea.celltypes %>% left_join(annotations) %>% 
    filter(cell.type %in% levels(annotations$cell.type)[i]) %>%
    ggplot(aes(y=NAME, x=nlog.FDRq)) + 
    geom_col(aes(fill = NES))+
    theme_minimal()+
    theme(legend.position="top")+
    labs(
      x="-log10(FDR q.val)",
      y="",
      fill = "NES"
    )+
    scale_x_continuous(breaks=c(0,2.5, 5))+
    scale_y_discrete(position="right")+
    scale_fill_continuous(trans="reverse")+
    facet_grid2(~cell.type+cluster, 
                strip = strip_nested(bleed=TRUE, 
                                     background_x = elem_list_rect(fill = c(cell.type.col[i], 
                                                                            "white",
                                                                            "white"), 
                                                                   colour=c(cell.type.col[i], 
                                                                            "white",
                                                                            "white")),
                                     text_x = elem_list_text(face = "bold", size=6),
                                     by_layer_x = TRUE))+
    theme(panel.spacing.x = unit(4, "mm"))
  p.titles[i] <- paste0("figs_ge_barplot_genesetc8-",levels(annotations$cell.type)[i],".tiff")
}

size <- 5
ggsave(p.titles[[1]], plot=ind.plots[[1]], units="in", dpi=300, compression = 'lzw',
       width=size*1.5, height=size*0.9)
ggsave(p.titles[[2]], plot=ind.plots[[2]], units="in", dpi=300, compression = 'lzw',
       width=size*1.2, height=size*0.6)
ggsave(p.titles[[3]], plot=ind.plots[[3]], units="in", dpi=300, compression = 'lzw',
       width=size*1.5, height=size*0.5)
ggsave(p.titles[[4]], plot=ind.plots[[4]], units="in", dpi=300, compression = 'lzw',
       width=size*1.3, height=size*0.6)
ggsave(p.titles[[5]], plot=ind.plots[[5]], units="in", dpi=300, compression = 'lzw',
       width=size*1, height=size*0.9)
ggsave(p.titles[[6]], plot=ind.plots[[6]], units="in", dpi=300, compression = 'lzw',
       width=size*1, height=size*0.9)
ggsave(p.titles[[7]], plot=ind.plots[[7]], units="in", dpi=300, compression = 'lzw',
       width=size*1.3, height=size*0.5)
ggsave(p.titles[[8]], plot=ind.plots[[8]], units="in", dpi=300, compression = 'lzw',
       width=size*1.3, height=size*0.5)
ggsave(p.titles[[9]], plot=ind.plots[[9]], units="in", dpi=300, compression = 'lzw',
       width=size*1.2, height=size*0.5)
ggsave(p.titles[[10]], plot=ind.plots[[10]], units="in", dpi=300, compression = 'lzw',
       width=size*1.15, height=size*0.5)

sessionInfo()
# R version 4.2.2 (2022-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# locale:
#     [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
# [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# attached base packages:
#     [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#     [1] conflicted_1.1.0 dplyr_1.0.10     ggplot2_3.3.6   
# 
# loaded via a namespace (and not attached):
#     [1] rstudioapi_0.14  magrittr_2.0.3   tidyselect_1.2.0 munsell_0.5.0   
# [5] colorspace_2.0-3 R6_2.5.1         rlang_1.0.6      fastmap_1.1.0   
# [9] fansi_1.0.3      tools_4.2.2      grid_4.2.2       gtable_0.3.1    
# [13] utf8_1.2.2       cli_3.4.1        DBI_1.2.3        withr_2.5.0     
# [17] assertthat_0.2.1 tibble_3.1.8     lifecycle_1.0.3  crayon_1.5.2    
# [21] vctrs_0.5.0      cachem_1.0.6     memoise_2.0.1    glue_1.6.2      
# [25] compiler_4.2.2   pillar_1.8.1     generics_0.1.3   scales_1.2.1    
# [29] pkgconfig_2.0.3
