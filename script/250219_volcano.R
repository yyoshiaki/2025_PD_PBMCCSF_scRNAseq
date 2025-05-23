library(tidyverse)
library(ggrepel)
library(ggrastr)
library(ggh4x)
library(EnhancedVolcano)

th_padj = 1e-2
th_padj_plot = 1e-200
th_lfc = 0.2
th_mean = 0.5
n_top = 50

res <- read_csv('./scanpy/241028_CSF_downstream/graph/deg_CSF Mac_group_RBD.csv')
res <- res %>% filter(mean > th_mean)

genes <- c("TNFRSF1A", 'TNFRSF1B', "CIITA", "HLA-DQA2", 
           "HLA-DQB2", 'HIF1A', 'SRGAP2', 'CSF2RA',
           'HLA-DQA1', 'LYZ', 'NEAT1', 'TLR2', 'LILRA6',
           'MBNL1', 'RUNX1', 'IL17RA', 'TGFBR2', 'PTPRC', 
           'CEBPD', 'IL10RA', 'JAK2', 'CD86', 'CSF1R',
           'ITGA4', 'LRRK2', 'FOS', 'JUNB', 'MAF', 'HLA-B', 'ACTB', 'IFITM3', 'SPP1', 'AIF1',
           'B2M')

EnhancedVolcano(res,
                lab = res$names,
                x = 'logfoldchanges',
                y = 'pvals_adj',
                selectLab = genes,
                title = 'CSF Mac RBD vs HC',
                subtitle = "",
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = th_padj,
                FCcutoff = th_lfc,
                pointSize = 2.0,
                labSize = 10.0,
                # ylim = c(0, -log10(10e-200)),
                labCol = 'black',
                col=c('black', 'black', 'black', 'red3'),
                boxedLabels = FALSE,
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                arrowheads = FALSE,
                widthConnectors = 1.0,
                colConnectors = 'black')
ggsave(filename = './scanpy/241028_CSF_downstream/graph/enhvolcano_CSFMac_RBD.pdf', width = 10, height = 10)
