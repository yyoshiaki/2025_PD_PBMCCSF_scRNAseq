library("Startrac")
library("tidyverse")
library("data.table")
library("tictoc")
library("ggpubr")
library("ggplot2")
library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")
library(ggpubr)

# setwd('/home/yyasumizu/media32TB/bioinformatics/drSaito/')

# in.dat <- read.csv('./data/startrac.input.v2.csv')

in.dat <- read.csv('./scanpy/241219_startrac/startrac_input_T_crosstissue.csv')
head(in.dat)

out <- Startrac.run(in.dat, proj="PD",verbose=F)

Startrac::plot(out, index.type="cluster.all",byPatient=T)

p <- Startrac::plot(out,index.type="cluster.all",byPatient=T)
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename = './scanpy/241219_startrac/graph/crosstissue_startrac_cluster_all.pdf', width = 8, height = 6)

p <- Startrac::plot(out,index.type="pairwise.migr",byPatient=T)
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename = './scanpy/241219_startrac/graph/crosstissue_startrac_cluster_all_migr.pdf', width = 8, height = 6)

pdf('./scanpy/241219_startrac/graph/crosstissue_startrac_cluster_all_trans.pdf', width = 7, height = 6)
Startrac::plot(out,index.type="pairwise.tran",byPatient=T)
dev.off()

obj <- out
df <- as.data.table(obj@cluster.sig.data)[aid!=obj@proj,][order(majorCluster),]
in.dat_unique <- unique(in.dat[, c("patient", "disease", "disease_rbd")])

merged_df <- merge(df, 
                   in.dat_unique, 
                   by.x = "aid", 
                   by.y = "patient", 
                   all.x = TRUE)
merged_df$disease <- factor(merged_df$disease, levels=c("HC", 'RBD', 'PD', 'PD-RBD'))
head(merged_df)

ggboxplot(df,
          x="majorCluster",y="value",palette = "npg",
          color = "index", add = "point", outlier.colour=NULL) +
  facet_wrap(~index,ncol=1,scales = "free_y") +
  theme(axis.text.x=element_text(angle = 60,hjust = 1))


ggboxplot(merged_df,
          x = "majorCluster",
          y = "value",
          palette = "npg",
          color = "disease",
          add = "point",
          outlier.colour = NULL) +
  facet_wrap(~ index, ncol = 1, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(filename = './scanpy/241219_startrac/graph/crosstissue_startrac_cluster_all_strat_disease.pdf', width = 8, height = 8)


ggboxplot(merged_df %>% filter(index == 'expa'),
          x="disease", y="value", palette = "npg",
          color = "disease", add = "point", outlier.colour=NULL, facet.by = "majorCluster") +
  # facet_wrap(~inde,ncol=1,scales = "free_y") +
  labs(y = 'Expansion') +
  theme(axis.text.x=element_text(angle = 60,hjust = 1)) +
  stat_compare_means(comparisons = list(c("HC", "PD"), c("HC", "RBD"), c("HC", "PD-RBD")), 
                     method = "t.test",
                     label = "p.signif",
                     hide.ns = FALSE)
ggsave(filename = './scanpy/241219_startrac/graph/crosstissue_startrac_cluster_all_strat_disease.stat.expa.pdf', width = 8, height = 8)

ggboxplot(merged_df %>% filter(index == 'tran'),
          x="disease", y="value", palette = "npg",
          color = "disease", add = "point", outlier.colour=NULL, facet.by = "majorCluster") +
  # facet_wrap(~inde,ncol=1,scales = "free_y") +
  labs(y = 'Transition') +
  theme(axis.text.x=element_text(angle = 60,hjust = 1)) +
  stat_compare_means(comparisons = list(c("HC", "PD"), c("HC", "RBD"), c("HC", "PD-RBD")), 
                     method = "t.test",
                     label = "p.signif",
                     hide.ns = FALSE)
ggsave(filename = './scanpy/241219_startrac/graph/crosstissue_startrac_cluster_all_strat_disease.stat.tran.pdf', width = 8, height = 8)


ggboxplot(merged_df %>% filter(index == 'gini'),
          x="disease", y="value", palette = "npg",
          color = "disease", add = "point", outlier.colour=NULL, facet.by = "majorCluster") +
  # facet_wrap(~inde,ncol=1,scales = "free_y") +
  labs(y = 'Gini') +
  theme(axis.text.x=element_text(angle = 60,hjust = 1)) +
  stat_compare_means(comparisons = list(c("HC", "PD"), c("HC", "RBD"), c("HC", "PD-RBD")), 
                     method = "t.test",
                     label = "p.signif",
                     hide.ns = FALSE)
ggsave(filename = './scanpy/241219_startrac/graph/crosstissue_startrac_cluster_all_strat_disease.stat.gini.pdf', width = 8, height = 8)

ggboxplot(merged_df %>% filter(index == 'migr'),
          x="disease", y="value", palette = "npg",
          color = "disease", add = "point", outlier.colour=NULL, facet.by = "majorCluster") +
  # facet_wrap(~inde,ncol=1,scales = "free_y") +
  labs(y = 'Migration') +
  theme(axis.text.x=element_text(angle = 60,hjust = 1)) +
  stat_compare_means(comparisons = list(c("HC", "PD"), c("HC", "RBD"), c("HC", "PD-RBD")), 
                     method = "t.test",
                     label = "p.signif",
                     hide.ns = FALSE) 
ggsave(filename = './scanpy/241219_startrac/graph/crosstissue_startrac_cluster_all_strat_disease.stat.migr.pdf', width = 8, height = 8)


obj <- out
df <- as.data.table(obj@pIndex.migr)[aid!=obj@proj,][order(majorCluster),]
in.dat_unique <- unique(in.dat[, c("patient", "disease", "disease_rbd")])

merged_df <- merge(df, 
                   in.dat_unique, 
                   by.x = "aid", 
                   by.y = "patient", 
                   all.x = TRUE)

head(merged_df)
merged_df$disease <- factor(merged_df$disease, levels=c("HC", 'RBD', 'PD', 'PD-RBD'))
ggboxplot(merged_df,
          x="disease",y="`BLD-CSF`",palette = "npg",
          color = "disease", add = "point", outlier.colour=NULL, facet.by = "majorCluster") +
  # facet_wrap(~inde,ncol=1,scales = "free_y") +
  theme(axis.text.x=element_text(angle = 60,hjust = 1)) +
  stat_compare_means(comparisons = list(c("HC", "PD"), c("HC", "RBD"), c("HC", "PD-RBD")), 
                     method = "t.test",
                     label = "p.signif",
                     hide.ns = FALSE) 

ggsave(filename = './scanpy/241219_startrac/graph/crosstissue_startrac_cluster_all_migr_strat_disease.pdf', width = 8, height = 6)
