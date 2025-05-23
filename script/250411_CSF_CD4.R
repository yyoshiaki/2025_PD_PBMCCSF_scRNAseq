library(tidyverse)
library(readxl)
library(data.table)
library(ggbeeswarm)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(extrafont)
library(patchwork)

loadfonts()

setwd("~/pi_hafler/ASAP")

prefix.output <- "./output/250411_CSF_CD4T"

df.obs <- read_csv('./scanpy/250407_CSF_downstream/obs.csv')
df.gex <- read_csv('./data/241024_cellranger_files_GEX.csv') %>% filter(SampleType == 'CSF')
dis.order.row <- c("diseaseRBD-LowIntm", "diseaseRBD-High", "diseasePD", "diseasePD-RBD")

unique_rows <- df.obs %>% distinct(PD_Probability, DonorID, Age, Gender)

list_df <- list()
for (i in 1:nrow(unique_rows)) {
  row <- unique_rows[i, ]
  cat("Processing Sample:", row$DonorID, "\n")
  
  f.CD4 <- df.gex %>% 
    filter(DonorID == row$DonorID) %>% 
    mutate(file_CD4 = paste0('/home/yy693/pi_hafler/ASAP/output/CD4T_screfmapping/', directory, '/', cellranger_sample_id, '/', cellranger_sample_id, '_Reference_Mapping.csv')) %>%
    select(file_CD4) %>%
    as.vector()
  
  l_df <- list()
  for (f in f.CD4$file_CD4) {
    cat("Checking file:", f, "\n")
    if (file.exists(f)) {
      cat("File exists. Reading file...\n")
      d <- read_csv(f, col_types = cols()) %>%
        mutate(age = row$Age, sex = row$Gender, disease = row$PD_Probability, sample = row$DonorID)
      l_df[[length(l_df) + 1]] <- d
    } else {
      cat("skip", f, "File not found.\n")
    }
  }
  list_df[[length(list_df) + 1]] <- bind_rows(l_df)
}
df.output <- bind_rows(list_df)

print(df.output)
write_csv(df.output, paste0(prefix.output, '_tidy_queryL2_output.full.csv'))

### count
df.output <- 
  df.output %>% 
  filter(clusterL2 != 'nan') %>%
  filter(age != 'NaN')

df.output <- df.output %>%
  group_by(age, sex, disease, sample, clusterL2) %>%
  summarise(n = n(), .groups = "drop")

df.output <- df.output %>%
  group_by(age, sex, disease, sample) %>%
  mutate(freq = n / sum(n)) %>%
  mutate(age = age / 25)

df.output$sex <- str_replace(df.output$sex, "Female", "female")
df.output$sex <- str_replace(df.output$sex, "Male", "male")
# df.output$disease <- str_replace(df.output$disease, "Control", "HC")
# df.output$disease <- str_replace(df.output$disease, "PD\\+RBD", "PDRBD")
# df.output <- df.output %>% mutate(disease = relevel(factor(disease), ref = "HC"))
df.output <- df.output %>%
  mutate(disease = factor(disease, levels = c("HC", "RBD-LowIntm", "RBD-High", "PD", "PD-RBD"))) %>%
  mutate(disease = relevel(disease, ref = "HC"))

write_csv(df.output, paste0(prefix.output, '_tidy_queryL2_output.csv'))


p <- ggplot(df.output, aes(x=clusterL2, y=freq, color=disease), scale = "width") + 
  geom_quasirandom(dodge.width=1, size=0.5) +
  theme(axis.text.x = element_text(angle = 90))
p
ggsave(paste0(prefix.output, '_queryL2_ggbeeswarm.pdf'), width = 10, height = 4)

df.ntotal <- df.output %>% group_by(sample) %>% summarise(n_total = sum(n))

df.output.wide <- df.output %>% 
  pivot_wider(values_from = n, names_from = clusterL2,
              id_cols = c(sample, disease, age, sex), values_fill = 0) %>%
  left_join(df.ntotal, by='sample') %>%
  drop_na(age)

df.output.wide <- as.data.frame(df.output.wide)
df.output.wide <- within(df.output.wide, disease <- relevel(factor(disease), ref = "HC"))
df.output.wide <- within(df.output.wide, sex <- relevel(factor(sex), ref = "male"))

cells <- c("Tnaive", "TnaiveAct", "TnaiveMX1", "TnaiveSOX4","TcmTh0","TcmTh0Act", 
           "TcmTfh", "TcmPHLDA3", "TcmTh17", "TcmTh2",
           "TemTh1pre", "TemTh1", "TemTh117", "TemTph", "TemraTh1",
           "TregNaive", "TregAct", "TregEff")

list.res <- list()
for (cell in cells) {
  res <- glm(as.formula(paste0("cbind(", cell, ", n_total-", cell,") ~ disease + age + sex")),
             data = df.output.wide, family = binomial)
  # res <- glm(as.formula(paste0("cbind(", cell, ", n_total-", cell,") ~ disease + age + sex")),
  #            data = df.output.wide, family = binomial)
  list.res[[cell]] <- as_tibble(summary(res)$coefficients, rownames = "var") %>%
    mutate(cell = cell)
}

df.res <- rbindlist(list.res, fill = TRUE) %>% 
  as_tibble() %>%
  mutate(padj = p.adjust(`Pr(>|z|)`))

df.plot <- df.res %>%
  filter(str_detect(var, "^disease|^age|^sex")) %>%
  mutate(cell = factor(cell, levels = cells))

lim.est <- 1.8
lim.nest <- -1.8
lim.lpa <- 10^-20
lim.lpa.lower <- 0.05

df.plot <- df.plot %>% 
  mutate(Estimate=if_else(Estimate>lim.est, lim.est , Estimate)) %>% 
  mutate(Estimate=if_else(Estimate< lim.nest, -lim.est, Estimate)) %>%
  mutate(padj=if_else(padj<lim.lpa, lim.lpa, padj)) %>%
  mutate(var=fct_relevel(var, "sexfemale", "age"))

# order.row <- df.plot %>% group_by(var) %>% filter(`z value` == max(`z value`)) %>% arrange(cell) %>% select(var)
# order.row <- as.character(order.row$var)
# order.row <- order.row[!order.row %in% c("age", "sexfemale")]
# dis.order.row
order.row <- append(dis.order.row, c('age', 'sexfemale'))
# df.plot$var <- factor(df.plot$var, levels = rev(order.row))
df.plot$var <- str_replace(df.plot$var, 'disease', 'dis.')
df.plot$var <- factor(df.plot$var, levels = rev(order.row) %>% str_replace('disease', 'dis.'))
df.plot <- df.plot %>% mutate(padj = replace(padj, padj > lim.lpa.lower, 1))
ggplot(df.plot, aes(x= cell, y=var, size=-log10(padj), color=Estimate)) + 
  geom_point() + 
  # scale_color_gradient2(low = "blue", mid = "white",  high = "red", space = "Lab", limit = c(-2,2)) +
  # scale_fill_brewer(palette = "RdBu", limit = c(-2,2)) +
  # scale_color_gradient(palette = "RdBu", space = "Lab", limit = c(-2,2)) +
  scale_color_gradientn(colours = rev(brewer.pal(10, "RdBu")), limit = c(lim.nest,lim.est)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_size(limits = c(-log10(lim.lpa.lower),-log10(lim.lpa))) +
  theme_bw()
ggsave(paste0(prefix.output, '_queryL2_glm.pdf'), width = 6.6, height = 3)


write_csv(df.res, paste0(prefix.output, '_queryL2_glm.csv'))
write_csv(df.output.wide, paste0(prefix.output, '_queryL2_output.csv'))

list.res <- list()
for (cell in cells) {
  res <- glm(as.formula(paste0("cbind(", cell, ", n_total-", cell,") ~ disease + sex")),
             data = df.output.wide, family = binomial)
  # res <- glm(as.formula(paste0("cbind(", cell, ", n_total-", cell,") ~ disease + age + sex")),
  #            data = df.output.wide, family = binomial)
  list.res[[cell]] <- as_tibble(summary(res)$coefficients, rownames = "var") %>%
    mutate(cell = cell)
}

df.res <- rbindlist(list.res, fill = TRUE) %>% 
  as_tibble() %>%
  mutate(padj = p.adjust(`Pr(>|z|)`))

df.plot <- df.res %>%
  filter(str_detect(var, "^disease|^sex")) %>%
  mutate(cell = factor(cell, levels = cells))

lim.est <- 1.8
lim.nest <- -1.8
lim.lpa <- 10^-20
lim.lpa.lower <- 0.05

df.plot <- df.plot %>% 
  mutate(Estimate=if_else(Estimate>lim.est, lim.est , Estimate)) %>% 
  mutate(Estimate=if_else(Estimate< lim.nest, -lim.est, Estimate)) %>%
  mutate(padj=if_else(padj<lim.lpa, lim.lpa, padj)) %>%
  mutate(var=fct_relevel(var, "sexfemale"))

# order.row <- df.plot %>% group_by(var) %>% filter(`z value` == max(`z value`)) %>% arrange(cell) %>% select(var)
# order.row <- as.character(order.row$var)
# order.row <- order.row[!order.row %in% c("sexfemale")]
# order.row <- append(order.row, c('sexfemale'))
order.row <- append(dis.order.row, c('sexfemale'))
# df.plot$var <- factor(df.plot$var, levels = rev(order.row))
df.plot$var <- str_replace(df.plot$var, 'disease', 'dis.')
df.plot$var <- factor(df.plot$var, levels = rev(order.row) %>% str_replace('disease', 'dis.'))
df.plot <- df.plot %>% mutate(padj = replace(padj, padj > lim.lpa.lower, 1))
ggplot(df.plot, aes(x= cell, y=var, size=-log10(padj), color=Estimate)) + 
  geom_point() + 
  # scale_color_gradient2(low = "blue", mid = "white",  high = "red", space = "Lab", limit = c(-2,2)) +
  # scale_fill_brewer(palette = "RdBu", limit = c(-2,2)) +
  # scale_color_gradient(palette = "RdBu", space = "Lab", limit = c(-2,2)) +
  scale_color_gradientn(colours = rev(brewer.pal(10, "RdBu")), limit = c(lim.nest,lim.est)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_size(limits = c(-log10(lim.lpa.lower),-log10(lim.lpa))) +
  theme_bw()
ggsave(paste0(prefix.output, '_queryL2_glm_without_age.pdf'), width = 6.6, height = 3)


write_csv(df.res, paste0(prefix.output, '_queryL2_glm_without_age.csv'))
write_csv(df.output.wide, paste0(prefix.output, '_queryL2_output.csv'))

#### without gender, sex
list.res <- list()
for (cell in cells) {
  res <- glm(as.formula(paste0("cbind(", cell, ", n_total-", cell,") ~ disease")),
             data = df.output.wide, family = binomial)
  # res <- glm(as.formula(paste0("cbind(", cell, ", n_total-", cell,") ~ disease + age + sex")),
  #            data = df.output.wide, family = binomial)
  list.res[[cell]] <- as_tibble(summary(res)$coefficients, rownames = "var") %>%
    mutate(cell = cell)
}

df.res <- rbindlist(list.res, fill = TRUE) %>% 
  as_tibble() %>%
  mutate(padj = p.adjust(`Pr(>|z|)`))

df.plot <- df.res %>%
  filter(str_detect(var, "^disease|^sex")) %>%
  mutate(cell = factor(cell, levels = cells))

lim.est <- 1.8
lim.nest <- -1.8
lim.lpa <- 10^-20
lim.lpa.lower <- 0.05

df.plot <- df.plot %>% 
  mutate(Estimate=if_else(Estimate>lim.est, lim.est , Estimate)) %>% 
  mutate(Estimate=if_else(Estimate< lim.nest, -lim.est, Estimate)) %>%
  mutate(padj=if_else(padj<lim.lpa, lim.lpa, padj)) 
# %>%
#   mutate(var=fct_relevel(var, "sexfemale"))

# order.row <- df.plot %>% group_by(var) %>% filter(`z value` == max(`z value`)) %>% arrange(cell) %>% select(var)
# order.row <- as.character(order.row$var)
# order.row <- order.row[!order.row %in% c("sexfemale")]
# order.row <- append(order.row, c('sexfemale'))
order.row <- dis.order.row
# df.plot$var <- factor(df.plot$var, levels = rev(order.row))
df.plot$var <- str_replace(df.plot$var, 'disease', 'dis.')
df.plot$var <- factor(df.plot$var, levels = rev(order.row) %>% str_replace('disease', 'dis.'))
df.plot <- df.plot %>% mutate(padj = replace(padj, padj > lim.lpa.lower, 1))
ggplot(df.plot, aes(x= cell, y=var, size=-log10(padj), color=Estimate)) + 
  geom_point() + 
  # scale_color_gradient2(low = "blue", mid = "white",  high = "red", space = "Lab", limit = c(-2,2)) +
  # scale_fill_brewer(palette = "RdBu", limit = c(-2,2)) +
  # scale_color_gradient(palette = "RdBu", space = "Lab", limit = c(-2,2)) +
  scale_color_gradientn(colours = rev(brewer.pal(10, "RdBu")), limit = c(lim.nest,lim.est)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_size(limits = c(-log10(lim.lpa.lower),-log10(lim.lpa))) +
  theme_bw()
ggsave(paste0(prefix.output, '_queryL2_glm_without_agegender.pdf'), width = 6.6, height = 3)


write_csv(df.res, paste0(prefix.output, '_queryL2_glm_without_agegender.csv'))
write_csv(df.output.wide, paste0(prefix.output, '_queryL2_output.csv'))

df.plot <- df.output %>%
  filter(clusterL2 == "TnaiveAct") %>%
  mutate(age = age * 25)

ggplot(df.plot, aes(x = age, y = freq, color = disease, size = n)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, aes(group = disease)) + 
  scale_size_continuous(range = c(3, 12)) +
  labs(x = "Age", y = "Frequency", 
       color = "Disease", size = "Count n") +
  theme_minimal() +
  theme(legend.position = "right") 

ggsave(paste0(prefix.output, '_queryL2_TnaiveAct_age.pdf'), width = 7, height = 7)

for (cell in cells) {
  df.plot <- df.output %>%
    filter(clusterL2 == cell) %>%
    mutate(age = age * 25)
  
  plot <- ggplot(df.plot, aes(x = age, y = freq, color = disease, size = n)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, aes(group = disease)) + 
    scale_size_continuous(range = c(1, 8)) +
    labs(x = "Age", y = "Frequency", 
         color = "Disease", size = "Count n") +
    ggtitle(cell) +
    theme_minimal() +
    theme(legend.position = "right")
  
  ggsave(paste0(prefix.output, '_queryL2_', cell, '_age.pdf'), plot = plot, width = 5, height = 5)
}


### NMF
cells <- c("Tnaive", "TnaiveAct", "TnaiveMX1", "TnaiveSOX4","TcmTh0","TcmTh0Act", 
           "TcmTfh", "TcmPHLDA3", "TcmTh17", "TcmTh2",
           "TemTh1pre", "TemTh1", "TemTh117", "TemTph", "TemraTh1",
           "TregNaive", "TregAct", "TregEff")

list_df <- list()
for (i in 1:nrow(unique_rows)) {
  row <- unique_rows[i, ]
  cat("Processing Sample:", row$DonorID, "\n")
  
  files_NMF <- df.gex %>% 
    filter(DonorID == row$DonorID) %>% 
    mutate(file_CD4 = paste0('/home/yy693/pi_hafler/ASAP/output/CD4T_screfmapping/', directory, '/', cellranger_sample_id, '/', cellranger_sample_id, '_CD4T_AssayData_projection.csv')) %>%
    select(file_CD4) %>%
    as.vector()
  
  files_refmap <- df.gex %>% 
    filter(DonorID == row$DonorID) %>% 
    mutate(file_CD4 = paste0('/home/yy693/pi_hafler/ASAP/output/CD4T_screfmapping/', directory, '/', cellranger_sample_id, '/', cellranger_sample_id, '_Reference_Mapping.csv')) %>%
    select(file_CD4) %>%
    as.vector()
  
  l_df <- list()
  for (i in length(files_NMF$file_CD4)) {
    f_NMF <- files_NMF$file_CD4[i]
    f_refmap <- files_refmap$file_CD4[i]
    cat("Checking file:", f_NMF, "\n")
    if (file.exists(f_NMF)) {
      cat("File exists. Reading file...\n")
      d_NMF <- read_csv(f_NMF, col_names = TRUE) %>%
        column_to_rownames(var = "...1")
      d_NMF <- as.data.frame(t(d_NMF)) %>%
        rownames_to_column(var = "Cell")
      d_refmap <- read_csv(f_refmap, col_types = cols()) %>%
        column_to_rownames(var = "...1") %>%
        rownames_to_column(var = "Cell")
      d_NMF <- d_NMF %>% left_join(d_refmap, by= join_by(Cell))
      d_NMF <- d_NMF %>%
        mutate(age = row$Age, sex = row$Gender, disease = row$PD_Probability, sample = row$DonorID)
      
      l_df[[length(l_df) + 1]] <- d_NMF
    } else {
      cat("skip", f, "File not found.\n")
    }
  }
  list_df[[length(list_df) + 1]] <- bind_rows(l_df)
}

df.output <- bind_rows(list_df)

print(df.output)

df.output <- 
  df.output %>% 
  filter(clusterL2 != 'nan') %>%
  filter(age != 'NaN')

df.output <- df.output %>%
  mutate(age = age / 25)

df.output$sex <- str_replace(df.output$sex, "Female", "female")
df.output$sex <- str_replace(df.output$sex, "Male", "male")
# df.output$disease <- str_replace(df.output$disease, "Control", "HC")
# df.output$disease <- str_replace(df.output$disease, "PD\\+RBD", "PDRBD")
# df.output <- df.output %>% mutate(disease = relevel(factor(disease), ref = "HC"))
df.output <- df.output %>%
  mutate(disease = factor(disease, levels = c("HC", "RBD-LowIntm", "RBD-High", "PD", "PD-RBD"))) %>%
  mutate(disease = relevel(disease, ref = "HC"))
df.output <- within(df.output, sex <- relevel(factor(sex), ref = "male"))

write_csv(df.output, paste0(prefix.output, '_tidy_NMFL2_output.csv'))

df.output %>% filter(clusterL2=="TnaiveMX1") %>% ggplot(aes(NMF_7)) + geom_histogram()
df.output %>% filter(clusterL2=="Tnaive") %>% ggplot(aes(NMF_7)) + geom_histogram()

df.plot <- df.output %>% group_by(clusterL1, disease) %>%
  summarise(NMF_0=mean(NMF_0), NMF_1=mean(NMF_1), NMF_2=mean(NMF_2),
            NMF_3=mean(NMF_3), NMF_4=mean(NMF_4), NMF_5=mean(NMF_5),
            NMF_6=mean(NMF_6), NMF_7=mean(NMF_7), NMF_8=mean(NMF_8),
            NMF_9=mean(NMF_9), NMF_10=mean(NMF_10), NMF_11=mean(NMF_11)) %>%
  mutate(description = paste0(clusterL1,'_',disease))

col.labels <- c('NMF0 Cytotoxic-F', 'NMF1 Treg-F', 'NMF2 Th17-F', 'NMF3 Naive-F', 
                'NMF4 Act-F', 'NMF5 TregEff/Th2-F', 'NMF6 Tfh-F', 'NMF7 IFN-F', 'NMF8 Cent.Mem.-F',
                'NMF9 Thy.Emi.-F', 'NMF10 Tissue-F', 'NMF11 Th1-F')
col.factors <- c('NMF_0', 'NMF_1', 'NMF_2',
                 'NMF_3', 'NMF_4', 'NMF_5',
                 'NMF_6', 'NMF_7', 'NMF_8',
                 'NMF_9', 'NMF_10', 'NMF_11')

p <- pheatmap(df.plot[col.factors],
              # color = cividis(100),
              cluster_cols = FALSE, cluster_rows = FALSE,
              labels_row = df.plot$description, labels_col = col.labels,
              scale = "column")
p

pdf(paste0(prefix.output, '_NMF_clusterL1_heatmap.pdf'), width = 6, height = 10)
p
dev.off()

df.plot <- df.output %>% group_by(clusterL2, disease) %>%
  summarise(NMF_0=mean(NMF_0), NMF_1=mean(NMF_1), NMF_2=mean(NMF_2),
            NMF_3=mean(NMF_3), NMF_4=mean(NMF_4), NMF_5=mean(NMF_5),
            NMF_6=mean(NMF_6), NMF_7=mean(NMF_7), NMF_8=mean(NMF_8),
            NMF_9=mean(NMF_9), NMF_10=mean(NMF_10), NMF_11=mean(NMF_11)) %>%
  mutate(description = paste0(clusterL2,'_',disease))

p <- pheatmap(df.plot[c('NMF_0', 'NMF_1', 'NMF_2',
                        'NMF_3', 'NMF_4', 'NMF_5',
                        'NMF_6', 'NMF_7', 'NMF_8',
                        'NMF_9', 'NMF_10', 'NMF_11')],
              # color = cividis(100),
              cluster_cols = FALSE, cluster_rows = FALSE,
              labels_row = df.plot$description, labels_col = col.labels,
              scale = "column")
p

pdf(paste0(prefix.output, '_NMF_clusterL2_heatmap.pdf'), width = 6, height = 32)
p
dev.off()


list.res <- list()
for (cell in cells) {
  for (comp in colnames(df.output) %>% str_subset("NMF_")){
    # d <- df.output %>% filter(clusterL2 == cell) %>% 
    #   group_by(sample) %>% 
    #   summarise(NMF = mean(.data[[comp]]), disease=first(disease), 
    #             age=first(age), sex=first(sex), project=first(project))
    # res <- glm(as.formula("NMF ~ disease + age + sex + project"),
    #            data = d, family = Gamma)
    # res <- glm(as.formula(paste0(comp, " ~ disease + age + sex + project")),
    #            data = df.output %>% filter(clusterL2 == cell))
    res <- glm(as.formula(paste0(comp, " ~ disease")),
               data = df.output %>% filter(clusterL2 == cell))
    list.res[[paste(cell, comp)]] <- as_tibble(summary(res)$coefficients, rownames = "var") %>%
      mutate(cell = cell, component = comp)
  }
}

df.res <- rbindlist(list.res, fill = TRUE) %>% 
  as_tibble() %>%
  mutate(padj = p.adjust(`Pr(>|t|)`))


lim.est <- 0.3
lim.nest <- -0.3
lim.lpa <- 10^-100
lim.scaledEstimate <- 3

for (c in cells) {
  df.plot <- df.res %>%
    filter(cell == c) %>%
    filter(str_detect(var, "^disease|^age|^sex")) %>%
    mutate(cell = factor(cell, levels = cells))
  
  df.plot <- df.plot %>% 
    mutate(Estimate=if_else(Estimate>lim.est, lim.est , Estimate)) %>% 
    mutate(Estimate=if_else(Estimate< lim.nest, -lim.est, Estimate)) %>%
    mutate(padj=if_else(padj<lim.lpa, lim.lpa, padj)) %>% 
    mutate(padj=if_else(padj>0.1, 0, padj))
  
  df.plot <- df.plot %>% mutate(component=as_factor(component)) %>% 
    mutate(component=fct_relevel(component, col.factors)) %>%
    mutate(component=fct_recode(component, 
                                'NMF0 Cytotoxic-F'='NMF_0', 'NMF1 Treg-F'='NMF_1', 
                                'NMF2 Th17-F'='NMF_2', 'NMF3 Naive-F'='NMF_3', 
                                'NMF4 Act-F'='NMF_4', 'NMF5 TregEff/Th2-F'='NMF_5', 
                                'NMF6 Tfh-F'='NMF_6', 'NMF7 IFN-F'='NMF_7', 
                                'NMF8 Cent.Mem.-F'='NMF_8', 
                                'NMF9 Thy.Emi.-F'='NMF_9', 
                                'NMF10 Tissue-F'='NMF_10',
                                'NMF11 Th1-F'='NMF_11')) %>%
    mutate(var=fct_relevel(var, "sexfemale", "age"))
  
  # order.row <- c('diseaseRBD_Hyposmia_neg', 'diseaseRBD_Hyposmia_pos', 'diseasePD', 'diseasePD-RBD', 'sexfemale', 'age')
  order.row <- dis.order.row
  df.plot$var <- str_replace(df.plot$var, 'disease', 'dis.')
  df.plot$var <- factor(df.plot$var, levels = rev(order.row) %>% str_replace('disease', 'dis.'))
  
  df.plot.int <- df.res %>%
    filter(var == "(Intercept)") %>%
    mutate(component=as_factor(component)) %>% 
    mutate(component=fct_relevel(component, col.factors)) %>%
    mutate(component=fct_recode(component, 
                                'NMF0 Cytotoxic-F'='NMF_0', 'NMF1 Treg-F'='NMF_1', 
                                'NMF2 Th17-F'='NMF_2', 'NMF3 Naive-F'='NMF_3', 
                                'NMF4 Act-F'='NMF_4', 'NMF5 TregEff/Th2-F'='NMF_5', 
                                'NMF6 Tfh-F'='NMF_6', 'NMF7 IFN-F'='NMF_7', 
                                'NMF8 Cent.Mem.-F'='NMF_8', 
                                'NMF9 Thy.Emi.-F'='NMF_9', 
                                'NMF10 Tissue-F'='NMF_10',
                                'NMF11 Th1-F'='NMF_11')) %>%
    group_by(component) %>%
    mutate(z_scaled_Estimate = as.numeric(scale(Estimate))) %>%
    ungroup() %>%
    mutate(z_scaled_Estimate=if_else(z_scaled_Estimate>lim.scaledEstimate, lim.scaledEstimate, z_scaled_Estimate)) %>% 
    # mutate(z_scaled_Estimate=if_else(z_scaled_Estimate<-lim.scaledEstimate, -lim.scaledEstimate, z_scaled_Estimate)) %>% 
    filter(cell == c) %>%
    mutate(cell = factor(cell, levels = cells))
  
  if (dim(df.plot)[1] > 0){
    # Create the dotplot
    dotplot <- ggplot(df.plot, aes(x= component, y=var, size=-log10(padj), color=Estimate)) + 
      geom_point() +
      scale_color_gradientn(colours = rev(brewer.pal(10, "RdBu")), limit = c(lim.nest,lim.est)) +
      scale_x_discrete(guide = guide_axis(angle = 90)) +
      scale_size(range = c(0, 4), limits = c(0,100)) + 
      xlab("") + ylab("") +
      labs(color = "Coef.") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, colour = "black"),
            axis.text = element_text(color = "black"),
            axis.title = element_text(color = "black"),
            legend.text = element_text(color = "black"),
            legend.title = element_text(color = "black"))
    
    # # Save the dotplot
    # ggsave(paste0(prefix.output, '_NMF_queryL2', c, '_glm_withproj.pdf'), width = 6, height = 5)
    
    # Create the heatmap
    heatmap <- ggplot(df.plot.int, aes(x = component, y = var, fill = z_scaled_Estimate)) +
      geom_tile() +
      ggtitle(c) +
      scale_fill_gradientn(colours = rev(brewer.pal(10, "RdGy")), limit = c(-lim.scaledEstimate,lim.scaledEstimate)) +
      theme_minimal() +
      labs(fill = "Intercept") +
      theme(
        axis.title.x = element_blank(), # Remove x-axis title
        axis.title.y = element_blank(), # Remove y-axis title
        axis.text.x = element_blank(),  # Remove x-axis text
        axis.text.y = element_blank(),  # Remove y-axis text
        panel.grid.major = element_blank(), # Remove major grid
        panel.grid.minor = element_blank(), # Remove minor grid
        legend.position = "right",
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black")
      )
    
    # Combine the plots
    heatmap / dotplot + plot_layout(heights = c(1, 8))
    ggsave(paste0(prefix.output, '_NMF_queryL2', c, '_glm_withscaledintercept.pdf'), width = 6, height = 4)
  }
}

write_csv(df.res, paste0(prefix.output, '_NMF_queryL2_glm.csv'))


