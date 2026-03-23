library("readxl")
library("reshape2")
library("irlba")
library("factoextra")
library("ggplot2")
library("cowplot")
library("Rtsne")
library("pheatmap")
library("tidyverse")
library("car")
library("ggpubr")
library("limma")
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")

### FUNCTIONS ###

## FUNCTIONS: Explorative analysis

pca_graph <- function(pca_data, panel, metadata = NULL, group = NULL) { #Section (I).2: PCA
  pca_scores <- as.data.frame(pca_data$x[, 1:2])
  pca_scores$Sample <- rownames(pca_scores)
  rownames(pca_scores) <- NULL
  # Adding groups
  pca_scores <- pca_scores %>% left_join(metadata[, c("Sample", group)], by = "Sample")  #Sorting samples and merging with metadata
  
  # Graph: PC1 vs PC2
  p<- ggplot(pca_scores, aes(x = PC1, y = PC2, color = .data[[group]])) +
        geom_point(size = 3) +                 
        stat_ellipse(aes(fill = .data[[group]]), alpha = 0.2, geom = "polygon") +
        labs(title = sprintf("PCA - %s panel for %s", panel, group),
             x = paste0("PC1 (", round(100 * summary(pca_data)$importance[2,1], 1), "% of variance)"),
             y = paste0("PC2 (", round(100 * summary(pca_data)$importance[2,2], 1), "% of variance)")) +
        theme_minimal() +
        theme(legend.title = element_blank(),
             plot.title = element_text(hjust = 0.5))
  return(p)
}


tsne_graph <- function (tsne_data, matrix, panel, metadata, group) { #Section (I).3: t-SNE
  tsne_scores <- as.data.frame(tsne_data$Y)
  colnames(tsne_scores) <- c("TSNE1", "TSNE2")
  
  tsne_scores$Sample <- rownames(matrix) # Reattaching Sample index
  
  tsne_scores <- tsne_scores %>% left_join(metadata[, c("Sample", group)], by = "Sample")  #Sorting samples and merging with metadata
  #Graph:
  gg<- ggplot(tsne_scores, aes(TSNE1, TSNE2, color = .data[[group]])) +
    geom_point(size = 3) +
    stat_ellipse(aes(fill = .data[[group]]), alpha = 0.2, geom = "polygon") +
    theme_minimal() +
    labs(title = sprintf("t-SNE for %s - %s panel", group, panel),
         x = "t-SNE 1",
         y = "t-SNE 2")
  return(gg)
} 


heatmap_cluster <- function (heatmap_matrix, panel, metadata) { #Section (I).4: Heatmap
  Conditions <- data.frame(
    Group  = metadata$Group,
    Gender = metadata$Gender,
    Age    = metadata$Age
  )
  rownames(Conditions) <- metadata$Sample
  
  pheat <- pheatmap(heatmap_matrix,
                    cluster_rows = FALSE,
                    cluster_cols = TRUE,
                    annotation_col = Conditions,
                    show_rownames = FALSE,
                    show_colnames = TRUE,
                    main = sprintf("Heatmap - %s panel", panel))
}


boxplot_distr <- function(NPX_long, variable, panel) { #Section (I).5: boxplot
  
  title_str <- sprintf("NPX distribution for %s - %s panel", variable, panel)
  gg<- ggplot(NPX_long, aes(x = .data[[variable]], y = NPX)) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = title_str)
  
  return(gg)
}


## FUNCTIONS: Differential Expression analysis ##


plot_violin <- function(NPX_data_full, protein_name, sig_pvalues_df, test_type) { #Section (II).3: violin plots
  # Trasforming data
  NPX_data <- NPX_data_full %>% filter(Assay == protein_name)
  sig_df_prot <- sig_pvalues_df %>% filter(Assay == protein_name & p_value < 0.05) %>% mutate(p_label = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    TRUE ~ "*"
  ))
  # If there are multiple significant comparisons
  if(nrow(sig_df_prot) > 0) {
    sig_df_prot <- sig_df_prot %>% 
      mutate(y.position = max(NPX_data$NPX) * 1.05 + (row_number() - 1) * 0.1 * max(NPX_data$NPX))
  }
  
  # Plotting:
  p <- ggplot(NPX_data, aes(x = Group, y = NPX, fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.7) +
    theme_minimal() +
    labs(y = "NPX", x = "Group",
         title = sprintf("Violin plot: %s - %s test", protein_name, test_type))
  # Adding significance bars
  if(nrow(sig_df_prot) > 0) {
    p <- p + stat_pvalue_manual(
      sig_df_prot,
      label = "p_label",
      tip.length = 0.02,
      label.size = 4
    )
  }
  return(p)
}


grid_graph <- function(plot_list, plots_per_page, row, col) { #Section (II).3: violin plots
  n_pages <- ceiling(length(plot_list) / plots_per_page)
  for (page in 1:n_pages) {
    # Batch indexes
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, length(plot_list))
    plots_batch <- plot_list[start_idx:end_idx]
    
    # Creating grid and printing:
    grid <- plot_grid(plotlist = plots_batch, nrow = row, ncol = col)
    print(grid)
  }
}


volcano_graph_limma <- function(tot_results, g1, g2) {
  
  volcano_df <- tot_results %>%
    filter(group1 == !!g1, group2 == !!g2) %>%
    mutate(
      neg_log10_p = -log10(p_value),
      signif = p_value < 0.05
    )
  
  p <- ggplot(volcano_df, aes(x = logFC, y = neg_log10_p)) +
    geom_point(aes(color = signif), alpha = 0.8) +
    scale_color_manual(values = c("grey60", "red")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    theme_minimal() +
    labs(title = sprintf("Volcano plot (limma): %s vs %s", g1, g2),
         x = "log2 Fold Change", y = "-log10(p-value)", color = "Significant")
  return(p)
}


### (I) CODE START: Explorative analysis ###

# All outputs (graphs) will be printed on this file.
pdf("Output Graphs.pdf", width = 20, height = 15)


## 1. Importing tha data, changing the format, filtering by LOD ##

# This is the path pointing the folder where the data is
Path = "/your_data_path/" 

#Cardiometabolic (Cm)
Cm_NPX_wide <- read_excel(paste0(Path, "Cardiometabolic_NPX values.xlsx"))
Cm_NPX_wide <- Cm_NPX_wide[, -ncol(Cm_NPX_wide)]
#Cm LODs
LOD_Cm <- read_excel(paste0(Path, "Cardiometabolic_NPX values.xlsx"), sheet = 2)
LOD_Cm <- as.data.frame(t(LOD_Cm))
colnames(LOD_Cm) <- LOD_Cm[1, ]
LOD_Cm <- LOD_Cm[-1,  ,drop = FALSE]


#Cardiovascular II (Cv)
Cv_NPX_wide <- read_excel(paste0(Path, "Cardiovascular II_NPX values.xlsx"))
Cv_NPX_wide <- Cv_NPX_wide[, -ncol(Cv_NPX_wide)]
#Cv LODs
LOD_Cv <- read_excel(paste0(Path, "Cardiovascular II_NPX values.xlsx"), sheet = 2)
LOD_Cv <- as.data.frame(t(LOD_Cv))
colnames(LOD_Cv) <- LOD_Cv[1, ]
LOD_Cv <- LOD_Cv[-1,  ,drop = FALSE]


#Immuno Oncology (Io)
#!There is an additional sample here, C2!
Io_NPX_wide <- read_excel(paste0(Path, "Immuno Oncology_NPX values.xlsx"))
Io_NPX_wide <- Io_NPX_wide[, -ncol(Io_NPX_wide)]
#Io LODs
LOD_Io <- read_excel(paste0(Path, "Immuno Oncology_NPX values.xlsx"), sheet = 2)
LOD_Io <- as.data.frame(t(LOD_Io))
colnames(LOD_Io) <- LOD_Io[1, ]
LOD_Io <- LOD_Io[-1,  ,drop = FALSE]

#Converting formats into "long"
Cm_NPX_long_full <- melt(Cm_NPX_wide, id="Sample",  variable.name = "Assay", value.name = "NPX")
Cv_NPX_long_full <- melt(Cv_NPX_wide, id="Sample",  variable.name = "Assay", value.name = "NPX")
Io_NPX_long_full <- melt(Io_NPX_wide, id="Sample",  variable.name = "Assay", value.name = "NPX")

#Filtering proteins with 80% if values below LOD
#a.Cardiometabolic
LOD_Cm$Assay <- rownames(LOD_Cm)
Cm_NPX_long <- merge(Cm_NPX_long_full, LOD_Cm, by = "Assay")
det_rate <- Cm_NPX_long %>% group_by(Assay) %>% summarise(detection_rate = mean(NPX > LOD, na.rm = TRUE))
good_assays_Cm <- det_rate %>% filter(detection_rate >= 0.8) %>% pull(Assay)
Cm_NPX_long <- Cm_NPX_long %>% filter(Assay %in% good_assays_Cm)
#b.Cardiovascular
LOD_Cv$Assay <- rownames(LOD_Cv)
Cv_NPX_long <- merge(Cv_NPX_long_full, LOD_Cv, by = "Assay")
det_rate <- Cv_NPX_long %>% group_by(Assay) %>% summarise(detection_rate = mean(NPX > LOD, na.rm = TRUE))
good_assays_Cv <- det_rate %>% filter(detection_rate >= 0.8) %>% pull(Assay)
Cv_NPX_long <- Cv_NPX_long %>% filter(Assay %in% good_assays_Cv)
#c.Immuno Oncology
LOD_Io$Assay <- rownames(LOD_Io)
Io_NPX_long <- merge(Io_NPX_long_full, LOD_Io, by = "Assay")
det_rate <- Io_NPX_long %>% group_by(Assay) %>% summarise(detection_rate = mean(NPX > LOD, na.rm = TRUE))
good_assays_Io <- det_rate %>% filter(detection_rate >= 0.8) %>% pull(Assay)
Io_NPX_long <- Io_NPX_long %>% filter(Assay %in% good_assays_Io)




## 2. Explorative Analysis: PCA ##

#Importing metadata
metadata <- read_excel(paste0(Path, "Metadata.xlsx"))

## Normalization (Using wide formats) ##
#Checking there aren't na values
colSums(is.na(Cm_NPX_wide))
colSums(is.na(Cv_NPX_wide))
colSums(is.na(Io_NPX_wide)) #IL4 ne ha 1 

#For a single value it's ok to substitute the na with the other values' mean:
Io_NPX_wide_noNA <- Io_NPX_wide
Io_NPX_wide_noNA[is.na(Io_NPX_wide_noNA[, "IL4"]), "IL4"] <- mean(Io_NPX_wide_noNA$IL4, na.rm = TRUE)
colSums(is.na(Io_NPX_wide_noNA))

#Normalizing
Cm_matrix <- scale(Cm_NPX_wide[,-1])
rownames(Cm_matrix) <- Cm_NPX_wide$Sample
Cv_matrix <- scale(Cv_NPX_wide[,-1])
rownames(Cv_matrix) <- Cv_NPX_wide$Sample
Io_matrix <- scale(Io_NPX_wide_noNA[, -1][, apply(Io_NPX_wide_noNA[, -1], 2, sd) > 0]) #There is a protein, IL33, with a variance of 0 (all values are identical). Removing it
rownames(Io_matrix) <- Io_NPX_wide_noNA$Sample
sum(colSums(is.na(Io_matrix)))

## PCA:
pca_Cm <- prcomp(Cm_matrix)
pca_Cv <- prcomp(Cv_matrix)
pca_Io <- prcomp(Io_matrix)
summary(pca_Cm)
summary(pca_Cv)
summary(pca_Io)

#Scree plot: Visualizing the explained variances
scree_Cm <- fviz_eig(pca_Cm, addlabels = TRUE, main = "PCA Scree plot - Cardiometabolic panel") #library("factoextra")
scree_Cv <- fviz_eig(pca_Cv, addlabels = TRUE, main = "PCA Scree plot - Cardiovascular panel")
scree_Io <- fviz_eig(pca_Io, addlabels = TRUE, main = "PCA Scree plot - Immuno Oncology panel")

plot_grid(scree_Cm, scree_Cv, scree_Io, nrow = 2, ncol = 2) #library(cowplot)

#Plot PCA1 vs PCA2 vs disease condition, Gender, Age
#a.Cardiometabolic
pca_Cm_group_graph <- pca_graph(pca_Cm, "Cardiometabolic", metadata, "Group")
pca_Cm_gender_graph <- pca_graph(pca_Cm, "Cardiometabolic", metadata, "Gender")
pca_Cm_age_graph <- pca_graph(pca_Cm, "Cardiometabolic", metadata, "Age")

#b.Cardiovascular
pca_Cv_group_graph <- pca_graph(pca_Cv, "Cardiovascular", metadata, "Group")
pca_Cv_gender_graph <- pca_graph(pca_Cv, "Cardiovascular", metadata, "Gender")
pca_Cv_age_graph <- pca_graph(pca_Cv, "Cardiovascular", metadata, "Age")

#c.Immuno Oncology
pca_Io_group_graph <- pca_graph(pca_Io, "Immuno Oncology", metadata, "Group")
pca_Io_gender_graph <- pca_graph(pca_Io, "Immuno Oncology", metadata, "Gender")
pca_Io_age_graph <- pca_graph(pca_Io, "Immuno Oncology", metadata, "Age")

#Visualizing all graphs together:
plot_grid(pca_Cm_group_graph, pca_Cv_group_graph, pca_Io_group_graph, pca_Cm_gender_graph, pca_Cv_gender_graph, pca_Io_gender_graph, pca_Cm_age_graph, pca_Cv_age_graph, pca_Io_age_graph, nrow = 3, ncol = 3)   #library(cowplot)

### Conclusions PCA: ###
# No extremely distinct pattern. The 'severe' group tends to always cluster closer and more centrally compared to the other groups: more similar characteristics among the severe group?
# The explained variances are very heterogeneous for Cv and Io, less so for Cm




## 3. Explorative analysis: t-SNE ##

#   Using scaled data, "matrix"
#   Arguments:
#   perplexity -> how many "neighbors" are considered important (Choose perplexity based on the rule: perplexity < (n_samples − 1) / 3)
#   dims -> in how many dimensions to visualize the results
#   verbose -> whether to display iterations in the output (not very important, used for debugging)
#   max_iter -> number of iterations to perform. Modify to obtain slightly different plots
tsne_Cm <- Rtsne(Cm_matrix, dims = 2, perplexity = 7, verbose = TRUE, max_iter = 1000)
tsne_Cv <- Rtsne(Cv_matrix, dims = 2, perplexity = 7, verbose = TRUE, max_iter = 1000)
tsne_Io <- Rtsne(Io_matrix, dims = 2, perplexity = 7, verbose = TRUE, max_iter = 1000)

#Visualizing tSNE1 vs tSNE2
#Graphs tSNE per "Group"
tsne_Cm_group_graph <- tsne_graph(tsne_Cm, Cm_matrix, "Cardiometabolic", metadata, "Group")
tsne_Cv_group_graph <- tsne_graph(tsne_Cv, Cv_matrix, "Cardiovascular", metadata, "Group")
tsne_Io_group_graph <- tsne_graph(tsne_Io, Io_matrix, "Immuno Oncology", metadata, "Group")
#Graphs tSNE per "Gender"
tsne_Cm_gender_graph <- tsne_graph(tsne_Cm, Cm_matrix, "Cardiometabolic", metadata, "Gender")
tsne_Cv_gender_graph <- tsne_graph(tsne_Cv, Cv_matrix, "Cardiovascular", metadata, "Gender")
tsne_Io_gender_graph <- tsne_graph(tsne_Io, Io_matrix, "Immuno Oncology", metadata, "Gender")
#Graphs tSNE per "Age"
tsne_Cm_age_graph <- tsne_graph(tsne_Cm, Cm_matrix, "Cardiometabolic", metadata, "Age")
tsne_Cv_age_graph <- tsne_graph(tsne_Cv, Cv_matrix, "Cardiovascular", metadata, "Age")
tsne_Io_age_graph <- tsne_graph(tsne_Io, Io_matrix, "Immuno Oncology", metadata, "Age")
 
#printing the graphs together:
plot_grid(tsne_Cm_group_graph, tsne_Cv_group_graph, tsne_Io_group_graph, tsne_Cm_gender_graph, tsne_Cv_gender_graph, tsne_Io_gender_graph, tsne_Cm_age_graph, tsne_Cv_age_graph, tsne_Io_age_graph, nrow = 3, ncol = 3)

### t-SNE Conclusions: ###
# Tests were conducted with max_iter: 1000 and 1500, and perplexity: 5 or 7.
# Not much of interest is observable at first glance. Similar conditions appear to have some nearby points, but not in a particularly generalizable pattern.
# Perhaps the "Cardiovascular" panel is slightly more clustered.
# In the case of (Cm, perplexity = 5, max_iter = 1000), there is an outlier in the control group – is this interesting or not??
  



## 4. Analisi esplorativa: Clustering (heatmap) ##

## Using scaled data, "matrix"
#The pheatmap function needs "samples" on the column, so reansposing:
heatmap_Cm_matrix <- t(Cm_matrix)  # colonne = campioni, righe = proteine
heatmap_Cv_matrix <- t(Cv_matrix) 
heatmap_Io_matrix <- t(Io_matrix)

#Heatmaps:
heatmap_cluster(heatmap_Cm_matrix, "Cardiometabolic", metadata)
heatmap_cluster(heatmap_Cv_matrix, "Cardiovascular", metadata)
heatmap_cluster(heatmap_Io_matrix, "Immuno Oncology", metadata)

## Heatmap Conclusions:
# There is no particularly evident clustering in general. Some observations include:
#  - In Cm, there is a clear separation in the expression of a central block, highly expressed --> cluster of samples with high expression (only Control and Mild).
#  - In Cv, there is a bit more clustering among samples, but the various blocks are quite mixed ("Severe" samples are particularly mixed, and again show lower expression levels).
#  - In Io, expressions are generally very low. Only 3 samples show more red, and again, they are not "Severe".




## 5. Explorative analysis: proteins and patients distribution (general Boxplot) ##

# The data is needed in 'long' format. Already obtained above.

#Now the boxplots:
#a. Cardiometabolic
box_Cm_assay <- boxplot_distr(Cm_NPX_long, "Assay", "Cardiometabolic")
box_Cm_sample <- boxplot_distr(Cm_NPX_long, "Sample", "Cardiometabolic")
#b. Cardiovascular
box_Cv_assay <- boxplot_distr(Cv_NPX_long, "Assay", "Cardiovascular")
box_Cv_sample <- boxplot_distr(Cv_NPX_long, "Sample", "Cardiovascular")
#c. Immuno Oncology
box_Io_assay <- boxplot_distr(Io_NPX_long, "Assay", "Immuno Oncology")
box_Io_sample <- boxplot_distr(Io_NPX_long, "Sample", "Immuno Oncology")

plot_grid(box_Cm_assay, box_Cm_sample, nrow = 2, ncol = 1)
plot_grid(box_Cv_assay, box_Cv_sample, nrow = 2, ncol = 1)
plot_grid(box_Io_assay, box_Io_sample, nrow = 2, ncol = 1)



### Explorative analysis conclusions ###
# To summarize: we have seen these data have little clustering.
# The "Severe" case stands out because it shows little variability in the PCA (points more centered at the origin), is poorly clusterable, and has lower expression levels in the heatmap.
# This suggests that this condition generally exhibits protein underexpression, but it is generalized, meaning possibly there are no specific protein clusters.
# PCA and t-SNE exploratory analysis showed no evident separation of samples by Age or Gender. So there doesn't seem to be any confunding role for "Age" and "Gender": they will not be included in the DE statistical analysis






### (II) CODE START: Differential Expression ###

## 1. Obtaining the work dataframe ##

#Merging all panels into one
all_NPX_long <- bind_rows(Cm_NPX_long, Cv_NPX_long, Io_NPX_long)
all_NPX_long$Assay <- as.character(all_NPX_long$Assay)  #altrimenti mantiene "levels" fantasma


#Removing the sample C2, that is only present in the Io panel. Needed in order to uniform the data
common_samples <- Reduce(intersect, list(Cm_NPX_long$Sample, Cv_NPX_long$Sample, Io_NPX_long$Sample))
NPX_long <- all_NPX_long[all_NPX_long$Sample %in% common_samples, ]  

#Transferring protein names on a variable (assay):
assay <- unique(NPX_long$Assay)

#Finally, obtaining the work dataframe:
NPX_long <- merge(NPX_long, metadata[, c("Sample", "Group")], by = "Sample")
NPX_long$Group <- factor(NPX_long$Group)




## 2. Differential Expression (DE) Statistical Analysis with limma ##

# limma is the preferred method here because there are few samples and many proteins.

#For limma, the format must have 1 sample per column. Obtaining the merged data in this format, starting from the wide format:
#Deleting the additional sample in "Io" and sorting:
Io_NPX_wide <- Io_NPX_wide[Io_NPX_wide$Sample %in% common_samples, ]  #common_samples è stato ottenuto alla fine di #1

temp <- Cm_NPX_wide %>% inner_join(Cv_NPX_wide, by = "Sample") %>% inner_join(Io_NPX_wide, by = "Sample")
temp <- as.data.frame(t(temp))
colnames(temp) <- temp[1, ]  
all_NPX_wide_full <- temp[-1, ]

#Filtering by LODs (on the wide dataset)
#From these, filtering proteins with 80% of values below the LOD. I already did this for the long format, so just retrieving the "good_assays" here.
all_NPX_wide_full <- all_NPX_wide_full %>% rownames_to_column("Assay")
good_assays_all <- c(good_assays_Cm, good_assays_Cv, good_assays_Io)  #The good-assays were already identified during the filtering of the long format. Merging them.

all_NPX_wide <- all_NPX_wide_full %>% filter(Assay %in% good_assays_all) 
all_NPX_wide <- as.data.frame(all_NPX_wide)
rownames(all_NPX_wide) <- all_NPX_wide$Assay
all_NPX_wide$Assay <- NULL


#Finally, limma needs a numerical matrix:
all_NPX_wide_mat <- as.matrix(sapply(all_NPX_wide, as.numeric))
rownames(all_NPX_wide_mat) <- rownames(all_NPX_wide)

#Samples must be alphabetically ordered (groups are A, B, C)
all_NPX_wide_mat <- all_NPX_wide_mat[, order(colnames(all_NPX_wide_mat))]
# The metadata are already sorted alphabetically. However, to make the code more formally correct, I sort the metadata anyway:
metadata_ord <- metadata %>% filter(Sample %in% common_samples) %>% arrange(Sample)
group <- factor(metadata_ord$Group)

#Test design
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
design

#fit
fit <- lmFit(all_NPX_wide_mat, design)

#Statistic test

contrasts <- makeContrasts(
  Mild_vs_Control   = Mild - Control,
  Severe_vs_Control = Severe - Control,
  Severe_vs_Mild = Severe - Mild,
  levels = design
)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

#Extracting fit informations and inserting them in a dataframe
contrast_names <- colnames(fit2$coefficients)

## This study has an exploratory nature and a limited sample size, so no multiple testing correction was applied. 
## Results with p < 0.05 (nominal) are reported as hypothesis-generating findings, and not as definitive results.
## DE results are expected to be moderate, therefore this is the best set-up for this particular case

#List for complete results
tot_res_list <- list()
for (cname in contrast_names) {
  #topTable for a specific contrast
  tt_c <- topTable(fit2, coef = cname, number = Inf, sort.by = "none", adjust.method = "none") %>%
    rownames_to_column("Assay") %>%
    dplyr::select(Assay, P.Value, logFC)  %>%
    mutate(group1 = strsplit(cname, "_vs_")[[1]][1],
           group2 = strsplit(cname, "_vs_")[[1]][2]) %>%
    dplyr::select(Assay, group1, group2, P.Value, logFC) %>%
    dplyr::rename(p_value = P.Value)
  tot_res_list[[cname]] <- tt_c
}
tot_p_values_limma <- bind_rows(tot_res_list)

#List for significant results
res_list <- list()
for (cname in contrast_names) {
  #topTable for a specific contrast
  tt_c <- topTable(fit2, coef = cname, number = Inf, sort.by = "none", adjust.method = "none") %>%
    rownames_to_column("Assay") %>%
    dplyr::select(Assay, P.Value, logFC) %>%
    filter(P.Value <= 0.05) %>%       # threshold
    mutate(group1 = strsplit(cname, "_vs_")[[1]][1],
           group2 = strsplit(cname, "_vs_")[[1]][2]) %>%
    dplyr::select(Assay, group1, group2, P.Value, logFC) %>%
    dplyr::rename(p_value = P.Value)
  res_list[[cname]] <- tt_c
}

#Merging contrasts
sign_p_values_limma <- bind_rows(res_list)


#Violin plots
sign_proteins <- unique(sign_p_values_limma$Assay)
vio_plot_list_limma <- list()
for (i in sign_proteins) {
  vio_plot_list_limma[[i]] <- plot_violin(NPX_long, i, sign_p_values_limma, "Limma")
}
#grid_graph function for visualization
grid_graph(vio_plot_list_limma, 9, 3, 3)




## 3. Volcano plot per data visualization: ##

## Volcano plot
vol_plot_list <- list()
vol_plot_list[[1]] <- volcano_graph_limma(tot_p_values_limma, "Mild", "Control")
vol_plot_list[[2]] <- volcano_graph_limma(tot_p_values_limma, "Severe", "Control")
vol_plot_list[[3]] <- volcano_graph_limma(tot_p_values_limma, "Severe", "Mild")

do.call(plot_grid, c(vol_plot_list, nrow = 2, ncol = 2))




## 4. Outlier detection ##
# To conclude, it's important to verify that the significance obtained does not depend on extreme outliers, through visual inspection of the violin plots (performed on limma results).







### (III) CODE START: Functional Enrichment analysis ###

## Renaming significant genes for clarity ##
sign_results <- sign_p_values_limma
tot_results <- tot_p_values_limma

all_id <- row.names(all_NPX_wide)
all_id <- as.data.frame(all_id)
colnames(all_id) <- "Assay"


## 1. ORA with HUGO SYMBOLS ##

#Preparing data
background_sym <- all_id$Assay[all_id$Assay %in% tot_results$Assay]   #Keeping the ones filtered for LOD (204)


#Checking if these are present in the database
all_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
mapped_symbols <- intersect(background_sym, all_symbols)
length(mapped_symbols)          # number of mapped symbols
# Only 158/204 are present in the database. This is because some HUGO names utilized are slightly wrong. 

#After checking databases (and using AI to identify the 46 missing genes in the list), 
#this is the conversion needed for the missing genes:
olink_to_hugo <- c(
  "BMP-6"                = "BMP6",
  "ADAM-TS13"            = "ADAMTS13",
  "IL-4RA"               = "IL4R",
  "IL-1ra"               = "IL1RN",
  "MCP-1"                = "CCL2",
  "MCP-2"                = "CCL8",
  "MCP-3"                = "CCL7",
  "MCP-4"                = "CCL13",
  "LOX-1"                = "OLR1",
  "HSP 27"               = "HSPB1",
  "PD-L1"                = "CD274",
  "IFN-gamma"            = "IFNG",
  "TRAIL"                = "TNFSF10",
  "TRAIL-R2"             = "TNFRSF10B",
  "TWEAK"                = "TNFSF12",
  "PAR-1"                = "F2R",
  "PSGL-1"               = "SELPLG",
  "VEGFR-2"              = "KDR",
  "CASP-8"               = "CASP8",
  "Gal-1"                = "LGALS1",
  "Dkk-1"                = "DKK1",
  "SCF"                  = "KITLG",
  "GH"                   = "GH1",
  "IL-27"                = "IL27",
  "IL8"                  = "CXCL8",
  "FGF-21"               = "FGF21",
  "RAGE"                 = "AGER",
  "KIM1"                 = "HAVCR1",
  "NEMO"                 = "IKBKG",
  "HB-EGF"               = "HBEGF",
  "GDF-2"                = "GDF2",
  "CTSL1"                = "CTSL",
  "CSF-1"                = "CSF1",
  "IL4RA"                = "IL4R",
  "LAP TGF-beta-1"       = "TGFB1",
  "MIC-A/B"              = "MICA",
  "IL12"                 = "IL12A",
  "GIF"                  = "CBLIF",
  "PIgR"                 = "PIGR",
  "FS"                   = "FST",
  "TM"                   = "THBD",
  "IgG Fc receptor II-b" = "FCGR2B",
  "GT"                   = "GCNT2",
  "hOSCAR"               = "OSCAR",
  "HAOX1"                = "HAO1",
  "CAIX"                 = "CA9",
  "MUC-16"               = "MUC16"
)

#Applying the conversion:
background_sym_fixed <- background_sym
names_to_fix <- names(olink_to_hugo)
background_sym_fixed <- ifelse(
  background_sym_fixed %in% names_to_fix,
  olink_to_hugo[background_sym_fixed],
  background_sym_fixed
)

#Repeating the mapping
mapped_symbols <- intersect(background_sym_fixed, all_symbols)
length(mapped_symbols) 
#This time 204/204 gene names were found. Proceeding

background_sym <- background_sym_fixed
background_sym <- background_sym[background_sym %in% mapped_symbols]

#Before proceeding it is important to check the same on the significant genes (sign_results).
#Checking that all the significant HUGO symbols are correctly named (now that some background names are different)
mapped_symbols_sign <- intersect(sign_results$Assay, background_sym_fixed)
length(mapped_symbols_sign)
#14/28. Only half names are common.
#Like before:
sign_assays_fixed <- ifelse(
  unique(sign_results$Assay) %in% names(olink_to_hugo),
  olink_to_hugo[unique(sign_results$Assay)],
  unique(sign_results$Assay)
)
#Checking:
mapped_symbols_sign <- intersect(sign_assays_fixed, background_sym_fixed)
length(mapped_symbols_sign)
#22/28, but the 6 missing were duplicates. So it is correct: 22/22


#ORA with HUGO SYMBOL: BP
ego_BP_sym <- enrichGO(
  gene          = sign_assays_fixed,   
  universe      = background_sym,            
  OrgDb         = org.Hs.eg.db,             
  keyType       = "SYMBOL",                  
  ont           = "BP",                      
  pvalueCutoff =  0.05,
  qvalueCutoff =  0.2
)

#ORA with HUGO SYMBOL: MF
ego_MF_sym <- enrichGO(
  gene          = sign_assays_fixed,        
  universe      = background_sym,            
  OrgDb         = org.Hs.eg.db,              
  keyType       = "SYMBOL",               
  ont           = "MF",                      
  #pAdjustMethod=  "BH",                     
  pvalueCutoff =  0.05,
  qvalueCutoff =  0.2
)

nrow(ego_BP_sym) #0
nrow(ego_MF_sym) #0
## ! No enrichment term were found with the cutoff at 0.05 for p-value and 0.2 for q-value. 
## This is not completely unexpected, as both the significant genes and the background list are very short. 


dev.off()