#!/usr/bin/env Rscript

# Bradley Feb 2024
# Checking the output of trans-eQTL results
# conda activate rplots

# Load packages
library(ggplot2)
library(dplyr)
library(readr)
library(tidyverse)
library(ggplot2)
library(ggrepel)

# Define script options
alpha=0.05
basedir = "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/eqtl/results/"
setwd(basedir)

# Define wrappers:
read_eqtls <- function(base_dir){
    sumstat.files <- list.files(base_dir, pattern="*trans-by-cis_bonf_fdr.tsv",
                                full.names=FALSE, recursive=TRUE)

    list.of.sumstat.dfs <- lapply(sumstat.files, function(file.name) {
    # Get info from filename
    split.file.name <- str_split(file.name,'/')[[1]]
    name <- split.file.name[2] # Diff to Tobi
    n.pcs <- split.file.name[1] # Diff to Tobi

    # Get info from name
    split.name <- str_split(name,'-')[[1]]
    annotation.type <- split.name[1]
    annotation <- split.name[2]
    method <- split.name[3]

    # Read file
    df <- read_tsv(paste0(base_dir, file.name))

    # Add new columns
    mutate(df, num_PCs = n.pcs,
            name = name,
            annotation_type = annotation.type,
            annotation = annotation,
            method = method)
    })

    #sumstat.df <- bind_rows(list.of.sumstat.dfs) %>% 
    #  left_join(dplyr::select(grch38, ensgene, symbol),
    #            by = c("phenotype_id" = "ensgene"),multiple='first')

    sumstat.df <- bind_rows(list.of.sumstat.dfs)

    geno_pheno <- read_tsv(paste0(base_dir,'../genotype_phenotype_mapping.tsv'))

    # How many individuals per annotation
    n.individuals.df <- geno_pheno %>% 
        group_by(Sample_Category) %>% 
        tally(name = 'n_indiv') %>% 
        ungroup() %>% 
        separate(Sample_Category, sep = '-', into = c('annotation_type',
                                                        'annotation',
                                                        'method'))

    # How many cells per individual-cell type pseudobulk
    # Filter out n < 5 because those are dropped before eQTL mapping is performed
    cell.metadata <- read_csv(paste0(base_dir,'../eqtl_processed.obs.csv'))
    cells.per.indiv.annot <- tibble()

    for (annot in unique(sumstat.df$annotation_type)){
    print(annot)
    tmp.df <- cell.metadata %>% 
        group_by(Corrected_genotyping_ID, get(annot)) %>% 
        mutate(n_cells=n()) %>% 
        filter(n_cells>=5) %>% 
        slice(1) %>% 
        ungroup() %>% 
        mutate(annotation_type = annot) %>% 
        rename(annotation = `get(annot)`)
        cells.per.indiv.annot <- rbind(cells.per.indiv.annot, tmp.df)
        } #Corrected_genotyping_ID not PostQC_genotyping_ID

    # # Get total number of CDs and healthies for each category
    # disease.per.annot <- cells.per.indiv.annot %>% 
    #   group_by(annotation, annotation_type) %>% 
    #   count(disease_status) %>% 
    #   ungroup() %>% 
    #   pivot_wider(names_from = disease_status, values_from = n) %>% 
    #   mutate(total = cd + healthy + other_inflammation) %>% 
    #   mutate(prop_cd = cd / total,
    #          prop_non_ibd = (healthy + other_inflammation) / total)

    # Get the average number of cells per annotation
    average.cells.per.annot <- cells.per.indiv.annot %>% 
        group_by(annotation, annotation_type) %>% 
        summarise(mean_n_cells = mean(n_cells)) %>% 
        ungroup()

    # Count number of genes expressed in each annot-indiv combination
    phenotype_file <- read_tsv(paste0(base_dir,'../phenotype_file.tsv'))

    count.genes.df <- phenotype_file %>%
        summarize_all(~sum(. != 0)) %>% 
        rename('All_genes' = ...1)

    average.genes.per.annot <- count.genes.df %>% 
    dplyr::select(-All_genes) %>% 
        pivot_longer(everything(), names_to = "annotation_sample", values_to = "n_genes") %>% 
        separate(annotation_sample,into = c('annotation', 'sample'), sep='-dMean_') %>% 
        group_by(annotation) %>% 
        summarise(mean_n_genes = mean(n_genes)) %>% 
        ungroup() %>% 
        separate(annotation, into = c('annotation_type', 'annotation'), sep='-')

    sumstat.df <- sumstat.df %>%
      group_by(name) %>%
      mutate(ntests = n())

    sumstat.all = sumstat.df

    sumstat.df <- sumstat.df %>%
        filter(pval_bonf_fdr < alpha)
    
    sumstat.df$name = gsub("-dMean", "", sumstat.df$name)
    sumstat.df$label__machine_old = unlist(strsplit(sumstat.df$name, "\\-"))[c(F,T)]


    sumstat.df <- sumstat.df %>% left_join(n.individuals.df, by = join_by(annotation_type, annotation, method)) %>%
        left_join(average.cells.per.annot, by = join_by(annotation_type, annotation)) %>%
        left_join(average.genes.per.annot, by = join_by(annotation_type, annotation)) %>%
        # left_join(disease.per.annot, by = join_by(annotation_type, annotation)) %>% 
        left_join(annot.mapping, by=c('label__machine_old')) %>%
        # Add a new column 'num_eGenes' with the number of rows
        group_by(annotation_type, annotation) %>%
        mutate(num_trans_eGenes = n())

    sumstat.df$category__machine = ifelse(sumstat.df$annotation_type == "category__machine", sumstat.df$category_new, "")
    sumstat.df$category__machine = ifelse(sumstat.df$annotation_type == "everything__machine", "Everything", sumstat.df$category_new)
    sumstat.df$class = ifelse(sumstat.df$annotation_type == "label__machine", "Cell type", "Major category")

  return(sumstat.df)
}

plot_negenes_power <- function(df, x.var){
  if (x.var == 'mean_n_cells'){
  #  x.scale.lower <- 5
  #  x.scale.higher <- 1000
    x.lab <- 'Mean number of cells per donor'
  }
  if (x.var == 'mean_n_genes'){
  #  x.scale.lower <- 3500
  #  x.scale.higher <- 22000
    x.lab <- 'Mean number of genes expressed '
  }
  if (x.var == 'n_indiv'){
  #  x.scale.lower <- 30
  #  x.scale.higher <- 240
    x.lab <- 'Number of individuals'
  }
  if (x.var == "ntests"){
    x.lab <- "Total number of tests"
  }
  ggplot(df, aes(y = num_trans_eGenes, x = .data[[x.var]])) +
    geom_point(aes(fill = category__machine, size = class), shape = 21, stroke = 0.7) +
    # scatterpie::geom_scatterpie(aes(fill = category_new, size = Annot_class_label),
    #                             shape = 21, stroke = 0.7) +
    # geom_text(aes(colour = category_new, label = label_new), size=4) +
    # geom_text(data=filter(df, h2.enriched == TRUE),
    #           aes(size=Annot_class_label),
    #           fontface = "bold", colour='white',vjust=0.65, label = '*') +
    geom_text_repel(data=df,
                    (aes(label = label_new)),box.padding = 0.5,point.padding = 1,
                    bg.color = "white") +
    scale_size_manual(values = c("Major category" = 7, "Cell type" = 3),
                      name = 'Annotation pseudobulk\nresolution') +
    #scale_x_log10(limits=c(x.scale.lower,x.scale.higher), n.breaks=6) +
    scale_x_log10(n.breaks=6) +
    scale_y_log10() +
    scale_fill_manual(values = umap.category.palette,
                      name = 'Major category\ncolour',
                      guide='none') +
    theme_classic() + 
    guides(size = guide_legend(order = 1),
           fill = 'none'
           # fill = guide_legend(override.aes = list(size=7), order = 0)
    ) +
    labs(x = x.lab, y = 'Number of trans-eGenes') +
    theme(legend.position = c(0.2, 0.8),
          legend.background = element_rect(fill = "white", color = "black"))
}

# Set up outdir
outdir="summary_plots/trans"
if(file.exists(outdir) == F){
    dir.create(outdir, recursive=T)
}

# Load the palettes and mapping files
umap.palette.df <- read_tsv('/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/proc_data/ti_umap_palette.tsv')
umap.category.palette <- deframe(dplyr::select(umap.palette.df, category, palette_category))
umap.celltype.palette <- deframe(dplyr::select(umap.palette.df, label, palette_label))

# Import upto date TI annotations
old_annots <- read.csv('/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/proc_data/data-cluster_labels.csv')
new_annots <- read.csv('/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/proc_data/data-clean_annotation.csv') 
annot.mapping <- left_join(old_annots, new_annots, by='cluster', suffix=c('_old', '_new')) %>% 
  dplyr::select(label__machine_old, category__machine_new, category_new, label_new)


# Generate the the metadata (eqtl_processed.obs.csv)
# adata.obs.to_csv(f"{basedir}/../eqtl_processed.obs.csv", index=False)
# Get the phenotype file
# find . -type f -name "phenotype_file.tsv" -exec cp {} ./ \;


# Load the eqtls
sumstat.df <- read_eqtls(basedir)
sumred = as.data.frame(sumstat.df[,-c(1:9)])
sumred$label_new = ifelse(sumred$class == "Major category", gsub("\\_", " ", sumred$annotation), "")
sumred = sumred[!duplicated(sumred),]
sumred$category__machine = ifelse(is.na(sumred$category__machine), sumred$label_new, sumred$category__machine)

# Make the plots
metrics = c("n_indiv", "mean_n_cells", "mean_n_genes", "ntests")

# Apply the function to each column and save plots
plot_list = vector("list", length = length(metrics))
for(m in metrics){
    print(m)
    # Make plot
    index = which(metrics == m)
    plot_list[[index]] = plot_negenes_power(sumred, m)
    # Save
    ggsave(paste0(outdir, "/num_trans_eGenes_vs_", m, ".png"), plot_list[[index]], width = 7, height = 5, dpi = 300)
}

# General plots of results
# 1. Number of trans-eQTL target genes per variant
sumstat.df$var_labal = paste0(sumstat.df$variant_id, "_", sumstat.df$name)
ntransvar = as.data.frame(table(sumstat.df$var_labal))
ntransvarplot = ggplot(ntransvar, aes(x = Freq)) +
                        geom_histogram(binwidth = 0.2, fill = "grey", color = "black", alpha = 0.7) +
                        labs(title = "Number of trans-eGenes per trans-eQTL", x = "Number of trans-eGenes", y = "Frequency") + 
                        theme_bw() + 
                        theme(plot.title = element_text(face = "bold"))

ggsave(paste0(outdir, "/num_trans_eGenes_per_eQTL.png"), ntransvarplot, width = 6.5, height = 5, dpi = 300)

# If we remove everything
notall = ntransvar[-grep("Everything", ntransvar$Var1),]
ntransvarplotnotall = ggplot(notall, aes(x = Freq)) +
                        geom_histogram(binwidth = 0.2, fill = "grey", color = "black", alpha = 0.7) +
                        labs(title = "Number of trans-eGenes per trans-eQTL (Not inc. everything)", x = "Number of trans-eGenes", y = "Frequency") + 
                        theme_bw() + 
                        theme(plot.title = element_text(face = "bold")) + 
                        scale_x_continuous(breaks = seq(min(notall$Freq), max(notall$Freq), by = 1))


ggsave(paste0(outdir, "/num_trans_eGenes_per_eQTL_not_all.png"), ntransvarplotnotall, width = 6.5, height = 5, dpi = 300)

# Make a pie chart for the number of tests and number of significant results
nall = length(unique(sumstat.all$variant_id))
nsig = length(unique(sumstat.df$variant_id))
df = data.frame(category = c("FDR>0.05", "FDR<0.05"), value = c(nall-nsig, nsig))
df$percentage <- df$value / sum(df$value) * 100

# Create a pie chart
pie_chart <- ggplot(df, aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("orange", "navy"))

# Add labels, connectors, background, and percentage
pie_chart +
  geom_text(
    aes(label = paste0(category, "\n", sprintf("%.1f%%", percentage))),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  geom_segment(
    aes(x = 0.5, y = 0, xend = 0.5, yend = 0.5),
    color = "black", size = 0.5, arrow = arrow(length = unit(0.2, "cm"))
  ) +
  theme(legend.position = "none") +
  labs(title = "Pie Chart with FDR Categories") +
  annotate(
    "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
    fill = "white", alpha = 0.8
  )

ggsave(paste0(outdir, "/num_trans_eGenes_per_tests_pie.png"), pie_chart, width = 6.5, height = 5, dpi = 300)

# Distance for cis-eQTL to trans-eGene
phenotype_file = '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_eqtl/eqtl/assets/gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt'
gene_pos = read.delim(phenotype_file)
colnames(gene_pos)[which(colnames(gene_pos)== "feature_id")] = "phenotype_id"
sumstat.df = merge(sumstat.df, gene_pos, by="phenotype_id")
sumstat.df$variant_chr = unlist(strsplit(sumstat.df$variant_id, "\\:"))[c(T,F,F,F)]
sumstat.df$variant_chr = gsub("chr", "", sumstat.df$variant_chr)
sumstat.df$variant_pos = unlist(strsplit(sumstat.df$variant_id, "\\:"))[c(F,T,F,F)]
sumstat.df$distance = ifelse(sumstat.df$strand == "+", as.numeric(sumstat.df$variant_pos) - sumstat.df$start, as.numeric(sumstat.df$variant_pos) - sumstat.df$end)
sumstat.df$distance = ifelse(sumstat.df$variant_chr != sumstat.df$chromosome, NA, sumstat.df$distance)
sumstat.df$distance = sumstat.df$distance/1e6
samec = sum(!is.na(sumstat.df$distance))
perc = 100*samec / nrow(sumstat.df)

distance_plot = ggplot(sumstat.df[!is.na(sumstat.df$distance),], aes(x = distance)) +
                        geom_histogram(binwidth = 0.2, fill = "grey", color = "black", alpha = 0.7) +
                        labs(title = paste0("Distance of trans-eQTL to eGene (same chromosome - ", samec, " [", signif(perc, 4), "%])"), x = "log10(Distance (Mb))", y = "Frequency") + 
                        theme_bw() + 
                        theme(plot.title = element_text(face = "bold", size=10)) + 
                        geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red", size = 1)

ggsave(paste0(outdir, "/distance_eQTL_target.png"), distance_plot, width = 6.5, height = 5, dpi = 300)

# QQ and histogram p-values
# Using everything
everything_raw = read.delim(gzfile(paste0(base_dir, "TensorQTL_eQTLS/everything__machine-Everything-dMean/trans-by-cis_all.tsv.gz")))
phist  <- ggplot(everything_raw, aes(x = pval)) +
                        geom_histogram(binwidth = 0.01, fill = "grey", color = "black", alpha = 0.7) +
                        labs(title = "P-value histogram - Everything", x = "raw p-value", y = "Frequency") + 
                        theme_bw() + 
                        theme(plot.title = element_text(face = "bold"))

ggsave(paste0(outdir, "/everything_pvalue_hist.png"), phist, width = 6.5, height = 5, dpi = 300)

# Check this
png(paste0(outdir, "/everything_qq.png"),  width = 600, height = 400, units = "px", res = 100)
qqplot(-log10(everything_raw$pval), qnorm(ppoints(length(everything_raw$pval))), main = "QQ Plot for Raw P-Values - Everything")
abline(0, 1, col = "red", lty = 2)  # Add a reference line
xlabel <- expression(-log10(p-values))
ylabel <- expression("Theoretical Quantiles (" * Z * ")")
# Close the PDF file
dev.off()


# Have a look at an interesting example
test_annot = sumred[sumred$label_new != "Everything",]
notall = as.data.frame(sumstat.df[sumstat.df$annotation != "Everything",])
notall = notall[order(-notall$num_trans_eGenes),]