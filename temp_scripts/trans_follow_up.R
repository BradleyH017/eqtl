# Bradley Feb 2024
# Checking the output of trans-eQTL results
# conda activate rplots

setwd("/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/eqtl")
library(ggplot2)
outdir="summary_plots/trans"
if(file.exists(outdir) == F){
    dir.create(outdir, recursive=T)
}

# list the input dirs
annots = list.files("results/TensorQTL_eQTLS")

# Define alpha
alpha = 0.1

# Load in the results per annot
res = vector("list", length = length(annots))
for(level in annots){
    print(level)
    index = which(annots == level)
    file=paste0("results/TensorQTL_eQTLS/", level, "/OPTIM_pcs/base_output__base/trans-by-cis_bonf_fdr.tsv")
    if(file.exists(file)){
        temp = read.delim(file, sep = "\t")
        temp = temp[temp$pval_bonf_fdr < alpha,]
        if(nrow(temp) > 0){
            temp$label = level
            res[[index]] = temp
        }
    }
}
res = do.call(rbind, res)

# Get the number of trans-effects per level
trans_per_annot = as.data.frame(table(res$label))
colnames(trans_per_annot) = c("label", "n_trans")
for(level in annots){
    if(level %in% trans_per_annot$label != T){
        temp = data.frame("label" = level, "n_trans" = 0)
        trans_per_annot = rbind(trans_per_annot, temp)
    }
}

# Make a histogram of these results
hist_plot <- ggplot(trans_per_annot, aes(x = n_trans)) +
  geom_histogram(binwidth = 0.5, fill = "grey", color = "black", alpha = 0.7) +
  labs(x = "Number of trans-eQTL effects - All annotations", y = "Frequency")

# Save the plot using ggsave
ggsave(paste0(outdir, "/histogram_ntrans.png"), plot = hist_plot, width = 6, height = 4, units = "in", dpi=300)