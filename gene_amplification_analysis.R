library(openxlsx)
library(data.table)
library(ggplot2)

##################################
## Load and calculate frequency ##
##################################

WGDpos <- read.xlsx('amplification_results.xlsx', sheet = 'WGD+amp+')
WGDneg <- read.xlsx('amplification_results.xlsx', sheet = 'WGD-amp+')
dataset <- rbind(WGDpos, WGDneg)
length(unique(dataset$number))
summary(dataset$copy_number_LOG2FC)
setDT(dataset)
dataset[!is.na(ploidy), WGD_label := fifelse(ploidy >= 2.7, 'positive', 'negative')]
# dataset[!is.na(copynumber), WGD_label := fifelse(copynumber >= 3, 'positive', 'negative')]


# Unique sample counts
total_counts <- unique(dataset[, .(number, WGD_label)])[, .N, by = WGD_label]
# total_counts <- dataset[, .(number, WGD_label)][, .N, by = WGD_label]
setkey(total_counts, WGD_label)
n_pos <- total_counts["positive", N]
n_neg <- total_counts["negative", N]
total_counts

# Amplification frequencies
amp_data <- dataset[!is.na(copy_number_GENE) & copy_number_SVTYPE == "<DUP>"]
gene_counts <- amp_data[, .N, by = .(copy_number_GENE, WGD_label)]
gene_counts <- dcast(gene_counts, copy_number_GENE ~ WGD_label, value.var = "N", fill = 0)
dim(gene_counts)
if (!"positive" %in% names(gene_counts)) gene_counts[, positive := 0]
if (!"negative" %in% names(gene_counts)) gene_counts[, negative := 0]
dim(gene_counts)

gene_counts[, `:=`(
  freq_pos = 100 * positive / n_pos,
  freq_neg = 100 * negative / n_neg
)]

# OR, p-value, FDR
gene_counts[, c("pval", "OR") := {
  test <- mapply(function(pos, neg, pos_total, neg_total) {
    mat <- matrix(c(pos, pos_total - pos, neg, neg_total - neg), nrow = 2)
    ft <- fisher.test(mat)
    odds <- (pos / (pos_total - pos + 1e-5)) / (neg / (neg_total - neg + 1e-5))
    list(ft$p.value, odds)
  }, positive, negative, n_pos, n_neg, SIMPLIFY = FALSE)
  list(sapply(test, `[[`, 1), sapply(test, `[[`, 2))
}]
gene_counts[, FDR := p.adjust(pval, method = "fdr")]

gene_counts[, `:=`(
  p_stars = fifelse(pval < 0.001, "***", fifelse(pval < 0.01, "**", fifelse(pval < 0.05, "*", ""))),
  fdr_stars = fifelse(FDR < 0.001, "***", fifelse(FDR < 0.01, "**", fifelse(FDR < 0.05, "*", "")))
)]
gene_counts[, stat_label := paste0("OR=", sprintf("%.2f", OR),
                                   ", p=", formatC(pval, format = "e", digits = 2), p_stars,
                                   ", FDR=", formatC(FDR, format = "e", digits = 2), fdr_stars)]
#######################
## Top 20 by p-value ## 
#######################
top_genes <- gene_counts[order(pval)][1:20]
plot_data1 <- melt(top_genes,
                   id.vars = c("copy_number_GENE", "freq_pos"),
                   measure.vars = c("freq_pos", "freq_neg"),
                   variable.name = "WGD_label", value.name = "frequency")
plot_data1[WGD_label == "freq_pos", WGD_label := "positive"]
plot_data1[WGD_label == "freq_neg", WGD_label := "negative"]
plot_data1[WGD_label == "positive", frequency := -frequency]
plot_data1[, gene_order := factor(copy_number_GENE, levels = rev(top_genes[order(pval)]$copy_number_GENE))]
plot_data1 <- merge(plot_data1, unique(top_genes[, .(copy_number_GENE, stat_label)]),
                    by = "copy_number_GENE", all.x = TRUE)


setDT(plot_data1)

plot_data1[, `:=`(
  is_zero   = frequency == 0,
  y_lab     = fifelse(frequency == 0,
                      fifelse(WGD_label == "positive", -0.2, 0.2),
                      frequency),
  hjust_lab = fcase(
    frequency < 0,  1.05,
    frequency > 0, -0.05,
    frequency == 0 & WGD_label == "positive", 1.05,
    default = -0.05
  )
)]

ggplot(plot_data1, aes(x = gene_order, y = frequency, fill = WGD_label)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = y_lab,
                label = paste0(sprintf("%.1f", abs(frequency)), "%"),
                hjust = hjust_lab),
            size = 3, check_overlap = TRUE) +
  geom_text(data = subset(plot_data1, WGD_label == "negative"),
            aes(x = gene_order, y = 20, label = stat_label),
            inherit.aes = FALSE,
            hjust = 0, size = 3, color = "black") +
  coord_flip() +
  labs(title = "Top 20 Differentially Amplified Genes by WGD Status",
       x = "Gene", y = "Amplification Frequency (%)", fill = "WGD Status") +
  scale_fill_manual(values = c("positive" = "blue", "negative" = "red")) +
  scale_y_continuous(labels = abs, expand = expansion(mult = c(0.05, 0.6))) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

######################################
## Top 10 by WGD-positive frequency ##
######################################
# top_freq <- gene_counts[order(-freq_pos)][1:10]
# plot_data2 <- melt(top_freq,
#                    id.vars = c("copy_number_GENE", "freq_pos"),
#                    measure.vars = c("freq_pos", "freq_neg"),
#                    variable.name = "WGD_label", value.name = "frequency")
# plot_data2[WGD_label == "freq_pos", WGD_label := "positive"]
# plot_data2[WGD_label == "freq_neg", WGD_label := "negative"]
# plot_data2[WGD_label == "positive", frequency := -frequency]
# plot_data2[, gene_order := factor(copy_number_GENE, levels = rev(top_freq[order(-freq_pos)]$copy_number_GENE))]
# plot_data2 <- merge(plot_data2, unique(top_freq[, .(copy_number_GENE, stat_label)]),
#                     by = "copy_number_GENE", all.x = TRUE)
# 
# ggplot(plot_data2, aes(x = gene_order, y = frequency, fill = WGD_label)) +
#   geom_bar(stat = "identity") +
#   geom_text(aes(label = paste0(sprintf("%.1f", abs(frequency)), "%"),
#                 hjust = ifelse(frequency < 0, 1.05, -0.05)),
#             size = 3) +
#   geom_text(data = subset(plot_data2, WGD_label == "negative"),
#             aes(x = gene_order, y = 20, label = stat_label),
#             inherit.aes = FALSE,
#             hjust = 0, size = 3, color = "black") +
#   coord_flip() +
#   labs(title = "Top 10 Genes with Highest Amplification Frequency in WGD-positive Samples",
#        x = "Gene", y = "Amplification Frequency (%)", fill = "WGD Status") +
#   scale_fill_manual(values = c("positive" = "blue", "negative" = "red")) +
#   scale_y_continuous(labels = abs, expand = expansion(mult = c(0.05, 0.6))) +
#   theme_minimal() +
#   theme(axis.text.y = element_text(size = 8))

