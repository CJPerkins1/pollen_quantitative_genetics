# Introduction ------------------------------------------------------------
# The Manhattan plot code from kgwasflow is not working with the tomato genome, 
# so I'll make a Manhattan plot here using the aligned kmer output.

library(dplyr)
library(Rsamtools)
library(stringr)
library(ggplot2)


# Formatting the data -----------------------------------------------------
# Read the BAM file
kgwasflow_output <- scanBam(file = file.path(getwd(), "R_data", "flower_test_aligned_kmers.bam"))

# Formatting
rname_factor <- kgwasflow_output[[1]]$rname
rname_vector <- levels(rname_factor)[rname_factor]

df <- data.frame(
  kmer_name = kgwasflow_output[[1]]$qname,
  chr = rname_vector,
  chr_coord = kgwasflow_output[[1]]$pos,
  p_value = kgwasflow_output[[1]]$qname
)

df <- df %>%
  mutate(
    kmer_name = str_extract(kmer_name, str_extract(kmer_name, "^[^_]*")),
    chr = str_replace(chr, "SL4.0", ""),
    p_value = as.numeric(str_extract(p_value, "(?<=_)[0-9\\.e\\-]+")),
    log10_p_value = -log10(p_value)
  ) 

# Converting the coordinates to genome coordinates from chromosome, and getting a 
# list of the chromosome boundaries for plotting.
chrom_lengths <- data.frame(
  chr = sprintf("ch%02d", 0:12),
  length = c(9643250, 90863682, 53473368, 65298490, 64459972, 65269487, 
             47258699, 67883646, 63995357, 68513564, 64792705, 54379777, 66688036)
)

chromosome_info <- chrom_lengths %>%
  arrange(chr) %>% 
  filter(chr != "ch00") %>%
  mutate(end = cumsum(length),
         start = lag(end, default = 0) + 1) %>%
  select(chr, length, start, end)

chromosome_info <- chromosome_info %>%
  mutate(midpoint = (start + end) / 2)

df <- df %>%
  left_join(chromosome_info, by = c("chr")) %>%
  mutate(genomic_coord = start + chr_coord)


# Making the plot ---------------------------------------------------------
df %>% 
  ggplot(aes(x = genomic_coord, y = log10_p_value)) + 
  geom_point(size = 2) +
  geom_vline(data = chromosome_info, aes(xintercept = start), linewidth = 0.2, color = "gray", linetype = "dashed") +
  geom_vline(data = chromosome_info, aes(xintercept = end), linewidth = 0.2, color = "gray", linetype = "dashed") +
  scale_x_continuous(breaks = chromosome_info$midpoint, 
                     labels = chromosome_info$chr,
                     limits = c(1, 772876783),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 12, 2),
                     labels = seq(0, 12, 2),
                     limits = c(0, 13)) +
  labs(title = "Anther/pistil ratio kmers-GWAS", x = "Chromosome", y = "-log10 p-value") +
  theme_bw() +
  theme(axis.title = element_text(size = 26, face = 'bold'),
        axis.text = element_text(size = 18, face = 'bold', color = 'black'),
        axis.text.x = element_text(size = 18, face = 'bold', color = 'black'),
        plot.title = element_text(size = 28, face = 'bold', margin = margin(0, 0, 10, 0)),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 1, color = 'black'),
        axis.ticks = element_line(linewidth = 1, color = 'black'), axis.ticks.length = unit(8, 'pt'),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm'),
        panel.grid = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 18, face = 'bold', color = 'black'),
        legend.text = element_text(size = 14, face = 'bold', color = 'black'),
        legend.key.width = unit(1.5, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside")

ggsave(filename = file.path(getwd(), "R_plots", "pistil_anther_ratio_kmers_gwas.png"),
       device = 'png',
       width = 14,
       height = 8,
       dpi = 400,
       units = 'in')


# Making locule number plot -----------------------------------------------
# Reading the locule bam file
kgwasflow_locule_output <- scanBam(file = file.path(getwd(), "R_data", "locule_number_kmers_alignment.bam"))

# Formatting
locule_rname_factor <- kgwasflow_locule_output[[1]]$rname
locule_rname_vector <- levels(locule_rname_factor)[locule_rname_factor]

locule_df <- data.frame(
  kmer_name = kgwasflow_locule_output[[1]]$qname,
  chr = locule_rname_vector,
  chr_coord = kgwasflow_locule_output[[1]]$pos, 
  p_value = kgwasflow_locule_output[[1]]$qname
)

locule_df <- locule_df %>%
  mutate(
    kmer_name = str_extract(kmer_name, str_extract(kmer_name, "^[^_]*")),
    chr = str_replace(chr, "SL4.0", ""),
    p_value = as.numeric(str_extract(p_value, "(?<=_)[0-9\\.e\\-]+")),
    log10_p_value = -log10(p_value)
  ) 

locule_df <- locule_df %>%
  left_join(chromosome_info, by = c("chr")) %>%
  mutate(genomic_coord = start + chr_coord)

locule_df %>% 
  ggplot(aes(x = genomic_coord, y = log10_p_value)) + 
  geom_point(size = 2, alpha = 0.2) +
  geom_vline(data = chromosome_info, aes(xintercept = start), linewidth = 0.2, color = "gray", linetype = "dashed") +
  geom_vline(data = chromosome_info, aes(xintercept = end), linewidth = 0.2, color = "gray", linetype = "dashed") +
  scale_x_continuous(breaks = chromosome_info$midpoint, 
                     labels = chromosome_info$chr,
                     limits = c(1, 772876783),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 18, 2),
                     labels = seq(0, 18, 2),
                     limits = c(0, 18)) +
  labs(title = "Locule number kmers-GWAS", x = "Chromosome", y = "-log10 p-value") +
  theme_bw() +
  theme(axis.title = element_text(size = 26, face = 'bold'),
        axis.text = element_text(size = 18, face = 'bold', color = 'black'),
        axis.text.x = element_text(size = 18, face = 'bold', color = 'black'),
        plot.title = element_text(size = 28, face = 'bold', margin = margin(0, 0, 10, 0)),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 1, color = 'black'),
        axis.ticks = element_line(linewidth = 1, color = 'black'), axis.ticks.length = unit(8, 'pt'),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm'),
        panel.grid = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 18, face = 'bold', color = 'black'),
        legend.text = element_text(size = 14, face = 'bold', color = 'black'),
        legend.key.width = unit(1.5, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside")

ggsave(filename = file.path(getwd(), "R_plots", "locule_number_gwas.png"),
       device = 'png',
       width = 14,
       height = 8,
       dpi = 400,
       units = 'in')

locule_df %>% 
  filter(chr == "ch02") %>%
  ggplot(aes(x = chr_coord, y = log10_p_value)) + 
  geom_point(size = 2, alpha = 0.2) +
  geom_vline(data = chromosome_info, aes(xintercept = start), linewidth = 0.2, color = "gray", linetype = "dashed") +
  geom_vline(data = chromosome_info, aes(xintercept = end), linewidth = 0.2, color = "gray", linetype = "dashed") +
  # scale_x_continuous(breaks = chromosome_info$midpoint, 
  #                    labels = chromosome_info$chr,
  #                    limits = c(135000000, 139000000),
  scale_x_continuous(limits = c(44700000, 46000000),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 18, 2),
                     labels = seq(0, 18, 2),
                     limits = c(0, 18)) +
  labs(title = "Locule number kmers-GWAS", x = "Chromosome 2", y = "-log10 p-value") +
  theme_bw() +
  theme(axis.title = element_text(size = 26, face = 'bold'),
        axis.text = element_text(size = 18, face = 'bold', color = 'black'),
        axis.text.x = element_text(size = 18, face = 'bold', color = 'black'),
        plot.title = element_text(size = 28, face = 'bold', margin = margin(0, 0, 10, 0)),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 1, color = 'black'),
        axis.ticks = element_line(linewidth = 1, color = 'black'), axis.ticks.length = unit(8, 'pt'),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm'),
        panel.grid = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 18, face = 'bold', color = 'black'),
        legend.text = element_text(size = 14, face = 'bold', color = 'black'),
        legend.key.width = unit(1.5, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside")

ggsave(filename = file.path(getwd(), "R_plots", "locule_number_gwas_ch02.png"),
       device = 'png',
       width = 14,
       height = 8,
       dpi = 400,
       units = 'in')
