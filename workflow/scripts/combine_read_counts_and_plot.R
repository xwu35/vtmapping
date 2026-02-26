#!/usr/bin/env Rscript

# this script is used to combine read counts from different steps and plot the composition of reads

suppressPackageStartupMessages(library("tidyverse"))

# define input
after_human <- snakemake@input$after_human
after_host <- snakemake@input$after_host

# after human
counts_after_human <- read.table(after_human, 
                                   sep = "\t", header = T) %>% 
  mutate(sample = sub("_R1.*", "", basename(file)),
         total_reads_after_human = num_seqs * 2) %>%
  select(sample, total_reads_after_human)

# after _host
counts_after_host <- read.table(after_host, 
                                sep = " ", header = T) %>%
  rename(unmapped = "read_counts") 

# merge all tables
df_lists = list(counts_after_human, counts_after_host)

combined_table <- Reduce(function(x, y) merge(x, y, by = "sample", all = TRUE), df_lists)

# calculate number of removed mapped reads
combined_table <- combined_table %>% 
  mutate(mapped = total_reads_after_human - unmapped)

# barplot of removed reads at each step
# define colors
colors <- c("#56B4E9", "#E69F00", "#009E73", "#999999", "#CC79A7")

# plot
p_composition_reads <- combined_table %>%
  select(sample, mapped, unmapped) %>%
  pivot_longer(!sample, names_to = "type", values_to = "count") %>% 
  ggplot(aes(x=sample, y=count, fill=type)) +
  geom_bar(position="fill", stat="identity") +
  theme_classic() +
  theme(axis.title = element_text(size = 11),
        axis.line = element_line(colour = "black"), 
        axis.text.x=element_text(size=8, colour = "black", angle=90, vjust=0.5),
        axis.text.y = element_text(size=10, colour = "black"),
        axis.ticks.x=element_blank(), 
        legend.title=element_blank(),
        legend.text=element_text(size=8),
        legend.position = "bottom") + 
  guides(fill=guide_legend(nrow=1)) +
  labs(x="", y="Proportion of reads") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = colors) # + coord_flip()  

# select few columns to write out
output_table <- combined_table %>%
  mutate(pct_unmapped = (unmapped/total_reads_after_human)*100,
         pct_mapped = (mapped/total_reads_after_human)*100)
  
# save table
write.table(output_table, snakemake@output$table,
            sep="\t", row.names = F, col.names=T, quote = F)

# save figure
ggsave(file=snakemake@output$figure,
       plot=p_composition_reads, 
       width=9, height=5)
