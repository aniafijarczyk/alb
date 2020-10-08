setwd("/media/anna/Volume/ALB/2020_10_snpcalling/01_test")

library(ggplot2)
library(dplyr)


dcov <- read.table("pool1_cov.tab", sep="\t", head=TRUE)
head(dcov)




ggplot(dcov) + aes(x = pos_cum, y = DP4, colour = as.factor(chrom)) +
  geom_hline(yintercept = c(12,70), colour = "black") +
  geom_point() + 
  scale_y_continuous(trans = 'log2') +
  scale_colour_manual(values = c("grey")) +
  labs(x = 'scaffold position (bp)', y = 'n reads') +
  #facet_wrap(~species,ncol=1) +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.title = element_text(size=14),
        #axis.text = element_text(size=10),
        legend.position = "bottom",
        strip.text = element_text(size=12),
        strip.background = element_blank())

