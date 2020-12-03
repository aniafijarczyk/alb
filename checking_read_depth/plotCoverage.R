library(ggplot2)
library(dplyr)


dcov <- read.table("snp_bcftools_cov.tab", sep="\t", head=TRUE)
head(dcov)

meancov = sum(dcov$DP4)/length(dcov$DP4)
cov_hi = meancov*2
cov_lo = meancov/10

cov_hi = 1000
cov_lo = 40

p1 <- ggplot(dcov) + aes(x = pos_cum, y = DP4, colour = as.factor(chrom)) +
  geom_point() + 
  scale_y_continuous(trans = 'log10') +
  geom_hline(yintercept = c(cov_lo,meancov,cov_hi), colour = "black") +
  labs(x = 'scaffold position (bp)', y = 'n reads') +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.title = element_text(size=14),
        legend.position = "none",
        strip.text = element_text(size=12),
        strip.background = element_blank())
p1

png('plotCoverage.png',width = 1000,height = 700,res=150)
p1
dev.off()
