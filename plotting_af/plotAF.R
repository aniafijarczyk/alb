library(ggplot2)
library(dplyr)

### Reading and merging tables with pool info

pools <- c('Pool10','pool1','pool2','pool3','Pool4',
           'Pool5','Pool6','Pool7','Pool8','Pool9')
pool_start <- pools[1]
dp <- read.table(paste0('snp_bcftools_annotated.f2.recode_',pool_start,'.tab'),sep='\t',header=FALSE,
                 col.names=c('chrom','pos','dp','adRef','adAlt','gt'))
dp_$af <- dp_$adAlt/(dp_$adAlt+dp_$adRef)
dp_$pool <- pool_start
dp_samp <- sample_n(data.frame(dp_), 10000)
DP <- dp_samp

for (i in 2:10) {
  dp <- read.table(paste0('snp_bcftools_annotated.f2.recode_',pools[i],'.tab'),sep='\t',header=FALSE,
                   col.names=c('chrom','pos','dp','adRef','adAlt','gt'))
  dp_$af <- dp_$adAlt/(dp_$adAlt+dp_$adRef)
  dp_$pool <- pools[i]
  dp_samp <- sample_n(data.frame(dp_), 10000)
  DP <- rbind(DP,dp_samp)
}

dim(DP)
head(DP)


### Plotting AFS per pool

DP_vars <- DP %>% filter(af>0.0 & af<1.0 & gt != './.')

p2 <- ggplot(DP_vars) + aes(x = af,colour = gt) +
  geom_freqpoly(alpha=0.75,bins=75) +
  scale_colour_manual(values=c('#F1BB7B','#FD6467','#5B1A18')) +
  facet_wrap(~pool,ncol=5) +
  theme(panel.background = element_blank())
p2

png('plotAF_af_perPool.png',width = 2000,height = 800,res=150)
p2
dev.off()


### Plotting AFS from all pools


p3 <- ggplot(DP_vars) + aes(x = af) +
  geom_freqpoly(alpha=0.75,bins=75,colour='#5B1A18') +
  #facet_wrap(~pool,ncol=5) +
  theme(panel.background = element_blank())
p3

png('plotAF_af_pool_together.png',width = 1200,height = 800,res=150)
p3
dev.off()


