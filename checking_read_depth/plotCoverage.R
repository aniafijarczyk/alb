library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

dcov <- read.table("snp_bcftools_cov.tab", sep="\t", head=TRUE)
head(dcov)

meancov = sum(dcov$DP4)/length(dcov$DP4)
meancov

summary(dcov) 
dg <- dcov %>% gather(key='pool',value='cov',-chrom,-pos,-pos_cum,-length,-DP4)
head(dg)
dgg <- dg %>% group_by(pool) %>% summarise(mean = mean(cov))
dgg %>% summarise(mean(mean),mode=getmode(mean))

pool_hi = mean(dgg$mean) + (4*sqrt(mean(dgg$mean)))
pool_hi*10


p1 <- ggplot(dcov) + aes(x = pos_cum, y = DP4, colour = as.factor(chrom)) +
  geom_point(pch=21) + 
  scale_y_continuous(trans = 'log10') +
  geom_hline(yintercept = c(meancov,pool_hi*10), lty=c('dashed','solid'), colour = c("grey40","black")) +
  labs(x = 'scaffold position (bp)', y = 'n reads') +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.title = element_text(size=14),
        #axis.text = element_text(size=10),
        legend.position = "none",
        strip.text = element_text(size=12),
        strip.background = element_blank())
p1


#png('plotCoverage_distribution.png',width = 1000,height = 700,res=150)
#p1
#dev.off()


# Plotting histogram of coverage across pools

text1 <- data.frame(x=c(350,800), y=c(1900,1900),
                    label=c(paste0('Mean = ',round(meancov,2)),
                            paste0('Higher limit = ',round(pool_hi*10,2),' \n[mean(pool)+4*sqrt(mean(pool))]*10')))
text1
                    
p2 <- ggplot(dcov) + aes(x = DP4) +
  geom_histogram(bins=40) + 
  scale_x_continuous(limits=c(0,1100)) +
  #geom_hline(yintercept = c(cov_lo,meancov,cov_hi), colour = "black") +
  geom_vline(xintercept = c(meancov,pool_hi*10),linetype=c('dashed','solid'),colour=c("grey80","black")) +
  geom_text(data=text1,aes(x=x,y=y,label=label),color='black') +
  labs(x = 'n reads', y = 'frequency') +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.title = element_text(size=14),
        #axis.text = element_text(size=10),
        legend.position = "none",
        strip.text = element_text(size=12),
        strip.background = element_blank())
p2

#png('plotCoverage_totaldensity.png',w=900,h=600,res=150)
#p2
#dev.off()


# Plotting histograms of coverage for each pool

text2 = data.frame(x = 65, y = 0.037, label = paste('Mean = ',round(mean(dgg$mean),2)))
text2

p3 <- ggplot(dg) + aes(x = cov, colour = pool) +
  geom_vline(xintercept=44,linetype='dashed',colour='grey80') +
  geom_density() +
  scale_color_viridis(discrete=TRUE,option='viridis') +
  geom_text(data=text2,aes(x=x,y=y,label=label),colour="black") +
  #scale_color_gradientn(colors=c('goldenrod','skyblue','firebrick'))
  labs(x = 'n reads', y = 'frequency') +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.title = element_text(size=14))
p3


#png('plotCoverage_density.png',w=900,h=600,res=150)
#p3
#dev.off()


png('plotCoverage_combined.png',w=1800,h=1000,res=150)
bottom_row <- plot_grid(p2,p3,ncol=2)
plot_grid(p1,bottom_row,ncol=1)
dev.off()
