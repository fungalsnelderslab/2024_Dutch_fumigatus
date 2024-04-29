library(ggplot2)
library(cowplot)

clin_env <- read.table("/Users/ben/Desktop/baseclear/fst/fst_clin_env.windowed.weir.fst",header=T)
clin_env$POS <- (clin_env$BIN_START+clin_env$BIN_END)/2

clin_env$CHROM[clin_env$CHROM == "NC_007194.1"] <- "1"
clin_env$CHROM[clin_env$CHROM == "NC_007195.1"] <- "2"
clin_env$CHROM[clin_env$CHROM == "NC_007196.1"] <- "3"
clin_env$CHROM[clin_env$CHROM == "NC_007197.1"] <- "4"
clin_env$CHROM[clin_env$CHROM == "NC_007198.1"] <- "5"
clin_env$CHROM[clin_env$CHROM == "NC_007199.1"] <- "6"
clin_env$CHROM[clin_env$CHROM == "NC_007200.1"] <- "7"
clin_env$CHROM[clin_env$CHROM == "NC_007201.1"] <- "8"

res_susc <- read.table("/Users/ben/Desktop/baseclear/fst/fst_res_susc.windowed.weir.fst",header=T)
res_susc$POS <- (res_susc$BIN_START+res_susc$BIN_END)/2

res_susc$CHROM[res_susc$CHROM == "NC_007194.1"] <- "1"
res_susc$CHROM[res_susc$CHROM == "NC_007195.1"] <- "2"
res_susc$CHROM[res_susc$CHROM == "NC_007196.1"] <- "3"
res_susc$CHROM[res_susc$CHROM == "NC_007197.1"] <- "4"
res_susc$CHROM[res_susc$CHROM == "NC_007198.1"] <- "5"
res_susc$CHROM[res_susc$CHROM == "NC_007199.1"] <- "6"
res_susc$CHROM[res_susc$CHROM == "NC_007200.1"] <- "7"
res_susc$CHROM[res_susc$CHROM == "NC_007201.1"] <- "8"

cyp51A <- data.frame(CHROM="4",
                     POS=1780204,
                     text="cyp51A")

top <- ggplot(clin_env) + 
  geom_segment(data=cyp51A,aes(x=POS,xend=POS,y=0,yend=0.55),inherit.aes = FALSE,color="grey30")+
  geom_text(data=cyp51A,parse=T,aes(x=POS,y=0.58,label='paste(italic("cyp51"),"A")'))+
  geom_point(aes(x=POS,y=WEIGHTED_FST,color=CHROM),size=0.3) + 
  facet_grid(cols=vars(CHROM),scales="free",space="free") +
  theme_classic(base_size=14) + ylim(-0.01,0.61) +
  geom_hline(aes(yintercept=0.3),lty=2,lwd=0.2)+
  theme(legend.position="none",
        panel.spacing=unit(0.01,"cm"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank()) +
  scale_x_continuous(breaks=c(1e6,2e6,3e6,4e6,5e6),labels=c(1,2,3,4,5)) +
  labs(x="Position (Mb)",y=expression(atop(Clin~vs.~Env,(Weighted~F[st]))))

bottom <- ggplot(res_susc) + 
  geom_segment(data=cyp51A,aes(x=POS,xend=POS,y=0,yend=0.55),inherit.aes = FALSE,color="grey30")+
  geom_point(aes(x=POS,y=WEIGHTED_FST,color=CHROM),size=0.3) + 
  geom_text(data=cyp51A,parse=T,aes(x=POS,y=0.58,label='paste(italic("cyp51"),"A")'))+
  facet_grid(cols=vars(CHROM),scales="free",space="free") +
  theme_classic(base_size=14) + ylim(-0.01,0.61) +
  geom_hline(aes(yintercept=0.3),lty=2,lwd=0.2)+
  theme(legend.position="none",
        panel.spacing = unit(0.01,"cm"),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_x_continuous(breaks=c(1e6,2e6,3e6,4e6,5e6),labels=c(1,2,3,4,5))+
  labs(x="Position (Mb)",y=expression(atop(Res~vs.~Sens,(Weighted~F[st]))))


plot_grid(top,bottom,nrow=2,rel_heights = c(0.5,0.53))
ggsave("/Users/ben/Desktop/baseclear/fst/dutch_fst.svg",width=10,height=4.5)



