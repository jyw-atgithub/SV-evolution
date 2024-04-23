library("ggplot2")
library("tidyverse")
library(RColorBrewer)
library(khroma)
setwd("/Users/Oscar/Desktop/Emerson_Lab_Work/SV-project")
#it is tricky that "-nan" is interpreted as "NaN"
#Discard the windows containing too few variable sites
#Negative Fst values are turned into zero

t1=read.table("fst_snp_AF-EU.windowed.weir.fst",header=TRUE) %>% filter(N_VARIANTS >5) %>%
  mutate( WEIGHTED_FST= replace(WEIGHTED_FST, WEIGHTED_FST < 0 ,0)) %>%
  mutate(pops = "AF-EU")
t2=read.table("fst_snp_EU-AM.windowed.weir.fst",header=TRUE) %>% filter(N_VARIANTS >5)%>%
  mutate( WEIGHTED_FST= replace(WEIGHTED_FST, WEIGHTED_FST < 0 ,0))%>%
  mutate(pops = "EU-AM")

t1.2R=t1%>%filter(CHROM=="2R") #%>%filter(BIN_END > 5000000 & BIN_END < 9000000 )
#p1.2R=ggplot(data=t1.2R, aes( x= BIN_END , y = WEIGHTED_FST ))
#p1.2R=p1.2R+ geom_point(color="lightblue")+  geom_smooth(method = "lm",data=t1.2R,aes(x=BIN_END, y=WEIGHTED_FST), color="#774150" ) 
  #+geom_smooth(method = "gam",data=t1.2R,aes(x=BIN_END, y=WEIGHTED_FST), color="#774150" )
#p1.2R
t2.2R=t2%>%filter(CHROM=="2R") #%>%filter(BIN_END > 5000000 & BIN_END < 9000000 )


p.a=ggplot(data=t1.2R, aes( x= BIN_END , y = WEIGHTED_FST )) +
  geom_point(aes(color = pops)) +
  geom_smooth(method = "lm",data=t1.2R,aes(x=BIN_END, y=WEIGHTED_FST, color = pops) ) +
  geom_point(data = t2.2R, aes(x = BIN_END, y = WEIGHTED_FST, color = pops )) +
  geom_smooth(method = "lm", data=t2.2R, aes(x=BIN_END, y=WEIGHTED_FST, color = pops) ) +
  labs(title="2R, Fst") +theme_light()
p.a

t3=read.table("fst_SV_AF-EU_2.windowed.weir.fst",header=TRUE) %>% filter(N_VARIANTS >5) %>%
  mutate( WEIGHTED_FST= replace(WEIGHTED_FST, WEIGHTED_FST < 0 ,0))%>%
  mutate(pops = "AF-EU")
t4=read.table("fst_SV_EU-AM_2.windowed.weir.fst",header=TRUE) %>% filter(N_VARIANTS >5) %>%
  mutate( WEIGHTED_FST= replace(WEIGHTED_FST, WEIGHTED_FST < 0 ,0))%>%
  mutate(pops = "EU-AM")
t3.2R=t3%>%filter(CHROM=="2R") 
t4.2R=t4%>%filter(CHROM=="2R")
p.b=ggplot(data=t3.2R, aes( x= BIN_END , y = WEIGHTED_FST )) +
  geom_point(aes(color = pops)) +
  geom_smooth(method = "lm",data=t3.2R,aes(x=BIN_END, y=WEIGHTED_FST, color = pops) ) +
  geom_point(data = t4.2R, aes(x = BIN_END, y = WEIGHTED_FST, color = pops )) +
  geom_smooth(method = "lm", data=t4.2R, aes(x=BIN_END, y=WEIGHTED_FST, color = pops) ) +
  labs(title="2R, Weighted Fst from SV") +theme_light()
p.b
