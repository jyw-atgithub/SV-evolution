library("ggplot2")
library("tidyverse")
library(RColorBrewer)
library(khroma)
setwd("/Users/Oscar/Desktop/Emerson_Lab_Work/SV-project")
#it is tricky that "-nan" is interpreted as "NaN"
t1=read.table("fst_snp_AF-EU.windowed.weir.fst",header=TRUE) #%>% filter(N_VARIANTS >5)

t1.2R=t1%>%filter(CHROM=="2R") #%>%filter(BIN_END > 5000000 & BIN_END < 9000000 )
p1.2R=ggplot(data=t1.2R,aes(x=BIN_END, y=MEAN_FST))
p1.2R=p1.2R+ geom_point(color="lightblue")+
  geom_smooth(method = "lm",data=t1.2R,aes(x=BIN_END, y=MEAN_FST), color="#774150" ) +
  geom_smooth(method = "gam",data=t1.2R,aes(x=BIN_END, y=MEAN_FST), color="#774150" )
p1.2R

t2=read.table("fst_snp_EU-AM.windowed.weir.fst",header=TRUE) %>% filter(N_VARIANTS >5)
t2.2R=t2%>%filter(CHROM=="2R") #%>%filter(BIN_END > 5000000 & BIN_END < 9000000 )
p2.2R=ggplot(data=t2.2R,aes(x=BIN_END, y=MEAN_FST))
p2.2R=p2.2R+ geom_point(color="red3") + 
  geom_smooth(method = lm,data=t2.2R,aes(x=BIN_END, y=MEAN_FST), color="#00FF36" ) 
p2.2R