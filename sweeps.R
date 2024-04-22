library("ggplot2")
library("tidyverse")
library(RColorBrewer)
library(khroma)
setwd("/Users/Oscar/Desktop/Emerson_Lab_Work/SV-project")
t.sv=read.table("2R.SV.asm.output.txt",header = TRUE)
t.snp=read.table("2R.snp.output.txt",header = TRUE)

p.sv = ggplot(data=t.sv ,mapping = aes(x=location,y=LR))
p.sv = p.sv + geom_point()+annotate(geom = "text", x=13399000, y=15, label="Cyp6g1" )+
  geom_vline(aes(xintercept=12187000))+
  labs(title = "sweepfinder2, 2R, SV-SV")
p.sv

p.snp = ggplot(data=t.snp ,mapping = aes(x=location,y=LR))
p.snp = p.snp + geom_point()+annotate(geom = "text", x=13399000, y=100, label="Cyp6g1" )+
geom_vline(aes(xintercept=12187000))+
  labs(title = "sweepfinder2, 2R, snp-snp")
p.snp

p = ggplot(data=t.sv ,mapping = aes(x=location,y=LR))+
  geom_point(color="blue")+
  geom_point(data = t.snp,color="red")
p