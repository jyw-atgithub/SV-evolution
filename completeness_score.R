library("ggplot2")
library("tidyverse")
library("dplyr")
library("gridExtra")
setwd("/Users/Oscar/Desktop/Emerson_Lab_Work/SV-project")
t1=read.table("completeness.tsv",header=FALSE)
colnames(t1)<-c("score","strain","assembler","polisher","program")
levels(factor(t1$polisher))
levels(factor(t1$strain))

t2=t1 %>% filter(program=="busco") %>% filter(polisher!="racontmp") %>% filter(polisher!="polished") %>% filter(assembler!="nextdenovo-30") %>% select(-c(strain))

t2=t1 %>% filter(program=="busco") %>% filter(polisher!="racontmp") %>% filter(polisher!="polished") %>% filter(assembler!="nextdenovo-30") %>% filter(strain!="A1"&strain!="A2"&strain!="A3"&strain!="A4"&strain!="A5"&strain!="A6"&strain!="A7"&strain!="AB8"&strain!="B1"&strain!="B2"&strain!="B3"&strain!="B4"&strain!="B6"&strain!="ORE")

t2$polisher <- factor(t2$polisher , levels=c("original", "ori-purged", "racon", "nextpolish", "polished-purged"))
p1<-ggplot(data=t2, aes(x=polisher,y=score)) +geom_boxplot(aes(fill=assembler, color=assembler))
p1

t3=t1 %>% filter(program=="compleasm") %>% filter(polisher!="racontmp") %>% select(-c(strain))
p3<-ggplot(data=t3, aes(x=polisher,y=score)) +geom_boxplot(aes(fill=assembler, color=assembler))
p3


png(filename="score.png",width=2000, height=1000)
g1 <-grid.arrange(p1,p3,ncol=2)
dev.off()