library("ggplot2")
library("tidyverse")
library("dplyr")
library("gridExtra")
setwd("/home/jenyu/Desktop/SV-project")
t1=read.table("completeness.tsv",header=FALSE)
colnames(t1)<-c("score","strain","assembler","polisher","program")

t2=t1 %>% filter(program=="busco") %>% filter(polisher!="racontmp") %>% select(-c(strain))
p1<-ggplot(data=t2, aes(x=polisher,y=score)) +geom_boxplot(aes(fill=assembler, color=assembler))
p1

t3=t1 %>% filter(program=="compleasm") %>% filter(polisher!="racontmp") %>% select(-c(strain))
p3<-ggplot(data=t3, aes(x=polisher,y=score)) +geom_boxplot(aes(fill=assembler, color=assembler))
p3


png(filename="score.png",width=2000, height=1000)
g1 <-grid.arrange(p1,p3,ncol=2)
dev.off()