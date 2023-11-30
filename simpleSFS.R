library("ggplot2")
library("tidyverse")
library("dplyr")

setwd("/Users/Oscar/Desktop/Emerson_Lab_Work/SV-project")
t1=read.table("extraction2.tsv",header=FALSE)
#colnames(t1) <- c("CHROM","POS", "SVTYPE","SVLEN", S[1:55])
#t1=t1 %>% mutate(V5 = str_replace(V5, './.', '0'))
#t1 %>%  mutate(conf = recode(conf, 'East' = 'E', 'West' = 'W', 'North' = 'N'))
t1=t1 %>%mutate(across(V5:V65, ~ recode(.x, './.' = "0", '0/1' = "1", '1/0' = "1", '1/1' = "1")))
t1$V66 <- apply(t1, 1, function(x) length(which(x==1)))
levels(as.factor(t1$V3))

##polarizing failed. so we keep going on
tdel = t1 %>% filter(V3=="DEL") %>% select(V1,V2,V3,V4,V66) 
sfsdel=as.data.frame(table(tdel$V66))
colnames(sfsdel)<-c("frequency","occurrence")
sfsdel$percentage=(sfsdel$occurrence*100/nrow(tdel))
sfsdel$type="DEL"
#sfsdel= sfsdel %>% mutate(frequency=as.character(frequency))
pdel <- ggplot(sfsdel, aes(x=frequency, y=percentage)) + 
  geom_bar(stat = "identity") + 
  ggtitle("DELETION, %")+ scale_color_grey() + theme_classic()
pdel


t.ins = t1 %>% filter(V3=="INS") %>% select(V1,V2,V3,V4,V66) 
sfs.ins=as.data.frame(table(t.ins$V66))
colnames(sfs.ins)<-c("frequency","occurrence")
sfs.ins$percentage=(sfs.ins$occurrence*100/nrow(t.ins))
sfs.ins$type="INS"

t.dup = t1 %>% filter(V3=="DUP") %>% select(V1,V2,V3,V4,V66) 
sfs.dup=as.data.frame(table(t.dup$V66))
colnames(sfs.dup)<-c("frequency","occurrence")
sfs.dup$percentage=(sfs.dup$occurrence*100/nrow(t.dup))
sfs.dup$type="DUP"

t.inv = t1 %>% filter(V3=="INV") %>% select(V1,V2,V3,V4,V66) 
sfs.inv=as.data.frame(table(t.inv$V66))
colnames(sfs.inv)<-c("frequency","occurrence")
sfs.inv$percentage=(sfs.inv$occurrence*100/nrow(t.inv))
sfs.inv$type="INV"

t.bnd = t1 %>% filter(V3=="BND") %>% select(V1,V2,V3,V4,V66) 
sfs.bnd=as.data.frame(table(t.bnd$V66))
colnames(sfs.bnd)<-c("frequency","occurrence")
sfs.bnd$percentage=(sfs.bnd$occurrence*100/nrow(t.bnd))
sfs.bnd$type="BND"


sfs.all=bind_rows(sfsdel,sfs.ins,sfs.dup,sfs.inv,sfs.bnd)

pall <- ggplot(sfs.all, aes(x=frequency, y=percentage, fill=type)) + 
  geom_bar(stat = "identity",position="dodge") +
  ggtitle("all, %")+ scale_color_grey() + theme_classic()
pall
