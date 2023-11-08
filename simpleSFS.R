library("ggplot2")
library("tidyverse")
library("dplyr")

setwd("/Users/Oscar/Desktop/Emerson_Lab_Work/SV-project")
t1=read.table("extraction.tsv",header=FALSE)
#colnames(t1) <- c("CHROM","POS", "SVTYPE","SVLEN", S[1:55])
#t1=t1 %>% mutate(V5 = str_replace(V5, './.', '0'))
#t1 %>%  mutate(conf = recode(conf, 'East' = 'E', 'West' = 'W', 'North' = 'N'))
t1=t1 %>%mutate(across(V5:V59, ~ recode(.x, './.' = "0", '0/1' = "1", '1/0' = "1", '1/1' = "1")))
t1$V60 <- apply(t1, 1, function(x) length(which(x==1)))
levels(as.factor(t1$V3))

##polarizing failed. so we keep going on
tdel = t1 %>% filter(V3=="DEL") %>% select(V1,V2,V3,V4,V60) 
sfsdel=as.data.frame(table(tdel$V60))
colnames(sfsdel)<-c("frequency","occurrence")
#sfsdel= sfsdel %>% mutate(frequency=as.character(frequency))
pdel <- ggplot(sfsdel, aes(x=frequency, y=occurrence)) + 
  geom_bar(stat = "identity")+ scale_y_log10()  + 
  ggtitle("DELETION, log")+ scale_color_grey() + theme_classic()
pdel


t.ins = t1 %>% filter(V3=="INS") %>% select(V1,V2,V3,V4,V60) 
sfs.ins=as.data.frame(table(t.ins$V60))
colnames(sfs.ins)<-c("frequency","occurrence")
p.ins <- ggplot(sfs.ins, aes(x=frequency, y=occurrence)) + 
  geom_bar(stat = "identity")+ scale_y_log10() + 
  ggtitle("INSERTION, log")+ scale_color_grey() + theme_classic()
p.ins


t.dup = t1 %>% filter(V3=="DUP") %>% select(V1,V2,V3,V4,V60) 
sfs.dup=as.data.frame(table(t.dup$V60))
colnames(sfs.dup)<-c("frequency","occurrence")
p.dup <- ggplot(sfs.dup, aes(x=frequency, y=occurrence)) + 
  geom_bar(stat = "identity")+
  ggtitle("DUPLICATION")+ scale_color_grey() + theme_classic()
p.dup


tinv = t1 %>% filter(V3=="INV") %>% select(V1,V2,V3,V4,V60) 
sfsinv=as.data.frame(table(tinv$V60))
colnames(sfsinv)<-c("frequency","occurrence")
pinv <- ggplot(sfsinv, aes(x=frequency, y=occurrence)) + 
  geom_bar(stat = "identity") + 
  ggtitle("INVERSION")+ scale_color_grey() + theme_classic()
pinv

t.bnd = t1 %>% filter(V3=="BND") %>% select(V1,V2,V3,V4,V60) 
sfs.bnd=as.data.frame(table(t.bnd$V60))
colnames(sfs.bnd)<-c("frequency","occurrence")
p.bnd <- ggplot(sfs.bnd, aes(x=frequency, y=occurrence)) + 
  geom_bar(stat = "identity")  + ggtitle("BND") + scale_color_grey() + theme_classic()
p.bnd
