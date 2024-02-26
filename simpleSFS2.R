library("ggplot2")
library("tidyverse")

setwd("/Users/Oscar/Desktop/Emerson_Lab_Work/SV-project")
t1=read.table("repeat_type_genotype.tsv", header=TRUE)%>% 
  mutate(across(V4:V65, ~ recode(.x, './.' = "0", '0/1' = "1", '1/0' = "1", '1/1' = "1")))
t1$AF <- apply(t1, 1, function(x) length(which(x==1)))

levels(as.factor(t1$SUPERFAMILY))

#t.Gypsy = t1 %>% filter(grepl("Gypsy", ORDER)) %>% select(CHROM, POS, ID, ORDER, SUPERFAMILY, AF)
t.Gypsy = t1 %>% filter(grepl("Gypsy", SUPERFAMILY)) %>% select(CHROM, POS, ID, ORDER, SUPERFAMILY, AF)
sfs.Gypsy=as.data.frame(table(t.Gypsy$AF))
colnames(sfs.Gypsy)<-c("frequency","occurrence")
sfs.Gypsy$ratio=(sfs.Gypsy$occurrence*10000/sum(sfs.Gypsy$occurrence))
p.Gypsy <- ggplot(sfs.Gypsy, aes(x=frequency, y=ratio)) + 
  geom_bar(stat = "identity") + 
  ggtitle("Gypsy, %")+ scale_color_grey() + theme_classic()+scale_y_continuous(trans = 'log2')
p.Gypsy

t.pao = t1 %>% filter(grepl("Pao", SUPERFAMILY)) %>% select(CHROM, POS, ID, ORDER, SUPERFAMILY, AF)
sfs.pao=as.data.frame(table(t.pao$AF))
colnames(sfs.pao)<-c("frequency","occurrence")
sfs.pao$ratio=(sfs.pao$occurrence*10000/sum(sfs.pao$occurrence))
p.pao <- ggplot(sfs.pao, aes(x=frequency, y=ratio)) + 
  geom_bar(stat = "identity") + 
  ggtitle("pao, %")+ scale_color_grey() + theme_classic()+scale_y_continuous(trans = 'log2')
p.pao

t.Helitron = t1 %>% filter(grepl("Helitron", SUPERFAMILY)) %>% select(CHROM, POS, ID, ORDER, SUPERFAMILY, AF)
sfs.Helitron=as.data.frame(table(t.Helitron$AF))
colnames(sfs.Helitron)<-c("frequency","occurrence")
sfs.Helitron$ratio=(sfs.Helitron$occurrence*10000/sum(sfs.Helitron$occurrence))
p.Helitron <- ggplot(sfs.Helitron, aes(x=frequency, y=ratio)) + 
  geom_bar(stat = "identity") + 
  ggtitle("Helitron, %")+ scale_color_grey() + theme_classic()+scale_y_continuous(trans = 'log2')
p.Helitron

t.tc1 = t1 %>% filter(grepl("Tc1", SUPERFAMILY)) %>% select(CHROM, POS, ID, ORDER, SUPERFAMILY, AF)
sfs.tc1=as.data.frame(table(t.tc1$AF))
colnames(sfs.tc1)<-c("frequency","occurrence")
sfs.tc1$ratio=(sfs.tc1$occurrence*10000/sum(sfs.tc1$occurrence))
p.tc1 <- ggplot(sfs.tc1, aes(x=frequency, y=ratio)) + 
  geom_bar(stat = "identity") + 
  ggtitle("tc1, %")+ scale_color_grey() + theme_classic()+scale_y_continuous(trans = 'log2')
p.tc1

t.r1 = t1 %>% filter(grepl("R1", SUPERFAMILY)) %>% select(CHROM, POS, ID, ORDER, SUPERFAMILY, AF)
sfs.r1=as.data.frame(table(t.r1$AF))
colnames(sfs.r1)<-c("frequency","occurrence")
sfs.r1$ratio=(sfs.r1$occurrence*10000/sum(sfs.r1$occurrence))
p.r1 <- ggplot(sfs.r1, aes(x=frequency, y=ratio)) + 
  geom_bar(stat = "identity") + 
  ggtitle("r1, %")+ scale_color_grey() + theme_classic()+scale_y_continuous(trans = 'log2')
p.r1

t_sr = t1 %>% filter(ORDER=="Simple_repeat") %>% select(CHROM, POS, ID, ORDER, SUPERFAMILY, AF)
sfs_sr=as.data.frame(table(t_sr$AF))
colnames(sfs_sr)<-c("frequency","occurrence")
sfs_sr$ratio=(sfs_sr$occurrence*10000/sum(sfs_sr$occurrence))
p_sr <- ggplot(sfs_sr, aes(x=frequency, y=ratio)) + 
  geom_bar(stat = "identity") + 
  ggtitle("Simple repeat, %")+ scale_color_grey() + theme_classic()+scale_y_continuous(trans = 'log2')
p_sr
