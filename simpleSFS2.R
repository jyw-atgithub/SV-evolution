library("ggplot2")
library("tidyverse")

setwd("/Users/Oscar/Desktop/Emerson_Lab_Work/SV-project")
t1=read.table("repeat_type_genotype.tsv", header=TRUE)
t1=t1 %>%mutate(across(colnames(t1)[6:ncol(t1)], ~ recode(.x, './.' = "0", '0/1' = "1", '1/0' = "1", '1/1' = "1")))
t1$AF <- apply(t1, 1, function(x) length(which(x==1)))

sfs.SNP=read.table("MAF_syn.tsv", header=FALSE)
colnames(sfs.SNP)=c("occurrence","MAF")
sfs.SNP$type="SNP"
sfs.SNP$frequency=round(sfs.SNP$MAF*62*2)
sfs.SNP=sfs.SNP %>% select(-c("MAF"))
sfs.SNP$ratio=(sfs.SNP$occurrence*10000/sum(sfs.SNP$occurrence))

#levels(as.factor(t1$SUPERFAMILY))

#t.Gypsy = t1 %>% filter(grepl("Gypsy", ORDER)) %>% select(CHROM, POS, ID, ORDER, SUPERFAMILY, AF)
t.Gypsy = t1 %>% filter(grepl("Gypsy", SUPERFAMILY)) %>% select(CHROM, POS, ID, ORDER, SUPERFAMILY, AF)
sfs.Gypsy=as.data.frame(table(t.Gypsy$AF))
colnames(sfs.Gypsy)<-c("frequency","occurrence")
sfs.Gypsy$ratio=(sfs.Gypsy$occurrence*10000/sum(sfs.Gypsy$occurrence))
sfs.Gypsy$type="Gypsy"
# p.Gypsy <- ggplot(sfs.Gypsy, aes(x=frequency, y=ratio)) + 
#   geom_bar(stat = "identity") + 
#   ggtitle("Gypsy, %")+ scale_color_grey() + theme_classic()+scale_y_continuous(trans = 'log2')
# p.Gypsy

t.pao = t1 %>% filter(grepl("Pao", SUPERFAMILY)) %>% select(CHROM, POS, ID, ORDER, SUPERFAMILY, AF)
sfs.pao=as.data.frame(table(t.pao$AF))
colnames(sfs.pao)<-c("frequency","occurrence")
sfs.pao$ratio=(sfs.pao$occurrence*10000/sum(sfs.pao$occurrence))
sfs.pao$type="pao"

t.Helitron = t1 %>% filter(grepl("Helitron", SUPERFAMILY)) %>% select(CHROM, POS, ID, ORDER, SUPERFAMILY, AF)
sfs.Helitron=as.data.frame(table(t.Helitron$AF))
colnames(sfs.Helitron)<-c("frequency","occurrence")
sfs.Helitron$ratio=(sfs.Helitron$occurrence*10000/sum(sfs.Helitron$occurrence))
sfs.Helitron$type="Helitron"

t.r1 = t1 %>% filter(grepl("R1", SUPERFAMILY)) %>% select(CHROM, POS, ID, ORDER, SUPERFAMILY, AF)
sfs.r1=as.data.frame(table(t.r1$AF))
colnames(sfs.r1)<-c("frequency","occurrence")
sfs.r1$ratio=(sfs.r1$occurrence*10000/sum(sfs.r1$occurrence))
sfs.r1$type="R1"

t.Jockey = t1 %>% filter(grepl("Jockey", SUPERFAMILY)) %>% select(CHROM, POS, ID, ORDER, SUPERFAMILY, AF)
sfs.Jockey=as.data.frame(table(t.Jockey$AF))
colnames(sfs.Jockey)<-c("frequency","occurrence")
sfs.Jockey$ratio=(sfs.Jockey$occurrence*10000/sum(sfs.Jockey$occurrence))
sfs.Jockey$type="Jockey"

sfs.all=rbind(sfs.SNP,sfs.Gypsy,sfs.pao,sfs.Helitron,sfs.r1,sfs.Jockey)
sfs.all$frequency=as.character(sfs.all$frequency)
sfs.all$frequency=factor(sfs.all$frequency,levels=unique(sfs.all$frequency))
p.all.log <- ggplot(sfs.all, aes(x=frequency, y=ratio, fill=type)) + 
  geom_bar(stat = "identity",position="dodge") +
  ggtitle("all, %")+ scale_color_grey() + theme_classic()+scale_y_continuous(trans = 'log2')+
  theme(legend.title = element_text(size=30), legend.text = element_text(size=30))+
  theme(legend.key.size = unit(1, 'cm'))
p.all.log

p.all <- ggplot(sfs.all, aes(x=frequency, y=ratio, fill=type)) + 
  geom_bar(stat = "identity",position="dodge") +
  ggtitle("Relative frequency by TE type")+ 
  scale_color_grey() + theme_classic()+
  theme(legend.title = element_text(size=30), legend.text = element_text(size=30))+
  theme(legend.key.size = unit(1, 'cm'))
p.all


t_sr = t1 %>% filter(ORDER=="Simple_repeat") %>% select(CHROM, POS, ID, ORDER, SUPERFAMILY, AF)
sfs_sr=as.data.frame(table(t_sr$AF))
colnames(sfs_sr)<-c("frequency","occurrence")
sfs_sr$ratio=(sfs_sr$occurrence*10000/sum(sfs_sr$occurrence))
sfs_sr$type="Simple_repeat"

t.nr = t1 %>% filter(grepl("not_repeat", SUPERFAMILY)) %>% select(CHROM, POS, ID, ORDER, SUPERFAMILY, AF)
sfs.nr=as.data.frame(table(t.nr$AF))
colnames(sfs.nr)<-c("frequency","occurrence")
sfs.nr$ratio=(sfs.nr$occurrence*10000/sum(sfs.nr$occurrence))
sfs.nr$type="not_repeat"

t.lc = t1 %>% filter(grepl("Low_complexity", SUPERFAMILY)) %>% select(CHROM, POS, ID, ORDER, SUPERFAMILY, AF)
sfs.lc=as.data.frame(table(t.lc$AF))
colnames(sfs.lc)<-c("frequency","occurrence")
sfs.lc$ratio=(sfs.lc$occurrence*10000/sum(sfs.lc$occurrence))
sfs.lc$type="Low_complexity"

sfs.nonTE=rbind(sfs.SNP,sfs_sr,sfs.nr,sfs.lc)
sfs.nonTE$frequency=as.character(sfs.nonTE$frequency)
sfs.nonTE$frequency=factor(sfs.nonTE$frequency,levels=unique(sfs.nonTE$frequency))

p.nonTE <- ggplot(sfs.nonTE, aes(x=frequency, y=ratio, fill=type)) + 
  geom_bar(stat = "identity",position="dodge") +
  ggtitle("Relative frequency of non-TE part")+ 
  scale_color_grey() + theme_classic()+scale_y_continuous(trans = 'log2')+
  theme(legend.title = element_text(size=30), legend.text = element_text(size=30))+
  theme(legend.key.size = unit(1, 'cm'))
p.nonTE
