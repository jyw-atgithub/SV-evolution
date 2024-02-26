library("ggplot2")
library("tidyverse")

setwd("/Users/Oscar/Desktop/Emerson_Lab_Work/SV-project")
t1=read.table("extraction.polarized.asm.tsv",header=FALSE)
t2=read.table("extraction.syn-snp.tsv",header=FALSE)
t3=read.table("extraction.syn-snp.outside-1000-svimasm.tsv",header=FALSE)


#colnames(t1) <- c("CHROM","POS", "SVTYPE","SVLEN", S[1:55])
#t1=t1 %>% mutate(V5 = str_replace(V5, './.', '0'))
#t1 %>%  mutate(conf = recode(conf, 'East' = 'E', 'West' = 'W', 'North' = 'N'))
nc=ncol(t1)
t1=t1 %>%mutate(across(V5:V65, ~ recode(.x, './.' = "0", '0/1' = "1", '1/0' = "1", '1/1' = "1")))
t2=t2 %>%mutate(across(V5:V65, ~ recode(.x, './.' = "0", '0/1' = "1", '1/0' = "1", '1/1' = "1")))
t3=t3 %>%mutate(across(V5:V65, ~ recode(.x, './.' = "0", '0/1' = "1", '1/0' = "1", '1/1' = "1")))
t1$V66 <- apply(t1, 1, function(x) length(which(x==1)))
t2$V66 <- apply(t2, 1, function(x) length(which(x==1)))
t3$V66 <- apply(t3, 1, function(x) length(which(x==1)))
levels(as.factor(t1$V3))


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


#sfs.all=bind_rows(sfsdel,sfs.ins,sfs.dup,sfs.inv,sfs.bnd)
sfs.all=bind_rows(sfsdel,sfs.ins,sfs.inv)
pall <- ggplot(sfs.all, aes(x=frequency, y=percentage, fill=type)) + 
  geom_bar(stat = "identity",position="dodge")  +
  ggtitle("all, %")+ scale_color_grey() + theme_classic()
pall

sfs.rare=bind_rows(sfsdel,sfs.dup,sfs.inv,sfs.bnd)
p.rare <- ggplot(sfs.rare, aes(x=frequency, y=percentage, fill=type)) + 
  geom_bar(stat = "identity",position="dodge") + 
  ggtitle("all, %")+ scale_color_grey() + theme_classic()
p.rare

t.synsnp = t2 %>% select(V1,V2,V3,V4,V66) 
sfs.synsnp=as.data.frame(table(t.synsnp$V66))
colnames(sfs.synsnp)<-c("frequency","occurrence")
sfs.synsnp$percentage=(sfs.synsnp$occurrence*100/nrow(t.synsnp))
sfs.synsnp$type="SNP"

t.1000snp = t3 %>% select(V1,V2,V3,V4,V66) 
sfs.1000snp=as.data.frame(table(t.1000snp$V66))
colnames(sfs.1000snp)<-c("frequency","occurrence")
sfs.1000snp$percentage=(sfs.1000snp$occurrence*100/nrow(t.1000snp))
sfs.1000snp$type="1000SNP"

sfs.allplus=bind_rows(sfsdel,sfs.ins,sfs.inv, sfs.synsnp, sfs.1000snp)
order_of_variants=c("1000SNP","SNP","INS","DEL","INV")
p.allplus <- ggplot(sfs.allplus, aes(x=frequency, y=percentage, fill=type)) + 
  geom_bar(stat = "identity",position = position_dodge2(preserve = "single")) +
  ggtitle("all, %")+ scale_color_grey() + 
  theme_classic() + 
  theme(legend.title = element_text(size=30), legend.text = element_text(size=30))+
  theme(legend.key.size = unit(1, 'cm')) +
  scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442"))
p.allplus

sfs.rare=bind_rows(sfsdel,sfs.dup,sfs.inv,sfs.bnd, sfs.synsnp)
p.rare <- ggplot(sfs.rare, aes(x=frequency, y=percentage, fill=type)) + 
  geom_bar(stat = "identity",position="dodge") +
  ggtitle("all, %")+ scale_color_grey() + theme_classic()
p.rare
