library("ggplot2")
library("tidyverse")

setwd("/Users/Oscar/Desktop/Emerson_Lab_Work/SV-project")
t1=read.table("repeat_type_genotype.tsv", header=TRUE)
t1=t1 %>%mutate(across(colnames(t1)[6:ncol(t1)], ~ recode(.x, './.' = "0", '0/1' = "1", '1/0' = "1", '1/1' = "1")))
t1$AF <- apply(t1, 1, function(x) length(which(x==1)))

t_del=t1 %>%filter(SVTYPE=="DEL") %>% count(ORDER,sort=TRUE) %>%filter(n>100)
temp=t1 %>%filter(SVTYPE=="DEL") %>% count(ORDER,sort=TRUE) %>%filter(n<100)
t_del=t_del%>% add_row(ORDER = "other", n =sum(temp$n) )
p=ggplot(data=t_del,aes(x="", y=n, fill=ORDER))+
  geom_bar(stat="identity", width=3, color="white")+
  coord_polar("y", start=0)
p

t_ins=t1 %>%filter(SVTYPE=="INS") %>% count(ORDER,sort=TRUE) %>%filter(n>100)
temp=t1 %>%filter(SVTYPE=="INS") %>% count(ORDER,sort=TRUE) %>%filter(n<100)
t_ins=t_ins%>% add_row(ORDER = "other", n =sum(temp$n) )
p_ins=ggplot(data=t_ins,aes(x="", y=n, fill=ORDER))+
  geom_bar(stat="identity", width=3, color="white")+
  coord_polar("y", start=0)
p_ins

t_inv=t1%>%filter(grepl("INV", ID)) %>% count(ORDER,sort=TRUE)
p_inv=ggplot(data=t_inv,aes(x="", y=n, fill=ORDER))+
  geom_bar(stat="identity", width=3, color="white")+
  coord_polar("y", start=0)
p_inv

t_dup=t1%>%filter(grepl("DUP", ID)) %>% count(ORDER,sort=TRUE)%>%filter(n>25)
temp=t1 %>%filter(grepl("DUP", ID)) %>% count(ORDER,sort=TRUE) %>%filter(n<25)
t_dup=t_dup%>% add_row(ORDER = "other", n =sum(temp$n) )
p_dup=ggplot(data=t_dup,aes(x="", y=n, fill=ORDER))+
  geom_bar(stat="identity", width=3, color="white")+
  coord_polar("y", start=0)
p_dup
