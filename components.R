library("ggplot2")
library("tidyverse")
library(RColorBrewer)
library(khroma)

setwd("/Users/Oscar/Desktop/Emerson_Lab_Work/SV-project")
t1=read.table("repeat_type_genotype.tsv", header=TRUE)
t1=t1 %>%mutate(across(colnames(t1)[6:ncol(t1)], ~ recode(.x, './.' = "0", '0/1' = "1", '1/0' = "1", '1/1' = "1")))
t1$AF <- apply(t1, 1, function(x) length(which(x==1)))

t_del=t1 %>%filter(SVTYPE=="DEL") %>% count(ORDER,sort=TRUE) %>%filter(n>100)
temp=t1 %>%filter(SVTYPE=="DEL") %>% count(ORDER,sort=TRUE) %>%filter(n<100)
t_del=t_del%>% add_row(ORDER = "other", n =sum(temp$n) ) %>% add_column(TYPE="DEL") 
t_del=t_del%>% mutate(PERCENT=n/sum(t_del$n))

t_ins=t1 %>%filter(SVTYPE=="INS") %>% count(ORDER,sort=TRUE) %>%filter(n>200)
temp=t1 %>%filter(SVTYPE=="INS") %>% count(ORDER,sort=TRUE) %>%filter(n<200)
t_ins=t_ins%>% add_row(ORDER = "other", n =sum(temp$n) )%>% add_column(TYPE="INS")
t_ins=t_ins %>% mutate(PERCENT=n/sum(t_ins$n))

t_dup=t1 %>%filter(SVTYPE=="DUP") %>% count(ORDER,sort=TRUE) %>%filter(n>80)
temp=t1 %>%filter(SVTYPE=="DUP") %>% count(ORDER,sort=TRUE) %>%filter(n<80)
t_dup=t_dup%>% add_row(ORDER = "other", n =sum(temp$n) )%>% add_column(TYPE="DUP")
t_dup=t_dup %>% mutate(PERCENT=n/sum(t_dup$n))

t_all=bind_rows(t_del, t_ins, t_dup) 
t_all= t_all %>% mutate(ORDER= case_match(ORDER, "Simple_repeat,Simple_repeat" ~ "Simple_repeat", "LTR,Simple_repeat,Simple_repeat" ~ "LTR,Simple_repeat",.default = ORDER)) %>%
  group_by(ORDER, TYPE) %>% summarise(n=sum(n), PERCENT=sum(PERCENT)) %>%
  mutate(ORDER=factor(ORDER, levels=c("LINE","LINE,LTR","LINE,Simple_repeat","LTR","LTR,Simple_repeat","DNA","DNA,Simple_repeat","RC","Satellite","Satellite,Unknown","Simple_repeat","Low_complexity","Low_complexity,Simple_repeat","not_repeat","other") ))

levels(as.factor(t_all$ORDER))
p0=ggplot(data=t_all,aes(x=ORDER, y=PERCENT*100))
#Each position adjustment can be recast as a function with manual width and height arguments:
p0=p0+geom_col(aes(fill=TYPE), position=position_dodge2(preserve = "single"), width=0.5) +
  labs(x="Type of repeat", y="percentage, %") + guides(x=guide_axis(angle=45) ) + 
  theme_minimal(base_size = 24 ) +scale_fill_vibrant(name="SV types") +
  theme(legend.text = element_text(size = 40),legend.title = element_text(size = 36))
p0

# 
# t_del=t1 %>%filter(SVTYPE=="DEL") %>% count(ORDER,sort=TRUE) %>%filter(n>100)
# temp=t1 %>%filter(SVTYPE=="DEL") %>% count(ORDER,sort=TRUE) %>%filter(n<100)
# t_del=t_del%>% add_row(ORDER = "other", n =sum(temp$n) )
# #Making donut chart
# #reference: https://r-graph-gallery.com/128-ring-or-donut-plot.html
# levels(as.factor( t_del$ORDER))
#                                     
# t_del$ORDER <- factor(t_del$ORDER, levels=c("DNA", "LINE","LINE,LTR","LINE,Simple_repeat","LTR","LTR,Simple_repeat","RC","Simple_repeat","Simple_repeat,Simple_repeat","Satellite" ,"Satellite,Unknown","Low_complexity","other" ))
# t_del$fraction=t_del$n/sum(t_del$n)
# t_del$ymax = cumsum(t_del$fraction)
# t_del$ymin = c(0, head(t_del$ymax, n=-1))
# t_del$labelPosition <- (t_del$ymax + t_del$ymin) / 2
# #t_del$label<-c("LTR","simple_repeat","LINE","Satellite","DNA","Helitron","","","","","","","")
# cbPalette <- c("#999999","#b8e186", "#c51b7d", "#E69F00","#de77ae", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#F0E442", "#8e0152","#7fbc41","#000000", "#9999CC","#e6f5d0", "#4d9221","#f1b6da","#276419", "#f7f7f7")
# 
# p2=ggplot(t_del, aes(y=fraction, ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=ORDER)) +
#   geom_rect(color="black",size=2) +
#   coord_polar(theta="y")+
#   geom_text( x=3.5, aes(y=labelPosition, label=paste(round(fraction,4)*100, "%")), 
#              size=10, color="white") +
#   xlim(c(1.5, 4))+ scale_fill_manual(values=cbPalette)+
#   theme_void(base_size = 45)+
#   theme(legend.text = element_text(size = 40),legend.title = element_text(size = 36))
# p2
# 
# 
# t_ins=t1 %>%filter(SVTYPE=="INS") %>% count(ORDER,sort=TRUE) %>%filter(n>240)
# temp=t1 %>%filter(SVTYPE=="INS") %>% count(ORDER,sort=TRUE) %>%filter(n<240)
# t_ins=t_ins%>% add_row(ORDER = "other", n =sum(temp$n) )
# t_ins$fraction=t_ins$n/sum(t_ins$n)
# t_ins$ymax = cumsum(t_ins$fraction)
# t_ins$ymin = c(0, head(t_ins$ymax, n=-1))
# t_ins$labelPosition <- (t_ins$ymax + t_ins$ymin) / 2
# levels(as.factor( t_ins$ORDER))
# # [1] "DNA" "#999999"                          "DNA,Simple_repeat" "#f1b6da"             
# # [3] "LINE" "#b8e186"                       "LINE,Simple_repeat""#E69F00"             
# # [5] "Low_complexity"  "#7fbc41"         "Low_complexity,Simple_repeat"   "#276419"
# # [7] "LTR"  "#de77ae"                   "LTR,Simple_repeat"  "#56B4E9"            
# # [9] "LTR,Simple_repeat,Simple_repeat""#9999CC" "other" "#000000"                         
# # [11] "Satellite"     "#F0E442"          "Satellite,Unknown"  "#8e0152"            
# # [13] "Simple_repeat"   "#D55E00"
# t_ins$ORDER <- factor(t_ins$ORDER, levels=c("DNA","DNA,Simple_repeat","LINE","LINE,Simple_repeat","Low_complexity","Low_complexity,Simple_repeat","LTR","LTR,Simple_repeat","LTR,Simple_repeat,Simple_repeat","other","Satellite","Satellite,Unknown","Simple_repeat"))
# cbPalette <- c("#999999","#f1b6da","#b8e186","#E69F00","#7fbc41","#276419","#de77ae","#56B4E9","#9999CC", "#000000","#F0E442","#8e0152","#D55E00")
# 
# p_ins=ggplot(t_ins, aes(y=fraction, ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=ORDER)) +
#   geom_rect(color="black",size=2) +
#   coord_polar(theta="y")+
#   geom_text( x=3.5, aes(y=labelPosition, label=paste(round(fraction,4)*100, "%")), 
#              size=10, color="white") +
#   xlim(c(1.5, 4))+ scale_fill_manual(values=cbPalette)+theme_void(base_size = 40)+
#   theme(legend.text = element_text(size = 30),legend.title = element_text(size = 30))
# p_ins
# 
# t_inv=t1%>%filter(grepl("INV", ID)) %>% count(ORDER,sort=TRUE)
# p_inv=ggplot(data=t_inv,aes(x="", y=n, fill=ORDER))+
#   geom_bar(stat="identity", width=3, color="white")+
#   coord_polar("y", start=0)
# p_inv
# 
# t_dup=t1%>%filter(grepl("DUP", ID)) %>% count(ORDER,sort=TRUE)%>%filter(n>25)
# temp=t1 %>%filter(grepl("DUP", ID)) %>% count(ORDER,sort=TRUE) %>%filter(n<25)
# t_dup=t_dup%>% add_row(ORDER = "other", n =sum(temp$n) )
# p_dup=ggplot(data=t_dup,aes(x="", y=n, fill=ORDER))+
#   geom_bar(stat="identity", width=3, color="white")+
#   coord_polar("y", start=0)
# p_dup
