##library("ggplot2")
library("tidyverse")
library(RColorBrewer)
library("khroma")
setwd("/Users/Oscar/Desktop/Emerson_Lab_Work/SV-project/iso1_temp")

reg.40k.winn=read.table("regular_40k.winn.read_info.tsv", header = FALSE)
colnames(reg.40k.winn) = c("reference_name","flag", "leftmost", "rightmost", "mapQ", "read_length")
adaptive.8k.winn=read.table("adaptive_8k.winn.read_info.tsv")
colnames(adaptive.8k.winn) = c("reference_name","flag", "leftmost", "rightmost", "mapQ", "read_length")

reg.40k.ONT_19.winn=read.table("regular_40k.ONT_19.winn.read_info.tsv")
colnames(reg.40k.ONT_19.winn)== c("reference_name","flag", "leftmost", "rightmost", "mapQ", "read_length")
adaptive.8k.ONT_19.winn=read.table("adaptive_8k.ONT_19.winn.read_info.tsv")
colnames(adaptive.8k.ONT_19.winn)=c("reference_name","flag", "leftmost", "rightmost", "mapQ", "read_length")


t1=reg.40k.winn %>% filter(reference_name=="chromosome_Y") %>%
  pivot_longer(cols = c(leftmost, rightmost), names_to = "side", values_to = "POS")
  
p1=ggplot(data= t1)+
  geom_histogram(bins = round(nrow(t1)/25), aes(x=POS, fill=side),position="dodge")+
  theme_bw() +
  scale_fill_highcontrast()+
  scale_y_continuous(
    trans = "log2",
    labels = scales::math_format(2^.y, format = log2)
  )

p1


