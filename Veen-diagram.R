library(tidyverse)
library(ggplot2)
library(VennDiagram)

setwd("/Users/Oscar/Desktop/Emerson_Lab_Work/SV-project")
tab <- read.table("stats-nv107.tsv", header = FALSE)
colnames(tab) <- c("prog", "chr", "start", "type", "length")
levels(as.factor(tab$prog))
levels(as.factor(tab$type))

tab = tab %>% mutate(abs_length=abs(as.numeric(tab$length))) 
# %>% mutate(new_type=recode(type,"DUP:INT"="DUP", "DUP:TANDEM"="DUP"))
#not required

png("2_1-sniffles_count.png", width = 1200, height = 1200)
p1 <- ggplot(data=tab%>%filter(prog=="sniffles")) + geom_bar(mapping=aes(x=type), fill = 'orange3')+
  scale_y_sqrt()+labs(title="sniffles", x ="SV type", y = "Count")+
  theme(text = element_text(size = 30)) + scale_y_continuous(limits = c(0,6000)) 
p1 
dev.off()

png("2_2-cuteSV_count.png", width = 1200, height = 1200)
p2 <- ggplot(data=tab%>%filter(prog=="cuteSV")) + geom_bar(mapping=aes(x=type), fill = 'green4')+
  scale_y_sqrt() +labs(title="cuteSV", x ="SV type", y = "Count")+
  theme(text = element_text(size = 30)) + scale_y_continuous(limits = c(0,6000)) 
p2
dev.off()

library(gridExtra)
library(grid)
lay <- rbind(c(1,2),c(1,2))
grid.arrange(p1, p2, ncol=2)
grid.arrange(p1,p2, layout_matrix = lay)


tab2 = tab %>% na.omit()
#%>% filter(abs_length >49)
tab3 = tab %>% filter(prog=="sniffles")
levels(as.factor(tab3$prog))
levels(as.factor(tab3$type))
summary(tab3$abs_length)
sum(tab3$abs_length)


png("1-mapping-length-hist.png", width = 6000, height = 2500)
p<-ggplot(data=tab2) +
  geom_histogram(mapping = (aes(x=abs_length, fill=prog)) , bins = 500)+
  scale_x_log10(breaks=c(1,10,30,100,300, 1000,3000, 10000, 1e+5, 1e+6, 5e+6))+
  scale_y_log10()+
  scale_fill_brewer(palette="Accent")+
  facet_grid(prog ~ .)+
  theme(text = element_text(size = 50)) 
p
dev.off()

t005 <- read.table("sample_merged_overlapp.005.txt")
venn.diagram(list(cuteSV=which(t005[,1]==1), sniffles=which(t005[,2]==1), SVIM=which(t005[,3]==1)) , 
             euler.d=TRUE, scaled=TRUE, 
             fill = c("gray", "orange" ,"blue") , alpha = c(0.5, 0.5, 0.5), 
             cex = 2, lty =0,
             imagetype = "png", filename = "my_sample_overlapp.005.png", disable.logging = TRUE)

t01 <- read.table("sample_merged_overlapp.01.txt")
venn.diagram(list(cuteSV=which(t01[,1]==1), sniffles=which(t01[,2]==1), SVIM=which(t01[,3]==1)) , 
             euler.d=TRUE, scaled=TRUE, 
             fill = c("gray", "orange" ,"blue") , alpha = c(0.5, 0.5, 0.5), 
             cex = 2, lty =0,
             imagetype = "png", filename = "my_sample_overlapp.01.png", disable.logging = TRUE)

t1000 <- read.table("sample_merged_overlapp.1000.txt")
venn.diagram(list(cuteSV=which(t1000[,1]==1), sniffles=which(t1000[,2]==1), SVIM=which(t1000[,3]==1)) , 
             euler.d=TRUE, scaled=TRUE, 
             fill = c("gray", "orange" ,"blue") , alpha = c(0.5, 0.5, 0.5), 
             cex = 2, lty =0,
             imagetype = "png", filename = "my_sample_overlapp.1000.png", disable.logging = TRUE)

