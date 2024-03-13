##Credit
#from "https://gist.github.com/eacooper400/c954c757db02eb6c4ccfae1aa090658c"
#from "https://eacooper400.github.io/gen8900/exercises/tajd.html"

## vk tajima & vcftools --TajimaD only works on SNPs
## vcftools --TajimaD Expected at least 2 parts in FORMAT entry: ID=DR,Number=2,Type=Integer

##packages
library("tidyverse")
library(devtools)
setwd("/Users/Oscar/Desktop/Emerson_Lab_Work/SV-project")

##Functions
read.vcf <- function(file, special.char="##", ...) {
  my.search.term=paste0(special.char, ".*")
  all.lines=readLines(file)
  clean.lines=gsub(my.search.term, "",  all.lines)
  clean.lines=gsub("#CHROM", "CHROM", clean.lines)
  read.table(..., text=paste(clean.lines, collapse="\n"))
}

get.field <- function(samples, format, fieldName) {
  x=strsplit(samples, split=":")
  fields=unlist(strsplit(format, split=":")) 
  i=which(fields==fieldName)
  if (!(fieldName %in% fields)) i=0
  #if (!(fieldName %in% fields)) stop('fieldName not found in format fields')
  return(sapply(x, `[[`, i)) 
}

count.genotypes <- function(genotypes) {
  genotypes = gsub("(\\||/)", "", genotypes) 
  gen.patterns = c("00", "01", "10", "11", "..") 
  my.counts=table(factor(genotypes, levels=gen.patterns)) 
  final.counts = c(my.counts[1], (my.counts[2]+my.counts[3]), my.counts[4:5]) 
  names(final.counts) = c("AA", "Aa", "aa", "NN") 
  return(final.counts)
}
derivedCount <- function(row) {
  row=as.vector(row, mode="character")
  x=count.genotypes(get.field(row[10:length(row)], row[9], "GT"))
  dc=(2*x["aa"])+x["Aa"]
  return(unname(dc))
}
variance.d <- function(n,S) {
  a1=sum(1/(seq(from=1, to=(n-1), by=1)))
  a2=sum(1/((seq(from=1, to=(n-1), by=1))**2))
  b1=(n+1)/(3*(n-1))
  b2=(2*((n**2)+n+3))/((9*n)*(n-1))
  c1=b1 - (1/a1)
  c2=b2-((n+2)/(a1*n)) + (a2/(a1**2))
  e1=c1/a1
  e2=c2/((a1**2)+a2)
  var=(e1*S) + (e2*S*(S-1))
  return(var)
}


#read the file
vcf <- read.vcf("3corrected.polarized.asm.vcf", header=TRUE, stringsAsFactors=FALSE)
## if the quality score is ".", this function will report error.
## sed s@"\t.\tPASS"@"\t10\tPASS"@g truvari.svimASM.vcf >truvari.svimASM.rn.vcf
#snp.vcf <- read.vcf("all.snps.vcf", header=TRUE, stringsAsFactors=FALSE)
#snp.outside.vcf <- read.vcf("syn.outside-1000-svimasm.snps.vcf", header=TRUE, stringsAsFactors=FALSE)
syn.snp.vcf <- read.vcf("synSNPs.vcf", header=TRUE, stringsAsFactors=FALSE)
#bed <- read.table("pure-outside-svimasm.bed", header = FALSE)

##variables
win.size=1000000
syn.win.size=3000
total_chr <- levels(as.factor(vcf$CHROM))
num.samples=ncol(vcf)-9
container <- data.frame(matrix(ncol=13, nrow=0))
colnames(container) <- c("mut.type", "SV.type", "Chromosome","Start","End","Count","Nchr", "ThetaW","Pi" ,"varD" ,"TajimaD")


#This works
#vcf.part = vcf %>% filter(CHROM == "2R")
ncol(vcf)
ncol(container)

for (chr in total_chr){
  for (what.to.filter in c("INS","DEL","DUP","INV", "TRA")){
    #print(chr)
    vcf.part = vcf %>% filter(CHROM == chr) #calculate the window first 
    my.win <- seq(min(vcf.part$POS), max(vcf.part$POS), by=win.size)
    #print(my.win)
    ##Use the INFO Tag, not the ID. There are few incongruence though
    vcf.part = vcf %>% filter(CHROM == chr)%>% filter(grepl(paste("SVTYPE=",what.to.filter,sep=""),INFO)) 
    my.results=data.frame(mut.type="SV", SV.type=what.to.filter, Chromosome=chr, Start=my.win, End=(my.win+win.size), Count=rep(0, length(my.win)))
    print(nrow(my.results))
    
    for (i in 1:nrow(my.results)) {
      d=subset(vcf.part, (vcf.part$POS>=my.results$Start[i] & vcf.part$POS<my.results$End[i]))
      my.results$Count[i]=nrow(d)
    }
    my.results$Nchr=2*(num.samples)
    a=seq(from=1, to=((2*num.samples)-1), by=1) #the sequence of numbers from 1 to 2N-1,
    my.results$ThetaW = my.results$Count/(sum(1/a)) #calculate the Theta-W
    my.results$ThetaW.persite = my.results$ThetaW/win.size #calculate the Theta-W, per bp
    
    my.results$Pi=rep(0, nrow(my.results)) # Set up an empty column to hold the results of the function
    for (i in 1:nrow(my.results)) {
      d=subset(vcf.part, (vcf.part$POS>=my.results$Start[i] & vcf.part$POS<my.results$End[i])) # get the subset in the window
      if (nrow(d)==0){ #to prevent error if it is empty
        my.results$Pi[i]=0
        my.results$Pi.persite[i]=0
      }else{
        j=apply(d, 1, FUN=derivedCount)
        c=rep(my.results$Nchr[i], length(j)) # a vector of my 2N values (same as c in Pi equation) that is the same length as my vector of j values
        my.results$Pi[i]=sum((2*j*(c-j))/(c*(c-1)))
        my.results$Pi.persite[i]=my.results$Pi[i]/win.size
      }
    }#calculate the pi
    
    my.results$varD=rep(0, nrow(my.results))
    for (i in 1:nrow(my.results)) { # Loop through the windows one more time
      my.results$varD[i] = variance.d(n=my.results$Nchr[i], S=my.results$Count[i])
    }
    my.results$TajimaD = (my.results$Pi - my.results$ThetaW)/(sqrt(my.results$varD))#calculate the Tajima's D
    container <- rbind(container, my.results)
  }
}

##Let's work on snps
## all SYN snps

for (chr in total_chr){
    #print(chr)
    vcf.part = syn.snp.vcf %>% filter(CHROM == chr) #calculate the window first 
    my.win <- seq(min(vcf.part$POS), max(vcf.part$POS), by=syn.win.size)
    #print(my.win)
    my.results=data.frame(mut.type="SNP", SV.type="all.syn.SNP", Chromosome=chr, Start=my.win, End=(my.win+syn.win.size), Count=rep(0, length(my.win)))
    print(nrow(my.results))
    
    for (i in 1:nrow(my.results)) {
      d=subset(vcf.part, (vcf.part$POS>=my.results$Start[i] & vcf.part$POS<my.results$End[i]))
      my.results$Count[i]=nrow(d)
    }
    my.results$Nchr=2*(num.samples)
    a=seq(from=1, to=((2*num.samples)-1), by=1) #the sequence of numbers from 1 to 2N-1,
    my.results$ThetaW = my.results$Count/(sum(1/a)) #calculate the Theta-W
    my.results$ThetaW.persite = my.results$ThetaW/syn.win.size #calculate the Theta-W
    
    my.results$Pi=rep(0, nrow(my.results)) # Set up an empty column to hold the results of the function
    for (i in 1:nrow(my.results)) {
      d=subset(vcf.part, (vcf.part$POS>=my.results$Start[i] & vcf.part$POS<my.results$End[i])) # get the subset in the window
      if (nrow(d)==0){ #to prevent error if it is empty
        my.results$Pi[i]=0
        my.results$Pi.persite[i]=0
      }else{
        j=apply(d, 1, FUN=derivedCount)
        c=rep(my.results$Nchr[i], length(j)) # a vector of my 2N values (same as c in Pi equation) that is the same length as my vector of j values
        my.results$Pi[i]=sum((2*j*(c-j))/(c*(c-1)))
        my.results$Pi.persite[i]=my.results$Pi[i]/syn.win.size
      }
    }#calculate the pi
    
    my.results$varD=rep(0, nrow(my.results))
    for (i in 1:nrow(my.results)) { # Loop through the windows one more time
      my.results$varD[i] = variance.d(n=my.results$Nchr[i], S=my.results$Count[i])
    }
    my.results$TajimaD = (my.results$Pi - my.results$ThetaW)/(sqrt(my.results$varD))#calculate the Tajima's D
    container <- rbind(container, my.results)
}

##Let's work on snps
## neighboring outside snps, in 1M windows
# for (chr in total_chr){
#   #print(chr)
#   vcf.part = snp.outside.vcf %>% filter(CHROM == chr) #calculate the window first 
#   my.win <- seq(min(vcf.part$POS), max(vcf.part$POS), by=win.size)
#   #print(my.win)
#   my.results=data.frame(mut.type="SNP", SV.type="outside.SNP", Chromosome=chr, Start=my.win, End=(my.win+win.size), Count=rep(0, length(my.win)))
#   print(nrow(my.results))
#   
#   for (i in 1:nrow(my.results)) {
#     d=subset(vcf.part, (vcf.part$POS>=my.results$Start[i] & vcf.part$POS<my.results$End[i]))
#     my.results$Count[i]=nrow(d)
#   }
#   my.results$Nchr=2*(num.samples)
#   a=seq(from=1, to=((2*num.samples)-1), by=1) #the sequence of numbers from 1 to 2N-1,
#   my.results$ThetaW = my.results$Count/(sum(1/a)) #calculate the Theta-W
#   
#   
#   my.results$Pi=rep(0, nrow(my.results)) # Set up an empty column to hold the results of the function
#   for (i in 1:nrow(my.results)) {
#     d=subset(vcf.part, (vcf.part$POS>=my.results$Start[i] & vcf.part$POS<my.results$End[i])) # get the subset in the window
#     if (nrow(d)==0){ #to prevent error if it is empty
#       my.results$Pi[i]=0
#     }else{
#       j=apply(d, 1, FUN=derivedCount)
#       c=rep(my.results$Nchr[i], length(j)) # a vector of my 2N values (same as c in Pi equation) that is the same length as my vector of j values
#       my.results$Pi[i]=sum((2*j*(c-j))/(c*(c-1)))
#     }
#   }#calculate the pi
#   
#   my.results$varD=rep(0, nrow(my.results))
#   for (i in 1:nrow(my.results)) { # Loop through the windows one more time
#     my.results$varD[i] = variance.d(n=my.results$Nchr[i], S=my.results$Count[i])
#   }
#   my.results$TajimaD = (my.results$Pi - my.results$ThetaW)/(sqrt(my.results$varD))#calculate the Tajima's D
#   container <- rbind(container, my.results)
# }

##Let's work on snps
## neighboring outside snps, NO windows !!!
# for (i in 1:nrow(bed)){
#   vcf.part = snp.outside.vcf %>% filter(CHROM == bed[i, 1])
#   my.results=data.frame(mut.type="SNP", SV.type="1k.SNP", Chromosome=bed[i, 1], Start=bed[i, 2], End=bed[i, 3], Count=0)
#   d=subset(vcf.part, (vcf.part$POS>=bed[i, 2] & vcf.part$POS<bed[i, 3]))
#   my.results$Count=nrow(d)
#   my.results$Nchr=2*(num.samples)
#   a=seq(from=1, to=((2*num.samples)-1), by=1)
#   #calculate the Theta-W
#   my.results$ThetaW = my.results$Count/(sum(1/a)) 
#   my.results$ThetaW.persite = my.results$ThetaW/(bed[i, 3]-bed[i, 2])
#   
#   my.results$Pi=0 # Set up an empty column to hold the results of the function
#   my.results$Pi.persite=0
#   if (nrow(d)==0){ #to prevent error if it is empty
#     my.results$Pi=0
#   }else{
#     j=apply(d, 1, FUN=derivedCount)
#     c=rep(my.results$Nchr, length(j)) # a vector of my 2N values (same as c in Pi equation) that is the same length as my vector of j values
#     my.results$Pi=sum((2*j*(c-j))/(c*(c-1)))
#     my.results$Pi.persite=my.results$Pi/(bed[i, 3]-bed[i, 2])
#   }
#   
#   my.results$varD=0
#   my.results$varD = variance.d(n=my.results$Nchr, S=my.results$Count)
#   my.results$TajimaD = (my.results$Pi - my.results$ThetaW)/(sqrt(my.results$varD))#calculate the Tajima's D
#   container <- rbind(container, my.results)
#   print(paste(round(100*i/nrow(bed),2),"%", sep = "", collapse=""))
# }

##Let's work on snps
## neighboring outside SYN snps, NO windows !!!
for (i in 1:nrow(bed)){
  vcf.part = syn.snp.vcf %>% filter(CHROM == bed[i, 1])
  my.results=data.frame(mut.type="SNP", SV.type="1k.syn.SNP", Chromosome=bed[i, 1], Start=bed[i, 2], End=bed[i, 3], Count=0)
  d=subset(vcf.part, (vcf.part$POS>=bed[i, 2] & vcf.part$POS<bed[i, 3]))
  my.results$Count=nrow(d)
  my.results$Nchr=2*(num.samples)
  a=seq(from=1, to=((2*num.samples)-1), by=1)
  #calculate the Theta-W
  my.results$ThetaW = my.results$Count/(sum(1/a)) 
  my.results$ThetaW.persite = my.results$ThetaW/(bed[i, 3]-bed[i, 2])
  
  my.results$Pi=0 # Set up an empty column to hold the results of the function
  if (nrow(d)==0){ #to prevent error if it is empty
    my.results$Pi=0
    my.results$Pi.persite=0
  }else{
    j=apply(d, 1, FUN=derivedCount)
    c=rep(my.results$Nchr, length(j)) # a vector of my 2N values (same as c in Pi equation) that is the same length as my vector of j values
    my.results$Pi=sum((2*j*(c-j))/(c*(c-1)))
    my.results$Pi.persite=my.results$Pi/(bed[i, 3]-bed[i, 2])
  }
  
  my.results$varD=0
  my.results$varD = variance.d(n=my.results$Nchr, S=my.results$Count)
  my.results$TajimaD = (my.results$Pi - my.results$ThetaW)/(sqrt(my.results$varD))#calculate the Tajima's D
  container <- rbind(container, my.results)
  print(paste(round(100*i/nrow(bed),2),"%", sep = "", collapse=""))
}

library(RColorBrewer)
library(khroma)

container$SV.type =gsub("all.syn.SNP","SNP",container$SV.type)
# levels(as.factor(container$SV.type))
# change the order of display
container$SV.type <- factor(container$SV.type , levels=c("SNP", "INS", "DEL", "DUP", "INV", "TRA"))
#x=reorder(SV.type, Pi.persite, FUN = mean, decreasing=TRUE)

##output dimensions=2400*2000
p1 <-ggplot(data= container)
p1 + geom_violin(mapping=aes(x= SV.type, 
                              y= Pi.persite, color=SV.type), lwd=2.8)  +
  geom_boxplot(mapping=aes(x= SV.type, y= Pi.persite), width=0.12)+scale_colour_bright() +
  scale_y_log10() + theme_light(base_size = 50) +
  #xlab("Type of variant") + 
  ylab("log(Pi-per-site)")+scale_colour_vibrant()+
  theme(legend.title = element_text(size=40), legend.text = element_text(size=48))

# Rename the column and the values in the factor
# http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/
#levels(container$SV.type)[levels(container$SV.type)=="all.syn.SNP"] <- "SYN*"
#levels(container$SV.type)[levels(container$SV.type)=="1k.syn.SNP"] <- "1k**"

p3 <-ggplot(data= container)
p3 + geom_violin(mapping=aes(x= SV.type, y= TajimaD, color=SV.type), lwd=2.8, trim=FALSE)+
  theme_bw(base_size = 50) +
  xlab("Type of variant") + ylab("Tajima's D") +
  geom_boxplot(mapping=aes(x= SV.type, y= TajimaD), width=0.12)+scale_colour_vibrant()+
  theme(legend.title = element_text(size=40), legend.text = element_text(size=48))

#write.csv(container, file="container.0304.csv")
#container=read.csv("container.1212.csv")
a = container %>% filter(SV.type=="all.syn.SNP")
sum(a$Count)
summary(a$ThetaW.persite)
summary(a$Pi.persite)
b = container %>% filter(SV.type=="INS")
res <- wilcox.test(na.omit(a$TajimaD), na.omit(b$TajimaD))
res











