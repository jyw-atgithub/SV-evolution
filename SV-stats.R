##Credit
#from"https://gist.github.com/eacooper400/c954c757db02eb6c4ccfae1aa090658c"
#from 

##global settings
setwd("/Users/Oscar/Desktop/plotting/SV-stats")
##packages
library("tidyverse")
library(devtools)

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
vcf <- read.vcf("all.consensus.vcf", header=TRUE, stringsAsFactors=FALSE)
##variables
win.size=1000000
total_chr <- levels(as.factor(vcf$CHROM))
num.samples=ncol(vcf)-9
container <- data.frame(matrix(ncol=11, nrow=0))
colnames(container) <- c("mut.type", "SV.type", "Chromosome","Start","End","Count","Nchr", "ThetaW","Pi" ,"varD" ,"TajimaD")
ncol(vcf)
#This works
#vcf.part = vcf %>% filter(CHROM == "2R")


container <- data.frame(matrix(ncol=11, nrow=0))
for (chr in total_chr){
  for (what.to.filter in c("INS","DEL","DUP","INV", "BND", "*")){
    #print(chr)
    ##Use the INFO Tag, not the ID. There are few incongruence though
    vcf.part = vcf %>% filter(CHROM == chr) #calculate the window first 
    my.win <- seq(min(vcf.part$POS), max(vcf.part$POS), by=win.size)
    #print(my.win)
    vcf.part = vcf %>% filter(CHROM == chr)%>% filter(grepl(what.to.filter,INFO)) 
    my.results=data.frame(mut.type="SV", SV.type=what.to.filter, Chromosome=chr, Start=my.win, End=(my.win+win.size), Count=rep(0, length(my.win)))
    print(nrow(my.results))
    
    for (i in 1:nrow(my.results)) {
      d=subset(vcf.part, (vcf.part$POS>=my.results$Start[i] & vcf.part$POS<my.results$End[i]))
      my.results$Count[i]=nrow(d)
    }
    my.results$Nchr=2*(num.samples)
    a=seq(from=1, to=((2*num.samples)-1), by=1) #the sequence of numbers from 1 to 2N-1,
    my.results$ThetaW = my.results$Count/(sum(1/a)) #calculate the Theta-W

    
    my.results$Pi=rep(0, nrow(my.results)) # Set up an empty column to hold the results of the function
    for (i in 1:nrow(my.results)) {
      d=subset(vcf.part, (vcf.part$POS>=my.results$Start[i] & vcf.part$POS<my.results$End[i])) # get the subset in the window
      if (nrow(d)==0){ #to prevent error if it is empty
        my.results$Pi[i]=0
      }else{
        j=apply(d, 1, FUN=derivedCount)
        c=rep(my.results$Nchr[i], length(j)) # a vector of my 2N values (same as c in Pi equation) that is the same length as my vector of j values
        my.results$Pi[i]=sum((2*j*(c-j))/(c*(c-1)))
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



p1 <-ggplot(data= container)
p1 + geom_boxplot(mapping=aes(x= SV.type, y= Pi/win.size, color=SV.type))

p2 <-ggplot(data= container)
p2 + geom_boxplot(mapping=aes(x= SV.type, y= Pi/win.size, color=Chromosome))


p3 <-ggplot(data= container)
p3 + geom_boxplot(mapping=aes(x= SV.type, y= TajimaD, color=SV.type))

p4 <-ggplot(data= container)
p4 + geom_boxplot(mapping=aes(x= SV.type, y= TajimaD, color=Chromosome))



